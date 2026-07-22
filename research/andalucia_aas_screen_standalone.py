#!/usr/bin/env python3
"""Whole-mitoproteome screen for Aas-specific modules in Andalucia godoyi.

The script is self-contained. It extracts all proteins from the published Table
S1 workbook, builds eight bacterial class profiles, and asks whether any
mitochondrial protein prefers an Aas-derived module to closely related enzyme
families. Provisional candidates are tested with exact-composition shuffles.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import random
import re
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np
import openpyxl
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "andalucia_aas_specific_screen"
WORK = OUT / "work"
MODELS = OUT / "models"
for directory in (OUT, WORK, MODELS):
    directory.mkdir(parents=True, exist_ok=True)

XLSX = ROOT / "inputs" / "Andalucia_TableS1.xlsx"
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO*]{40,}$", re.I)
ID_RE = re.compile(r"(?:ANDGO_)?_?(\d{5})\.mRNA\.1", re.I)
HX4D_RE = re.compile(r"H....D")
N_CLASSES = ["AasN", "LpaT", "PlsC", "PlsB"]
C_CLASSES = ["AasC", "split_AasC", "FadD", "Acs"]
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "proto-mito-Aas-specific-screen/1.0"})


def run(command: list[str], *, stdout: Path | None = None, quiet: bool = False) -> None:
    if stdout is None:
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.DEVNULL if quiet else None,
            stderr=subprocess.DEVNULL if quiet else None,
        )
    else:
        with stdout.open("w", encoding="utf-8") as handle:
            subprocess.run(command, check=True, stdout=handle, stderr=subprocess.DEVNULL if quiet else None)


def clean_sequence(value: object) -> str:
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(value or "").upper())


def get_text(url: str, params: dict[str, Any], attempts: int = 8) -> str:
    delay = 1.5
    last = ""
    for _ in range(attempts):
        try:
            response = SESSION.get(url, params=params, timeout=300)
            last = f"{response.status_code} {response.text[:300]}"
            if response.status_code == 200:
                return response.text
            if response.status_code not in {429, 500, 502, 503, 504}:
                break
        except Exception as exc:  # noqa: BLE001
            last = repr(exc)
        time.sleep(delay)
        delay = min(delay * 1.8, 25)
    raise RuntimeError(f"UniProt request failed: {last}")


def fetch_uniprot(query: str, size: int) -> pd.DataFrame:
    fields = "accession,id,protein_name,gene_names,organism_name,organism_id,length,sequence,reviewed,proteome,lineage"
    text = get_text(
        "https://rest.uniprot.org/uniprotkb/search",
        {"query": query, "format": "tsv", "fields": fields, "size": str(size)},
    )
    frame = pd.read_csv(io.StringIO(text), sep="\t")
    rename = {
        "Entry": "accession", "Entry Name": "entry_name", "Protein names": "protein_name",
        "Gene Names": "gene", "Organism": "organism", "Organism (ID)": "taxid",
        "Length": "length", "Sequence": "sequence", "Reviewed": "reviewed",
        "Proteomes": "proteomes", "Taxonomic lineage": "lineage",
    }
    frame = frame.rename(columns=rename)
    for column in rename.values():
        if column not in frame:
            frame[column] = ""
    frame["sequence"] = frame["sequence"].map(clean_sequence)
    frame = frame[frame["sequence"].str.len() >= 40].drop_duplicates("accession").reset_index(drop=True)
    frame["genus"] = frame["organism"].astype(str).str.replace(r"^\[|\]$", "", regex=True).str.split().str[0]
    return frame


def first_hx4d(sequence: str) -> int | None:
    match = HX4D_RE.search(sequence)
    return match.start() if match else None


def motif_window(sequence: str, before: int = 35, after: int = 155) -> str:
    position = first_hx4d(sequence)
    if position is None:
        return sequence[:220]
    return sequence[max(0, position - before):min(len(sequence), position + after)]


def one_per_genus(frame: pd.DataFrame, maximum: int) -> pd.DataFrame:
    if frame.empty:
        return frame
    frame = frame.sort_values(["genus", "accession"]).drop_duplicates("genus")
    if len(frame) <= maximum:
        return frame.reset_index(drop=True)
    indices = np.linspace(0, len(frame) - 1, maximum, dtype=int)
    return frame.iloc[indices].reset_index(drop=True)


def class_panels() -> dict[str, pd.DataFrame]:
    fused = fetch_uniprot('protein_name:"Bifunctional protein Aas" AND taxonomy_id:2', 300)
    fused = fused[fused["sequence"].str.len().between(650, 820)].copy()
    fused_n = fused.copy()
    fused_n["domain_sequence"] = fused_n["sequence"].str.slice(0, 180)
    fused_c = fused.copy()
    fused_c["domain_sequence"] = fused_c["sequence"].str.slice(190)
    panels: dict[str, pd.DataFrame] = {
        "AasN": one_per_genus(fused_n, 70),
        "AasC": one_per_genus(fused_c, 70),
    }
    queries = {
        "LpaT": '(gene_exact:lpaT OR protein_name:"lysophospholipid acyltransferase") AND taxonomy_id:2',
        "PlsC": 'gene_exact:plsC AND taxonomy_id:2',
        "PlsB": 'gene_exact:plsB AND taxonomy_id:2',
        "split_AasC": '(gene_exact:aasC OR protein_name:"acyl-acyl carrier protein synthetase") AND taxonomy_id:2',
        "FadD": 'gene_exact:fadD AND taxonomy_id:2',
        "Acs": '(gene_exact:acs OR protein_name:"acetyl-CoA synthetase") AND taxonomy_id:2',
    }
    for label, query in queries.items():
        frame = fetch_uniprot(query, 1000)
        names = frame["protein_name"].astype(str).str.lower()
        if label == "LpaT":
            frame = frame[names.str.contains("lysophospholipid|acyltransferase", regex=True)]
            frame = frame[frame["sequence"].str.len().between(120, 450)]
        elif label == "PlsC":
            frame = frame[frame["sequence"].str.len().between(150, 450)]
        elif label == "PlsB":
            frame = frame[frame["sequence"].str.len().between(600, 1000)]
        elif label == "split_AasC":
            frame = frame[~names.str.contains("bifunctional protein aas", regex=False)]
            frame = frame[frame["sequence"].str.len().between(350, 700)]
        elif label == "FadD":
            frame = frame[frame["sequence"].str.len().between(450, 700)]
        elif label == "Acs":
            frame = frame[names.str.contains("acetyl-coa synthetase", regex=False)]
            frame = frame[frame["sequence"].str.len().between(500, 800)]
        frame = frame.copy()
        if label in N_CLASSES:
            frame["domain_sequence"] = frame["sequence"].map(motif_window)
            frame = frame[frame["domain_sequence"].str.contains(r"H....D", regex=True)]
        else:
            frame["domain_sequence"] = frame["sequence"]
        panels[label] = one_per_genus(frame, 70)
    sizes = {label: len(frame) for label, frame in panels.items()}
    (OUT / "class_panel_sizes.json").write_text(json.dumps(sizes, indent=2), encoding="utf-8")
    if any(size < 5 for size in sizes.values()):
        raise RuntimeError(f"Insufficient class panel: {sizes}")
    rows = []
    for label, frame in panels.items():
        for row in frame.itertuples():
            rows.append({
                "class": label, "accession": row.accession, "organism": row.organism,
                "protein_name": row.protein_name, "full_length": len(row.sequence),
                "domain_length": len(row.domain_sequence), "domain_sequence": row.domain_sequence,
            })
    pd.DataFrame(rows).to_csv(OUT / "class_sequence_panels.tsv.gz", sep="\t", index=False, compression="gzip")
    return panels


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 100):
                handle.write(sequence[start:start + 100] + "\n")


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def build_hmm(label: str, frame: pd.DataFrame) -> Path:
    prefix = MODELS / label
    fasta = Path(str(prefix) + ".faa")
    alignment = Path(str(prefix) + ".aln.faa")
    hmm = Path(str(prefix) + ".hmm")
    write_fasta(fasta, [(f"{label}|{row.accession}", str(row.domain_sequence)) for row in frame.itertuples()])
    run(["mafft", "--auto", "--quiet", str(fasta)], stdout=alignment)
    run(["hmmbuild", "--amino", "-n", label, str(hmm), str(alignment)], quiet=True)
    return hmm


def parse_tbl(path: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 6:
                result[fields[0]] = float(fields[5])
    return result


def score_hmm(hmm: Path, fasta: Path, tag: str) -> dict[str, float]:
    tbl = WORK / f"{tag}.tbl"
    run([
        "hmmsearch", "--max", "--noali", "-E", "1000000", "--tblout", str(tbl),
        str(hmm), str(fasta),
    ], quiet=True)
    return parse_tbl(tbl)


def extract_mitoproteome() -> pd.DataFrame:
    workbook = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for sheet in workbook.worksheets:
        title = sheet.title
        if title.lower().startswith("key") or "stat" in title.lower():
            continue
        category = title.strip()[:1]
        for row_number, values in enumerate(sheet.iter_rows(values_only=True), start=1):
            cells = ["" if value is None else str(value).strip() for value in values]
            sequences = []
            for cell in cells:
                compact = re.sub(r"\s+", "", cell).upper().rstrip("*")
                if AA_RE.fullmatch(compact) and not set(compact) <= {"A", "C", "G", "T", "N"}:
                    sequences.append(compact)
            if not sequences:
                continue
            sequence = max(sequences, key=len)
            identifier = ""
            for cell in cells:
                match = ID_RE.search(cell)
                if match:
                    identifier = f"ANDGO_{match.group(1)}"
                    break
            if not identifier:
                identifier = f"{category}_ROW_{row_number:04d}"
            if (identifier, sequence) in seen:
                continue
            seen.add((identifier, sequence))
            annotation = []
            for cell in cells:
                compact = re.sub(r"\s+", "", cell).upper().rstrip("*")
                if not cell or compact == sequence or len(cell) > 500 or ID_RE.search(cell):
                    continue
                if cell.isdigit():
                    continue
                annotation.append(cell)
            rows.append({
                "query": identifier, "sheet": title, "category": category,
                "source_row": row_number, "length": len(sequence),
                "annotation": " | ".join(annotation[:8]), "sequence": sequence,
            })
    frame = pd.DataFrame(rows).drop_duplicates(["query", "sequence"]).reset_index(drop=True)
    if len(frame) < 850:
        raise RuntimeError(f"Only {len(frame)} proteins extracted")
    frame.to_csv(OUT / "mitoproteome_inventory.tsv.gz", sep="\t", index=False, compression="gzip")
    write_fasta(OUT / "andalucia_mitoproteome.faa", [(row.query, row.sequence) for row in frame.itertuples()])
    return frame


def classify(scores: pd.DataFrame) -> pd.DataFrame:
    scores = scores.copy()
    scores["N_best_class"] = scores[N_CLASSES].idxmax(axis=1)
    scores["N_best_score"] = scores[N_CLASSES].max(axis=1)
    scores["N_Aas_like_score"] = scores[["AasN", "LpaT"]].max(axis=1)
    scores["N_Aas_specific_margin"] = scores["N_Aas_like_score"] - scores[["PlsC", "PlsB"]].max(axis=1)
    scores["C_best_class"] = scores[C_CLASSES].idxmax(axis=1)
    scores["C_best_score"] = scores[C_CLASSES].max(axis=1)
    scores["C_Aas_like_score"] = scores[["AasC", "split_AasC"]].max(axis=1)
    scores["C_Aas_specific_margin"] = scores["C_Aas_like_score"] - scores[["FadD", "Acs"]].max(axis=1)
    scores["N_provisional"] = scores["N_Aas_like_score"].ge(20) & scores["N_Aas_specific_margin"].ge(8)
    scores["C_provisional"] = scores["C_Aas_like_score"].ge(30) & scores["C_Aas_specific_margin"].ge(8)
    return scores


def stable_seed(text: str) -> int:
    return int(hashlib.sha256((text + "|20260722").encode()).hexdigest()[:16], 16)


def composition_test(row: pd.Series, hmms: dict[str, Path], side: str, replicates: int = 100) -> dict[str, Any]:
    sequence = str(row["sequence"])
    rng = random.Random(stable_seed(str(row["query"]) + side))
    records = [("natural", sequence)]
    for index in range(replicates):
        chars = list(sequence)
        rng.shuffle(chars)
        records.append((f"shuffle_{index:03d}", "".join(chars)))
    fasta = WORK / f"{row['query']}_{side}_nulls.faa"
    write_fasta(fasta, records)
    maps = {label: score_hmm(hmm, fasta, f"{row['query']}_{side}_{label}") for label, hmm in hmms.items()}
    margins = []
    intended_scores = []
    for name, _ in records:
        if side == "N":
            intended = max(maps[label].get(name, -50.0) for label in ["AasN", "LpaT"])
            alternative = max(maps[label].get(name, -50.0) for label in ["PlsC", "PlsB"])
        else:
            intended = max(maps[label].get(name, -50.0) for label in ["AasC", "split_AasC"])
            alternative = max(maps[label].get(name, -50.0) for label in ["FadD", "Acs"])
        intended_scores.append(float(intended))
        margins.append(float(intended - alternative))
    p_margin = (1 + sum(value >= margins[0] for value in margins[1:])) / (replicates + 1)
    p_score = (1 + sum(value >= intended_scores[0] for value in intended_scores[1:])) / (replicates + 1)
    return {
        "query": row["query"], "side": side,
        "natural_intended_score": intended_scores[0],
        "natural_specific_margin": margins[0],
        "null_margin_mean": float(np.mean(margins[1:])),
        "null_margin_q95": float(np.quantile(margins[1:], 0.95)),
        "null_margin_max": float(np.max(margins[1:])),
        "empirical_p_specific_margin": p_margin,
        "empirical_p_intended_score": p_score,
        "replicates": replicates,
        "null_margins_json": json.dumps(margins[1:]),
    }


def main() -> None:
    mitoproteome = extract_mitoproteome()
    panels = class_panels()
    hmms = {label: build_hmm(label, frame) for label, frame in panels.items()}
    fasta = OUT / "andalucia_mitoproteome.faa"
    identifiers = [name for name, _ in parse_fasta(fasta)]
    scores = pd.DataFrame({"query": identifiers})
    for label, hmm in hmms.items():
        score_map = score_hmm(hmm, fasta, f"mitoproteome_{label}")
        scores[label] = scores["query"].map(score_map).fillna(-50.0)
    scores = classify(scores).merge(mitoproteome, on="query", how="left")
    scores.to_csv(OUT / "all_mitoprotein_class_scores.tsv", sep="\t", index=False)

    selected: dict[tuple[str, str], pd.Series] = {}
    for side, margin, flag in [
        ("N", "N_Aas_specific_margin", "N_provisional"),
        ("C", "C_Aas_specific_margin", "C_provisional"),
    ]:
        subset = pd.concat([scores[scores[flag]], scores.nlargest(8, margin)], ignore_index=True).drop_duplicates("query")
        for _, row in subset.iterrows():
            selected[(str(row["query"]), side)] = row
    nulls = pd.DataFrame([composition_test(row, hmms, side) for (_, side), row in selected.items()])
    nulls.to_csv(OUT / "candidate_exact_composition_nulls.tsv", sep="\t", index=False)
    if len(nulls):
        pivot = nulls.pivot(index="query", columns="side", values="empirical_p_specific_margin").rename(columns={"N": "N_null_p", "C": "C_null_p"})
        scores = scores.merge(pivot, on="query", how="left")
    else:
        scores["N_null_p"] = np.nan
        scores["C_null_p"] = np.nan
    scores["N_final_Aas_specific"] = scores["N_provisional"] & scores["N_null_p"].le(0.05)
    scores["C_final_Aas_specific"] = scores["C_provisional"] & scores["C_null_p"].le(0.05)
    scores.to_csv(OUT / "final_mitoprotein_Aas_specificity.tsv", sep="\t", index=False)

    summary = {
        "mitoproteins": int(len(scores)),
        "N_provisional": int(scores["N_provisional"].sum()),
        "C_provisional": int(scores["C_provisional"].sum()),
        "N_final_Aas_specific": int(scores["N_final_Aas_specific"].sum()),
        "C_final_Aas_specific": int(scores["C_final_Aas_specific"].sum()),
        "top_N_candidates": scores.nlargest(12, "N_Aas_specific_margin")[[
            "query", "annotation", "N_Aas_like_score", "N_best_class", "N_best_score", "N_Aas_specific_margin", "N_null_p"
        ]].to_dict("records"),
        "top_C_candidates": scores.nlargest(12, "C_Aas_specific_margin")[[
            "query", "annotation", "C_Aas_like_score", "C_best_class", "C_best_score", "C_Aas_specific_margin", "C_null_p"
        ]].to_dict("records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
