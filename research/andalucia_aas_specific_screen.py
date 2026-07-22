#!/usr/bin/env python3
"""Proteome-wide search for genuinely Aas-specific modules in Andalucia.

The published Table S1 mitoproteome is screened with eight independently
validated bacterial profile HMMs. A candidate must prefer the Aas-specific
profile to all related PlsC/LpaT/PlsB or FadD/Acs profiles and must retain that
margin against exact-composition shuffled versions of the candidate sequence.
"""
from __future__ import annotations

import csv
import hashlib
import json
import random
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import openpyxl
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import aas_module_strengthening as core

OUT = ROOT / "results" / "andalucia_aas_specific_screen"
MODELS = OUT / "models"
WORK = OUT / "work"
for path in (OUT, MODELS, WORK):
    path.mkdir(parents=True, exist_ok=True)

XLSX = ROOT / "inputs" / "Andalucia_TableS1.xlsx"
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO*]{40,}$", re.I)
ID_RE = re.compile(r"(?:ANDGO_)?_?(\d{5})\.mRNA\.1", re.I)
N_CLASSES = ["AasN", "LpaT", "PlsC", "PlsB"]
C_CLASSES = ["AasC", "split_AasC", "FadD", "Acs"]


def stable_seed(text: str, replicate: int = 0) -> int:
    digest = hashlib.sha256(f"{text}|{replicate}|20260722".encode()).hexdigest()
    return int(digest[:16], 16)


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
            sequence_candidates = []
            for cell in cells:
                compact = re.sub(r"\s+", "", cell).upper().rstrip("*")
                if AA_RE.fullmatch(compact) and not set(compact) <= {"A", "C", "G", "T", "N"}:
                    sequence_candidates.append(compact)
            if not sequence_candidates:
                continue
            sequence = max(sequence_candidates, key=len)
            identifier = ""
            for cell in cells:
                match = ID_RE.search(cell)
                if match:
                    identifier = f"ANDGO_{match.group(1)}"
                    break
            if not identifier:
                identifier = f"{category}_ROW_{row_number:04d}"
            key = (identifier, sequence)
            if key in seen:
                continue
            seen.add(key)
            annotation_cells = []
            for cell in cells:
                compact = re.sub(r"\s+", "", cell).upper().rstrip("*")
                if not cell or compact == sequence or len(cell) > 500:
                    continue
                if ID_RE.search(cell):
                    continue
                if cell.isdigit():
                    continue
                annotation_cells.append(cell)
            rows.append({
                "query": identifier,
                "sheet": title,
                "category": category,
                "source_row": row_number,
                "length": len(sequence),
                "annotation": " | ".join(annotation_cells[:8]),
                "sequence": sequence,
            })
    frame = pd.DataFrame(rows).drop_duplicates(["query", "sequence"]).reset_index(drop=True)
    if len(frame) < 850:
        raise RuntimeError(f"Only {len(frame)} mitoproteins were extracted from Table S1")
    frame.to_csv(OUT / "mitoproteome_inventory.tsv.gz", sep="\t", index=False, compression="gzip")
    core.write_fasta(OUT / "andalucia_mitoproteome.faa", [(row.query, row.sequence) for row in frame.itertuples()])
    return frame


def build_profiles() -> dict[str, Path]:
    panels_n, panels_c = core.class_panels()
    panels = {**panels_n, **panels_c}
    profiles: dict[str, Path] = {}
    for label in N_CLASSES + C_CLASSES:
        frame = panels[label]
        profiles[label] = core.build_hmm(
            [(f"{label}|{row.accession}", str(row.domain_sequence)) for row in frame.itertuples()],
            MODELS / label,
            label,
        )
    return profiles


def score_all(profiles: dict[str, Path], fasta: Path) -> pd.DataFrame:
    all_ids = [name for name, _ in core.parse_fasta(fasta)]
    table = pd.DataFrame({"query": all_ids})
    for label, hmm in profiles.items():
        scores = core.score_hmm(hmm, fasta, f"mitoproteome_{label}", WORK)
        table[label] = table["query"].map(scores).fillna(-50.0)
    return table


def classify(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    table["N_best_class"] = table[N_CLASSES].idxmax(axis=1)
    table["N_best_score"] = table[N_CLASSES].max(axis=1)
    table["N_Aas_specific_margin"] = table["AasN"] - table[["LpaT", "PlsC", "PlsB"]].max(axis=1)
    table["C_best_class"] = table[C_CLASSES].idxmax(axis=1)
    table["C_best_score"] = table[C_CLASSES].max(axis=1)
    table["C_Aas_specific_score"] = table[["AasC", "split_AasC"]].max(axis=1)
    table["C_Aas_specific_margin"] = table["C_Aas_specific_score"] - table[["FadD", "Acs"]].max(axis=1)
    table["N_provisional"] = table["AasN"].ge(20) & table["N_Aas_specific_margin"].ge(8)
    table["C_provisional"] = table["C_Aas_specific_score"].ge(30) & table["C_Aas_specific_margin"].ge(8)
    return table


def exact_composition_nulls(row: pd.Series, profiles: dict[str, Path], side: str, replicates: int = 100) -> dict[str, Any]:
    sequence = str(row["sequence"])
    rng = random.Random(stable_seed(str(row["query"])))
    records = [("natural", sequence)]
    for replicate in range(replicates):
        chars = list(sequence)
        rng.shuffle(chars)
        records.append((f"shuffle_{replicate:03d}", "".join(chars)))
    fasta = WORK / f"{row['query']}_{side}_nulls.faa"
    core.write_fasta(fasta, records)
    score_maps = {label: core.score_hmm(hmm, fasta, f"{row['query']}_{side}_{label}", WORK) for label, hmm in profiles.items()}
    margins = []
    intended = []
    for name, _ in records:
        if side == "N":
            target = score_maps["AasN"].get(name, -50.0)
            margin = target - max(score_maps[label].get(name, -50.0) for label in ["LpaT", "PlsC", "PlsB"])
        else:
            target = max(score_maps[label].get(name, -50.0) for label in ["AasC", "split_AasC"])
            margin = target - max(score_maps[label].get(name, -50.0) for label in ["FadD", "Acs"])
        intended.append(float(target))
        margins.append(float(margin))
    natural_margin = margins[0]
    natural_score = intended[0]
    shuffled = margins[1:]
    p_margin = (1 + sum(value >= natural_margin for value in shuffled)) / (1 + len(shuffled))
    p_score = (1 + sum(value >= natural_score for value in intended[1:])) / (1 + len(shuffled))
    return {
        "query": row["query"],
        "side": side,
        "natural_intended_score": natural_score,
        "natural_specific_margin": natural_margin,
        "shuffle_margin_mean": float(np.mean(shuffled)),
        "shuffle_margin_q95": float(np.quantile(shuffled, 0.95)),
        "shuffle_margin_max": float(np.max(shuffled)),
        "empirical_p_specific_margin": p_margin,
        "empirical_p_intended_score": p_score,
        "replicates": replicates,
        "null_margins_json": json.dumps(shuffled),
    }


def main() -> None:
    mitoproteome = extract_mitoproteome()
    profiles = build_profiles()
    scores = classify(score_all(profiles, OUT / "andalucia_mitoproteome.faa"))
    scores = scores.merge(mitoproteome, on="query", how="left")
    scores.to_csv(OUT / "all_mitoprotein_class_scores.tsv", sep="\t", index=False)

    # Test all provisional candidates and the top five margins on each side.
    candidates: dict[tuple[str, str], pd.Series] = {}
    for side, margin_col, flag_col in [
        ("N", "N_Aas_specific_margin", "N_provisional"),
        ("C", "C_Aas_specific_margin", "C_provisional"),
    ]:
        subset = scores[scores[flag_col]].copy()
        subset = pd.concat([subset, scores.nlargest(5, margin_col)], ignore_index=True).drop_duplicates("query")
        for _, row in subset.iterrows():
            candidates[(str(row["query"]), side)] = row
    null_rows = [exact_composition_nulls(row, profiles, side) for (query, side), row in candidates.items()]
    nulls = pd.DataFrame(null_rows)
    nulls.to_csv(OUT / "candidate_exact_composition_nulls.tsv", sep="\t", index=False)

    if len(nulls):
        scores = scores.merge(
            nulls[["query", "side", "empirical_p_specific_margin"]].pivot(index="query", columns="side", values="empirical_p_specific_margin").rename(columns={"N": "N_null_p", "C": "C_null_p"}),
            on="query", how="left",
        )
    scores["N_final_Aas_specific"] = scores["N_provisional"] & scores.get("N_null_p", pd.Series(index=scores.index, dtype=float)).le(0.05)
    scores["C_final_Aas_specific"] = scores["C_provisional"] & scores.get("C_null_p", pd.Series(index=scores.index, dtype=float)).le(0.05)
    scores.to_csv(OUT / "final_mitoprotein_Aas_specificity.tsv", sep="\t", index=False)

    summary = {
        "mitoproteins": int(len(scores)),
        "N_provisional": int(scores["N_provisional"].sum()),
        "C_provisional": int(scores["C_provisional"].sum()),
        "N_final_Aas_specific": int(scores["N_final_Aas_specific"].sum()),
        "C_final_Aas_specific": int(scores["C_final_Aas_specific"].sum()),
        "top_N_candidates": scores.nlargest(10, "N_Aas_specific_margin")[["query", "annotation", "AasN", "N_best_class", "N_best_score", "N_Aas_specific_margin"]].to_dict("records"),
        "top_C_candidates": scores.nlargest(10, "C_Aas_specific_margin")[["query", "annotation", "C_Aas_specific_score", "C_best_class", "C_best_score", "C_Aas_specific_margin"]].to_dict("records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
