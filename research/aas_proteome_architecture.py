#!/usr/bin/env python3
"""Representative-proteome survey of fused and split Aas-like architectures.

This is an orthogonal donor-specificity test.  It searches complete proteins from
one reference proteome per genus in Alphaproteobacteria and all available
Chlamydiota reference proteomes with four independently trained profiles:
fused-Aas N terminus, LpaT, fused-Aas C terminus, and split AasC.
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import aas_module_strengthening as core

OUT = ROOT / "results" / "aas_proteome_architecture"
WORK = OUT / "work"
MODELS = OUT / "models"
for path in (OUT, WORK, MODELS):
    path.mkdir(parents=True, exist_ok=True)

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "proto-mito-Aas-architecture/1.0"})


def proteome_id(value: object) -> str:
    text = str(value or "")
    match = re.search(r"UP\d{9}", text)
    return match.group(0) if match else ""


def select_proteomes(taxon: str, maximum: int) -> pd.DataFrame:
    seed = core.fetch_uniprot(f'taxonomy_name:"{taxon}" AND reviewed:true', 100000)
    seed["proteome_id"] = seed["proteomes"].map(proteome_id)
    seed = seed[seed["proteome_id"].ne("")].copy()
    seed = seed.sort_values(["genus", "proteome_id"]).drop_duplicates("genus")
    if len(seed) > maximum:
        indices = np.linspace(0, len(seed) - 1, maximum, dtype=int)
        seed = seed.iloc[indices]
    return seed[["proteome_id", "genus", "organism", "taxid"]].drop_duplicates("proteome_id").reset_index(drop=True)


def get_bytes(url: str, params: dict[str, str], attempts: int = 7) -> bytes:
    delay = 1.5
    last = ""
    for _ in range(attempts):
        try:
            response = SESSION.get(url, params=params, timeout=600)
            last = f"{response.status_code} {response.text[:200] if response.headers.get('content-type','').startswith('text') else ''}"
            if response.status_code == 200:
                return response.content
            if response.status_code not in {429, 500, 502, 503, 504}:
                break
        except Exception as exc:  # noqa: BLE001
            last = repr(exc)
        import time
        time.sleep(delay)
        delay = min(delay * 1.8, 25)
    raise RuntimeError(f"UniProt proteome download failed: {last}")


def download_one(pid: str) -> tuple[str, str, int]:
    path = WORK / f"{pid}.faa"
    if path.exists() and path.stat().st_size > 1000:
        text = path.read_text(encoding="utf-8")
        return pid, text, text.count("\n>") + int(text.startswith(">"))
    content = get_bytes(
        "https://rest.uniprot.org/uniprotkb/stream",
        {"query": f"proteome:{pid}", "format": "fasta", "compressed": "true"},
    )
    if content[:2] == b"\x1f\x8b":
        text = gzip.decompress(content).decode("utf-8")
    else:
        text = content.decode("utf-8")
    path.write_text(text, encoding="utf-8")
    return pid, text, text.count("\n>") + int(text.startswith(">"))


def combine_proteomes(taxon: str, table: pd.DataFrame) -> tuple[Path, pd.DataFrame]:
    output = OUT / f"{taxon}_representative_proteomes.faa"
    inventory = []
    results: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(download_one, pid): pid for pid in table["proteome_id"]}
        for future in as_completed(futures):
            pid, text, count = future.result()
            results[pid] = text
            row = table[table["proteome_id"].eq(pid)].iloc[0]
            inventory.append({
                "taxon": taxon,
                "proteome_id": pid,
                "genus": row["genus"],
                "organism": row["organism"],
                "taxid": row["taxid"],
                "proteins": count,
            })
    with output.open("w", encoding="utf-8") as dest:
        for pid in sorted(results):
            current_header = None
            chunks: list[str] = []
            for raw in results[pid].splitlines():
                if raw.startswith(">"):
                    if current_header is not None:
                        dest.write(f">{pid}|{current_header}\n{''.join(chunks)}\n")
                    token = raw[1:].split()[0]
                    accession = token.split("|")[1] if token.count("|") >= 2 else token
                    current_header = accession
                    chunks = []
                elif current_header is not None:
                    chunks.append(raw.strip())
            if current_header is not None:
                dest.write(f">{pid}|{current_header}\n{''.join(chunks)}\n")
    frame = pd.DataFrame(inventory).sort_values(["taxon", "genus", "proteome_id"])
    frame.to_csv(OUT / f"{taxon}_proteome_inventory.tsv", sep="\t", index=False)
    return output, frame


def build_profiles() -> dict[str, Path]:
    panels_n, panels_c = core.class_panels()
    profiles: dict[str, Path] = {}
    for label, frame in {
        "AasN": panels_n["AasN"],
        "LpaT": panels_n["LpaT"],
        "AasC": panels_c["AasC"],
        "split_AasC": panels_c["split_AasC"],
    }.items():
        profiles[label] = core.build_hmm(
            [(f"{label}|{row.accession}", str(row.domain_sequence)) for row in frame.itertuples()],
            MODELS / label,
            label,
        )
    return profiles


def search_profile(label: str, hmm: Path, fasta: Path, taxon: str) -> pd.DataFrame:
    domtbl = OUT / f"{taxon}_{label}.domtbl"
    core.run([
        "hmmsearch", "--max", "--noali", "-E", "100", "--domE", "100",
        "--domtblout", str(domtbl), str(hmm), str(fasta),
    ], quiet=True)
    frame = core.parse_domtbl(domtbl)
    if frame.empty:
        return frame
    frame["profile"] = label
    frame["proteome_id"] = frame["target"].astype(str).str.split("|").str[0]
    frame["accession"] = frame["target"].astype(str).str.split("|").str[1]
    frame["side"] = "N" if label in {"AasN", "LpaT"} else "C"
    threshold = 1e-3 if frame["side"].iloc[0] == "N" else 1e-5
    score_threshold = 15 if frame["side"].iloc[0] == "N" else 20
    frame["passes"] = (
        frame["i_evalue"].le(threshold)
        & frame["score"].ge(score_threshold)
        & frame["hmm_coverage"].ge(0.35)
    )
    return frame


def classify_taxon(taxon: str, fasta: Path, inventory: pd.DataFrame, profiles: dict[str, Path]) -> dict[str, object]:
    frames = [search_profile(label, hmm, fasta, taxon) for label, hmm in profiles.items()]
    frames = [frame for frame in frames if not frame.empty]
    all_hits = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    all_hits.to_csv(OUT / f"{taxon}_all_profile_hits.tsv.gz", sep="\t", index=False, compression="gzip")
    passed = all_hits[all_hits["passes"]].copy() if len(all_hits) else pd.DataFrame()
    passed.to_csv(OUT / f"{taxon}_passing_profile_hits.tsv", sep="\t", index=False)
    if passed.empty:
        pd.DataFrame().to_csv(OUT / f"{taxon}_split_architectures.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(OUT / f"{taxon}_fused_architectures.tsv", sep="\t", index=False)
        return {"taxon": taxon, "proteomes": len(inventory), "proteins": int(inventory["proteins"].sum()), "N_hits": 0, "C_hits": 0, "split_proteomes": 0, "fused_proteins": 0}

    # Keep the best profile within each side for every protein.
    best = passed.sort_values(["target", "side", "i_evalue", "score"], ascending=[True, True, True, False]).drop_duplicates(["target", "side"])
    n = best[best["side"].eq("N")].copy()
    c = best[best["side"].eq("C")].copy()

    fused_rows = []
    for target in sorted(set(n["target"]) & set(c["target"])):
        nr = n[n["target"].eq(target)].iloc[0]
        cr = c[c["target"].eq(target)].iloc[0]
        overlap = max(0, min(nr.ali_to, cr.ali_to) - max(nr.ali_from, cr.ali_from) + 1)
        shorter = min(nr.ali_to - nr.ali_from + 1, cr.ali_to - cr.ali_from + 1)
        fused_rows.append({
            "taxon": taxon, "proteome_id": nr.proteome_id, "target": target,
            "N_profile": nr.profile, "C_profile": cr.profile,
            "N_score": nr.score, "C_score": cr.score,
            "N_i_evalue": nr.i_evalue, "C_i_evalue": cr.i_evalue,
            "N_start": nr.ali_from, "N_end": nr.ali_to,
            "C_start": cr.ali_from, "C_end": cr.ali_to,
            "overlap_fraction_shorter": overlap / max(shorter, 1),
            "nonoverlapping": overlap / max(shorter, 1) <= 0.25,
        })
    fused = pd.DataFrame(fused_rows)
    fused.to_csv(OUT / f"{taxon}_fused_architectures.tsv", sep="\t", index=False)

    split_rows = []
    for pid in sorted(set(n["proteome_id"]) & set(c["proteome_id"])):
        ns = n[n["proteome_id"].eq(pid)].sort_values(["i_evalue", "score"], ascending=[True, False])
        cs = c[c["proteome_id"].eq(pid)].sort_values(["i_evalue", "score"], ascending=[True, False])
        for nr in ns.head(5).itertuples():
            for cr in cs.head(5).itertuples():
                if nr.target == cr.target:
                    continue
                split_rows.append({
                    "taxon": taxon, "proteome_id": pid,
                    "N_target": nr.target, "N_accession": nr.accession, "N_profile": nr.profile,
                    "N_score": nr.score, "N_i_evalue": nr.i_evalue,
                    "C_target": cr.target, "C_accession": cr.accession, "C_profile": cr.profile,
                    "C_score": cr.score, "C_i_evalue": cr.i_evalue,
                })
    split = pd.DataFrame(split_rows)
    split.to_csv(OUT / f"{taxon}_split_architectures.tsv", sep="\t", index=False)
    return {
        "taxon": taxon,
        "proteomes": len(inventory),
        "proteins": int(inventory["proteins"].sum()),
        "N_hits": int(n["target"].nunique()),
        "C_hits": int(c["target"].nunique()),
        "split_proteomes": int(split["proteome_id"].nunique()) if len(split) else 0,
        "fused_proteins": int(fused.loc[fused["nonoverlapping"], "target"].nunique()) if len(fused) else 0,
    }


def targeted_chlamydia_controls(profiles: dict[str, Path]) -> None:
    query = (
        'taxonomy_name:"Chlamydiota" AND '
        '(gene:CT775 OR gene:CT776 OR gene:lpaT OR gene:aasC OR '
        'protein_name:"lysophospholipid acyltransferase" OR '
        'protein_name:"acyl-acyl carrier protein synthetase")'
    )
    controls = core.fetch_uniprot(query, 1000)
    controls.to_csv(OUT / "Chlamydiota_targeted_control_sequences.tsv", sep="\t", index=False)
    if controls.empty:
        return
    fasta = OUT / "Chlamydiota_targeted_controls.faa"
    core.write_fasta(fasta, [(f"CONTROL|{row.accession}", row.sequence) for row in controls.itertuples()])
    rows = []
    for label, hmm in profiles.items():
        score_map = core.score_hmm(hmm, fasta, f"control_{label}", OUT)
        for row in controls.itertuples():
            rows.append({
                "accession": row.accession, "protein_name": row.protein_name,
                "gene": row.gene, "organism": row.organism, "profile": label,
                "score": score_map.get(f"CONTROL|{row.accession}", -50.0),
            })
    pd.DataFrame(rows).to_csv(OUT / "Chlamydiota_targeted_control_profile_scores.tsv", sep="\t", index=False)


def main() -> None:
    profiles = build_profiles()
    controls = targeted_chlamydia_controls(profiles)
    summaries = []
    selections = {
        "Chlamydiota": select_proteomes("Chlamydiota", 100),
        "Alphaproteobacteria": select_proteomes("Alphaproteobacteria", 100),
    }
    for taxon, table in selections.items():
        fasta, inventory = combine_proteomes(taxon, table)
        summaries.append(classify_taxon(taxon, fasta, inventory, profiles))
    pd.DataFrame(summaries).to_csv(OUT / "proteome_architecture_summary.tsv", sep="\t", index=False)
    (OUT / "summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
