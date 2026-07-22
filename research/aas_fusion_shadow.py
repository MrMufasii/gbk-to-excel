#!/usr/bin/env python3
"""Test whether the tafazzin+fatty-synthetase Aas signal is a fusion shadow.

A fusion shadow occurs when two proteins from broad superfamilies tile a fused
bacterial protein even though neither is specifically descended from that fused
family. We build independent class profiles for AasN, LpaT, PlsC, PlsB and for
AasC, split AasC, FadD, Acs, search them against held-out full-length Aas
proteins, and quantify the domain coverage obtained by every N x C class pair.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import aas_module_strengthening as core

OUT = ROOT / "results" / "aas_fusion_shadow"
MODELS = OUT / "models"
OUT.mkdir(parents=True, exist_ok=True)
MODELS.mkdir(parents=True, exist_ok=True)

N_CLASSES = ["AasN", "LpaT", "PlsC", "PlsB"]
C_CLASSES = ["AasC", "split_AasC", "FadD", "Acs"]


def search_domains(label: str, hmm: Path, fasta: Path) -> pd.DataFrame:
    domtbl = OUT / f"{label}.domtbl"
    core.run([
        "hmmsearch", "--max", "--noali", "-E", "1000000", "--domE", "1000000",
        "--domtblout", str(domtbl), str(hmm), str(fasta),
    ], quiet=True)
    frame = core.parse_domtbl(domtbl)
    if frame.empty:
        return frame
    frame["profile"] = label
    return frame.sort_values(
        ["target", "i_evalue", "score"], ascending=[True, True, False]
    ).drop_duplicates("target")


def interval_union(a0: int, a1: int, b0: int, b1: int, length: int) -> dict[str, float | bool]:
    a0, a1 = sorted((int(a0), int(a1)))
    b0, b1 = sorted((int(b0), int(b1)))
    la = a1 - a0 + 1
    lb = b1 - b0 + 1
    overlap = max(0, min(a1, b1) - max(a0, b0) + 1)
    union = la + lb - overlap
    return {
        "N_start": a0,
        "N_end": a1,
        "C_start": b0,
        "C_end": b1,
        "N_length": la,
        "C_length": lb,
        "overlap": overlap,
        "overlap_fraction_shorter": overlap / max(1, min(la, lb)),
        "combined_coverage": union / max(1, int(length)),
        "correct_order": a0 < b0,
        "nonoverlapping": overlap / max(1, min(la, lb)) <= 0.25,
    }


def main() -> None:
    panels_n, panels_c = core.class_panels()
    panels = {**panels_n, **panels_c}

    # Alternate Aas accessions between profile training and architecture testing.
    aas = panels_n["AasN"].drop_duplicates("accession").reset_index(drop=True)
    indices = np.arange(len(aas))
    test_mask = indices % 2 == 1
    aas_test = aas.loc[test_mask].copy()
    aas_train_accessions = set(aas.loc[~test_mask, "accession"])

    aas_fasta = OUT / "heldout_full_length_Aas.faa"
    core.write_fasta(aas_fasta, [(f"AAS|{r.accession}", r.sequence) for r in aas_test.itertuples()])
    length_column = "full_length" if "full_length" in aas_test.columns else "length"
    inventory = aas_test[["accession", "organism", "taxid", length_column]].copy()
    inventory = inventory.rename(columns={length_column: "full_length"})
    inventory.to_csv(OUT / "heldout_Aas_inventory.tsv", sep="\t", index=False)

    hmms: dict[str, Path] = {}
    for label in N_CLASSES + C_CLASSES:
        frame = panels[label].copy()
        if label in {"AasN", "AasC"}:
            frame = frame[frame["accession"].isin(aas_train_accessions)]
        hmms[label] = core.build_hmm(
            [(f"{label}|{row.accession}", str(row.domain_sequence)) for row in frame.itertuples()],
            MODELS / label,
            label,
        )

    hits = {label: search_domains(label, hmm, aas_fasta) for label, hmm in hmms.items()}
    for label, frame in hits.items():
        frame.to_csv(OUT / f"{label}_heldout_Aas_domains.tsv", sep="\t", index=False)

    rows = []
    for n_class in N_CLASSES:
        n_hits = hits[n_class]
        for c_class in C_CLASSES:
            c_hits = hits[c_class]
            common = sorted(set(n_hits.get("target", [])) & set(c_hits.get("target", [])))
            for target in common:
                nr = n_hits[n_hits["target"].eq(target)].iloc[0]
                cr = c_hits[c_hits["target"].eq(target)].iloc[0]
                metrics = interval_union(
                    nr.ali_from, nr.ali_to, cr.ali_from, cr.ali_to,
                    max(nr.target_length, cr.target_length),
                )
                rows.append({
                    "N_class": n_class,
                    "C_class": c_class,
                    "target": target,
                    "N_score": float(nr.score),
                    "N_i_evalue": float(nr.i_evalue),
                    "N_hmm_coverage": float(nr.hmm_coverage),
                    "C_score": float(cr.score),
                    "C_i_evalue": float(cr.i_evalue),
                    "C_hmm_coverage": float(cr.hmm_coverage),
                    **metrics,
                })
    pair_hits = pd.DataFrame(rows)
    pair_hits.to_csv(OUT / "class_pair_target_architectures.tsv", sep="\t", index=False)

    summaries = []
    for (n_class, c_class), group in pair_hits.groupby(["N_class", "C_class"]):
        valid = group[
            group["correct_order"]
            & group["nonoverlapping"]
            & group["N_score"].ge(10)
            & group["C_score"].ge(20)
        ]
        summaries.append({
            "N_class": n_class,
            "C_class": c_class,
            "targets_with_both_domains": len(group),
            "valid_order_nonoverlap_targets": len(valid),
            "valid_fraction": len(valid) / max(1, len(group)),
            "median_combined_coverage_all": float(group["combined_coverage"].median()),
            "median_combined_coverage_valid": float(valid["combined_coverage"].median()) if len(valid) else math.nan,
            "q25_combined_coverage_valid": float(valid["combined_coverage"].quantile(0.25)) if len(valid) else math.nan,
            "minimum_combined_coverage_valid": float(valid["combined_coverage"].min()) if len(valid) else math.nan,
            "median_N_score": float(group["N_score"].median()),
            "median_C_score": float(group["C_score"].median()),
        })
    summary = pd.DataFrame(summaries).sort_values(
        ["valid_fraction", "median_combined_coverage_valid", "median_N_score", "median_C_score"],
        ascending=[False, False, False, False],
    )
    summary.to_csv(OUT / "class_pair_architecture_summary.tsv", sep="\t", index=False)

    observed_coverage = 0.901252
    observed: dict[str, object] = {
        "andalucia_median_combined_coverage": observed_coverage,
        "andalucia_N_class_by_profile": "PlsC",
        "andalucia_C_class_by_profile": "FadD",
        "class_pairs_with_median_valid_coverage_at_least_observed": int(
            summary["median_combined_coverage_valid"].ge(observed_coverage).fillna(False).sum()
        ),
        "total_class_pairs": int(len(summary)),
    }
    pls_fad = summary[(summary["N_class"].eq("PlsC")) & (summary["C_class"].eq("FadD"))]
    if len(pls_fad):
        observed["PlsC_FadD_matrix_result"] = pls_fad.iloc[0].to_dict()
    aas_pair = summary[(summary["N_class"].eq("AasN")) & (summary["C_class"].eq("AasC"))]
    if len(aas_pair):
        observed["AasN_AasC_positive_control"] = aas_pair.iloc[0].to_dict()
    (OUT / "fusion_shadow_summary.json").write_text(
        json.dumps(observed, indent=2, default=str), encoding="utf-8"
    )
    print(json.dumps(observed, indent=2, default=str))


if __name__ == "__main__":
    main()
