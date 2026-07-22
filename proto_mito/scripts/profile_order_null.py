from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run(command: list[str], *, stdout: Path | None = None) -> subprocess.CompletedProcess[str]:
    if stdout is not None:
        with stdout.open("w", encoding="utf-8") as handle:
            return subprocess.run(command, check=True, text=True, stdout=handle, stderr=subprocess.PIPE)
    return subprocess.run(command, check=True, text=True, capture_output=True)


def build_and_score(alignment: Path, targets: Path, work: Path, prefix: str) -> dict[str, float]:
    hmm = work / f"{prefix}.hmm"
    tbl = work / f"{prefix}.tbl"
    run(["hmmbuild", "--amino", "--cpu", "1", "-n", prefix, str(hmm), str(alignment)])
    run(
        [
            "hmmsearch", "--max", "--noali", "--cpu", "1",
            "-E", "1000000", "--domE", "1000000", "--incE", "1000000", "--incdomE", "1000000",
            "--tblout", str(tbl), str(hmm), str(targets),
        ]
    )
    scores: dict[str, float] = {}
    for line in tbl.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        scores[parts[0]] = float(parts[5])
    return scores


def clean_identifier(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", value)[:120]


def write_alignment(records: list[SeqRecord], path: Path) -> None:
    AlignIO.write(MultipleSeqAlignment(records), path, "fasta")


def global_column_shuffle(records: list[SeqRecord], rng: np.random.Generator) -> list[SeqRecord]:
    length = len(records[0].seq)
    order = rng.permutation(length)
    return [
        SeqRecord(Seq("".join(str(record.seq)[int(index)] for index in order)), id=record.id, description="")
        for record in records
    ]


def within_row_residue_shuffle(records: list[SeqRecord], rng: np.random.Generator) -> list[SeqRecord]:
    shuffled: list[SeqRecord] = []
    for record in records:
        chars = list(str(record.seq))
        positions = [index for index, residue in enumerate(chars) if residue not in {"-", "."}]
        residues = [chars[index] for index in positions]
        rng.shuffle(residues)
        for index, residue in zip(positions, residues):
            chars[index] = residue
        shuffled.append(SeqRecord(Seq("".join(chars)), id=record.id, description=""))
    return shuffled


def circular_column_shift(records: list[SeqRecord], rng: np.random.Generator) -> list[SeqRecord]:
    length = len(records[0].seq)
    shift = int(rng.integers(1, length))
    return [
        SeqRecord(Seq(str(record.seq)[shift:] + str(record.seq)[:shift]), id=record.id, description="")
        for record in records
    ]


def stats(observed: float, null: list[float]) -> dict[str, Any]:
    values = np.asarray(null, dtype=float)
    finite = values[np.isfinite(values)]
    if not len(finite):
        return {
            "observed": observed,
            "n": len(values),
            "finite_n": 0,
            "empirical_p_ge": 1.0,
        }
    return {
        "observed": observed,
        "n": len(values),
        "finite_n": len(finite),
        "null_mean": float(np.mean(finite)),
        "null_sd": float(np.std(finite, ddof=1)) if len(finite) > 1 else 0.0,
        "null_median": float(np.median(finite)),
        "null_q90": float(np.quantile(finite, 0.90)),
        "null_q95": float(np.quantile(finite, 0.95)),
        "null_q99": float(np.quantile(finite, 0.99)),
        "null_max": float(np.max(finite)),
        "empirical_p_ge": float((1 + np.sum(values >= observed)) / (1 + len(values))),
        "z_score": float((observed - np.mean(finite)) / np.std(finite, ddof=1)) if len(finite) > 1 and np.std(finite, ddof=1) > 0 else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--accessions", required=True, help="Comma-separated target accessions")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260722)
    parser.add_argument("--out-root", type=Path, default=Path("proto_mito_results/profile_nulls"))
    args = parser.parse_args()

    out = args.out_root / args.label
    out.mkdir(parents=True, exist_ok=True)
    target_accessions = [value.strip() for value in args.accessions.split(",") if value.strip()]

    alignment = AlignIO.read(args.alignment, "fasta")
    records = [SeqRecord(Seq(str(record.seq).upper()), id=clean_identifier(record.id), description="") for record in alignment]
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise RuntimeError("Input is not a rectangular alignment")
    if not records:
        raise RuntimeError("Empty alignment")

    all_targets = {record.id: str(record.seq) for record in SeqIO.parse(args.targets, "fasta")}
    missing = [accession for accession in target_accessions if accession not in all_targets]
    if missing:
        raise RuntimeError(f"Target sequences missing: {missing}; available={sorted(all_targets)}")
    selected_target_path = out / "selected_targets.fasta"
    SeqIO.write(
        [SeqRecord(Seq(all_targets[accession]), id=accession, description="") for accession in target_accessions],
        selected_target_path,
        "fasta",
    )
    natural_alignment = out / "natural_alignment.fasta"
    write_alignment(records, natural_alignment)

    rng = np.random.default_rng(args.seed)
    null_scores: dict[str, dict[str, list[float]]] = {
        accession: {"column_permutation": [], "row_composition": [], "circular_shift": []}
        for accession in target_accessions
    }
    detailed_rows: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="profile_null_") as tmp:
        work = Path(tmp)
        natural = build_and_score(natural_alignment, selected_target_path, work, f"{args.label}_natural")
        natural_scores = {accession: natural.get(accession, float("-inf")) for accession in target_accessions}

        generators = {
            "column_permutation": global_column_shuffle,
            "row_composition": within_row_residue_shuffle,
            "circular_shift": circular_column_shift,
        }
        for null_type, generator in generators.items():
            for iteration in range(args.permutations):
                randomized = generator(records, rng)
                alignment_path = work / f"{null_type}_{iteration:05d}.fasta"
                write_alignment(randomized, alignment_path)
                scores = build_and_score(alignment_path, selected_target_path, work, f"{args.label}_{null_type}_{iteration:05d}")
                for accession in target_accessions:
                    score = scores.get(accession, float("-inf"))
                    null_scores[accession][null_type].append(score)
                    detailed_rows.append(
                        {
                            "label": args.label,
                            "target_accession": accession,
                            "null_type": null_type,
                            "iteration": iteration + 1,
                            "bitscore": score,
                        }
                    )
                # HMMER products are no longer needed after parsing.
                for suffix in [".hmm", ".tbl"]:
                    path = work / f"{args.label}_{null_type}_{iteration:05d}{suffix}"
                    path.unlink(missing_ok=True)
                alignment_path.unlink(missing_ok=True)
                if (iteration + 1) % 100 == 0:
                    print(f"{args.label}: {null_type} {iteration+1}/{args.permutations}")

    summary_rows: list[dict[str, Any]] = []
    for accession in target_accessions:
        for null_type, values in null_scores[accession].items():
            summary_rows.append(
                {
                    "label": args.label,
                    "target_accession": accession,
                    "null_type": null_type,
                    **stats(natural_scores[accession], values),
                }
            )

    summary_fields = sorted({key for row in summary_rows for key in row})
    with (out / "profile_null_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary_rows)
    with (out / "profile_null_scores.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["label", "target_accession", "null_type", "iteration", "bitscore"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(detailed_rows)

    metadata = {
        "label": args.label,
        "alignment": str(args.alignment),
        "alignment_sequences": len(records),
        "alignment_columns": len(records[0].seq),
        "target_accessions": target_accessions,
        "permutations_per_null": args.permutations,
        "null_types": list(null_scores[target_accessions[0]]),
        "seed": args.seed,
        "natural_scores": natural_scores,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({"metadata": metadata, "summary": summary_rows}, indent=2))


if __name__ == "__main__":
    main()
