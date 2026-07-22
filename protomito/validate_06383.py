#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import random
import re
from collections import Counter
from pathlib import Path

import openpyxl

CANDIDATE = (
    "MDPKVLRQVLRVLRSSSISQLNNLSETLKNTVPAFVIPKERILKLLEIGSLEDLISKVQIHGLNAHPQLKILSDSLNPHDKELTIYCRGFLSSDRFQDWIITHDRLVQMHGWSKKAVGWTWPSGRVLPIPPLLPPLPVLYKPAVLTVATALMIGGQIWLQWVLSQQRAVERARDLAEQLTLLRPHFDRIRIVSHSLGCRHVVEACLLLEPNVRPDYLHLCAPAFHDDDMLVRWSEVARRKTIVYYTPKDMLLQTVYRTFNKARVPFGIAPPETFRTDVFQSNRYQSVKVCDHFDLWVHQAYQKRFANFAVNDLGE"
)
MATURE = CANDIDATE[30:]
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO*]{40,}$", re.I)
ID_RE = re.compile(r"(?:ANDGO_)?_?(\d{5})(?:\.mRNA\.1)?", re.I)
PROFILES = {
    "PF01738_DLH": "PF01738",
    "PF00561_Abhydrolase1": "PF00561",
    "PF07859_Abhydrolase3": "PF07859",
    "PF12697_Abhydrolase6": "PF12697",
    "PF01083_Cutinase": "PF01083",
    "PF00756_Esterase": "PF00756",
    "PF03959_FSH1": "PF03959",
    "PF12146_Hydrolase4": "PF12146",
    "PF07224_Chlorophyllase": "PF07224",
}
OUT = Path("results/ANDGO06383_DLH")


def norm(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    seq = re.sub(r"\s+", "", value).upper().rstrip("*")
    if not AA_RE.fullmatch(seq) or set(seq) <= {"A", "C", "G", "T", "N"}:
        return None
    return seq


def fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, seq in rows:
            handle.write(f">{name}\n{seq}\n")


def block_shuffle(seq: str, rng: random.Random, lo: int, hi: int) -> str:
    blocks, pos = [], 0
    while pos < len(seq):
        size = rng.randint(lo, hi)
        blocks.append(seq[pos:pos + size])
        pos += size
    rng.shuffle(blocks)
    return "".join(blocks)


def prepare(workbook: Path, n: int = 5000, seed: int = 6383) -> dict[str, int]:
    OUT.mkdir(parents=True, exist_ok=True)
    book = openpyxl.load_workbook(workbook, data_only=True, read_only=True)
    records: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    counts: Counter[str] = Counter()
    for ws in book.worksheets:
        for row_number, row in enumerate(ws.iter_rows(values_only=True), 1):
            values = list(row)
            seqs = [s for s in (norm(v) for v in values) if s]
            if not seqs:
                continue
            seq = max(seqs, key=len)
            identifier = ""
            for value in values:
                if isinstance(value, str):
                    match = ID_RE.search(value.strip())
                    if match:
                        identifier = f"ANDGO_{match.group(1)}"
                        break
            if not identifier:
                identifier = f"UNRESOLVED_{len(records)+1:05d}"
            key = (identifier, seq)
            if key in seen:
                continue
            seen.add(key)
            counts[identifier] += 1
            panel_id = identifier if counts[identifier] == 1 else f"{identifier}__dup{counts[identifier]}"
            records.append({"panel_id": panel_id, "base_id": identifier, "sheet": ws.title,
                            "row": row_number, "length": len(seq), "sequence": seq})
    matches = [r for r in records if r["base_id"] == "ANDGO_06383"]
    if not matches:
        exact = [r for r in records if r["sequence"] == CANDIDATE]
        if exact:
            exact[0]["base_id"] = "ANDGO_06383"
            exact[0]["panel_id"] = "ANDGO_06383"
            matches = exact
    if not matches:
        raise RuntimeError("ANDGO_06383 absent from workbook parse")
    if str(matches[0]["sequence"]) != CANDIDATE:
        raise RuntimeError("Workbook ANDGO_06383 sequence differs from frozen candidate")
    if len(records) < 700:
        raise RuntimeError(f"Only {len(records)} sequence-containing mitochondrial records parsed")
    fasta(OUT / "mitochondrial_panel.faa", [(str(r["panel_id"]), str(r["sequence"])) for r in records])
    with (OUT / "mitochondrial_panel_inventory.tsv").open("w", newline="") as handle:
        fields = ["panel_id", "base_id", "sheet", "row", "length"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows({k: r[k] for k in fields} for r in records)
    fasta(OUT / "candidate_queries.faa", [
        ("ANDGO_06383_full", CANDIDATE),
        ("ANDGO_06383_mature", MATURE),
        ("ANDGO_06383_structural_core", MATURE[36:283]),
        ("ANDGO_06383_reverse", MATURE[::-1]),
    ])
    rng = random.Random(seed)
    nulls = [("observed", MATURE)]
    for i in range(n):
        chars = list(MATURE); rng.shuffle(chars)
        nulls.append((f"iid_{i:05d}", "".join(chars)))
    for i in range(n):
        nulls.append((f"block3_7_{i:05d}", block_shuffle(MATURE, rng, 3, 7)))
    for i in range(n):
        nulls.append((f"block8_15_{i:05d}", block_shuffle(MATURE, rng, 8, 15)))
    for shift in range(1, len(MATURE)):
        nulls.append((f"circular_{shift:04d}", MATURE[shift:] + MATURE[:shift]))
    fasta(OUT / "candidate_nulls.faa", nulls)
    metadata = {"panel_size": len(records), "candidate_full_length": len(CANDIDATE),
                "candidate_mature_length": len(MATURE), "iid_n": n, "block3_7_n": n,
                "block8_15_n": n, "circular_n": len(MATURE)-1, "seed": seed}
    (OUT / "input_metadata.json").write_text(json.dumps(metadata, indent=2))
    return metadata


def tbl(path: Path) -> dict[str, dict[str, float | str]]:
    result: dict[str, dict[str, float | str]] = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split(maxsplit=18)
        if len(p) < 18:
            continue
        result[p[0]] = {"evalue": float(p[4]), "score": float(p[5]), "bias": float(p[6]),
                        "domain_evalue": float(p[7]), "domain_score": float(p[8]),
                        "description": p[18] if len(p) > 18 else ""}
    return result


def empirical(scores: list[float], observed: float) -> float:
    return (1 + sum(s >= observed for s in scores)) / (1 + len(scores))


def summarize() -> None:
    meta = json.loads((OUT / "input_metadata.json").read_text())
    inventory = list(csv.DictReader((OUT / "mitochondrial_panel_inventory.tsv").open(), delimiter="\t"))
    panel_ids = [r["panel_id"] for r in inventory]
    matched_ids = [r["panel_id"] for r in inventory if 0.70*len(MATURE) <= int(r["length"]) <= 1.30*len(MATURE)]
    schemes = {
        "iid": [f"iid_{i:05d}" for i in range(meta["iid_n"])],
        "block3_7": [f"block3_7_{i:05d}" for i in range(meta["block3_7_n"])],
        "block8_15": [f"block8_15_{i:05d}" for i in range(meta["block8_15_n"])],
        "circular": [f"circular_{i:04d}" for i in range(1, len(MATURE))],
    }
    rows, null_rows = [], []
    for label, accession in PROFILES.items():
        cand = tbl(OUT / "hmmer" / f"{label}.candidate.tbl")
        panel = tbl(OUT / "hmmer" / f"{label}.panel.tbl")
        nulls = tbl(OUT / "hmmer" / f"{label}.nulls.tbl")
        obs = cand.get("ANDGO_06383_mature", {"evalue": 1.0, "score": -math.inf, "domain_score": -math.inf})
        score = float(obs["score"])
        panel_scores = [float(panel.get(x, {"score": -math.inf})["score"]) for x in panel_ids]
        matched_scores = [float(panel.get(x, {"score": -math.inf})["score"]) for x in matched_ids]
        rank = 1 + sum(s > score for s in panel_scores)
        mrank = 1 + sum(s > score for s in matched_scores)
        rows.append({"profile": label, "pfam_accession": accession,
                     "candidate_evalue": float(obs["evalue"]), "candidate_score": score,
                     "candidate_domain_score": float(obs["domain_score"]),
                     "panel_rank": rank, "panel_size": len(panel_ids),
                     "panel_empirical_p": rank/(len(panel_ids)+1),
                     "length_matched_rank": mrank, "length_matched_size": len(matched_ids),
                     "length_matched_empirical_p": mrank/(len(matched_ids)+1)})
        for scheme, ids in schemes.items():
            vals = [float(nulls.get(x, {"score": -math.inf})["score"]) for x in ids]
            null_rows.append({"profile": label, "scheme": scheme, "n": len(vals),
                              "reported_n": sum(x in nulls for x in ids), "observed_score": score,
                              "null_mean_finite": sum(v for v in vals if math.isfinite(v))/max(1,sum(math.isfinite(v) for v in vals)),
                              "null_max": max(vals), "empirical_p": empirical(vals, score)})
    rows.sort(key=lambda r: (r["candidate_evalue"], -r["candidate_score"]))
    for i, row in enumerate(rows, 1): row["candidate_profile_rank"] = i
    for name, data in (("profile_specificity_summary.tsv", rows), ("sequence_order_nulls.tsv", null_rows)):
        with (OUT / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(data)
    dlh = next(r for r in rows if r["profile"] == "PF01738_DLH")
    dlh_nulls = [r for r in null_rows if r["profile"] == "PF01738_DLH"]
    gates = {"dlh_is_best_profile": rows[0]["profile"] == "PF01738_DLH",
             "dlh_panel_p_le_0_01": dlh["panel_empirical_p"] <= 0.01,
             "dlh_length_matched_p_le_0_01": dlh["length_matched_empirical_p"] <= 0.01,
             "all_order_null_p_le_0_01": all(r["empirical_p"] <= 0.01 for r in dlh_nulls)}
    result = {"candidate": "ANDGO_06383", "best_profile": rows[0], "dlh_profile": dlh,
              "dlh_nulls": dlh_nulls, "interpretation_gate": gates,
              "passes_family_specificity_gate": all(gates.values())}
    (OUT / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False))
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["prepare", "summarize"])
    parser.add_argument("--workbook", type=Path)
    parser.add_argument("--nulls", type=int, default=5000)
    args = parser.parse_args()
    if args.mode == "prepare":
        if not args.workbook: raise SystemExit("--workbook is required")
        prepare(args.workbook, args.nulls)
    else:
        summarize()
