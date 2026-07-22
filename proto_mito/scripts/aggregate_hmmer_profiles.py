from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

ROOT = Path("downloaded_hmmer")
OUT = Path("proto_mito_results/hmmer_aggregate")
OUT.mkdir(parents=True, exist_ok=True)

summaries: list[dict[str, Any]] = []
hits: list[dict[str, Any]] = []
failures: list[dict[str, str]] = []

for summary_path in sorted(ROOT.rglob("summary.json")):
    try:
        summary = json.loads(summary_path.read_text())
        summaries.append(summary)
        label = str(summary.get("label") or summary_path.parent.name)
        hit_path = summary_path.parent / "hits_flat.tsv"
        if hit_path.exists():
            with hit_path.open(encoding="utf-8") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    hits.append(
                        {
                            "label": label,
                            "KO": summary.get("KO", ""),
                            "DB": summary.get("DB", ""),
                            "family_id": summary.get("family_id", ""),
                            **row,
                        }
                    )
    except Exception as exc:
        failures.append({"path": str(summary_path), "error": repr(exc)})

expected_labels = {
    "MRPL35_A", "MRPL35_B", "MRPL35_C", "UQCRQ_A", "UQCRQ_C",
    "QCR7_A", "QCR7_B", "QCR7_C", "QCR6_A", "QCR9_A", "COX17_A",
    "NDUFA2_A", "ATP23_A", "MPC2_A", "ATP5H_d_A", "MRPS23_A",
    "MRPL41_A", "COX19_A",
}
observed_labels = {str(row.get("label")) for row in summaries}
for missing in sorted(expected_labels - observed_labels):
    failures.append({"path": missing, "error": "missing summary artifact"})

summary_fields = sorted({key for row in summaries for key in row})
with (OUT / "profile_summaries.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(summaries)

hit_fields = [
    "label", "KO", "DB", "family_id", "json_path", "accession", "name",
    "description", "evalue", "bitscore", "bias", "taxid", "species",
    "domain_count", "raw_json",
]
with (OUT / "all_hits_flat.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=hit_fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(hits)

with (OUT / "aggregation_failures.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["path", "error"], delimiter="\t")
    writer.writeheader()
    writer.writerows(failures)

run_summary = {
    "expected_profiles": len(expected_labels),
    "completed_profiles": len(observed_labels),
    "flattened_hit_rows": len(hits),
    "aggregation_failures": len(failures),
    "labels": sorted(observed_labels),
}
(OUT / "summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")
print(json.dumps(run_summary, indent=2))
if failures:
    raise RuntimeError(f"HMMER aggregate incomplete: {failures}")
