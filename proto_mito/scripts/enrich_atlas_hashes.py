from __future__ import annotations

import csv
import gzip
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

import requests

ROOT = Path("proto_mito")
OUT = Path("proto_mito_results/atlas_context")
RAW = OUT / "raw"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

HASHES = [
    line.strip()
    for line in (ROOT / "input/atlas_hit_hashes.txt").read_text().splitlines()
    if line.strip()
]
SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "proto-mito-atlas-enrichment/1.0 (public reproducible scientific workflow)",
        "Accept": "application/json",
    }
)


def get_json(url: str, params: dict[str, Any] | None = None, attempts: int = 12) -> dict[str, Any]:
    delay = 1.0
    last = ""
    for _ in range(attempts):
        try:
            response = SESSION.get(url, params=params, timeout=240)
            last = f"{response.status_code} {response.text[:500]}"
        except Exception as exc:
            last = repr(exc)
            time.sleep(delay)
            delay = min(delay * 1.8, 30)
            continue
        if response.status_code == 200:
            value = response.json()
            if isinstance(value, dict):
                return value
            last = f"non-object JSON: {type(value).__name__}"
        elif response.status_code not in {429, 500, 502, 503, 504}:
            raise RuntimeError(f"GET {response.url} failed: {last}")
        time.sleep(delay)
        delay = min(delay * 1.8, 30)
    raise RuntimeError(f"GET failed after retries: {url}; last={last}")


def compact(value: Any) -> str:
    return json.dumps(value if value is not None else {}, sort_keys=True, separators=(",", ":"))


def tax_class(name: str, rank: str, top_phyla: dict[str, Any]) -> str:
    n = (name or "").strip().lower()
    r = (rank or "").strip().lower()
    if n == "bacteria":
        return "Bacteria"
    if n == "archaea":
        return "Archaea"
    if n == "eukaryota":
        return "Eukaryota"
    euk_names = {
        "opisthokonta", "fungi", "metazoa", "holozoa", "amoebozoa",
        "archaeplastida", "viridiplantae", "sar", "alveolata", "stramenopiles",
        "rhizaria", "excavata", "discoba", "euglenozoa", "apusozoa",
    }
    if n in euk_names:
        return "Eukaryota"
    # Phylum summaries are useful context but not a valid replacement for the
    # cluster LCA; keep root/mixed clusters unresolved.
    return "Mixed_or_unresolved"


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


protein_rows: list[dict[str, Any]] = []
cluster_rows: list[dict[str, Any]] = []
failures: list[dict[str, str]] = []
cluster_cache: dict[str, dict[str, Any]] = {}

for index, hash_ in enumerate(HASHES, 1):
    try:
        protein = get_json(
            f"https://biohub.ai/esm/protein/api/v1alpha1/proteins/{hash_}",
            {"topk_features": 20, "fold_on_miss": "false"},
        )
        cluster_hash = str(protein.get("cluster_rep_protein_hash") or hash_)
        protein_rows.append(
            {
                "requested_hit_hash": hash_,
                "protein_hash": protein.get("protein_hash", hash_),
                "accession": protein.get("accession", ""),
                "source": protein.get("source", ""),
                "header": protein.get("header", ""),
                "sequence_length": protein.get("sequence_length", ""),
                "ptm": protein.get("ptm", ""),
                "mean_plddt": protein.get("mean_plddt", ""),
                "cluster_rep_hash": cluster_hash,
            }
        )
        sanitized = dict(protein)
        for key in ["pdb", "residues_plddt", "per_residue_activations", "protein_activations"]:
            sanitized.pop(key, None)
        with gzip.open(RAW / f"protein_{hash_}.json.gz", "wt", encoding="utf-8") as handle:
            json.dump(sanitized, handle, indent=2)

        cluster = cluster_cache.get(cluster_hash)
        if cluster is None:
            cluster = get_json(
                f"https://biohub.ai/esm/protein/api/v1alpha1/clusters/{cluster_hash}",
                {"topk_features": 20},
            )
            cluster_cache[cluster_hash] = cluster
            with gzip.open(RAW / f"cluster_{cluster_hash}.json.gz", "wt", encoding="utf-8") as handle:
                json.dump(cluster, handle, indent=2)

        tax = cluster.get("cluster_taxonomy_info") or {}
        tax_name = str(tax.get("name") or "") if isinstance(tax, dict) else ""
        tax_rank = str(tax.get("rank") or "") if isinstance(tax, dict) else ""
        top_phyla = cluster.get("top_phyla") or {}
        top_pfams = cluster.get("cluster_top_pfam_domains") or {}
        cluster_rows.append(
            {
                "requested_hit_hash": hash_,
                "cluster_rep_hash": cluster_hash,
                "cluster_rep_name": cluster.get("protein_name", ""),
                "cluster_rep_accession": cluster.get("accession", ""),
                "cluster_rep_source": cluster.get("source", ""),
                "cluster_size": cluster.get("cluster_size", ""),
                "cluster_pct_characterized": cluster.get("cluster_pct_characterized", ""),
                "cluster_mean_domain_coverage": cluster.get("cluster_mean_domain_coverage", ""),
                "taxonomy_rank": tax_rank,
                "taxonomy_name": tax_name,
                "taxonomy_class": tax_class(tax_name, tax_rank, top_phyla),
                "top_phyla": compact(top_phyla),
                "top_pfam_domains": compact(top_pfams),
                "cluster_member_count_returned": len(cluster.get("member_protein_hashes") or []),
            }
        )
    except Exception as exc:
        failures.append({"requested_hit_hash": hash_, "error": repr(exc)})
    if index % 10 == 0:
        print(f"resolved {index}/{len(HASHES)}")
    time.sleep(0.2)

protein_fields = [
    "requested_hit_hash", "protein_hash", "accession", "source", "header",
    "sequence_length", "ptm", "mean_plddt", "cluster_rep_hash",
]
cluster_fields = [
    "requested_hit_hash", "cluster_rep_hash", "cluster_rep_name",
    "cluster_rep_accession", "cluster_rep_source", "cluster_size",
    "cluster_pct_characterized", "cluster_mean_domain_coverage",
    "taxonomy_rank", "taxonomy_name", "taxonomy_class", "top_phyla",
    "top_pfam_domains", "cluster_member_count_returned",
]
write_tsv(OUT / "protein_context.tsv", protein_rows, protein_fields)
write_tsv(OUT / "cluster_context.tsv", cluster_rows, cluster_fields)
write_tsv(OUT / "failures.tsv", failures, ["requested_hit_hash", "error"])

summary = {
    "requested_hashes": len(HASHES),
    "proteins_resolved": len(protein_rows),
    "clusters_resolved": len(cluster_rows),
    "unique_cluster_representatives": len(cluster_cache),
    "failures": len(failures),
    "taxonomy_classes": dict(Counter(row["taxonomy_class"] for row in cluster_rows)),
    "bacterial_hit_hashes": [row["requested_hit_hash"] for row in cluster_rows if row["taxonomy_class"] == "Bacteria"],
    "archaeal_hit_hashes": [row["requested_hit_hash"] for row in cluster_rows if row["taxonomy_class"] == "Archaea"],
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
if failures:
    raise RuntimeError(f"{len(failures)} of {len(HASHES)} Atlas lookups failed")
