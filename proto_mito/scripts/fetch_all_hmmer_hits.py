from __future__ import annotations

import csv
import gzip
import json
import math
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import requests

OUT = Path("proto_mito_results/hmmer_all_hits")
RAW = OUT / "raw_pages"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)
BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"

PROFILES = [
    {"label":"ATP23_A","KO":"K18156","DB":"TOLDBA","family_id":"OGA0001337","job_id":"62dd04c7-00fd-4af8-885a-c9094e2826be"},
    {"label":"ATP5H_d_A","KO":"K02138","DB":"TOLDBA","family_id":"OGA0004353","job_id":"c6d4d192-1753-49f1-86c8-c6f646e784dc"},
    {"label":"COX17_A","KO":"K02260","DB":"TOLDBA","family_id":"OGA0004921_1","job_id":"db9e2a68-4afc-48c0-ab1a-3ca5d79f1f35"},
    {"label":"COX19_A","KO":"K18183","DB":"TOLDBA","family_id":"OGA0000444_1","job_id":"4c93c85c-3ea9-48c6-947c-35c2f5869403"},
    {"label":"MPC2_A","KO":"K22139","DB":"TOLDBA","family_id":"OGA0000778_1","job_id":"3481d062-8f6d-4489-8cd7-607deb0c4b8b"},
    {"label":"MRPL35_A","KO":"K17439","DB":"TOLDBA","family_id":"OGA0006231","job_id":"d04632ec-3274-483c-84b5-4e7e2372d7e0"},
    {"label":"MRPL35_B","KO":"K17439","DB":"TOLDBB","family_id":"OGB0007321","job_id":"dc013532-f70b-4c4f-a8ad-61b9cb60b97d"},
    {"label":"MRPL35_C","KO":"K17439","DB":"TOLDBC","family_id":"OGC0006364","job_id":"13977123-bec5-4525-9830-ed07ce7684a5"},
    {"label":"MRPL41_A","KO":"K17422","DB":"TOLDBA","family_id":"OGA0003251_1","job_id":"bb781ecf-e0c0-497f-b5ff-f280e17599cc"},
    {"label":"MRPS23_A","KO":"K17402","DB":"TOLDBA","family_id":"OGA0010739","job_id":"711e9624-2d72-4fd2-bf6a-2f36abb309d0"},
    {"label":"NDUFA2_A","KO":"K03946","DB":"TOLDBA","family_id":"OGA0002544","job_id":"7b051f08-818d-4c65-9bde-2f8514840554"},
    {"label":"QCR6_A","KO":"K00416","DB":"TOLDBA","family_id":"OGA0003627","job_id":"3aa1e0c5-50c5-4873-9397-1d1fc9b85f6a"},
    {"label":"QCR7_A","KO":"K00417","DB":"TOLDBA","family_id":"OGA0002602","job_id":"000153b8-a1ec-483d-920a-efc602427ecb"},
    {"label":"QCR7_B","KO":"K00417","DB":"TOLDBB","family_id":"OGB0002633","job_id":"edd3f44d-a210-401b-80da-c32bacbd4c21"},
    {"label":"QCR7_C","KO":"K00417","DB":"TOLDBC","family_id":"OGC0002401","job_id":"af7ec9af-5f33-453d-af0a-856a676307d9"},
    {"label":"QCR9_A","KO":"K00419","DB":"TOLDBA","family_id":"OGA0004703","job_id":"661feada-d130-4e7a-a8ac-c00a779c66bf"},
    {"label":"UQCRQ_A","KO":"K00418","DB":"TOLDBA","family_id":"OGA0012053","job_id":"eba9fb1a-1a27-4695-9d10-973bd57a4d99"},
    {"label":"UQCRQ_C","KO":"K00418","DB":"TOLDBC","family_id":"OGC0009468","job_id":"8295d257-b360-434f-b2ad-ee98ceac7e62"},
]


def session() -> requests.Session:
    value = requests.Session()
    value.headers.update({
        "User-Agent": "proto-mito-full-hmmer-hit-audit/1.0 (public reproducible scientific workflow)",
        "Accept": "application/json",
    })
    return value


def get_json(url: str, params: dict[str, Any] | None = None, attempts: int = 12) -> dict[str, Any]:
    s = session()
    delay = 1.0
    last = ""
    for _ in range(attempts):
        try:
            response = s.get(url, params=params, timeout=300)
            last = f"{response.status_code} {response.text[:500]}"
        except Exception as exc:
            last = repr(exc)
            time.sleep(delay)
            delay = min(delay * 1.7, 45)
            continue
        if response.status_code == 200:
            data = response.json()
            if isinstance(data, dict):
                return data
            last = f"non-object JSON {type(data).__name__}"
        elif response.status_code not in {429, 500, 502, 503, 504}:
            raise RuntimeError(f"GET {response.url} failed: {last}")
        time.sleep(delay)
        delay = min(delay * 1.7, 45)
    raise RuntimeError(f"GET failed after retries: {url}; last={last}")


def fetch_page(profile: dict[str, str], page: int) -> tuple[str, int, dict[str, Any] | None, str | None]:
    try:
        payload = get_json(f"{BASE}/result/{profile['job_id']}", {"page": page})
        return profile["label"], page, payload, None
    except Exception as exc:
        return profile["label"], page, None, repr(exc)


def flatten_hit(profile: dict[str, str], hit: dict[str, Any]) -> dict[str, Any]:
    metadata = hit.get("metadata") or {}
    return {
        **profile,
        "rank": (int(hit.get("index")) + 1) if hit.get("index") is not None else "",
        "evalue": hit.get("evalue", ""),
        "bitscore": hit.get("score", ""),
        "bias": hit.get("bias", ""),
        "is_included": bool(hit.get("is_included")),
        "is_reported": bool(hit.get("is_reported")),
        "nregions": hit.get("nregions", ""),
        "ndom": hit.get("ndom", ""),
        "accession": metadata.get("accession") or metadata.get("uniprot_accession") or "",
        "identifier": metadata.get("identifier") or metadata.get("uniprot_identifier") or "",
        "description": metadata.get("description") or "",
        "kingdom": metadata.get("kingdom") or "",
        "phylum": metadata.get("phylum") or "",
        "species": metadata.get("species") or "",
        "taxonomy_id": metadata.get("taxonomy_id") or "",
        "lineage": json.dumps(metadata.get("lineage") or [], separators=(",", ":")),
        "architecture": metadata.get("architecture") or "",
        "architecture_score": metadata.get("architecture_score") or "",
        "structures": json.dumps(metadata.get("structures") or [], separators=(",", ":")),
        "raw_hit_json": json.dumps(hit, sort_keys=True, separators=(",", ":")),
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# First page establishes immutable result size and number of pages.
profiles_by_label = {profile["label"]: profile for profile in PROFILES}
first_pages: dict[str, dict[str, Any]] = {}
failures: list[dict[str, Any]] = []
for profile in PROFILES:
    try:
        payload = get_json(f"{BASE}/result/{profile['job_id']}", {"page": 1})
        if payload.get("status") != "SUCCESS":
            raise RuntimeError(f"unexpected status {payload.get('status')}")
        first_pages[profile["label"]] = payload
    except Exception as exc:
        failures.append({**profile, "page": 1, "error": repr(exc)})

page_tasks: list[tuple[dict[str, str], int]] = []
for label, payload in first_pages.items():
    profile = profiles_by_label[label]
    page_count = int(payload.get("page_count") or 1)
    for page in range(2, page_count + 1):
        page_tasks.append((profile, page))

all_pages: dict[tuple[str, int], dict[str, Any]] = {
    (label, 1): payload for label, payload in first_pages.items()
}
with ThreadPoolExecutor(max_workers=8) as pool:
    futures = {pool.submit(fetch_page, profile, page): (profile, page) for profile, page in page_tasks}
    for index, future in enumerate(as_completed(futures), 1):
        label, page, payload, error = future.result()
        profile = profiles_by_label[label]
        if payload is not None:
            all_pages[(label, page)] = payload
        else:
            failures.append({**profile, "page": page, "error": error or "unknown"})
        if index % 50 == 0:
            print(f"pages fetched: {index}/{len(futures)}")

all_hits: list[dict[str, Any]] = []
profile_summaries: list[dict[str, Any]] = []
first_non_euk_rows: list[dict[str, Any]] = []
for profile in PROFILES:
    label = profile["label"]
    pages = sorted((page, payload) for (this_label, page), payload in all_pages.items() if this_label == label)
    profile_hits: list[dict[str, Any]] = []
    stats: dict[str, Any] = {}
    expected_page_count = int(first_pages.get(label, {}).get("page_count") or 0)
    profile_dir = RAW / label
    profile_dir.mkdir(parents=True, exist_ok=True)
    for page, payload in pages:
        with gzip.open(profile_dir / f"page_{page:04d}.json.gz", "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
        result = payload.get("result") or {}
        if not stats:
            stats = result.get("stats") or {}
        for hit in result.get("hits") or []:
            profile_hits.append(flatten_hit(profile, hit))
    profile_hits.sort(key=lambda row: int(row["rank"]) if row["rank"] != "" else 10**12)
    all_hits.extend(profile_hits)

    reported = [row for row in profile_hits if row["is_reported"]]
    included = [row for row in profile_hits if row["is_included"]]
    reported_counts = Counter(row["kingdom"] or "Unknown" for row in reported)
    included_counts = Counter(row["kingdom"] or "Unknown" for row in included)
    non_euk = [row for row in reported if row["kingdom"] not in {"Eukaryota", ""}]
    non_euk_included = [row for row in included if row["kingdom"] not in {"Eukaryota", ""}]
    bacterial_included = [row for row in included if row["kingdom"] == "Bacteria"]
    archaeal_included = [row for row in included if row["kingdom"] == "Archaea"]
    viral_included = [row for row in included if row["kingdom"] in {"Viruses", "Virus"}]

    for category, values in [
        ("first_non_euk_reported", non_euk),
        ("first_non_euk_included", non_euk_included),
        ("first_bacterial_included", bacterial_included),
        ("first_archaeal_included", archaeal_included),
    ]:
        if values:
            first_non_euk_rows.append({"category": category, **values[0]})

    profile_summaries.append({
        **profile,
        "expected_pages": expected_page_count,
        "pages_retrieved": len(pages),
        "reported_hits_stats": stats.get("nreported", ""),
        "included_hits_stats": stats.get("nincluded", ""),
        "reported_hits_retrieved": len(reported),
        "included_hits_retrieved": len(included),
        "reported_eukaryota": reported_counts.get("Eukaryota", 0),
        "reported_bacteria": reported_counts.get("Bacteria", 0),
        "reported_archaea": reported_counts.get("Archaea", 0),
        "reported_viruses": reported_counts.get("Viruses", 0) + reported_counts.get("Virus", 0),
        "reported_unknown": reported_counts.get("Unknown", 0),
        "included_eukaryota": included_counts.get("Eukaryota", 0),
        "included_bacteria": included_counts.get("Bacteria", 0),
        "included_archaea": included_counts.get("Archaea", 0),
        "included_viruses": included_counts.get("Viruses", 0) + included_counts.get("Virus", 0),
        "included_unknown": included_counts.get("Unknown", 0),
        "first_non_euk_reported_rank": non_euk[0]["rank"] if non_euk else "",
        "first_non_euk_reported_kingdom": non_euk[0]["kingdom"] if non_euk else "",
        "first_non_euk_reported_accession": non_euk[0]["accession"] if non_euk else "",
        "first_non_euk_reported_description": non_euk[0]["description"] if non_euk else "",
        "first_non_euk_reported_evalue": non_euk[0]["evalue"] if non_euk else "",
        "first_bacterial_included_rank": bacterial_included[0]["rank"] if bacterial_included else "",
        "first_bacterial_included_accession": bacterial_included[0]["accession"] if bacterial_included else "",
        "first_bacterial_included_description": bacterial_included[0]["description"] if bacterial_included else "",
        "first_bacterial_included_evalue": bacterial_included[0]["evalue"] if bacterial_included else "",
        "first_archaeal_included_rank": archaeal_included[0]["rank"] if archaeal_included else "",
        "first_archaeal_included_accession": archaeal_included[0]["accession"] if archaeal_included else "",
        "first_archaeal_included_description": archaeal_included[0]["description"] if archaeal_included else "",
        "first_archaeal_included_evalue": archaeal_included[0]["evalue"] if archaeal_included else "",
        "database_nseqs": stats.get("nseqs", ""),
        "elapsed_seconds": stats.get("elapsed", ""),
    })

hit_fields = [
    "label", "KO", "DB", "family_id", "job_id", "rank", "evalue", "bitscore", "bias",
    "is_included", "is_reported", "nregions", "ndom", "accession", "identifier",
    "description", "kingdom", "phylum", "species", "taxonomy_id", "lineage",
    "architecture", "architecture_score", "structures", "raw_hit_json",
]
summary_fields = list(profile_summaries[0]) if profile_summaries else []
write_tsv(OUT / "all_reported_hits.tsv", all_hits, hit_fields)
write_tsv(OUT / "profile_taxonomy_summary.tsv", profile_summaries, summary_fields)
write_tsv(OUT / "first_non_eukaryotic_hits.tsv", first_non_euk_rows, ["category"] + hit_fields)
write_tsv(OUT / "failures.tsv", failures, list(PROFILES[0]) + ["page", "error"])

# A compact, directly testable list of all significant non-eukaryotic hits.
sig_non_euk = [
    row for row in all_hits
    if row["is_included"] and row["kingdom"] not in {"Eukaryota", ""}
]
write_tsv(OUT / "significant_non_eukaryotic_hits.tsv", sig_non_euk, hit_fields)

run_summary = {
    "profiles": len(PROFILES),
    "pages_expected": sum(int(first_pages.get(p["label"], {}).get("page_count") or 0) for p in PROFILES),
    "pages_retrieved": len(all_pages),
    "reported_hit_rows": len(all_hits),
    "significant_hit_rows": sum(bool(row["is_included"]) for row in all_hits),
    "significant_non_eukaryotic_hits": len(sig_non_euk),
    "profiles_with_significant_bacterial_hits": sum(int(row["included_bacteria"] or 0) > 0 for row in profile_summaries),
    "profiles_with_significant_archaeal_hits": sum(int(row["included_archaea"] or 0) > 0 for row in profile_summaries),
    "failures": len(failures),
}
(OUT / "summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")
print(json.dumps(run_summary, indent=2))
if failures:
    raise RuntimeError(f"Pagination incomplete: {len(failures)} page failures")
