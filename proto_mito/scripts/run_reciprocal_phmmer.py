from __future__ import annotations

import csv
import gzip
import json
import re
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable

import requests

OUT = Path("proto_mito_results/reciprocal_phmmer")
RAW = OUT / "raw"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)
HMMER = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
UNIPROT = "https://rest.uniprot.org"

QUERIES = [
    {
        "label": "NDUFA2_bacterial_record",
        "source": "uniprot",
        "identifier": "A0A5B2VBN4",
        "expected_pattern": r"NDUFA2|NADH dehydrogenase.*B8|complex I.*B8",
    },
    {
        "label": "ATP23_candidate_1",
        "source": "uniparc",
        "identifier": "UPI000B6B2022",
        "historical_accession": "A0A202DQE6",
        "expected_pattern": r"ATP23|ATP synthase subunit 6 assembly|peptidase M76|metalloprotease",
    },
    {
        "label": "ATP23_candidate_2",
        "source": "uniparc",
        "identifier": "UPI000B6C1C34",
        "historical_accession": "A0A202DNV3",
        "expected_pattern": r"ATP23|ATP synthase subunit 6 assembly|peptidase M76|metalloprotease",
    },
    {
        "label": "ATP23_candidate_3",
        "source": "uniparc",
        "identifier": "UPI000B75C84A",
        "historical_accession": "A0A202DMW5",
        "expected_pattern": r"ATP23|ATP synthase subunit 6 assembly|peptidase M76|metalloprotease",
    },
    {
        "label": "ATP23_candidate_4",
        "source": "uniparc",
        "identifier": "UPI000B6C0926",
        "historical_accession": "A0A202DQA8",
        "expected_pattern": r"ATP23|ATP synthase subunit 6 assembly|peptidase M76|metalloprotease",
    },
    {
        "label": "ATP5H_repeat_control_1",
        "source": "uniprot",
        "identifier": "A0A656D6G0",
        "expected_pattern": r"ATP synthase.*subunit d|ATP5H",
    },
    {
        "label": "ATP5H_repeat_control_2",
        "source": "uniprot",
        "identifier": "A0A2Z5IPW5",
        "expected_pattern": r"ATP synthase.*subunit d|ATP5H",
    },
    {
        "label": "UQCRQ_Acinetobacter_record",
        "source": "literal",
        "identifier": "MFX6572341",
        "sequence": "MGKQPVRMKAVVYALSPFQQKVMPGLWKDITSKVSHKITDNWISTTLLLAPLVGTYSYVQWYLEKEKMEHRY",
        "expected_pattern": r"UQCRQ|cytochrome b-c1.*subunit 8|Qcr8|UcrQ|ubiquinone-binding",
    },
    {
        "label": "UQCRQ_Ecoli_record",
        "source": "literal",
        "identifier": "MQL41719",
        "sequence": "EMAKATVPVKSVIYALSPFQQKIMSGLWKDLPSKLHHKVSENWISATLLLAPLVGTYAYVQNYLEKEKLAHRY",
        "expected_pattern": r"UQCRQ|cytochrome b-c1.*subunit 8|Qcr8|UcrQ|ubiquinone-binding",
    },
    {
        "label": "UQCRQ_archaeal_record",
        "source": "literal",
        "identifier": "RYH14666",
        "sequence": "MDGQVSQHLSPFEQKIVGPLFKDVPLKVFKRVKDFVSEAGLGLGLGIIVFYWGDAKHKELAFHHRA",
        "expected_pattern": r"UQCRQ|cytochrome b-c1.*subunit 8|Qcr8|UcrQ|ubiquinone-binding",
    },
]


def new_session() -> requests.Session:
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "proto-mito-reciprocal-phmmer/1.0 (public reproducible scientific workflow)",
            "Accept": "application/json",
        }
    )
    return session


def request(
    method: str,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    payload: dict[str, Any] | None = None,
    attempts: int = 12,
    timeout: int = 600,
    allow_pending: bool = False,
) -> requests.Response:
    session = new_session()
    delay = 1.0
    last = ""
    for _ in range(attempts):
        try:
            response = session.request(method, url, params=params, json=payload, timeout=timeout)
            last = f"{response.status_code} {response.text[:500]}"
        except Exception as exc:
            last = repr(exc)
            time.sleep(delay)
            delay = min(delay * 1.8, 45)
            continue
        if response.status_code in {200, 201, 202}:
            return response
        if allow_pending and response.status_code in {404, 409, 425}:
            return response
        if response.status_code in {429, 500, 502, 503, 504}:
            time.sleep(delay)
            delay = min(delay * 1.8, 45)
            continue
        raise RuntimeError(f"{method} {response.url} failed: {last}")
    raise RuntimeError(f"{method} {url} failed after retries: {last}")


def json_object(response: requests.Response) -> dict[str, Any]:
    value = response.json()
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected object, received {type(value).__name__}")
    return value


def recursive_dicts(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from recursive_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from recursive_dicts(child)


def sequence_from_payload(payload: dict[str, Any]) -> str:
    sequence = payload.get("sequence")
    if isinstance(sequence, dict):
        value = sequence.get("value") or sequence.get("sequence")
        if isinstance(value, str):
            cleaned = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", value.upper())
            if cleaned:
                return cleaned
    if isinstance(sequence, str):
        cleaned = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", sequence.upper())
        if cleaned:
            return cleaned
    for item in recursive_dicts(payload):
        for key in ["sequence", "value"]:
            value = item.get(key)
            if isinstance(value, str) and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{20,}", value.upper()):
                return value.upper()
    raise RuntimeError("No amino-acid sequence found")


def query_sequence(spec: dict[str, Any]) -> str:
    if spec["source"] == "literal":
        return str(spec["sequence"])
    if spec["source"] == "uniprot":
        payload = json_object(request("GET", f"{UNIPROT}/uniprotkb/{spec['identifier']}.json"))
    elif spec["source"] == "uniparc":
        payload = json_object(request("GET", f"{UNIPROT}/uniparc/{spec['identifier']}"))
    else:
        raise RuntimeError(f"Unsupported source {spec['source']}")
    (RAW / f"source_{spec['label']}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return sequence_from_payload(payload)


def find_job_id(payload: Any) -> str:
    for item in recursive_dicts(payload):
        for key in ["id", "job_id", "jobid", "uuid", "search_id", "searchid"]:
            value = item.get(key)
            if isinstance(value, str) and re.fullmatch(r"[0-9a-fA-F-]{20,64}", value):
                return value
    raise RuntimeError(f"No job ID found: {payload}")


def submit(spec: dict[str, Any], sequence: str) -> str:
    response = request(
        "POST",
        f"{HMMER}/search/phmmer",
        payload={"database": "refprot", "input": sequence},
    )
    payload = json_object(response)
    (RAW / f"submission_{spec['label']}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return find_job_id(payload)


def poll(label: str, job_id: str) -> dict[str, Any]:
    history: list[dict[str, Any]] = []
    for attempt in range(240):
        response = request("GET", f"{HMMER}/result/{job_id}", params={"page": 1}, attempts=4, allow_pending=True)
        if response.status_code in {200, 201, 202}:
            payload = json_object(response)
            history.append({"attempt": attempt + 1, "status": payload.get("status"), "bytes": len(response.content)})
            if payload.get("status") == "SUCCESS" and isinstance(payload.get("result"), dict):
                (RAW / f"poll_{label}.json").write_text(json.dumps(history, indent=2), encoding="utf-8")
                return payload
            if payload.get("status") in {"FAILURE", "ERROR", "CANCELLED"}:
                raise RuntimeError(f"phmmer job {job_id} failed: {payload}")
        else:
            history.append({"attempt": attempt + 1, "http_status": response.status_code})
        time.sleep(8)
    raise TimeoutError(f"phmmer job {job_id} did not finish")


def fetch_page(label: str, job_id: str, page: int) -> tuple[str, int, dict[str, Any] | None, str | None]:
    try:
        payload = json_object(request("GET", f"{HMMER}/result/{job_id}", params={"page": page}))
        return label, page, payload, None
    except Exception as exc:
        return label, page, None, repr(exc)


def flatten(spec: dict[str, Any], job_id: str, hit: dict[str, Any]) -> dict[str, Any]:
    metadata = hit.get("metadata") or {}
    description = str(metadata.get("description") or "")
    expected = bool(re.search(str(spec["expected_pattern"]), description, re.I))
    return {
        "query_label": spec["label"],
        "query_source": spec["source"],
        "query_identifier": spec["identifier"],
        "historical_accession": spec.get("historical_accession", ""),
        "job_id": job_id,
        "rank": int(hit.get("index")) + 1 if hit.get("index") is not None else "",
        "evalue": hit.get("evalue", ""),
        "bitscore": hit.get("score", ""),
        "bias": hit.get("bias", ""),
        "is_reported": bool(hit.get("is_reported")),
        "is_included": bool(hit.get("is_included")),
        "accession": metadata.get("accession") or "",
        "identifier": metadata.get("identifier") or "",
        "description": description,
        "kingdom": metadata.get("kingdom") or "",
        "phylum": metadata.get("phylum") or "",
        "species": metadata.get("species") or "",
        "taxonomy_id": metadata.get("taxonomy_id") or "",
        "lineage": json.dumps(metadata.get("lineage") or [], separators=(",", ":")),
        "expected_family_name": expected,
        "raw_hit_json": json.dumps(hit, sort_keys=True, separators=(",", ":")),
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


sequences: dict[str, str] = {}
jobs: dict[str, str] = {}
first_pages: dict[str, dict[str, Any]] = {}
failures: list[dict[str, str]] = []
for spec in QUERIES:
    try:
        sequence = query_sequence(spec)
        sequences[spec["label"]] = sequence
        job_id = submit(spec, sequence)
        jobs[spec["label"]] = job_id
    except Exception as exc:
        failures.append({"query_label": spec["label"], "stage": "prepare_or_submit", "error": repr(exc)})

# Poll simultaneously to avoid serial server latency.
with ThreadPoolExecutor(max_workers=5) as pool:
    futures = {pool.submit(poll, label, job_id): label for label, job_id in jobs.items()}
    for future in as_completed(futures):
        label = futures[future]
        try:
            first_pages[label] = future.result()
        except Exception as exc:
            failures.append({"query_label": label, "stage": "poll", "error": repr(exc)})

page_tasks: list[tuple[str, str, int]] = []
for label, first_page in first_pages.items():
    for page in range(2, int(first_page.get("page_count") or 1) + 1):
        page_tasks.append((label, jobs[label], page))

pages: dict[tuple[str, int], dict[str, Any]] = {(label, 1): payload for label, payload in first_pages.items()}
with ThreadPoolExecutor(max_workers=8) as pool:
    futures = {pool.submit(fetch_page, *task): task for task in page_tasks}
    for index, future in enumerate(as_completed(futures), 1):
        label, page, payload, error = future.result()
        if payload is not None:
            pages[(label, page)] = payload
        else:
            failures.append({"query_label": label, "stage": f"page_{page}", "error": error or "unknown"})
        if index % 50 == 0:
            print(f"pages fetched {index}/{len(futures)}")

spec_by_label = {spec["label"]: spec for spec in QUERIES}
hit_rows: list[dict[str, Any]] = []
summary_rows: list[dict[str, Any]] = []
for label, spec in spec_by_label.items():
    label_pages = sorted((page, payload) for (this_label, page), payload in pages.items() if this_label == label)
    local: list[dict[str, Any]] = []
    stats: dict[str, Any] = {}
    label_raw = RAW / label
    label_raw.mkdir(parents=True, exist_ok=True)
    for page, payload in label_pages:
        with gzip.open(label_raw / f"page_{page:04d}.json.gz", "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
        result = payload.get("result") or {}
        if not stats:
            stats = result.get("stats") or {}
        for hit in result.get("hits") or []:
            local.append(flatten(spec, jobs.get(label, ""), hit))
    local.sort(key=lambda row: int(row["rank"]) if row["rank"] != "" else 10**12)
    hit_rows.extend(local)
    included = [row for row in local if row["is_included"]]
    expected = [row for row in included if row["expected_family_name"]]
    by_kingdom = Counter(row["kingdom"] or "Unknown" for row in included)
    top = included[0] if included else {}
    top_expected = expected[0] if expected else {}
    summary_rows.append(
        {
            "query_label": label,
            "query_identifier": spec["identifier"],
            "historical_accession": spec.get("historical_accession", ""),
            "sequence_length": len(sequences.get(label, "")),
            "job_id": jobs.get(label, ""),
            "pages_expected": int(first_pages.get(label, {}).get("page_count") or 0),
            "pages_retrieved": len(label_pages),
            "reported_hits": stats.get("nreported", len(local)),
            "included_hits": stats.get("nincluded", len(included)),
            "included_eukaryota": by_kingdom.get("Eukaryota", 0),
            "included_bacteria": by_kingdom.get("Bacteria", 0),
            "included_archaea": by_kingdom.get("Archaea", 0),
            "included_viruses": by_kingdom.get("Viruses", 0) + by_kingdom.get("Virus", 0),
            "included_unknown": by_kingdom.get("Unknown", 0),
            "top_hit_rank": top.get("rank", ""),
            "top_hit_accession": top.get("accession", ""),
            "top_hit_description": top.get("description", ""),
            "top_hit_kingdom": top.get("kingdom", ""),
            "top_hit_evalue": top.get("evalue", ""),
            "top_expected_rank": top_expected.get("rank", ""),
            "top_expected_accession": top_expected.get("accession", ""),
            "top_expected_description": top_expected.get("description", ""),
            "top_expected_kingdom": top_expected.get("kingdom", ""),
            "top_expected_evalue": top_expected.get("evalue", ""),
        }
    )

hit_fields = [
    "query_label", "query_source", "query_identifier", "historical_accession", "job_id",
    "rank", "evalue", "bitscore", "bias", "is_reported", "is_included", "accession",
    "identifier", "description", "kingdom", "phylum", "species", "taxonomy_id",
    "lineage", "expected_family_name", "raw_hit_json",
]
summary_fields = list(summary_rows[0]) if summary_rows else []
write_tsv(OUT / "all_hits.tsv", hit_rows, hit_fields)
write_tsv(OUT / "query_summary.tsv", summary_rows, summary_fields)
write_tsv(OUT / "failures.tsv", failures, ["query_label", "stage", "error"])
with (OUT / "query_sequences.fasta").open("w", encoding="utf-8") as handle:
    for label, sequence in sequences.items():
        handle.write(f">{label}|{spec_by_label[label]['identifier']}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start+80] + "\n")

run_summary = {
    "queries_requested": len(QUERIES),
    "queries_submitted": len(jobs),
    "queries_completed": len(first_pages),
    "pages_retrieved": len(pages),
    "hit_rows": len(hit_rows),
    "failures": len(failures),
}
(OUT / "summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")
print(json.dumps(run_summary, indent=2))
print(json.dumps(summary_rows, indent=2))
if failures:
    raise RuntimeError(f"Reciprocal phmmer incomplete: {failures}")
