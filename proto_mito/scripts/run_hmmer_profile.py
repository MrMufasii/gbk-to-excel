from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import Any, Iterable

import requests

BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "proto-mito-hmmer-validation/1.0 (public reproducible scientific workflow)",
        "Accept": "application/json",
    }
)


def request(
    method: str,
    url: str,
    *,
    json_payload: dict[str, Any] | None = None,
    timeout: int = 300,
    attempts: int = 10,
    allow_pending: bool = False,
) -> requests.Response:
    delay = 2.0
    last: Exception | None = None
    for _ in range(attempts):
        try:
            response = SESSION.request(method, url, json=json_payload, timeout=timeout)
            if response.status_code in {200, 201, 202}:
                return response
            if allow_pending and response.status_code in {404, 409, 425}:
                return response
            if response.status_code in {429, 500, 502, 503, 504}:
                time.sleep(delay)
                delay = min(delay * 1.7, 45)
                continue
            response.raise_for_status()
        except Exception as exc:
            last = exc
            time.sleep(delay)
            delay = min(delay * 1.7, 45)
    raise RuntimeError(f"{method} {url} failed after retries: {last}")


def safe_json(response: requests.Response) -> Any:
    try:
        return response.json()
    except Exception:
        return {
            "_http_status": response.status_code,
            "_content_type": response.headers.get("content-type", ""),
            "_text": response.text[:500000],
        }


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def recursive_values(value: Any, wanted: set[str]) -> Iterable[Any]:
    if isinstance(value, dict):
        for key, child in value.items():
            if key.lower() in wanted:
                yield child
            yield from recursive_values(child, wanted)
    elif isinstance(value, list):
        for child in value:
            yield from recursive_values(child, wanted)


def find_job_id(payload: Any) -> str:
    patterns = {"id", "job_id", "jobid", "uuid", "search_id", "searchid"}
    for value in recursive_values(payload, patterns):
        if isinstance(value, str) and re.fullmatch(r"[0-9a-fA-F-]{20,64}", value):
            return value
    raise RuntimeError(f"No HMMER job identifier found in submission response: {payload}")


def all_statuses(value: Any) -> list[str]:
    values = []
    for item in recursive_values(value, {"status", "state", "job_status"}):
        if isinstance(item, str):
            values.append(item.lower().strip())
    return values


def looks_complete(value: Any) -> bool:
    statuses = all_statuses(value)
    pending = {"pending", "queued", "running", "started", "submitted", "processing", "waiting"}
    failed = {"failed", "error", "cancelled", "canceled", "expired"}
    if any(status in failed for status in statuses):
        raise RuntimeError(f"HMMER job reported failure: {statuses}")
    if any(status in pending for status in statuses):
        return False
    if any(status in {"done", "finished", "complete", "completed", "success", "successful"} for status in statuses):
        return True
    if isinstance(value, dict):
        signal_keys = {"hits", "results", "domains", "query", "stats", "search"}
        if signal_keys.intersection({str(key).lower() for key in value}):
            return True
    return False


def fetch_optional(url: str) -> Any:
    try:
        response = request("GET", url, attempts=5, allow_pending=True)
        if response.status_code not in {200, 201, 202}:
            return {"_http_status": response.status_code, "_text": response.text[:10000]}
        return safe_json(response)
    except Exception as exc:
        return {"_error": repr(exc)}


def scalar(value: Any, names: list[str]) -> Any:
    lowered = {name.lower() for name in names}
    if not isinstance(value, dict):
        return ""
    for key, child in value.items():
        if str(key).lower() in lowered and isinstance(child, (str, int, float, bool)):
            return child
    return ""


def candidate_hit_dicts(value: Any, path: str = "root") -> Iterable[tuple[str, dict[str, Any]]]:
    if isinstance(value, dict):
        keys = {str(key).lower() for key in value}
        identity = {
            "acc", "accession", "target", "target_acc", "name", "description",
            "hit", "hit_acc", "sequence_accession", "taxid", "species",
        }
        score = {"evalue", "e_value", "score", "bitscore", "bit_score", "bias"}
        if keys.intersection(identity) and keys.intersection(score):
            yield path, value
        for key, child in value.items():
            yield from candidate_hit_dicts(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from candidate_hit_dicts(child, f"{path}[{index}]")


def flatten_hits(result: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for path, hit in candidate_hit_dicts(result):
        row = {
            "json_path": path,
            "accession": scalar(hit, ["acc", "accession", "target_acc", "hit_acc", "sequence_accession"]),
            "name": scalar(hit, ["name", "target", "hit", "target_name"]),
            "description": scalar(hit, ["description", "desc", "protein_name"]),
            "evalue": scalar(hit, ["evalue", "e_value", "eval", "full_evalue"]),
            "bitscore": scalar(hit, ["bitscore", "bit_score", "score", "full_score"]),
            "bias": scalar(hit, ["bias"]),
            "taxid": scalar(hit, ["taxid", "tax_id", "taxonomy_id"]),
            "species": scalar(hit, ["species", "organism", "taxname"]),
            "domain_count": scalar(hit, ["ndom", "n_domains", "domain_count"]),
            "raw_json": json.dumps(hit, sort_keys=True, separators=(",", ":")),
        }
        key = json.dumps(row, sort_keys=True)
        if key not in seen:
            rows.append(row)
            seen.add(key)
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--ko", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--family", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--out-root", default="proto_mito_results/hmmer")
    args = parser.parse_args()

    out = Path(args.out_root) / args.label
    out.mkdir(parents=True, exist_ok=True)
    alignment_path = Path(args.alignment)
    alignment = alignment_path.read_text(encoding="utf-8")
    metadata = {
        "label": args.label,
        "KO": args.ko,
        "DB": args.db,
        "family_id": args.family,
        "alignment_path": str(alignment_path),
        "alignment_bytes": len(alignment.encode()),
        "alignment_sequence_count": alignment.count(">"),
        "algorithm": "hmmsearch",
        "target_database": "refprot",
        "submission_endpoint": f"{BASE}/search/hmmsearch",
    }
    write_json(out / "metadata.json", metadata)
    (out / "query_alignment.fasta").write_text(alignment, encoding="utf-8")

    submission_response = request(
        "POST",
        f"{BASE}/search/hmmsearch",
        json_payload={"database": "refprot", "input": alignment},
        timeout=600,
    )
    submission = safe_json(submission_response)
    write_json(out / "submission.json", submission)
    job_id = find_job_id(submission)
    (out / "job_id.txt").write_text(job_id + "\n", encoding="utf-8")

    poll_history: list[dict[str, Any]] = []
    result: Any = None
    for attempt in range(240):
        result_response = request(
            "GET",
            f"{BASE}/result/{job_id}",
            attempts=3,
            timeout=300,
            allow_pending=True,
        )
        result = safe_json(result_response)
        statuses = all_statuses(result)
        poll_history.append(
            {
                "attempt": attempt + 1,
                "http_status": result_response.status_code,
                "statuses": statuses,
                "response_bytes": len(result_response.content),
            }
        )
        if result_response.status_code in {200, 201, 202} and looks_complete(result):
            break
        time.sleep(10)
    else:
        write_json(out / "poll_history.json", poll_history)
        write_json(out / "result_last.json", result)
        raise TimeoutError(f"HMMER job {job_id} did not complete within polling window")

    write_json(out / "poll_history.json", poll_history)
    write_json(out / "result.json", result)
    search_details = fetch_optional(f"{BASE}/search/{job_id}")
    query = fetch_optional(f"{BASE}/search/{job_id}/query")
    taxonomy = fetch_optional(f"{BASE}/taxonomy/{job_id}/tree")
    taxonomy_distribution = fetch_optional(f"{BASE}/taxonomy/{job_id}/distribution")
    domains = fetch_optional(f"{BASE}/result/{job_id}/domains")
    architecture = fetch_optional(f"{BASE}/architecture/{job_id}")
    downloads = fetch_optional(f"{BASE}/download/{job_id}")
    for filename, value in [
        ("search_details.json", search_details),
        ("query.json", query),
        ("taxonomy_tree.json", taxonomy),
        ("taxonomy_distribution.json", taxonomy_distribution),
        ("domains.json", domains),
        ("architecture.json", architecture),
        ("downloads.json", downloads),
    ]:
        write_json(out / filename, value)

    hit_rows = flatten_hits(result)
    hit_fields = [
        "json_path", "accession", "name", "description", "evalue", "bitscore",
        "bias", "taxid", "species", "domain_count", "raw_json",
    ]
    write_tsv(out / "hits_flat.tsv", hit_rows, hit_fields)
    summary = {
        **metadata,
        "job_id": job_id,
        "poll_count": len(poll_history),
        "result_statuses": all_statuses(result),
        "flattened_hit_rows": len(hit_rows),
        "result_top_level_keys": sorted(result.keys()) if isinstance(result, dict) else [],
        "taxonomy_available": not (isinstance(taxonomy, dict) and "_error" in taxonomy),
        "domains_available": not (isinstance(domains, dict) and "_error" in domains),
    }
    write_json(out / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
