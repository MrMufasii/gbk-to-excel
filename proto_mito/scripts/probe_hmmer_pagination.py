from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import requests

OUT = Path("proto_mito_results/hmmer_pagination_probe")
OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
JOB = "eba9fb1a-1a27-4695-9d10-973bd57a4d99"  # UQCRQ_A, immutable completed job
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "proto-mito-hmmer-page-probe/1.0", "Accept": "application/json"})


def get(url: str, params: dict[str, Any] | None = None) -> requests.Response:
    delay = 1.0
    last = None
    for _ in range(8):
        try:
            response = SESSION.get(url, params=params, timeout=300)
            if response.status_code in {200, 201, 202}:
                return response
            last = f"{response.status_code} {response.text[:500]}"
            if response.status_code not in {429, 500, 502, 503, 504}:
                return response
        except Exception as exc:
            last = repr(exc)
        time.sleep(delay)
        delay = min(delay * 1.8, 30)
    raise RuntimeError(last)


def metadata(payload: Any) -> dict[str, Any]:
    result = payload.get("result") if isinstance(payload, dict) else None
    hits = result.get("hits") if isinstance(result, dict) else None
    rows = []
    if isinstance(hits, list):
        for hit in hits[:3]:
            md = hit.get("metadata") or {}
            rows.append({
                "index": hit.get("index"),
                "name": hit.get("name"),
                "accession": md.get("accession"),
                "description": md.get("description"),
                "kingdom": md.get("kingdom"),
                "species": md.get("species"),
                "evalue": hit.get("evalue"),
                "score": hit.get("score"),
            })
    return {
        "status": payload.get("status") if isinstance(payload, dict) else None,
        "page_count": payload.get("page_count") if isinstance(payload, dict) else None,
        "hit_count_returned": len(hits) if isinstance(hits, list) else None,
        "first_hits": rows,
    }

# Preserve the live OpenAPI schema for exact parameter discovery.
openapi_response = get(f"{BASE}/openapi.json")
(OUT / "openapi.json").write_bytes(openapi_response.content)

probes = {
    "default": {},
    "page_1": {"page": 1},
    "page_2": {"page": 2},
    "page_0": {"page": 0},
    "page_2_size_50": {"page": 2, "size": 50},
    "page_2_page_size_50": {"page": 2, "page_size": 50},
    "page_2_per_page_50": {"page": 2, "per_page": 50},
    "offset_50": {"offset": 50},
    "start_50_size_50": {"start": 50, "size": 50},
    "range_50_99": {"range": "50,99"},
}
report: dict[str, Any] = {}
for label, params in probes.items():
    response = get(f"{BASE}/result/{JOB}", params=params)
    try:
        payload = response.json()
    except Exception:
        payload = {"_text": response.text[:10000]}
    (OUT / f"{label}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    report[label] = {
        "url": response.url,
        "http_status": response.status_code,
        **metadata(payload),
    }

# Taxonomy and architecture views were still materializing immediately after
# the search. Re-query them now and record both status and payload.
for name, suffix in {
    "taxonomy_tree": f"taxonomy/{JOB}/tree",
    "taxonomy_distribution": f"taxonomy/{JOB}/distribution",
    "architecture": f"architecture/{JOB}",
}.items():
    response = get(f"{BASE}/{suffix}")
    try:
        payload = response.json()
    except Exception:
        payload = {"_text": response.text[:10000]}
    (OUT / f"{name}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    report[name] = {"url": response.url, "http_status": response.status_code, "payload_keys": list(payload) if isinstance(payload, dict) else []}

(OUT / "probe_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
