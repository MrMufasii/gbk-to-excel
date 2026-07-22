from __future__ import annotations

import csv
import json
import re
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable

import requests

OUT = Path("proto_mito_results/uqcrq_cluster_members")
OUT.mkdir(parents=True, exist_ok=True)

CLUSTERS = {
    "UQCRQ_plant_mixed_63": "79d76a05d243053e6f72ccebe6eb036c",
    "UQCRQ_plant_mixed_189": "a0b6c89f5e017309e789c417ab20fd64",
    "UQCRQ_broad_mixed_347": "af38011d514e4d32e8b5547884d40b31",
}
ATLAS_BASE = "https://biohub.ai/esm/protein/api/v1alpha1"
UNIPROT_BASE = "https://rest.uniprot.org"


def new_session() -> requests.Session:
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "proto-mito-uqcrq-member-audit/1.0 (public reproducible scientific workflow)",
            "Accept": "application/json",
        }
    )
    return session


def get_json(url: str, params: dict[str, Any] | None = None, attempts: int = 12) -> dict[str, Any]:
    session = new_session()
    delay = 1.0
    last = ""
    for _ in range(attempts):
        try:
            response = session.get(url, params=params, timeout=300)
            last = f"{response.status_code} {response.text[:500]}"
        except Exception as exc:
            last = repr(exc)
            time.sleep(delay)
            delay = min(delay * 1.8, 45)
            continue
        if response.status_code == 200:
            payload = response.json()
            if isinstance(payload, dict):
                return payload
            last = f"non-object JSON {type(payload).__name__}"
        elif response.status_code not in {429, 500, 502, 503, 504}:
            raise RuntimeError(f"GET {response.url} failed: {last}")
        time.sleep(delay)
        delay = min(delay * 1.8, 45)
    raise RuntimeError(f"GET failed after retries: {url}; last={last}")


def scalar(value: Any, names: set[str]) -> str:
    if not isinstance(value, dict):
        return ""
    for key, child in value.items():
        if str(key).lower() in names and isinstance(child, (str, int, float, bool)):
            return str(child)
    return ""


def recursive_dicts(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from recursive_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from recursive_dicts(child)


def extract_taxon_pairs(value: Any) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for item in recursive_dicts(value):
        name = scalar(
            item,
            {
                "scientificname", "scientific_name", "organismname", "organism_name",
                "taxonname", "taxon_name", "organism",
            },
        )
        taxid = scalar(item, {"taxonid", "taxon_id", "taxonomyid", "taxonomy_id", "organismid", "organism_id"})
        if name or taxid:
            pairs.add((taxid, name))
    return pairs


def extract_crossrefs(value: dict[str, Any]) -> list[dict[str, str]]:
    refs: list[dict[str, str]] = []
    lists = []
    for key in ["uniParcCrossReferences", "crossReferences", "dbReferences"]:
        child = value.get(key)
        if isinstance(child, list):
            lists.extend(child)
    for item in lists:
        if not isinstance(item, dict):
            continue
        database = scalar(item, {"database", "databaseid", "database_id", "type"})
        accession = scalar(item, {"id", "accession", "primaryid", "primary_id"})
        active = scalar(item, {"active"})
        version = scalar(item, {"version", "sequenceversion", "sequence_version"})
        taxa = sorted(extract_taxon_pairs(item))
        refs.append(
            {
                "database": database,
                "accession": accession,
                "active": active,
                "version": version,
                "taxa": json.dumps(taxa, separators=(",", ":")),
            }
        )
    return refs


def lineage_names(taxonomy: dict[str, Any]) -> list[str]:
    names: list[str] = []
    scientific = taxonomy.get("scientificName")
    if isinstance(scientific, str):
        names.append(scientific)
    lineage = taxonomy.get("lineage")
    if isinstance(lineage, list):
        for item in lineage:
            if isinstance(item, dict):
                name = item.get("scientificName") or item.get("name")
                if isinstance(name, str):
                    names.append(name)
            elif isinstance(item, str):
                names.append(item)
    return names


def classify_lineage(names: Iterable[str]) -> str:
    lowered = {name.lower() for name in names if name}
    if "bacteria" in lowered:
        return "Bacteria"
    if "archaea" in lowered:
        return "Archaea"
    if "eukaryota" in lowered:
        return "Eukaryota"
    if "viruses" in lowered or "virus" in lowered:
        return "Viruses"
    return "Unresolved"


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# Retrieve cluster membership exactly as represented by the current Atlas release.
cluster_payloads: dict[str, dict[str, Any]] = {}
membership: list[dict[str, str]] = []
for label, cluster_hash in CLUSTERS.items():
    payload = get_json(f"{ATLAS_BASE}/clusters/{cluster_hash}", {"topk_features": 10})
    cluster_payloads[cluster_hash] = payload
    members = payload.get("member_protein_hashes") or []
    for member_hash in members:
        membership.append({"cluster_label": label, "cluster_hash": cluster_hash, "member_hash": str(member_hash)})
    (OUT / f"cluster_{cluster_hash}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

write_tsv(OUT / "cluster_membership.tsv", membership, ["cluster_label", "cluster_hash", "member_hash"])
member_to_clusters: dict[str, list[dict[str, str]]] = defaultdict(list)
for row in membership:
    member_to_clusters[row["member_hash"]].append(row)


def fetch_atlas_member(member_hash: str) -> tuple[str, dict[str, Any] | None, str | None]:
    try:
        payload = get_json(f"{ATLAS_BASE}/proteins/{member_hash}", {"topk_features": 0, "fold_on_miss": "false"})
        return member_hash, payload, None
    except Exception as exc:
        return member_hash, None, repr(exc)


atlas_members: dict[str, dict[str, Any]] = {}
failures: list[dict[str, str]] = []
with ThreadPoolExecutor(max_workers=6) as pool:
    futures = {pool.submit(fetch_atlas_member, member_hash): member_hash for member_hash in member_to_clusters}
    for index, future in enumerate(as_completed(futures), 1):
        member_hash, payload, error = future.result()
        if payload is not None:
            atlas_members[member_hash] = payload
        else:
            failures.append({"stage": "atlas_protein", "identifier": member_hash, "error": error or "unknown"})
        if index % 50 == 0:
            print(f"Atlas members resolved: {index}/{len(futures)}")


def fetch_source_record(accession: str, source: str) -> tuple[str, dict[str, Any] | None, str | None]:
    try:
        if accession.startswith("UPI"):
            url = f"{UNIPROT_BASE}/uniparc/{accession}"
        elif accession and source.lower() in {"uniprot", "uniprotkb", "swissprot", "trembl"}:
            url = f"{UNIPROT_BASE}/uniprotkb/{accession}"
        elif accession and re.fullmatch(r"[A-Z0-9]{6,10}(?:-\d+)?", accession):
            url = f"{UNIPROT_BASE}/uniprotkb/{accession}"
        else:
            return accession, None, "unsupported source/accession"
        return accession, get_json(url), None
    except Exception as exc:
        return accession, None, repr(exc)


accession_source = {
    str(payload.get("accession") or ""): str(payload.get("source") or "")
    for payload in atlas_members.values()
    if payload.get("accession")
}
source_records: dict[str, dict[str, Any]] = {}
with ThreadPoolExecutor(max_workers=6) as pool:
    futures = {
        pool.submit(fetch_source_record, accession, source): accession
        for accession, source in accession_source.items()
    }
    for index, future in enumerate(as_completed(futures), 1):
        accession, payload, error = future.result()
        if payload is not None:
            source_records[accession] = payload
        elif error != "unsupported source/accession":
            failures.append({"stage": "source_record", "identifier": accession, "error": error or "unknown"})
        if index % 50 == 0:
            print(f"Source records resolved: {index}/{len(futures)}")

# Collect taxon IDs from source records, then retrieve authoritative lineages.
member_taxa: dict[str, set[tuple[str, str]]] = defaultdict(set)
member_refs: dict[str, list[dict[str, str]]] = defaultdict(list)
for member_hash, protein in atlas_members.items():
    accession = str(protein.get("accession") or "")
    record = source_records.get(accession, {})
    member_taxa[member_hash].update(extract_taxon_pairs(record))
    if accession.startswith("UPI"):
        member_refs[member_hash] = extract_crossrefs(record)
    else:
        organism = record.get("organism") if isinstance(record, dict) else None
        if isinstance(organism, dict):
            name = str(organism.get("scientificName") or "")
            taxid = str(organism.get("taxonId") or "")
            if name or taxid:
                member_taxa[member_hash].add((taxid, name))

unique_taxids = sorted({taxid for pairs in member_taxa.values() for taxid, _ in pairs if taxid.isdigit()})

def fetch_taxonomy(taxid: str) -> tuple[str, dict[str, Any] | None, str | None]:
    try:
        return taxid, get_json(f"{UNIPROT_BASE}/taxonomy/{taxid}"), None
    except Exception as exc:
        return taxid, None, repr(exc)


taxonomy_records: dict[str, dict[str, Any]] = {}
with ThreadPoolExecutor(max_workers=6) as pool:
    futures = {pool.submit(fetch_taxonomy, taxid): taxid for taxid in unique_taxids}
    for index, future in enumerate(as_completed(futures), 1):
        taxid, payload, error = future.result()
        if payload is not None:
            taxonomy_records[taxid] = payload
        else:
            failures.append({"stage": "taxonomy", "identifier": taxid, "error": error or "unknown"})
        if index % 50 == 0:
            print(f"Taxonomy records resolved: {index}/{len(futures)}")

member_rows: list[dict[str, Any]] = []
source_ref_rows: list[dict[str, Any]] = []
bacterial_fasta: list[tuple[str, str]] = []
for member_hash, memberships in member_to_clusters.items():
    protein = atlas_members.get(member_hash, {})
    accession = str(protein.get("accession") or "")
    taxa = sorted(member_taxa.get(member_hash, set()))
    lineage_by_taxid: dict[str, list[str]] = {}
    domain_classes: set[str] = set()
    for taxid, name in taxa:
        taxonomy = taxonomy_records.get(taxid, {})
        names = lineage_names(taxonomy)
        if not names and name:
            names = [name]
        lineage_by_taxid[taxid or name] = names
        domain_classes.add(classify_lineage(names))
    domain_classes.discard("Unresolved")
    if len(domain_classes) == 1:
        domain_class = next(iter(domain_classes))
    elif len(domain_classes) > 1:
        domain_class = "Cross-domain sequence"
    else:
        domain_class = "Unresolved"
    sequence = str(protein.get("sequence") or "")
    for membership_row in memberships:
        member_rows.append(
            {
                **membership_row,
                "accession": accession,
                "source": protein.get("source", ""),
                "sequence_length": protein.get("sequence_length", len(sequence)),
                "sequence": sequence,
                "domain_class": domain_class,
                "taxon_pairs": json.dumps(taxa, separators=(",", ":")),
                "lineages": json.dumps(lineage_by_taxid, sort_keys=True, separators=(",", ":")),
                "cross_reference_count": len(member_refs.get(member_hash, [])),
            }
        )
    if domain_class in {"Bacteria", "Archaea", "Cross-domain sequence"} and sequence:
        bacterial_fasta.append((f"{member_hash}|{accession}|{domain_class}", sequence))
    for ref in member_refs.get(member_hash, []):
        source_ref_rows.append({"member_hash": member_hash, "upi": accession, **ref})

member_fields = [
    "cluster_label", "cluster_hash", "member_hash", "accession", "source",
    "sequence_length", "sequence", "domain_class", "taxon_pairs", "lineages",
    "cross_reference_count",
]
write_tsv(OUT / "member_taxonomy.tsv", member_rows, member_fields)
write_tsv(
    OUT / "uniparc_cross_references.tsv",
    source_ref_rows,
    ["member_hash", "upi", "database", "accession", "active", "version", "taxa"],
)
write_tsv(OUT / "failures.tsv", failures, ["stage", "identifier", "error"])
with (OUT / "prokaryotic_or_cross_domain_members.fasta").open("w", encoding="utf-8") as handle:
    for header, sequence in bacterial_fasta:
        handle.write(f">{header}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start+80] + "\n")

cluster_summary: list[dict[str, Any]] = []
for label, cluster_hash in CLUSTERS.items():
    rows = [row for row in member_rows if row["cluster_hash"] == cluster_hash]
    counts = Counter(row["domain_class"] for row in rows)
    cluster_summary.append(
        {
            "cluster_label": label,
            "cluster_hash": cluster_hash,
            "member_count": len(rows),
            "bacteria": counts.get("Bacteria", 0),
            "archaea": counts.get("Archaea", 0),
            "eukaryota": counts.get("Eukaryota", 0),
            "cross_domain": counts.get("Cross-domain sequence", 0),
            "unresolved": counts.get("Unresolved", 0),
            "atlas_top_phyla": json.dumps(cluster_payloads[cluster_hash].get("top_phyla") or {}, sort_keys=True),
            "atlas_taxonomy": json.dumps(cluster_payloads[cluster_hash].get("cluster_taxonomy_info") or {}, sort_keys=True),
        }
    )
write_tsv(
    OUT / "cluster_domain_summary.tsv",
    cluster_summary,
    [
        "cluster_label", "cluster_hash", "member_count", "bacteria", "archaea",
        "eukaryota", "cross_domain", "unresolved", "atlas_top_phyla", "atlas_taxonomy",
    ],
)

summary = {
    "clusters": len(CLUSTERS),
    "membership_rows": len(membership),
    "unique_members": len(member_to_clusters),
    "atlas_members_resolved": len(atlas_members),
    "source_records_resolved": len(source_records),
    "taxon_ids_resolved": len(taxonomy_records),
    "domain_counts": dict(Counter(row["domain_class"] for row in member_rows)),
    "prokaryotic_or_cross_domain_fasta_records": len(bacterial_fasta),
    "failures": len(failures),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
if len(atlas_members) < len(member_to_clusters) - 3:
    raise RuntimeError("Too many Atlas member lookups failed")
