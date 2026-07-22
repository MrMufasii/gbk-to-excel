from __future__ import annotations

import csv
import json
import re
import time
from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Any

import requests
from Bio import SeqIO

OUT = Path("proto_mito_results/candidate_panels")
RAW = OUT / "raw"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

UNIPROT = "https://rest.uniprot.org"
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ATLAS = "https://biohub.ai/esm/protein/api/v1alpha1"
SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "proto-mito-candidate-panel/1.0 (public reproducible scientific workflow)",
        "Accept": "*/*",
    }
)

TARGET_ACCESSIONS = [
    "A0A5B2VBN4",  # NDUFA2/PF05047 bacterial candidate
    "A0A202DQE6", "A0A202DNV3", "A0A202DMW5", "A0A202DQA8",  # ATP23/M76 candidates
    "A0A656D6G0", "A0A2Z5IPW5",  # repeat-rich ATP5H-d negative controls
]

QUERIES = {
    "PF05047_bacteria": "(xref:pfam-PF05047) AND (taxonomy_id:2)",
    "PF05047_archaea": "(xref:pfam-PF05047) AND (taxonomy_id:2157)",
    "PF05047_eukaryota": "(xref:pfam-PF05047) AND (taxonomy_id:2759)",
    "NDUFA2_eukaryota": "(gene_exact:NDUFA2) AND (taxonomy_id:2759)",
    "PF09768_bacteria": "(xref:pfam-PF09768) AND (taxonomy_id:2)",
    "PF09768_archaea": "(xref:pfam-PF09768) AND (taxonomy_id:2157)",
    "PF09768_eukaryota": "(xref:pfam-PF09768) AND (taxonomy_id:2759)",
    "ATP23_eukaryota": "(gene_exact:ATP23) AND (taxonomy_id:2759)",
}

FIELDS = [
    "accession", "id", "reviewed", "protein_name", "gene_names", "organism_name",
    "organism_id", "lineage", "length", "sequence", "fragment", "xref_pfam",
    "xref_refseq", "xref_embl", "xref_proteomes",
]


def get(url: str, *, params: dict[str, Any] | None = None, attempts: int = 10, timeout: int = 300) -> requests.Response:
    delay = 1.0
    last = ""
    for _ in range(attempts):
        try:
            response = SESSION.get(url, params=params, timeout=timeout)
            last = f"{response.status_code} {response.text[:500]}"
        except Exception as exc:
            last = repr(exc)
            time.sleep(delay)
            delay = min(delay * 1.8, 45)
            continue
        if response.status_code == 200:
            return response
        if response.status_code in {429, 500, 502, 503, 504}:
            time.sleep(delay)
            delay = min(delay * 1.8, 45)
            continue
        raise RuntimeError(f"GET {response.url} failed: {last}")
    raise RuntimeError(f"GET failed after retries: {url}; last={last}")


def get_json(url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    payload = get(url, params=params).json()
    if not isinstance(payload, dict):
        raise RuntimeError(f"Expected JSON object from {url}, received {type(payload).__name__}")
    return payload


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str), encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def next_link(headers: requests.structures.CaseInsensitiveDict[str]) -> str | None:
    link = headers.get("Link", "")
    for item in link.split(","):
        if 'rel="next"' in item:
            match = re.search(r"<([^>]+)>", item)
            if match:
                return match.group(1)
    return None


def uniprot_search(label: str, query: str, max_records: int = 6000) -> list[dict[str, str]]:
    url = f"{UNIPROT}/uniprotkb/search"
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(FIELDS),
        "size": 500,
    }
    rows: list[dict[str, str]] = []
    page = 0
    while url and len(rows) < max_records:
        response = get(url, params=params if page == 0 else None)
        text = response.text
        (RAW / f"uniprot_{label}_page_{page+1:03d}.tsv").write_text(text, encoding="utf-8")
        parsed = list(csv.DictReader(text.splitlines(), delimiter="\t"))
        if page > 0 and parsed and rows and parsed[0].keys() != rows[0].keys():
            raise RuntimeError(f"UniProt page schema changed for {label}")
        rows.extend(parsed)
        url = next_link(response.headers)
        params = None
        page += 1
    rows = rows[:max_records]
    for row in rows:
        row["panel"] = label
        row["query"] = query
    return rows


def entry_sequence(entry: dict[str, Any]) -> str:
    sequence = entry.get("sequence") or {}
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(sequence.get("value") or "").upper())


def entry_crossrefs(entry: dict[str, Any]) -> list[dict[str, Any]]:
    refs = entry.get("uniProtKBCrossReferences") or []
    return [ref for ref in refs if isinstance(ref, dict)]


def properties(ref: dict[str, Any]) -> dict[str, str]:
    values: dict[str, str] = {}
    for prop in ref.get("properties") or []:
        if isinstance(prop, dict):
            values[str(prop.get("key") or "")] = str(prop.get("value") or "")
    return values


def ncbi_efetch(db: str, accession: str, rettype: str = "gb") -> str:
    response = get(
        f"{NCBI}/efetch.fcgi",
        params={"db": db, "id": accession, "rettype": rettype, "retmode": "text"},
    )
    time.sleep(0.36)
    return response.text


def parse_genbank(text: str):
    records = list(SeqIO.parse(StringIO(text), "genbank"))
    return records[0] if records else None


def coded_by(record) -> list[str]:
    if record is None:
        return []
    values: list[str] = []
    for feature in record.features:
        values.extend(str(x) for x in feature.qualifiers.get("coded_by", []))
    return values


def nucleotide_accessions(expressions: list[str]) -> list[str]:
    result: list[str] = []
    for expression in expressions:
        result.extend(re.findall(r"([A-Z]{1,6}_?[A-Z0-9]*\d+(?:\.\d+)?)\s*:", expression))
    return list(dict.fromkeys(result))


def feature_qualifier(feature, name: str) -> str:
    return ";".join(str(value) for value in feature.qualifiers.get(name, []))


def feature_bounds(feature) -> tuple[int, int, str]:
    return int(feature.location.start) + 1, int(feature.location.end), "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "?"


def source_metadata(record) -> dict[str, str]:
    if record is None:
        return {}
    for feature in record.features:
        if feature.type == "source":
            return {key: ";".join(map(str, values)) for key, values in feature.qualifiers.items()}
    return {}


def provenance(accession: str, entry: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    refs = entry_crossrefs(entry)
    candidates: list[str] = []
    for ref in refs:
        db = str(ref.get("database") or "")
        ref_id = str(ref.get("id") or "")
        props = properties(ref)
        if db in {"RefSeq", "EMBL", "EMBL-CDS"}:
            for value in [ref_id, props.get("ProteinId", "")]:
                if value and re.fullmatch(r"[A-Z]{1,6}_?[A-Z0-9]+(?:\.\d+)?", value):
                    candidates.append(value)
    candidates = list(dict.fromkeys(candidates))
    protein_records = []
    chosen_record = None
    chosen_accession = ""
    for candidate in candidates[:12]:
        try:
            text = ncbi_efetch("protein", candidate)
            (RAW / f"ncbi_protein_{accession}_{candidate.replace('.', '_')}.gb").write_text(text, encoding="utf-8")
            record = parse_genbank(text)
            protein_records.append({"accession": candidate, "retrieved": bool(record), "description": record.description if record else ""})
            if chosen_record is None and record is not None:
                chosen_record = record
                chosen_accession = candidate
        except Exception as exc:
            protein_records.append({"accession": candidate, "retrieved": False, "error": repr(exc)})

    coded = coded_by(chosen_record)
    nt_accessions = nucleotide_accessions(coded)
    contigs: list[dict[str, Any]] = []
    neighbourhoods: list[dict[str, Any]] = []
    for nt_accession in nt_accessions[:6]:
        try:
            text = ncbi_efetch("nuccore", nt_accession, "gbwithparts")
            (RAW / f"ncbi_nucleotide_{accession}_{nt_accession.replace('.', '_')}.gb").write_text(text, encoding="utf-8")
            record = parse_genbank(text)
            if record is None:
                continue
            source = source_metadata(record)
            target_features = []
            candidate_roots = {x.split(".")[0] for x in candidates}
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                protein_ids = {x.split(".")[0] for x in feature.qualifiers.get("protein_id", [])}
                if protein_ids & candidate_roots:
                    target_features.append(feature)
            target_start = target_end = None
            if target_features:
                target_start, target_end, _ = feature_bounds(target_features[0])
            for feature in record.features:
                if feature.type not in {"CDS", "tRNA", "rRNA", "ncRNA", "repeat_region", "misc_feature"}:
                    continue
                start, end, strand = feature_bounds(feature)
                if target_start is not None:
                    distance = 0 if not (end < target_start or start > target_end) else min(abs(end-target_start), abs(start-target_end))
                    if distance > 30000:
                        continue
                else:
                    distance = ""
                neighbourhoods.append(
                    {
                        "target_uniprot": accession,
                        "protein_record": chosen_accession,
                        "nucleotide_accession": record.id,
                        "target_start": target_start or "",
                        "target_end": target_end or "",
                        "feature_type": feature.type,
                        "feature_start": start,
                        "feature_end": end,
                        "feature_strand": strand,
                        "distance_to_target_bp": distance,
                        "protein_id": feature_qualifier(feature, "protein_id"),
                        "locus_tag": feature_qualifier(feature, "locus_tag"),
                        "gene": feature_qualifier(feature, "gene"),
                        "product": feature_qualifier(feature, "product"),
                        "note": feature_qualifier(feature, "note"),
                        "inference": feature_qualifier(feature, "inference"),
                    }
                )
            contigs.append(
                {
                    "accession": record.id,
                    "description": record.description,
                    "length": len(record.seq),
                    "annotations": record.annotations,
                    "source": source,
                    "target_start": target_start,
                    "target_end": target_end,
                    "nearby_feature_count": sum(row["nucleotide_accession"] == record.id for row in neighbourhoods),
                }
            )
        except Exception as exc:
            contigs.append({"accession": nt_accession, "error": repr(exc)})

    summary = {
        "target_uniprot": accession,
        "candidate_protein_accessions": candidates,
        "protein_record_attempts": protein_records,
        "chosen_protein_record": chosen_accession,
        "coded_by": coded,
        "nucleotide_accessions": nt_accessions,
        "contigs": contigs,
    }
    return summary, neighbourhoods


def atlas_context(accession: str, sequence: str) -> dict[str, Any]:
    search = get_json(
        f"{ATLAS}/similarity-search",
        {
            "sequence": sequence[:800],
            "topk_results": 30,
            "topk_features": 20,
            "include_cluster_info": "true",
        },
    )
    query_hash = str(search.get("protein_hash") or "")
    result: dict[str, Any] = {"accession": accession, "query_hash": query_hash, "search": search}
    if query_hash:
        protein = get_json(f"{ATLAS}/proteins/{query_hash}", {"topk_features": 20, "fold_on_miss": "true"})
        result["protein"] = protein
        cluster_hash = str(protein.get("cluster_rep_protein_hash") or query_hash)
        result["cluster_hash"] = cluster_hash
        try:
            result["cluster"] = get_json(f"{ATLAS}/clusters/{cluster_hash}", {"topk_features": 20})
        except Exception as exc:
            result["cluster_error"] = repr(exc)
        pdb = protein.get("pdb")
        if isinstance(pdb, str) and pdb.strip():
            (OUT / f"atlas_{accession}.pdb").write_text(pdb, encoding="utf-8")
        sanitized = dict(protein)
        for key in ["pdb", "residues_plddt", "per_residue_activations", "protein_activations"]:
            sanitized.pop(key, None)
        result["protein"] = sanitized
    return result


# Exact entries, provenance and Atlas context.
entry_rows: list[dict[str, Any]] = []
provenance_rows: list[dict[str, Any]] = []
neighbourhood_rows: list[dict[str, Any]] = []
atlas_rows: list[dict[str, Any]] = []
target_sequences: list[tuple[str, str]] = []
failures: list[dict[str, str]] = []
for accession in TARGET_ACCESSIONS:
    try:
        entry = get_json(f"{UNIPROT}/uniprotkb/{accession}.json")
        write_json(RAW / f"uniprot_{accession}.json", entry)
        sequence = entry_sequence(entry)
        organism = entry.get("organism") or {}
        description = entry.get("proteinDescription") or {}
        entry_rows.append(
            {
                "accession": accession,
                "entry_type": entry.get("entryType", ""),
                "organism": organism.get("scientificName", ""),
                "taxon_id": organism.get("taxonId", ""),
                "lineage": ";".join(organism.get("lineage") or []),
                "sequence_length": len(sequence),
                "sequence": sequence,
                "protein_description": json.dumps(description, sort_keys=True),
                "fragment": bool((entry.get("sequence") or {}).get("molWeight") is None) or "fragment" in json.dumps(entry).lower(),
                "annotation_score": entry.get("annotationScore", ""),
                "protein_existence": entry.get("proteinExistence", ""),
            }
        )
        target_sequences.append((accession, sequence))
        provenance_summary, neighbours = provenance(accession, entry)
        provenance_rows.append(provenance_summary)
        neighbourhood_rows.extend(neighbours)
        try:
            atlas = atlas_context(accession, sequence)
            write_json(RAW / f"atlas_{accession}.json", atlas)
            cluster = atlas.get("cluster") or {}
            tax = cluster.get("cluster_taxonomy_info") or {}
            atlas_rows.append(
                {
                    "accession": accession,
                    "query_hash": atlas.get("query_hash", ""),
                    "cluster_hash": atlas.get("cluster_hash", ""),
                    "cluster_size": cluster.get("cluster_size", ""),
                    "cluster_taxonomy_rank": tax.get("rank", "") if isinstance(tax, dict) else "",
                    "cluster_taxonomy_name": tax.get("name", "") if isinstance(tax, dict) else "",
                    "top_phyla": json.dumps(cluster.get("top_phyla") or {}, sort_keys=True),
                    "top_pfams": json.dumps(cluster.get("cluster_top_pfam_domains") or {}, sort_keys=True),
                    "pdb_available": (OUT / f"atlas_{accession}.pdb").exists(),
                }
            )
        except Exception as exc:
            failures.append({"target": accession, "stage": "atlas", "error": repr(exc)})
    except Exception as exc:
        failures.append({"target": accession, "stage": "uniprot_entry", "error": repr(exc)})

entry_fields = sorted({key for row in entry_rows for key in row})
write_tsv(OUT / "target_entries.tsv", entry_rows, entry_fields)
provenance_fields = sorted({key for row in provenance_rows for key in row})
write_tsv(OUT / "target_provenance.tsv", provenance_rows, provenance_fields)
write_tsv(
    OUT / "target_neighbourhoods_30kb.tsv",
    neighbourhood_rows,
    [
        "target_uniprot", "protein_record", "nucleotide_accession", "target_start", "target_end",
        "feature_type", "feature_start", "feature_end", "feature_strand", "distance_to_target_bp",
        "protein_id", "locus_tag", "gene", "product", "note", "inference",
    ],
)
write_tsv(
    OUT / "target_atlas_context.tsv",
    atlas_rows,
    [
        "accession", "query_hash", "cluster_hash", "cluster_size", "cluster_taxonomy_rank",
        "cluster_taxonomy_name", "top_phyla", "top_pfams", "pdb_available",
    ],
)
with (OUT / "target_candidates.fasta").open("w", encoding="utf-8") as handle:
    for accession, sequence in target_sequences:
        handle.write(f">{accession}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start+80] + "\n")

# Current UniProt family panels.
panel_rows: list[dict[str, str]] = []
query_summaries: list[dict[str, Any]] = []
for label, query in QUERIES.items():
    try:
        rows = uniprot_search(label, query)
        panel_rows.extend(rows)
        query_summaries.append({"panel": label, "query": query, "records": len(rows), "status": "success"})
    except Exception as exc:
        failures.append({"target": label, "stage": "uniprot_search", "error": repr(exc)})
        query_summaries.append({"panel": label, "query": query, "records": 0, "status": repr(exc)})

# Normalize TSV headers from UniProt into stable machine labels while retaining raw values.
normalized: list[dict[str, Any]] = []
for row in panel_rows:
    normalized.append(
        {
            "panel": row.get("panel", ""),
            "query": row.get("query", ""),
            "accession": row.get("Entry", ""),
            "entry_name": row.get("Entry Name", ""),
            "reviewed": row.get("Reviewed", ""),
            "protein_name": row.get("Protein names", ""),
            "gene_names": row.get("Gene Names", ""),
            "organism": row.get("Organism", ""),
            "taxon_id": row.get("Organism (ID)", ""),
            "lineage": row.get("Taxonomic lineage", ""),
            "length": row.get("Length", ""),
            "sequence": re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", row.get("Sequence", "").upper()),
            "fragment": row.get("Fragment", ""),
            "pfam": row.get("Pfam", ""),
            "refseq": row.get("RefSeq", ""),
            "embl": row.get("EMBL", ""),
            "proteomes": row.get("Proteomes", ""),
        }
    )

panel_fields = list(normalized[0]) if normalized else []
write_tsv(OUT / "uniprot_family_panels.tsv", normalized, panel_fields)
write_tsv(OUT / "query_summaries.tsv", query_summaries, ["panel", "query", "records", "status"])
for label in QUERIES:
    values = [row for row in normalized if row["panel"] == label and row["sequence"]]
    with (OUT / f"panel_{label}.fasta").open("w", encoding="utf-8") as handle:
        for row in values:
            safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", row["entry_name"] or row["accession"])
            handle.write(f">{label}|{row['accession']}|{safe_name}|taxid={row['taxon_id']}\n")
            sequence = row["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start+80] + "\n")

write_tsv(OUT / "failures.tsv", failures, ["target", "stage", "error"])
summary = {
    "exact_targets_requested": len(TARGET_ACCESSIONS),
    "exact_targets_resolved": len(entry_rows),
    "provenance_records": len(provenance_rows),
    "neighbourhood_features": len(neighbourhood_rows),
    "atlas_contexts": len(atlas_rows),
    "family_panel_records": len(normalized),
    "panel_counts": {label: sum(row["panel"] == label for row in normalized) for label in QUERIES},
    "failures": len(failures),
}
write_json(OUT / "summary.json", summary)
print(json.dumps(summary, indent=2))
if len(entry_rows) < len(TARGET_ACCESSIONS) - 1:
    raise RuntimeError("Too many exact candidate entries failed")
if any(row["status"] != "success" for row in query_summaries):
    raise RuntimeError("At least one UniProt family query failed")
