from __future__ import annotations

import csv
import json
import re
import time
from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Any, Iterable

import requests
from Bio import SeqIO

OUT = Path("proto_mito_results/atp23_uniparc")
RAW = OUT / "raw"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

UNIPROT = "https://rest.uniprot.org"
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ATLAS = "https://biohub.ai/esm/protein/api/v1alpha1"
SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "proto-mito-atp23-uniparc/1.0 (public reproducible scientific workflow)",
    "Accept": "*/*",
})

TARGETS = {
    "A0A202DQE6": "UPI000B6B2022",
    "A0A202DNV3": "UPI000B6C1C34",
    "A0A202DMW5": "UPI000B75C84A",
    "A0A202DQA8": "UPI000B6C0926",
}


def get(url: str, params: dict[str, Any] | None = None, attempts: int = 12, timeout: int = 300) -> requests.Response:
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
    value = get(url, params=params).json()
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected object from {url}, received {type(value).__name__}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str), encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def recursive_dicts(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from recursive_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from recursive_dicts(child)


def sequence_from_uniparc(payload: dict[str, Any]) -> str:
    for key in ["sequence", "sequenceObject"]:
        value = payload.get(key)
        if isinstance(value, dict):
            sequence = value.get("value") or value.get("sequence")
            if isinstance(sequence, str):
                sequence = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", sequence.upper())
                if sequence:
                    return sequence
        elif isinstance(value, str):
            sequence = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", value.upper())
            if sequence:
                return sequence
    for item in recursive_dicts(payload):
        sequence = item.get("sequence") or item.get("value")
        if isinstance(sequence, str) and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{20,}", sequence.upper()):
            return sequence.upper()
    raise RuntimeError("No amino-acid sequence found in UniParc payload")


def scalar(value: Any, names: set[str]) -> str:
    if not isinstance(value, dict):
        return ""
    for key, child in value.items():
        if str(key).lower() in names and isinstance(child, (str, int, float, bool)):
            return str(child)
    return ""


def extract_crossrefs(payload: dict[str, Any]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for key in ["uniParcCrossReferences", "crossReferences", "dbReferences"]:
        value = payload.get(key)
        if isinstance(value, list):
            result.extend(item for item in value if isinstance(item, dict))
    if result:
        return result
    # Fall back to recursively locating cross-reference-shaped objects.
    for item in recursive_dicts(payload):
        database = scalar(item, {"database", "databaseid", "database_id", "type"})
        accession = scalar(item, {"id", "accession", "primaryid", "primary_id"})
        if database and accession:
            result.append(item)
    return result


def taxon_pairs(value: Any) -> list[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for item in recursive_dicts(value):
        taxid = scalar(item, {"taxonid", "taxon_id", "taxonomyid", "taxonomy_id", "organismid", "organism_id"})
        name = scalar(item, {"scientificname", "scientific_name", "organismname", "organism_name", "taxonname", "taxon_name", "organism"})
        if taxid or name:
            pairs.add((taxid, name))
    return sorted(pairs)


def ncbi_efetch(db: str, accession: str, rettype: str = "gb") -> str:
    response = get(
        f"{NCBI}/efetch.fcgi",
        {"db": db, "id": accession, "rettype": rettype, "retmode": "text"},
    )
    time.sleep(0.36)
    return response.text


def parse_genbank(text: str):
    records = list(SeqIO.parse(StringIO(text), "genbank"))
    return records[0] if records else None


def feature_bounds(feature) -> tuple[int, int, str]:
    return int(feature.location.start) + 1, int(feature.location.end), "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "?"


def qualifier(feature, name: str) -> str:
    return ";".join(str(value) for value in feature.qualifiers.get(name, []))


def source_qualifiers(record) -> dict[str, str]:
    if record is None:
        return {}
    for feature in record.features:
        if feature.type == "source":
            return {key: ";".join(map(str, values)) for key, values in feature.qualifiers.items()}
    return {}


def candidate_source_accessions(crossrefs: list[dict[str, Any]]) -> list[str]:
    result: list[str] = []
    for ref in crossrefs:
        database = scalar(ref, {"database", "databaseid", "database_id", "type"}).lower()
        accession = scalar(ref, {"id", "accession", "primaryid", "primary_id"})
        if database in {"emblwgs", "embl", "refseq", "uniprotkb/trembl", "uniprotkb/swiss-prot", "ensemblbacteria"}:
            if accession and re.fullmatch(r"[A-Z]{1,8}_?[A-Z0-9]+(?:\.\d+)?", accession):
                result.append(accession)
        for item in recursive_dicts(ref):
            for key, value in item.items():
                if str(key).lower() in {"proteinid", "protein_id", "accession", "id"} and isinstance(value, str):
                    if re.fullmatch(r"[A-Z]{1,8}_?[A-Z0-9]+(?:\.\d+)?", value):
                        result.append(value)
    return list(dict.fromkeys(result))


def atlas_context(accession: str, sequence: str) -> dict[str, Any]:
    search = get_json(
        f"{ATLAS}/similarity-search",
        {
            "sequence": sequence[:800],
            "topk_results": 50,
            "topk_features": 20,
            "include_cluster_info": "true",
        },
    )
    query_hash = str(search.get("protein_hash") or "")
    result: dict[str, Any] = {"query_hash": query_hash, "search": search}
    if query_hash:
        protein = get_json(f"{ATLAS}/proteins/{query_hash}", {"topk_features": 20, "fold_on_miss": "true"})
        pdb = protein.get("pdb")
        if isinstance(pdb, str) and pdb.strip():
            (OUT / f"atlas_{accession}.pdb").write_text(pdb, encoding="utf-8")
        cluster_hash = str(protein.get("cluster_rep_protein_hash") or query_hash)
        result["cluster_hash"] = cluster_hash
        try:
            result["cluster"] = get_json(f"{ATLAS}/clusters/{cluster_hash}", {"topk_features": 20})
        except Exception as exc:
            result["cluster_error"] = repr(exc)
        sanitized = dict(protein)
        for key in ["pdb", "residues_plddt", "per_residue_activations", "protein_activations"]:
            sanitized.pop(key, None)
        result["protein"] = sanitized
    return result


entry_rows: list[dict[str, Any]] = []
crossref_rows: list[dict[str, Any]] = []
provenance_rows: list[dict[str, Any]] = []
neighbour_rows: list[dict[str, Any]] = []
atlas_rows: list[dict[str, Any]] = []
failures: list[dict[str, str]] = []
sequences: dict[str, str] = {}

for old_accession, upi in TARGETS.items():
    try:
        payload = get_json(f"{UNIPROT}/uniparc/{upi}")
        write_json(RAW / f"uniparc_{upi}.json", payload)
        sequence = sequence_from_uniparc(payload)
        sequences[old_accession] = sequence
        crossrefs = extract_crossrefs(payload)
        entry_rows.append({
            "old_uniprot_accession": old_accession,
            "uniparc_id": upi,
            "sequence_length": len(sequence),
            "sequence": sequence,
            "taxon_pairs": json.dumps(taxon_pairs(payload), separators=(",", ":")),
            "cross_reference_count": len(crossrefs),
        })
        for ref in crossrefs:
            crossref_rows.append({
                "old_uniprot_accession": old_accession,
                "uniparc_id": upi,
                "database": scalar(ref, {"database", "databaseid", "database_id", "type"}),
                "accession": scalar(ref, {"id", "accession", "primaryid", "primary_id"}),
                "active": scalar(ref, {"active"}),
                "version": scalar(ref, {"version", "sequenceversion", "sequence_version"}),
                "taxon_pairs": json.dumps(taxon_pairs(ref), separators=(",", ":")),
                "raw_json": json.dumps(ref, sort_keys=True, separators=(",", ":")),
            })

        source_accessions = candidate_source_accessions(crossrefs)
        protein_attempts: list[dict[str, Any]] = []
        parent_contigs: list[dict[str, Any]] = []
        for source_accession in source_accessions[:30]:
            # Old UniProt and WGS nucleotide accessions are not valid protein
            # queries; only preserve failures and continue to the next source.
            try:
                text = ncbi_efetch("protein", source_accession)
                (RAW / f"protein_{old_accession}_{source_accession.replace('.', '_')}.gb").write_text(text, encoding="utf-8")
                record = parse_genbank(text)
                if record is None:
                    continue
                protein_attempts.append({"accession": source_accession, "description": record.description, "length": len(record.seq)})
                coded: list[str] = []
                for feature in record.features:
                    coded.extend(str(value) for value in feature.qualifiers.get("coded_by", []))
                nucleotide_accessions = []
                for expression in coded:
                    nucleotide_accessions.extend(re.findall(r"([A-Z]{1,8}_?[A-Z0-9]*\d+(?:\.\d+)?)\s*:", expression))
                for nucleotide_accession in list(dict.fromkeys(nucleotide_accessions)):
                    try:
                        nt_text = ncbi_efetch("nuccore", nucleotide_accession, "gbwithparts")
                        (RAW / f"nucleotide_{old_accession}_{nucleotide_accession.replace('.', '_')}.gb").write_text(nt_text, encoding="utf-8")
                        nt_record = parse_genbank(nt_text)
                        if nt_record is None:
                            continue
                        source = source_qualifiers(nt_record)
                        target_features = []
                        for feature in nt_record.features:
                            if feature.type != "CDS":
                                continue
                            ids = {x.split(".")[0] for x in feature.qualifiers.get("protein_id", [])}
                            if source_accession.split(".")[0] in ids:
                                target_features.append(feature)
                        target_start = target_end = None
                        if target_features:
                            target_start, target_end, _ = feature_bounds(target_features[0])
                        for feature in nt_record.features:
                            if feature.type not in {"CDS", "tRNA", "rRNA", "ncRNA", "repeat_region", "misc_feature"}:
                                continue
                            start, end, strand = feature_bounds(feature)
                            if target_start is None:
                                distance: int | str = ""
                            elif end < target_start:
                                distance = target_start - end
                            elif start > target_end:
                                distance = start - target_end
                            else:
                                distance = 0
                            if distance != "" and int(distance) > 30000:
                                continue
                            neighbour_rows.append({
                                "old_uniprot_accession": old_accession,
                                "uniparc_id": upi,
                                "source_protein_accession": source_accession,
                                "nucleotide_accession": nt_record.id,
                                "target_start": target_start or "",
                                "target_end": target_end or "",
                                "feature_type": feature.type,
                                "feature_start": start,
                                "feature_end": end,
                                "feature_strand": strand,
                                "distance_to_target_bp": distance,
                                "protein_id": qualifier(feature, "protein_id"),
                                "locus_tag": qualifier(feature, "locus_tag"),
                                "gene": qualifier(feature, "gene"),
                                "product": qualifier(feature, "product"),
                                "note": qualifier(feature, "note"),
                                "inference": qualifier(feature, "inference"),
                            })
                        parent_contigs.append({
                            "accession": nt_record.id,
                            "description": nt_record.description,
                            "length": len(nt_record.seq),
                            "source": source,
                            "target_start": target_start,
                            "target_end": target_end,
                        })
                    except Exception as exc:
                        parent_contigs.append({"accession": nucleotide_accession, "error": repr(exc)})
            except Exception as exc:
                protein_attempts.append({"accession": source_accession, "error": repr(exc)})
        provenance_rows.append({
            "old_uniprot_accession": old_accession,
            "uniparc_id": upi,
            "source_accessions": json.dumps(source_accessions),
            "protein_attempts": json.dumps(protein_attempts, default=str),
            "parent_contigs": json.dumps(parent_contigs, default=str),
        })

        try:
            atlas = atlas_context(old_accession, sequence)
            write_json(RAW / f"atlas_{old_accession}.json", atlas)
            cluster = atlas.get("cluster") or {}
            tax = cluster.get("cluster_taxonomy_info") or {}
            atlas_rows.append({
                "old_uniprot_accession": old_accession,
                "uniparc_id": upi,
                "query_hash": atlas.get("query_hash", ""),
                "cluster_hash": atlas.get("cluster_hash", ""),
                "cluster_size": cluster.get("cluster_size", ""),
                "taxonomy_rank": tax.get("rank", "") if isinstance(tax, dict) else "",
                "taxonomy_name": tax.get("name", "") if isinstance(tax, dict) else "",
                "top_phyla": json.dumps(cluster.get("top_phyla") or {}, sort_keys=True),
                "top_pfams": json.dumps(cluster.get("cluster_top_pfam_domains") or {}, sort_keys=True),
                "pdb_available": (OUT / f"atlas_{old_accession}.pdb").exists(),
            })
        except Exception as exc:
            failures.append({"target": old_accession, "stage": "atlas", "error": repr(exc)})
    except Exception as exc:
        failures.append({"target": old_accession, "stage": "uniparc", "error": repr(exc)})

with (OUT / "atp23_targets.fasta").open("w", encoding="utf-8") as handle:
    for accession, sequence in sequences.items():
        handle.write(f">{accession}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start+80] + "\n")

entry_fields = list(entry_rows[0]) if entry_rows else []
write_tsv(OUT / "atp23_uniparc_entries.tsv", entry_rows, entry_fields)
write_tsv(
    OUT / "atp23_cross_references.tsv", crossref_rows,
    ["old_uniprot_accession", "uniparc_id", "database", "accession", "active", "version", "taxon_pairs", "raw_json"],
)
write_tsv(
    OUT / "atp23_provenance.tsv", provenance_rows,
    ["old_uniprot_accession", "uniparc_id", "source_accessions", "protein_attempts", "parent_contigs"],
)
write_tsv(
    OUT / "atp23_neighbourhoods_30kb.tsv", neighbour_rows,
    [
        "old_uniprot_accession", "uniparc_id", "source_protein_accession", "nucleotide_accession",
        "target_start", "target_end", "feature_type", "feature_start", "feature_end", "feature_strand",
        "distance_to_target_bp", "protein_id", "locus_tag", "gene", "product", "note", "inference",
    ],
)
write_tsv(
    OUT / "atp23_atlas_context.tsv", atlas_rows,
    [
        "old_uniprot_accession", "uniparc_id", "query_hash", "cluster_hash", "cluster_size",
        "taxonomy_rank", "taxonomy_name", "top_phyla", "top_pfams", "pdb_available",
    ],
)
write_tsv(OUT / "failures.tsv", failures, ["target", "stage", "error"])

# Pairwise sequence identity documents whether the four HMMER hits are distinct
# homologues or overlapping fragments from one larger protein family.
pairwise_rows: list[dict[str, Any]] = []
for a, seq_a in sequences.items():
    for b, seq_b in sequences.items():
        if a >= b:
            continue
        shorter = min(len(seq_a), len(seq_b))
        # Ungapped sliding identity is deliberately reported descriptively;
        # formal homology support comes from the profile/null and tree analyses.
        best_identity = 0.0
        best_overlap = 0
        best_offset = 0
        for offset in range(-len(seq_b) + 20, len(seq_a) - 19):
            matches = overlap = 0
            for i, aa in enumerate(seq_a):
                j = i - offset
                if 0 <= j < len(seq_b):
                    overlap += 1
                    matches += aa == seq_b[j]
            if overlap >= 20 and matches / overlap > best_identity:
                best_identity = matches / overlap
                best_overlap = overlap
                best_offset = offset
        pairwise_rows.append({
            "accession_a": a, "accession_b": b, "length_a": len(seq_a), "length_b": len(seq_b),
            "best_ungapped_identity": best_identity, "overlap": best_overlap, "offset_a_minus_b": best_offset,
        })
write_tsv(
    OUT / "atp23_pairwise_fragment_identity.tsv", pairwise_rows,
    ["accession_a", "accession_b", "length_a", "length_b", "best_ungapped_identity", "overlap", "offset_a_minus_b"],
)

summary = {
    "targets": len(TARGETS),
    "sequences_recovered": len(sequences),
    "cross_reference_rows": len(crossref_rows),
    "provenance_rows": len(provenance_rows),
    "neighbourhood_features": len(neighbour_rows),
    "atlas_contexts": len(atlas_rows),
    "failures": len(failures),
}
write_json(OUT / "summary.json", summary)
print(json.dumps(summary, indent=2))
if len(sequences) != len(TARGETS):
    raise RuntimeError(f"Recovered {len(sequences)} of {len(TARGETS)} inactive ATP23 candidates")
