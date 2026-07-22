from __future__ import annotations

import csv
import json
import math
import re
import time
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any

import requests
from Bio import SeqIO
from Bio.Seq import Seq

OUT = Path("proto_mito_results/uqcrq_prokaryotic_provenance")
RAW = OUT / "raw"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

TARGETS = [
    {
        "atlas_member_hash": "5d386212018f7b73cc57625d94617b28",
        "uniparc": "UPI003B7D02FA",
        "protein_accession": "MFX6572341",
        "reported_taxon": "Acinetobacter baumannii",
        "reported_taxid": "470",
        "sequence": "MGKQPVRMKAVVYALSPFQQKVMPGLWKDITSKVSHKITDNWISTTLLLAPLVGTYSYVQWYLEKEKMEHRY",
    },
    {
        "atlas_member_hash": "6362397839c833e0fc8324168c1a4af9",
        "uniparc": "UPI00128F8A4A",
        "protein_accession": "MQL41719",
        "secondary_protein_accession": "WP_152932624",
        "reported_taxon": "Escherichia coli",
        "reported_taxid": "562",
        "sequence": "EMAKATVPVKSVIYALSPFQQKIMSGLWKDLPSKLHHKVSENWISATLLLAPLVGTYAYVQNYLEKEKLAHRY",
    },
    {
        "atlas_member_hash": "33ad6577b2299848ce3489120dc9662b",
        "uniparc": "UPI0010AF5458",
        "protein_accession": "RYH14666",
        "reported_taxon": "archaeon",
        "reported_taxid": "1906665",
        "sequence": "MDGQVSQHLSPFEQKIVGPLFKDVPLKVFKRVKDFVSEAGLGLGLGIIVFYWGDAKHKELAFHHRA",
    },
]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA = "https://www.ebi.ac.uk/ena/browser/api"
SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "proto-mito-provenance-audit/1.0 (public reproducible scientific workflow; no API key)",
        "Accept": "*/*",
    }
)


def request_text(url: str, params: dict[str, Any] | None = None, attempts: int = 10) -> str:
    delay = 1.0
    last = ""
    for _ in range(attempts):
        try:
            response = SESSION.get(url, params=params, timeout=300)
            last = f"{response.status_code} {response.text[:500]}"
        except Exception as exc:
            last = repr(exc)
            time.sleep(delay)
            delay = min(delay * 1.8, 45)
            continue
        if response.status_code == 200:
            time.sleep(0.36)  # NCBI unauthenticated E-utilities rate limit
            return response.text
        if response.status_code in {404, 400}:
            raise FileNotFoundError(last)
        if response.status_code not in {429, 500, 502, 503, 504}:
            raise RuntimeError(last)
        time.sleep(delay)
        delay = min(delay * 1.8, 45)
    raise RuntimeError(f"GET failed after retries: {url}; last={last}")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str), encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def efetch(db: str, accession: str, rettype: str = "gb", retmode: str = "text") -> str:
    return request_text(
        f"{EUTILS}/efetch.fcgi",
        {"db": db, "id": accession, "rettype": rettype, "retmode": retmode},
    )


def esearch(db: str, term: str) -> list[str]:
    text = request_text(
        f"{EUTILS}/esearch.fcgi",
        {"db": db, "term": term, "retmode": "json", "retmax": 50},
    )
    return list((json.loads(text).get("esearchresult") or {}).get("idlist") or [])


def elink(dbfrom: str, db: str, ids: list[str]) -> list[str]:
    if not ids:
        return []
    text = request_text(
        f"{EUTILS}/elink.fcgi",
        {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "json"},
    )
    payload = json.loads(text)
    linked: list[str] = []
    for linkset in payload.get("linksets") or []:
        for dbset in linkset.get("linksetdbs") or []:
            linked.extend(str(x) for x in dbset.get("links") or [])
    return sorted(set(linked))


def esummary(db: str, ids: list[str]) -> dict[str, Any]:
    if not ids:
        return {}
    text = request_text(
        f"{EUTILS}/esummary.fcgi",
        {"db": db, "id": ",".join(ids), "retmode": "json", "version": "2.0"},
    )
    return json.loads(text)


def parse_record(text: str, fmt: str = "genbank"):
    from io import StringIO

    records = list(SeqIO.parse(StringIO(text), fmt))
    return records[0] if records else None


def coded_by_values(record) -> list[str]:
    values: list[str] = []
    if record is None:
        return values
    for feature in record.features:
        for value in feature.qualifiers.get("coded_by", []):
            values.append(str(value))
    return values


def nucleotide_accessions(coded_by: list[str]) -> list[str]:
    values: set[str] = set()
    for expression in coded_by:
        for accession in re.findall(r"([A-Z]{1,6}_?[A-Z0-9]*\d+(?:\.\d+)?)\s*:", expression):
            values.add(accession)
    return sorted(values)


def coordinate_span(coded_by: list[str]) -> tuple[int | None, int | None, str]:
    all_positions: list[int] = []
    strand = "-" if any("complement" in value.lower() for value in coded_by) else "+"
    for expression in coded_by:
        for start, end in re.findall(r"(\d+)\.\.(\d+)", expression):
            all_positions.extend([int(start), int(end)])
    if not all_positions:
        return None, None, strand
    return min(all_positions), max(all_positions), strand


def seq_gc(sequence: str) -> float | None:
    sequence = re.sub(r"[^ACGT]", "", sequence.upper())
    if not sequence:
        return None
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def feature_bounds(feature) -> tuple[int, int, str]:
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "?"
    return start, end, strand


def qualifier(feature, name: str) -> str:
    values = feature.qualifiers.get(name) or []
    return ";".join(str(value) for value in values)


def location_distance(start: int, end: int, target_start: int, target_end: int) -> int:
    if end < target_start:
        return target_start - end
    if start > target_end:
        return start - target_end
    return 0


def source_qualifiers(record) -> dict[str, str]:
    if record is None:
        return {}
    for feature in record.features:
        if feature.type == "source":
            return {key: ";".join(map(str, values)) for key, values in feature.qualifiers.items()}
    return {}


def parse_ena_text(accession: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for endpoint, suffix in [("text", "txt"), ("fasta", "fasta")]:
        try:
            text = request_text(f"{ENA}/{endpoint}/{accession}", {"download": "false"})
            write_text(RAW / f"ena_{accession}.{suffix}", text)
            result[f"ena_{endpoint}_bytes"] = len(text.encode())
            result[f"ena_{endpoint}_available"] = True
        except Exception as exc:
            result[f"ena_{endpoint}_available"] = False
            result[f"ena_{endpoint}_error"] = repr(exc)
    return result


summary_rows: list[dict[str, Any]] = []
neighbour_rows: list[dict[str, Any]] = []
link_payloads: dict[str, Any] = {}
failures: list[dict[str, str]] = []

for target in TARGETS:
    primary = target["protein_accession"]
    accessions = [primary]
    if target.get("secondary_protein_accession"):
        accessions.append(str(target["secondary_protein_accession"]))

    protein_record = None
    selected_accession = ""
    protein_text = ""
    accession_status: list[dict[str, Any]] = []
    for accession in accessions:
        try:
            text = efetch("protein", accession, "gb", "text")
            write_text(RAW / f"protein_{accession}.gb", text)
            record = parse_record(text)
            accession_status.append(
                {
                    "accession": accession,
                    "retrieved": True,
                    "record_id": record.id if record else "",
                    "description": record.description if record else "",
                    "annotations": record.annotations if record else {},
                }
            )
            if protein_record is None and record is not None:
                protein_record = record
                selected_accession = accession
                protein_text = text
        except Exception as exc:
            accession_status.append({"accession": accession, "retrieved": False, "error": repr(exc)})

    write_json(RAW / f"protein_accession_status_{primary}.json", accession_status)
    if protein_record is None:
        failures.append({"target": primary, "stage": "protein_efetch", "error": json.dumps(accession_status)})
        summary_rows.append({**target, "selected_accession": "", "record_available": False})
        continue

    coded = coded_by_values(protein_record)
    nt_accessions = nucleotide_accessions(coded)
    target_start, target_end, target_strand = coordinate_span(coded)
    translated_sequence = str(protein_record.seq).replace("*", "")
    sequence_identity = (
        sum(a == b for a, b in zip(translated_sequence, target["sequence"])) / max(len(translated_sequence), len(target["sequence"]))
        if translated_sequence and target["sequence"] else None
    )

    # Resolve Entrez links by numeric protein UID.
    protein_uids: list[str] = []
    for accession in accessions:
        try:
            protein_uids.extend(esearch("protein", f"{accession}[Accession]"))
        except Exception as exc:
            failures.append({"target": primary, "stage": f"esearch_{accession}", "error": repr(exc)})
    protein_uids = sorted(set(protein_uids))
    linked: dict[str, Any] = {"protein_uids": protein_uids}
    for db in ["nuccore", "assembly", "bioproject", "biosample"]:
        try:
            ids = elink("protein", db, protein_uids)
            linked[db] = {"ids": ids, "summary": esummary(db, ids)}
        except Exception as exc:
            linked[db] = {"ids": [], "error": repr(exc)}
    link_payloads[primary] = linked
    write_json(RAW / f"entrez_links_{primary}.json", linked)

    # Include linked nucleotide accessions if coded_by was absent or incomplete.
    linked_nuc_ids = list((linked.get("nuccore") or {}).get("ids") or [])
    linked_nuc_summary = (linked.get("nuccore") or {}).get("summary") or {}
    result = linked_nuc_summary.get("result") if isinstance(linked_nuc_summary, dict) else None
    if isinstance(result, dict):
        for uid in linked_nuc_ids:
            item = result.get(str(uid))
            if isinstance(item, dict):
                accession = str(item.get("accessionversion") or item.get("caption") or "")
                if accession:
                    nt_accessions.append(accession)
    nt_accessions = list(dict.fromkeys(nt_accessions))

    contig_summaries: list[dict[str, Any]] = []
    for nt_accession in nt_accessions[:10]:
        try:
            nt_text = efetch("nuccore", nt_accession, "gbwithparts", "text")
            write_text(RAW / f"nucleotide_{nt_accession.replace('.', '_')}.gb", nt_text)
            nt_record = parse_record(nt_text)
            if nt_record is None:
                continue
            source = source_qualifiers(nt_record)
            contig_length = len(nt_record.seq)
            contig_gc = seq_gc(str(nt_record.seq))

            # If coded_by refers to this contig, preserve exact coordinates.  If
            # not, locate the translated peptide on six-frame translations as a
            # fallback and mark this explicitly.
            this_start, this_end = target_start, target_end
            coordinate_method = "coded_by"
            if nt_accession.split(".")[0] not in " ".join(coded):
                this_start = this_end = None
            if this_start is None or this_end is None:
                coordinate_method = "feature_protein_id"
                for feature in nt_record.features:
                    if feature.type != "CDS":
                        continue
                    protein_ids = feature.qualifiers.get("protein_id") or []
                    if any(str(pid).split(".")[0] in {acc.split(".")[0] for acc in accessions} for pid in protein_ids):
                        this_start, this_end, _ = feature_bounds(feature)
                        break
            if this_start is None or this_end is None:
                coordinate_method = "unresolved"

            target_nt = ""
            target_gc = None
            target_gc3 = None
            if this_start is not None and this_end is not None:
                target_nt = str(nt_record.seq[this_start - 1:this_end])
                if target_strand == "-":
                    target_nt = str(Seq(target_nt).reverse_complement())
                target_gc = seq_gc(target_nt)
                third = target_nt[2::3]
                target_gc3 = seq_gc(third)

            feature_rows = []
            for feature in nt_record.features:
                if feature.type not in {"CDS", "rRNA", "tRNA", "ncRNA", "misc_feature", "repeat_region"}:
                    continue
                start, end, strand = feature_bounds(feature)
                distance = (
                    location_distance(start, end, int(this_start), int(this_end))
                    if this_start is not None and this_end is not None else 10**18
                )
                if distance > 30000:
                    continue
                row = {
                    "target_protein": primary,
                    "selected_protein_accession": selected_accession,
                    "nucleotide_accession": nt_record.id,
                    "contig_organism": source.get("organism", ""),
                    "target_start": this_start if this_start is not None else "",
                    "target_end": this_end if this_end is not None else "",
                    "target_strand": target_strand,
                    "feature_type": feature.type,
                    "feature_start": start,
                    "feature_end": end,
                    "feature_strand": strand,
                    "distance_to_target_bp": distance if distance < 10**18 else "",
                    "overlaps_target": bool(distance == 0),
                    "protein_id": qualifier(feature, "protein_id"),
                    "locus_tag": qualifier(feature, "locus_tag"),
                    "gene": qualifier(feature, "gene"),
                    "product": qualifier(feature, "product"),
                    "note": qualifier(feature, "note"),
                    "inference": qualifier(feature, "inference"),
                    "db_xref": qualifier(feature, "db_xref"),
                }
                feature_rows.append(row)
                neighbour_rows.append(row)

            contig_summaries.append(
                {
                    "nucleotide_accession": nt_record.id,
                    "description": nt_record.description,
                    "contig_length": contig_length,
                    "contig_gc": contig_gc,
                    "source_organism": source.get("organism", ""),
                    "source_taxon": source.get("db_xref", ""),
                    "source_isolation_source": source.get("isolation_source", ""),
                    "source_host": source.get("host", ""),
                    "source_country": source.get("country", ""),
                    "source_collection_date": source.get("collection_date", ""),
                    "source_metagenome": source.get("metagenome_source", ""),
                    "source_environmental_sample": source.get("environmental_sample", ""),
                    "source_strain": source.get("strain", ""),
                    "source_bio_material": source.get("bio_material", ""),
                    "source_culture_collection": source.get("culture_collection", ""),
                    "source_specimen_voucher": source.get("specimen_voucher", ""),
                    "coordinate_method": coordinate_method,
                    "target_start": this_start,
                    "target_end": this_end,
                    "target_strand": target_strand,
                    "target_nt_length": len(target_nt),
                    "target_gc": target_gc,
                    "target_gc3": target_gc3,
                    "distance_to_left_edge": (int(this_start) - 1) if this_start is not None else None,
                    "distance_to_right_edge": (contig_length - int(this_end)) if this_end is not None else None,
                    "nearby_feature_count_30kb": len(feature_rows),
                    "nearby_products": [row["product"] for row in feature_rows if row["product"]],
                }
            )
        except Exception as exc:
            failures.append({"target": primary, "stage": f"nuccore_{nt_accession}", "error": repr(exc)})

    ena_fields = parse_ena_text(selected_accession)
    for nt_accession in nt_accessions[:5]:
        ena_fields.update({f"{nt_accession}_{key}": value for key, value in parse_ena_text(nt_accession).items()})

    write_json(RAW / f"contig_summaries_{primary}.json", contig_summaries)
    summary_rows.append(
        {
            **target,
            "selected_accession": selected_accession,
            "record_available": True,
            "protein_record_id": protein_record.id,
            "protein_description": protein_record.description,
            "protein_length": len(protein_record.seq),
            "sequence_identity_to_atlas": sequence_identity,
            "coded_by": ";".join(coded),
            "nucleotide_accessions": ";".join(nt_accessions),
            "target_start": target_start,
            "target_end": target_end,
            "target_strand": target_strand,
            "protein_annotations": json.dumps(protein_record.annotations, sort_keys=True, default=str),
            "protein_uids": ";".join(protein_uids),
            "linked_assembly_ids": ";".join((linked.get("assembly") or {}).get("ids") or []),
            "linked_bioproject_ids": ";".join((linked.get("bioproject") or {}).get("ids") or []),
            "linked_biosample_ids": ";".join((linked.get("biosample") or {}).get("ids") or []),
            "contig_count_retrieved": len(contig_summaries),
            "contig_summaries": json.dumps(contig_summaries, sort_keys=True, default=str),
            "ena_status": json.dumps(ena_fields, sort_keys=True),
        }
    )

summary_fields = sorted({key for row in summary_rows for key in row})
write_tsv(OUT / "protein_and_contig_summary.tsv", summary_rows, summary_fields)
write_tsv(
    OUT / "genomic_neighbourhoods_30kb.tsv",
    neighbour_rows,
    [
        "target_protein", "selected_protein_accession", "nucleotide_accession",
        "contig_organism", "target_start", "target_end", "target_strand",
        "feature_type", "feature_start", "feature_end", "feature_strand",
        "distance_to_target_bp", "overlaps_target", "protein_id", "locus_tag",
        "gene", "product", "note", "inference", "db_xref",
    ],
)
write_tsv(OUT / "failures.tsv", failures, ["target", "stage", "error"])
write_json(OUT / "entrez_link_payloads.json", link_payloads)
run_summary = {
    "targets": len(TARGETS),
    "protein_records_retrieved": sum(bool(row.get("record_available")) for row in summary_rows),
    "contig_records_retrieved": sum(int(row.get("contig_count_retrieved") or 0) for row in summary_rows),
    "neighbourhood_features": len(neighbour_rows),
    "failures": len(failures),
}
write_json(OUT / "summary.json", run_summary)
print(json.dumps(run_summary, indent=2))
if sum(bool(row.get("record_available")) for row in summary_rows) < 2:
    raise RuntimeError("Too few source protein records recovered")
