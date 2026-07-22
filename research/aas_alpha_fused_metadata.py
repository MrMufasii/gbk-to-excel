#!/usr/bin/env python3
"""Resolve the 15 Alphaproteobacterial fused Aas-like candidates in UniProt."""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import requests

OUT = Path("results/aas_alpha_fused_metadata")
OUT.mkdir(parents=True, exist_ok=True)
ACCESSIONS = [
    "A0A0H3LVH7", "Q68WB7", "B2IEC0", "Q0BPJ5", "Q89SS6", "B1LZM2",
    "D7A3V3", "Q2W2Z9", "B1ZF95", "D4Z7D4", "V6F2V7", "A0A0P0IVY7",
    "A0A2D2D0R7", "A0A2D2D0S7", "A0A494TK52",
]
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "proto-mito-Aas-metadata/1.0"})


def protein_name(payload: dict) -> str:
    desc = payload.get("proteinDescription") or {}
    rec = desc.get("recommendedName") or {}
    value = (rec.get("fullName") or {}).get("value")
    if value:
        return str(value)
    for key in ["submissionNames", "alternativeNames"]:
        values = desc.get(key) or []
        if values:
            value = (values[0].get("fullName") or {}).get("value")
            if value:
                return str(value)
    return ""


def main() -> None:
    rows = []
    raw = {}
    for accession in ACCESSIONS:
        response = SESSION.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=180)
        response.raise_for_status()
        payload = response.json()
        raw[accession] = payload
        organism = payload.get("organism") or {}
        genes = payload.get("genes") or []
        gene = ""
        locus = ""
        if genes:
            gene = ((genes[0].get("geneName") or {}).get("value") or "")
            locus_values = genes[0].get("orderedLocusNames") or []
            locus = ";".join(str((item or {}).get("value") or "") for item in locus_values)
        cross = payload.get("uniProtKBCrossReferences") or []
        interpro = []
        pfam = []
        for ref in cross:
            if ref.get("database") == "InterPro":
                interpro.append(str(ref.get("id") or ""))
            if ref.get("database") == "Pfam":
                pfam.append(str(ref.get("id") or ""))
        comments = payload.get("comments") or []
        function_text = []
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                for text in comment.get("texts") or []:
                    if text.get("value"):
                        function_text.append(str(text["value"]))
        rows.append({
            "accession": accession,
            "entry_type": payload.get("entryType", ""),
            "protein_name": protein_name(payload),
            "gene": gene,
            "ordered_locus": locus,
            "organism": organism.get("scientificName", ""),
            "taxid": organism.get("taxonId", ""),
            "lineage": ";".join(organism.get("lineage") or []),
            "length": (payload.get("sequence") or {}).get("length", ""),
            "function": " | ".join(function_text),
            "interpro": ";".join(interpro),
            "pfam": ";".join(pfam),
        })
        time.sleep(0.15)
    with (OUT / "alpha_fused_Aas_like_metadata.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "uniprot_raw.json").write_text(json.dumps(raw, indent=2), encoding="utf-8")
    summary = {
        "proteins": len(rows),
        "annotated_as_Aas_or_acyl_ACP_synthetase": sum(
            any(term in row["protein_name"].lower() for term in ["aas", "acyl-[acyl carrier protein]", "acyl-acp"])
            for row in rows
        ),
        "rickettsiales_entries": sum("Rickettsiales" in row["lineage"] for row in rows),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
