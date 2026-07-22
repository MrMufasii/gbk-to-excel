#!/usr/bin/env python3
from pathlib import Path

p = Path("scripts/aas_module_strengthening.py")
text = p.read_text(encoding="utf-8")

old = '''def fetch_reference_taxon(query_name: str, output: Path) -> tuple[Path, pd.DataFrame]:
    query = f'taxonomy_name:"{query_name}" AND keyword:KW-1185'
    fields = "accession,id,organism_name,organism_id,xref_proteomes,length,sequence,lineage"
    url = "https://rest.uniprot.org/uniprotkb/stream"
    response = get(url, params={"query": query, "format": "tsv", "fields": fields, "compressed": "true"}, timeout=1200)
    content = response.content
    if content[:2] == b"\\x1f\\x8b":
        text = gzip.decompress(content).decode("utf-8")
    else:
        text = content.decode("utf-8")
    frame = pd.read_csv(io.StringIO(text), sep="\\t")
    frame = frame.rename(columns={
        "Entry": "accession", "Entry Name": "entry_name", "Organism": "organism",
        "Organism (ID)": "taxid", "Proteomes": "proteomes", "Length": "length",
        "Sequence": "sequence", "Taxonomic lineage": "lineage",
    })
    frame["sequence"] = frame["sequence"].map(clean_sequence)
    frame = frame[frame["sequence"].str.len() >= 40].drop_duplicates("accession")
    tsv = output / f"{query_name.replace(' ', '_')}_reference_proteins.tsv.gz"
    frame.to_csv(tsv, sep="\\t", index=False, compression="gzip")
    fasta = output / f"{query_name.replace(' ', '_')}_reference_proteins.faa"
    write_fasta(
        fasta,
        [
            (
                f"{row.accession}|{str(row.proteomes).split(';')[0]}|{row.taxid}",
                row.sequence,
            )
            for row in frame.itertuples()
        ],
    )
    return fasta, frame
'''
new = '''def fetch_reference_taxon(query_name: str, output: Path) -> tuple[Path, pd.DataFrame]:
    # Restrict to curated Swiss-Prot entries. This keeps the survey tractable,
    # retains the experimentally characterized Chlamydia split system, and
    # avoids fragile very-large compressed streams.
    query = f'taxonomy_name:"{query_name}" AND reviewed:true'
    frame = fetch_uniprot(query, 100000)
    if frame.empty:
        raise RuntimeError(f"No reviewed proteins returned for {query_name}")
    frame = frame.drop_duplicates("accession").copy()
    tsv = output / f"{query_name.replace(' ', '_')}_reviewed_proteins.tsv.gz"
    frame.to_csv(tsv, sep="\\t", index=False, compression="gzip")
    fasta = output / f"{query_name.replace(' ', '_')}_reviewed_proteins.faa"
    write_fasta(
        fasta,
        [
            (
                f"{row.accession}|{str(row.proteomes).split(';')[0]}|{row.taxid}",
                row.sequence,
            )
            for row in frame.itertuples()
        ],
    )
    return fasta, frame
'''
if old not in text:
    raise SystemExit("fetch_reference_taxon block not found")
text = text.replace(old, new)

old = '''    representatives = {
        "Andalucia_tafazzin": sequences["tafazzin"],
        "Andalucia_fatty_synthetase": sequences["fatty_synthetase"],
    }
    for label, frame in panels_n.items():
        row = frame.iloc[0]
        representatives[f"N_{label}_{row.accession}"] = row.sequence
    for label, frame in panels_c.items():
        row = frame.iloc[0]
        representatives[f"C_{label}_{row.accession}"] = row.sequence
'''
new = '''    representatives = {
        "Andalucia_tafazzin_domain": motif_window(sequences["tafazzin"]),
        "Andalucia_fatty_synthetase": sequences["fatty_synthetase"],
    }
    for label, frame in panels_n.items():
        row = frame.iloc[0]
        representatives[f"N_{label}_{row.accession}"] = row.domain_sequence
    for label, frame in panels_c.items():
        row = frame.iloc[0]
        representatives[f"C_{label}_{row.accession}"] = row.domain_sequence
'''
if old not in text:
    raise SystemExit("structure representative block not found")
text = text.replace(old, new)
text = text.replace(
    'best("Andalucia_tafazzin", ["Aas_N_domain", "N_LpaT_", "N_PlsC_", "N_PlsB_"])',
    'best("Andalucia_tafazzin_domain", ["N_AasN_", "N_LpaT_", "N_PlsC_", "N_PlsB_"])',
)
text = text.replace(
    'for query_pattern in ["Andalucia_tafazzin", "Andalucia_fatty_synthetase"]:',
    'for query_pattern in ["Andalucia_tafazzin_domain", "Andalucia_fatty_synthetase"]:',
)
start = text.find("    # Crop the fused Aas representative to the experimentally established module boundaries.")
end = text.find('    foldseek = shutil.which("foldseek")', start)
if start != -1 and end != -1:
    text = text[:start] + text[end:]

p.write_text(text, encoding="utf-8")
