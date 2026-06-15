#!/usr/bin/env python3
"""
gbk_to_excel.py - Extract GenBank (.gbk/.gbff) annotations into a searchable Excel workbook.

Designed for Bakta/Prokka-style bacterial annotation output, but works on any GenBank file.

For every annotated feature (CDS, tRNA, rRNA, ncRNA, tmRNA, regulatory, repeat_region,
rep_origin, oriT, ...) it pulls out the gene name, product description, coordinates,
every cross-reference / qualifier of interest, and the gene's nucleotide sequence, then
writes one row per gene to an Excel "Annotations" table.

The workbook has several sheets:
  * Lookup       - a single search box; type any text (gene, product keyword, EC/KEGG/COG
                   id, locus tag, sample ...) and matching genes appear below. Built with
                   AGGREGATE/INDEX so it works in any Excel (2016+) and LibreOffice - no
                   dynamic-array support required.
  * Annotations  - the master Excel Table with filter dropdowns on every column.
  * AMR_overview - a Sample x resistance-class matrix (how many AMR genes per isolate).
  * AMR_genes    - every gene flagged as an antimicrobial-resistance determinant.
  * Summary      - per-sample feature counts.
  * README       - in-workbook help.

Usage:
    python gbk_to_excel.py [INPUTS ...] [-o OUTPUT] [options]

INPUTS may be GenBank files, directories (searched recursively for *.gbk/*.gbff/*.gb),
or glob patterns. If no input is given the script searches the current directory and its
parent for GenBank files.

By default it writes ONE workbook per genome (named <sample>_annotations.xlsx). With -o
as a folder the files go there; otherwise they go next to the input. Use --combined to
instead merge every genome into a single workbook.

Options:
    -o, --output PATH     Output folder for per-genome files, or (with --combined) the
                          single .xlsx path. Default: next to the input.
    --combined            Merge all genomes into ONE workbook instead of one per genome.
    --no-sequence         Omit the per-gene nucleotide sequence column (smaller files).
    --translation         Also include the full amino-acid translation column.
    --no-genes            Drop orphan 'gene' features (kept by default).
    -h, --help            Show this help
"""

import argparse
import glob
import hashlib
import os
import re
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

GBK_EXTENSIONS = (".gbk", ".gbff", ".gb", ".genbank", ".gbk.txt")

# db_xref prefixes that get their own column; everything else -> "Other_Xref"
XREF_COLUMNS = ["EC", "COG", "KEGG", "GO", "UniRef", "RefSeq", "RFAM", "SO"]
# Sources whose accession is canonically written WITH the prefix (e.g. GO:0003677).
XREF_KEEP_PREFIX = {"GO"}

# Feature types we never emit as their own row.
SKIP_TYPES = {"source"}

# Final column order. Translation (if requested) is appended after this list so the
# Lookup sheet can return the whole "Sample:Protein_ID" block without the bulky sequence.
BASE_COLUMNS = [
    "Sample", "File", "Contig", "Type", "Locus_Tag", "Gene", "Product",
    "Start", "End", "Strand", "Length_bp", "Length_aa",
    "EC", "COG", "KEGG", "GO", "UniRef", "RefSeq", "RFAM", "SO",
    "ncRNA_class", "Pseudo", "Note", "Inference", "Other_Xref", "Protein_ID",
]

EXCEL_CELL_LIMIT = 32760  # Excel hard limit is 32767 chars per cell.

# Hidden helper column: a single concatenated string per row that the Lookup sheet
# searches against. Keeps the search formula fast and simple.
SEARCH_KEY_NAME = "Search_Key"
SEARCH_KEY_FIELDS = ["Sample", "Type", "Locus_Tag", "Gene", "Product",
                     "EC", "COG", "KEGG", "RefSeq", "Note"]

# Columns shown (in this order) in the Lookup search results - gene/product first so
# the most useful info is left-most.
LOOKUP_COLUMNS = ["Gene", "Product", "Type", "Locus_Tag", "Sample", "Contig",
                  "Start", "End", "Strand", "Length_bp", "Length_aa",
                  "EC", "COG", "KEGG", "GO", "UniRef", "RefSeq", "Pseudo", "Note"]
LOOKUP_MAX_ROWS = 400  # how many matching rows the Lookup sheet can display at once

# ---------------------------------------------------------------------------
# Antimicrobial-resistance (AMR) classification
# ---------------------------------------------------------------------------
# Each rule: (class label, gene-symbol regex, [product keywords]).
# A feature is flagged if its gene symbol matches the regex OR its product contains a
# keyword. Bakta writes genuine *acquired* resistance genes with parentheses or a
# number suffix - tet(B), mph(A), aac(6'), sul1, dfrA1 - while look-alike housekeeping
# genes do not - tetC, msrA, aphA, sulA. The patterns below exploit that to stay precise.
AMR_RULES = [
    ("Beta-lactam", re.compile(r"^bla\w", re.I),
     ["beta-lactamase", "carbapenem-hydrolyzing", "cephalosporinase",
      "metallo-beta-lactamase", "extended-spectrum beta"]),
    ("Aminoglycoside",
     re.compile(r"^(aac\(|aad[A-Z]|aph\(|ant\(|arm[A-Z]|rmt[A-Z]|npm[A-Z]|str[AB]\b|sat[\d-]|spc[A-Z\d])", re.I),
     ["aminoglycoside", "streptomycin", "kanamycin", "gentamicin", "spectinomycin",
      "streptothricin"]),
    ("Fluoroquinolone", re.compile(r"^(qnr|qep[A-Z]|crpp|smqnr|oqx)", re.I),
     ["fluoroquinolone-acetylating", "quinolone resistance"]),
    ("Tetracycline", re.compile(r"^tet\(", re.I), ["tetracycline resistance",
     "tetracycline efflux", "tetracycline-resistant"]),
    ("MLS (macrolide/lincosamide/streptogramin)",
     re.compile(r"^(mph\(|erm\(|ere\(|msr\(|mef\(|lnu\(|vga\(|vat\(|lsa\(|vgb\()", re.I),
     ["macrolide", "lincosamide", "streptogramin", "erythromycin resistance"]),
    ("Sulfonamide", re.compile(r"^sul\d", re.I), ["sulfonamide"]),
    ("Trimethoprim", re.compile(r"^dfr[A-Z]", re.I), ["trimethoprim"]),
    ("Phenicol", re.compile(r"^(cat[A-Z\d]|cml|flo[A-Z]|cmx|fex[A-Z]|cfr)", re.I),
     ["chloramphenicol", "florfenicol", "phenicol"]),
    ("Colistin/Polymyxin", re.compile(r"^mcr", re.I), ["colistin", "polymyxin"]),
    ("Fosfomycin", re.compile(r"^fos[A-Z]", re.I), ["fosfomycin"]),
    ("Rifamycin", re.compile(r"^arr", re.I), ["rifampin", "rifamycin"]),
    ("Glycopeptide (vancomycin)", re.compile(r"^van[A-Z]", re.I),
     ["vancomycin resistance", "glycopeptide resistance"]),
    ("Bleomycin", re.compile(r"^ble[A-Z]", re.I), ["bleomycin"]),
]
# AMR class display order (matrix columns)
AMR_CLASS_ORDER = [r[0] for r in AMR_RULES]

# If the product contains any of these, it is a regulator/sensor/biosynthetic gene, not a
# resistance determinant - skip it (kills false positives like "beta-lactamase regulator
# AmpE", "erythromycin resistance repressor", "cell division inhibitor SulA").
AMR_EXCLUDE_PRODUCT = ["regulator", "repressor", "sensor", "activator", "inhibitor",
                       "biosynthesis", "carbamoyl", "sulfoxide reductase", "two-component"]

# Core/intrinsic genes that slip through product matching but are not acquired resistance.
AMR_GENE_BLACKLIST = {"maca", "macb", "mdfa", "acra", "acrb", "acrd", "emra", "emrb"}


def amr_classify(gene, product):
    """Return (class, matched_on) if the gene looks like an acquired resistance gene."""
    g = (gene or "").strip()
    p = (product or "").lower()
    if g and g.lower() in AMR_GENE_BLACKLIST:
        return None, None
    if any(x in p for x in AMR_EXCLUDE_PRODUCT):
        return None, None
    for cls, gene_re, keywords in AMR_RULES:
        if g and gene_re.search(g):
            return cls, f"gene:{g}"
        for kw in keywords:
            if kw in p:
                return cls, f"product:{kw}"
    return None, None


# ---------------------------------------------------------------------------
# GenBank parsing
# ---------------------------------------------------------------------------

# A feature definition line: 5 leading spaces, the feature key, then the location.
_FEATURE_RE = re.compile(r"^ {5}(\S+)\s+(.*)$")
# A qualifier / continuation line: 6+ leading spaces.
_QUALIFIER_RE = re.compile(r"^ {6,}(.*)$")


def parse_location(loc):
    """Return (start, end, strand, joined) from a GenBank location string.

    Handles complement(), join()/order(), and partial markers (< >).
    start/end are 1-based; strand is '+' or '-'.
    """
    strand = "-" if "complement" in loc else "+"
    joined = "join" in loc or "order" in loc
    nums = [int(n) for n in re.findall(r"\d+", loc)]
    if nums:
        return min(nums), max(nums), strand, joined
    return None, None, strand, joined


_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                            "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
_SEGMENT_RE = re.compile(r"[<>]?(\d+)\.\.[<>]?(\d+)")


def revcomp(seq):
    """Reverse-complement a nucleotide string (IUPAC aware)."""
    return seq.translate(_COMPLEMENT)[::-1]


def extract_feature_seq(full_seq, location):
    """Slice a feature's nucleotide sequence out of its contig sequence.

    Coordinates are 1-based inclusive. Handles complement() and join()/order()
    (segments are concatenated in the listed order, then reverse-complemented if the
    whole location is on the minus strand). Returns "" if no sequence is available.
    """
    if not full_seq:
        return ""
    segs = _SEGMENT_RE.findall(location)
    if segs:
        sub = "".join(full_seq[int(a) - 1:int(b)] for a, b in segs if int(a) <= int(b))
    else:
        m = re.search(r"\d+", location)
        if not m:
            return ""
        pos = int(m.group())
        sub = full_seq[pos - 1:pos]
    if "complement" in location:
        sub = revcomp(sub)
    return sub


class Feature:
    __slots__ = ("type", "location", "quals")

    def __init__(self, ftype, location):
        self.type = ftype
        self.location = location
        self.quals = defaultdict(list)  # qualifier name -> list of values


def _clean(value):
    return value.strip()


def iter_records(path):
    """Yield (contig_name, [Feature, ...], sequence) for each LOCUS in a GenBank file.

    The FEATURES block is fully parsed; the ORIGIN block is read into a plain
    nucleotide string (digits/whitespace stripped) so per-gene sequences can be sliced.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        contig = None
        in_features = False
        in_origin = False
        features = []
        seq_parts = []
        cur = None
        cur_key = None
        cur_parts = []
        in_quote = False

        def commit_qualifier():
            nonlocal cur_key, cur_parts, in_quote
            if cur is not None and cur_key is not None:
                sep = "" if cur_key == "translation" else " "
                cur.quals[cur_key].append(_clean(sep.join(cur_parts)))
            cur_key, cur_parts, in_quote = None, [], False

        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")

            if line.startswith("LOCUS"):
                # flush any in-progress record from the previous LOCUS
                commit_qualifier()
                if contig is not None and features:
                    yield contig, features, "".join(seq_parts)
                parts = line.split()
                contig = parts[1] if len(parts) > 1 else "unknown"
                in_features = in_origin = False
                features, seq_parts, cur = [], [], None
                continue

            if line.startswith("FEATURES"):
                in_features = True
                continue

            if in_origin:
                if line.startswith("//"):
                    yield contig, features, "".join(seq_parts)
                    contig, features, seq_parts, cur = None, [], [], None
                    in_origin = False
                else:  # a sequence line: keep only the bases
                    seq_parts.append(re.sub(r"[^A-Za-z]", "", line))
                continue

            if not in_features:
                # handle "//" with no ORIGIN block, and ignore header lines otherwise
                if line.startswith("//") and contig is not None and features:
                    commit_qualifier()
                    yield contig, features, "".join(seq_parts)
                    contig, features, seq_parts, cur = None, [], [], None
                continue

            # --- inside the FEATURES block ---
            if line.startswith("ORIGIN"):
                commit_qualifier()
                in_features, in_origin, seq_parts = False, True, []
                continue
            if line.startswith("//"):
                commit_qualifier()
                if contig is not None and features:
                    yield contig, features, "".join(seq_parts)
                    contig, features, seq_parts, cur = None, [], [], None
                in_features = False
                continue
            if line and not line[0].isspace():  # some other top-level keyword
                commit_qualifier()
                in_features = False
                continue

            m = _FEATURE_RE.match(line)
            if m:
                # New feature definition.
                commit_qualifier()
                ftype, location = m.group(1), m.group(2).strip()
                cur = Feature(ftype, location)
                features.append(cur)
                continue

            q = _QUALIFIER_RE.match(line)
            if not q or cur is None:
                continue
            text = q.group(1).strip()

            if in_quote:
                # Continuation of a multi-line quoted value.
                if text.endswith('"'):
                    cur_parts.append(text[:-1])
                    in_quote = False
                else:
                    cur_parts.append(text)
                continue

            if text.startswith("/"):
                commit_qualifier()
                body = text[1:]
                if "=" in body:
                    key, val = body.split("=", 1)
                    cur_key = key
                    if val.startswith('"'):
                        val = val[1:]
                        if val.endswith('"'):  # opens and closes on one line
                            cur_parts = [val[:-1]]
                            in_quote = False
                        else:
                            cur_parts = [val]
                            in_quote = True
                    else:  # unquoted value, e.g. /rpt_type=direct
                        cur_parts = [val]
                        in_quote = False
                else:  # valueless flag, e.g. /pseudo
                    cur_key = body
                    cur_parts = [""]
                    in_quote = False
            else:
                # Continuation of an unquoted multi-line location or value.
                if cur_key is not None:
                    cur_parts.append(text)

        # End of file.
        commit_qualifier()
        if contig is not None and features:
            yield contig, features, "".join(seq_parts)


# ---------------------------------------------------------------------------
# Feature -> row conversion
# ---------------------------------------------------------------------------

def first(quals, key):
    vals = quals.get(key)
    return vals[0] if vals else ""


def bucket_xrefs(quals):
    """Split db_xref values into per-source buckets plus an 'Other' catch-all."""
    buckets = {col: [] for col in XREF_COLUMNS}
    other = []
    for ref in quals.get("db_xref", []):
        prefix, _, _ = ref.partition(":")
        if prefix in buckets:
            if prefix in XREF_KEEP_PREFIX:
                buckets[prefix].append(ref)
            else:
                buckets[prefix].append(ref.split(":", 1)[1] if ":" in ref else ref)
        else:
            other.append(ref)
    return buckets, other


def features_to_rows(sample, fname, contig, features, seq="", want_seq=False):
    """Convert one contig's features into annotation rows.

    Bakta emits a bare 'gene' feature plus a partner CDS/RNA feature sharing the same
    locus_tag. We emit one row per primary feature and fold the gene feature's gene name
    and /pseudogene flag onto it. Orphan gene features (pseudogenes with no CDS) are kept.
    """
    gene_info = {}            # locus_tag -> {"gene", "pseudo", "feature"}
    primary_tags = set()      # locus_tags that have a non-gene partner

    for feat in features:
        lt = first(feat.quals, "locus_tag")
        if feat.type == "gene":
            if lt:
                gene_info[lt] = {
                    "gene": first(feat.quals, "gene"),
                    "pseudo": bool(feat.quals.get("pseudo") or feat.quals.get("pseudogene")),
                    "feature": feat,
                }
        elif feat.type not in SKIP_TYPES:
            if lt:
                primary_tags.add(lt)

    rows = []

    def make_row(feat, force_gene="", force_pseudo=False):
        lt = first(feat.quals, "locus_tag")
        start, end, strand, joined = parse_location(feat.location)
        buckets, other = bucket_xrefs(feat.quals)

        gene = first(feat.quals, "gene") or force_gene
        translation = first(feat.quals, "translation")
        length_bp = (end - start + 1) if (start is not None and end is not None) else ""
        is_pseudo = (
            force_pseudo
            or bool(feat.quals.get("pseudo") or feat.quals.get("pseudogene"))
        )
        notes = list(feat.quals.get("note", []))
        if feat.quals.get("pseudogene"):
            notes.append("pseudogene=" + ";".join(feat.quals["pseudogene"]))

        row = {
            "Sample": sample,
            "File": fname,
            "Contig": contig,
            "Type": feat.type,
            "Locus_Tag": lt,
            "Gene": gene,
            "Product": first(feat.quals, "product"),
            "Start": start if start is not None else "",
            "End": end if end is not None else "",
            "Strand": strand,
            "Length_bp": length_bp,
            "Length_aa": len(translation) if translation else "",
            "ncRNA_class": first(feat.quals, "ncRNA_class")
            or first(feat.quals, "regulatory_class")
            or first(feat.quals, "rpt_family"),
            "Pseudo": "yes" if is_pseudo else "",
            "Note": " | ".join(notes),
            "Inference": "; ".join(feat.quals.get("inference", [])),
            "Other_Xref": "; ".join(other),
            "Protein_ID": first(feat.quals, "protein_id"),
            "Nt_Sequence": extract_feature_seq(seq, feat.location) if want_seq else "",
            "Translation": translation,
        }
        for col in XREF_COLUMNS:
            row[col] = "; ".join(buckets[col])
        return row

    for feat in features:
        if feat.type in SKIP_TYPES:
            continue
        if feat.type == "gene":
            continue  # handled via backfill / orphan pass below
        lt = first(feat.quals, "locus_tag")
        g = gene_info.get(lt, {})
        rows.append(make_row(feat, force_gene=g.get("gene", ""),
                              force_pseudo=g.get("pseudo", False)))

    return rows, gene_info, primary_tags


# ---------------------------------------------------------------------------
# Input discovery & dedup
# ---------------------------------------------------------------------------

def discover_inputs(inputs):
    """Expand files / dirs / globs into a sorted list of GenBank file paths."""
    found = []
    for item in inputs:
        if os.path.isdir(item):
            for root, _dirs, files in os.walk(item):
                for f in files:
                    if f.lower().endswith(GBK_EXTENSIONS):
                        found.append(os.path.join(root, f))
        elif os.path.isfile(item):
            found.append(item)
        else:
            found.extend(glob.glob(item, recursive=True))
    # unique, stable order
    seen = set()
    out = []
    for p in sorted(found):
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            out.append(p)
    return out


def content_key(path, nbytes=4 * 1024 * 1024):
    """Cheap content fingerprint: file size + sha1 of the first nbytes."""
    size = os.path.getsize(path)
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        h.update(fh.read(nbytes))
    return (size, h.hexdigest())


def dedupe(paths):
    """Drop byte-identical duplicates (e.g. the root .gbk vs annotation/.gbff).

    When two paths share content, prefer the .gbk, otherwise the shortest path.
    """
    by_key = defaultdict(list)
    for p in paths:
        by_key[content_key(p)].append(p)
    kept, dupes = [], []
    for _key, group in by_key.items():
        group.sort(key=lambda p: (0 if p.lower().endswith(".gbk") else 1, len(p), p))
        kept.append(group[0])
        dupes.extend(group[1:])
    kept.sort()
    return kept, dupes


def default_output(input_args, kept):
    """Pick a sensible output path when -o is not given.

    For a single file/dir argument the workbook is written right next to it so the
    right-click context menu drops the .xlsx in an obvious place.
    """
    if len(input_args) == 1:
        item = input_args[0]
        if os.path.isfile(item):
            stem = os.path.splitext(os.path.basename(item))[0]
            return os.path.join(os.path.dirname(os.path.abspath(item)),
                                stem + "_annotations.xlsx")
        if os.path.isdir(item):
            d = os.path.abspath(item)
            return os.path.join(d, os.path.basename(d.rstrip(os.sep)) + "_annotations.xlsx")
    return os.path.abspath("gbk_annotations.xlsx")


def safe_name(s):
    """Make a sample label safe to use as a file name."""
    return re.sub(r'[<>:"/\\|?*]+', "_", s).strip() or "genome"


def per_genome_out_dir(args):
    """Directory to drop the per-genome workbooks into."""
    if args.output:
        o = args.output
        if o.lower().endswith(".xlsx"):
            return os.path.dirname(os.path.abspath(o)) or "."
        return o  # treat as a directory
    if len(args.inputs) == 1:
        item = args.inputs[0]
        if os.path.isdir(item):
            return item
        if os.path.isfile(item):
            return os.path.dirname(os.path.abspath(item)) or "."
    return os.getcwd()


def open_path(target):
    """Open a file or folder with the OS default handler (Windows)."""
    try:
        os.startfile(os.path.abspath(target))  # noqa: S606 (Windows only)
    except Exception as exc:  # pragma: no cover
        print(f"[i] Could not auto-open {target}: {exc}")


def sample_name(path):
    """Derive a human-friendly sample label from the file's location."""
    parts = os.path.normpath(os.path.abspath(path)).split(os.sep)
    if "annotation" in parts:
        i = parts.index("annotation")
        if i >= 1:
            return parts[i - 1]
    # otherwise use the immediate parent folder name, falling back to the file stem
    parent = parts[-2] if len(parts) >= 2 else ""
    if parent and parent.lower() not in ("annotation", ""):
        return parent
    stem = os.path.splitext(os.path.basename(path))[0]
    return re.sub(r"(_\d+)?_reference$", "", stem) or stem


# ---------------------------------------------------------------------------
# Excel writing
# ---------------------------------------------------------------------------

def col_letter(idx):
    """1-based column index -> Excel column letters."""
    s = ""
    while idx:
        idx, r = divmod(idx - 1, 26)
        s = chr(65 + r) + s
    return s


def truncate(value):
    if isinstance(value, str) and len(value) > EXCEL_CELL_LIMIT:
        return value[:EXCEL_CELL_LIMIT] + " ...[truncated]"
    return value


WIDTHS = {
    "Sample": 32, "File": 28, "Contig": 11, "Type": 10, "Locus_Tag": 15,
    "Gene": 12, "Product": 50, "Start": 10, "End": 10, "Strand": 7,
    "Length_bp": 10, "Length_aa": 10, "EC": 12, "COG": 12, "KEGG": 12,
    "GO": 22, "UniRef": 26, "RefSeq": 18, "RFAM": 12, "SO": 14,
    "ncRNA_class": 16, "Pseudo": 8, "Note": 40, "Inference": 40,
    "Other_Xref": 24, "Protein_ID": 22, "Nt_Sequence": 40, "Translation": 40,
    "Search_Key": 30, "AMR_Class": 24, "Matched_On": 22,
}


def write_workbook(rows, columns, summary_rows, output, dupes):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill, Border, Side
    from openpyxl.worksheet.table import Table, TableStyleInfo
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    # Recalculate on open so the Lookup formulas evaluate immediately.
    try:
        wb.calculation.fullCalcOnLoad = True
    except Exception:
        pass

    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="2F5597")
    amr_fill = PatternFill("solid", fgColor="C0392B")
    title_font = Font(bold=True, size=14)
    search_fill = PatternFill("solid", fgColor="FFF2CC")
    search_border = Border(*[Side(style="medium", color="BF8F00")] * 4)

    def style_header(ws, row, headers, fill=header_fill):
        for j, h in enumerate(headers, start=1):
            c = ws.cell(row=row, column=j, value=h)
            c.font = header_font
            c.fill = fill

    def set_widths(ws, cols, extra=None):
        for i, c in enumerate(cols, start=1):
            ws.column_dimensions[get_column_letter(i)].width = (extra or {}).get(c, WIDTHS.get(c, 16))

    # ---- Annotations sheet (the master table) ----------------------------
    ws = wb.active
    ws.title = "Annotations"
    ws.append(columns)
    for r in rows:
        ws.append([truncate(r.get(c, "")) for c in columns])

    last_col = col_letter(len(columns))
    last_row = len(rows) + 1
    table = Table(displayName="Annotations", ref=f"A1:{last_col}{last_row}")
    table.tableStyleInfo = TableStyleInfo(
        name="TableStyleMedium2", showRowStripes=True, showColumnStripes=False)
    ws.add_table(table)
    ws.freeze_panes = "B2"
    set_widths(ws, columns)
    # Hide the internal search-key helper column.
    if SEARCH_KEY_NAME in columns:
        ws.column_dimensions[get_column_letter(columns.index(SEARCH_KEY_NAME) + 1)].hidden = True

    # ---- Lookup sheet ----------------------------------------------------
    # Pure AGGREGATE/INDEX extractor: works in Excel 2016+ and LibreOffice, with no
    # dynamic-array / implicit-intersection pitfalls.
    lk = wb.create_sheet("Lookup", 0)  # first sheet
    lk["A1"] = "Gene / annotation lookup"
    lk["A1"].font = title_font
    lk["A2"] = ("Type a gene, product keyword, EC/KEGG/COG id, locus tag or sample below, "
                "then press Enter.")
    lk["A2"].font = Font(italic=True, color="808080")
    lk["A3"] = "Search:"
    lk["A3"].font = Font(bold=True)
    lk["B3"].fill = search_fill
    lk["B3"].border = search_border
    lk["B3"].alignment = Alignment(horizontal="left")
    lk["A4"] = "Matches:"
    lk["A4"].font = Font(bold=True)
    lk["B4"].font = Font(bold=True, color="C0392B")

    disp = [c for c in LOOKUP_COLUMNS if c in columns]
    n_disp = len(disp)
    hdr_row = 6
    first_row = hdr_row + 1
    helper_col = n_disp + 2          # one blank column gap, then the hidden helper
    helper_letter = get_column_letter(helper_col)
    style_header(lk, hdr_row, disp)
    lk.cell(row=hdr_row, column=helper_col, value="_match_index")

    sk = f"Annotations[{SEARCH_KEY_NAME}]"
    sk_hdr = f"Annotations[[#Headers],[{SEARCH_KEY_NAME}]]"

    # Total match count (cheap, scans once).
    lk["B4"] = f'=IF($B$3="","",SUMPRODUCT(--ISNUMBER(SEARCH($B$3,{sk}))))'
    lk["D4"] = (f'=IF($B$3="","",IF(B4>{LOOKUP_MAX_ROWS},'
                f'"(showing first {LOOKUP_MAX_ROWS} - narrow your search to see the rest)",'
                f'""))')
    lk["D4"].font = Font(italic=True, color="808080")

    for k in range(LOOKUP_MAX_ROWS):
        r = first_row + k
        # Helper: the 1-based position (within the table body) of the k-th match.
        lk.cell(row=r, column=helper_col, value=(
            f'=IF($B$3="","",IFERROR(AGGREGATE(15,6,'
            f'(ROW({sk})-ROW({sk_hdr}))/(ISNUMBER(SEARCH($B$3,{sk}))),'
            f'ROW()-{hdr_row}),""))'))
        for j, cname in enumerate(disp, start=1):
            lk.cell(row=r, column=j, value=(
                f'=IF(${helper_letter}{r}="","",'
                f'IF(INDEX(Annotations[{cname}],${helper_letter}{r})="","",'
                f'INDEX(Annotations[{cname}],${helper_letter}{r})))'))

    lk.freeze_panes = "A7"
    set_widths(lk, disp)
    lk.column_dimensions[helper_letter].hidden = True
    lk.column_dimensions["B"].width = 30

    # ---- AMR sheets ------------------------------------------------------
    amr_rows = []
    sample_order = [s["sample"] for s in summary_rows]
    counts = {s: defaultdict(int) for s in sample_order}
    for r in rows:
        cls, matched = amr_classify(r.get("Gene", ""), r.get("Product", ""))
        if not cls:
            continue
        if r.get("Sample") in counts:
            counts[r["Sample"]][cls] += 1
        amr_rows.append({**r, "AMR_Class": cls, "Matched_On": matched})

    # AMR_overview: Sample x class matrix
    ov = wb.create_sheet("AMR_overview")
    ov["A1"] = "Antimicrobial-resistance genes per isolate"
    ov["A1"].font = title_font
    ov["A2"] = ("Counts of genes flagged as resistance determinants (by gene symbol or "
                "product). See the AMR_genes sheet for the full list.")
    ov["A2"].font = Font(italic=True, color="808080")
    classes_present = [c for c in AMR_CLASS_ORDER if any(counts[s].get(c) for s in sample_order)]
    ov_headers = ["Sample"] + classes_present + ["TOTAL"]
    ovh = 4
    style_header(ov, ovh, ov_headers)
    for i, s in enumerate(sample_order, start=ovh + 1):
        ov.cell(row=i, column=1, value=s)
        total = 0
        for j, c in enumerate(classes_present, start=2):
            v = counts[s].get(c, 0)
            total += v
            if v:
                ov.cell(row=i, column=j, value=v)
        tcell = ov.cell(row=i, column=len(ov_headers), value=total)
        tcell.font = Font(bold=True)
    # column totals row
    trow = ovh + 1 + len(sample_order)
    ov.cell(row=trow, column=1, value="TOTAL").font = Font(bold=True)
    for j, c in enumerate(classes_present, start=2):
        tot = sum(counts[s].get(c, 0) for s in sample_order)
        cell = ov.cell(row=trow, column=j, value=tot)
        cell.font = Font(bold=True)
    grand = ov.cell(row=trow, column=len(ov_headers),
                    value=sum(sum(counts[s].values()) for s in sample_order))
    grand.font = Font(bold=True)
    ov.freeze_panes = "B5"
    ov.column_dimensions["A"].width = 40
    for j in range(2, len(ov_headers) + 1):
        ov.column_dimensions[get_column_letter(j)].width = 14
    ov.row_dimensions[ovh].height = 90
    for j in range(2, len(ov_headers)):
        ov.cell(row=ovh, column=j).alignment = Alignment(textRotation=60, vertical="bottom")

    # AMR_genes: detailed, filterable table
    ag = wb.create_sheet("AMR_genes")
    amr_cols = ["Sample", "Contig", "AMR_Class", "Matched_On", "Gene", "Product",
                "Locus_Tag", "Start", "End", "Strand", "Length_aa",
                "EC", "KEGG", "RefSeq", "UniRef", "Other_Xref", "Note"]
    amr_cols = [c for c in amr_cols if c in (set(columns) | {"AMR_Class", "Matched_On"})]
    ag.append(amr_cols)
    amr_rows.sort(key=lambda r: (r.get("AMR_Class", ""), r.get("Sample", ""),
                                 r.get("Gene", "")))
    for r in amr_rows:
        ag.append([truncate(r.get(c, "")) for c in amr_cols])
    if amr_rows:
        atable = Table(displayName="AMR_Genes",
                       ref=f"A1:{col_letter(len(amr_cols))}{len(amr_rows) + 1}")
        atable.tableStyleInfo = TableStyleInfo(
            name="TableStyleMedium3", showRowStripes=True, showColumnStripes=False)
        ag.add_table(atable)
    ag.freeze_panes = "A2"
    set_widths(ag, amr_cols)

    # ---- Summary sheet ---------------------------------------------------
    sm = wb.create_sheet("Summary")
    sm["A1"] = "Per-sample summary"
    sm["A1"].font = title_font
    sum_headers = ["Sample", "File", "Contigs", "Total features"] + \
        sorted({k for row in summary_rows for k in row["types"]})
    hdr_row = 3
    for j, h in enumerate(sum_headers, start=1):
        cell = sm.cell(row=hdr_row, column=j, value=h)
        cell.font = header_font
        cell.fill = header_fill
    for i, row in enumerate(summary_rows, start=hdr_row + 1):
        sm.cell(row=i, column=1, value=row["sample"])
        sm.cell(row=i, column=2, value=row["file"])
        sm.cell(row=i, column=3, value=row["contigs"])
        sm.cell(row=i, column=4, value=row["total"])
        for j, h in enumerate(sum_headers[4:], start=5):
            sm.cell(row=i, column=j, value=row["types"].get(h, 0))
    sm.freeze_panes = "A4"
    sm.column_dimensions["A"].width = 34
    sm.column_dimensions["B"].width = 28

    # ---- README sheet ----------------------------------------------------
    rd = wb.create_sheet("README")
    bold = Font(bold=True)
    readme = [
        ("How to use this workbook", title_font),
        ("", None),
        ("Generated from GenBank (.gbk/.gbff) annotation files by gbk_to_excel.py.", None),
        ("", None),
        ("LOOKUP sheet - quick search (works in any Excel 2016+ and LibreOffice):", bold),
        ("   - Type a term in the yellow Search box (B3) and press Enter.", None),
        ("   - Matching genes appear below; searches gene, product, locus tag,", None),
        ("     sample, type, EC/KEGG/COG ids and notes (case-insensitive).", None),
        ("   - 'Matches' shows the total found; up to "
         f"{LOOKUP_MAX_ROWS} rows are listed - narrow the term to see more.", None),
        ("   - Examples: blaKPC, carbapenem, gyrA, 23S ribosomal, K02313, transposase.", None),
        ("", None),
        ("ANNOTATIONS sheet - the full table:", bold),
        ("   - Filter dropdown on every column (click the header arrow).", None),
        ("   - Filter/sort by Gene, Product, Sample, Type, EC, KEGG, position, length...", None),
        ("", None),
        ("AMR_OVERVIEW sheet - resistance genes per isolate:", bold),
        ("   - A Sample x antibiotic-class matrix: how many AMR genes each isolate carries.", None),
        ("AMR_GENES sheet - the detailed, filterable list of every flagged AMR gene.", bold),
        ("   - 'AMR_Class' = antibiotic class; 'Matched_On' = what triggered the flag", None),
        ("     (the gene symbol or a product keyword) so you can sanity-check each hit.", None),
        ("   - NOTE: this is a fast keyword/gene-symbol screen, NOT a validated AMR caller.", None),
        ("     For publication, confirm hits with AMRFinderPlus / ResFinder / CARD-RGI.", None),
        ("", None),
        ("COLUMNS", bold),
        ("   Sample/File/Contig - where the gene came from.", None),
        ("   Type - feature type (CDS, tRNA, rRNA, ncRNA, tmRNA, regulatory, ...).", None),
        ("   Gene / Product - gene symbol and functional description.", None),
        ("   Start/End/Strand/Length_bp/Length_aa - location and size.", None),
        ("   EC/COG/KEGG/GO/UniRef/RefSeq/RFAM/SO - database cross-references.", None),
        ("   Pseudo - 'yes' for pseudogenes.  Note - free-text notes / frameshifts.", None),
        ("   Inference - evidence for the annotation.  Other_Xref - any other db_xref.", None),
        ("   Nt_Sequence - the gene's nucleotide sequence (minus-strand genes are", None),
        ("     reverse-complemented, i.e. coding-strand 5'->3'), for variant checks.", None),
        ("     Sequences longer than ~32,760 chars are truncated (Excel cell limit).", None),
        ("   Translation - amino-acid sequence (only if built with --translation).", None),
        ("   (Search_Key is a hidden helper column used by the Lookup sheet.)", None),
        ("", None),
    ]
    for i, (text, font) in enumerate(readme, start=1):
        cell = rd.cell(row=i, column=1, value=text)
        if font:
            cell.font = font
    if dupes:
        base = len(readme) + 2
        rd.cell(row=base, column=1,
                value=f"Skipped {len(dupes)} duplicate file(s) with identical content:").font = Font(bold=True)
        for k, d in enumerate(dupes, start=1):
            rd.cell(row=base + k, column=1, value="   " + d)
    rd.column_dimensions["A"].width = 100

    wb.save(output)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Extract GenBank annotations into a searchable Excel workbook.")
    ap.add_argument("inputs", nargs="*",
                    help="GenBank files, directories, or glob patterns.")
    ap.add_argument("-o", "--output", default=None,
                    help="Output .xlsx path (default: next to the input file/folder).")
    ap.add_argument("--open", dest="open_after", action="store_true",
                    help="Open the workbook when finished (Windows).")
    ap.add_argument("--combined", action="store_true",
                    help="Write ALL genomes into a single workbook (default: one file per genome).")
    ap.add_argument("--no-sequence", dest="sequence", action="store_false",
                    help="Omit the per-gene nucleotide sequence column (smaller file).")
    ap.add_argument("--translation", action="store_true",
                    help="Also include the full amino-acid translation column.")
    ap.add_argument("--no-genes", dest="genes", action="store_false",
                    help="Drop orphan 'gene' features (kept by default).")
    ap.set_defaults(genes=True, sequence=True)
    args = ap.parse_args(argv)

    inputs = args.inputs
    if not inputs:
        # default search locations
        here = os.getcwd()
        inputs = [here, os.path.dirname(here)]
        print(f"[i] No input given; searching for GenBank files under:\n    {here}\n    {os.path.dirname(here)}")

    paths = discover_inputs(inputs)
    if not paths:
        print("[!] No GenBank files (*.gbk/*.gbff/*.gb) found.", file=sys.stderr)
        return 1

    kept, dupes = dedupe(paths)
    print(f"[i] Found {len(paths)} GenBank file(s); {len(kept)} unique after dedup.")
    if dupes:
        print(f"[i] Skipped {len(dupes)} duplicate(s) (identical content).")

    columns = list(BASE_COLUMNS)
    if args.sequence:
        columns.append("Nt_Sequence")
    if args.translation:
        columns.append("Translation")
    columns.append(SEARCH_KEY_NAME)  # hidden helper, kept last

    def sort_key(r):
        s = r.get("Start", "")
        return (r["Sample"], r["Contig"], s if isinstance(s, int) else 1 << 62)

    def finalize(rows):
        """Sort and attach the hidden Search_Key the Lookup sheet searches against."""
        rows.sort(key=sort_key)
        for r in rows:
            r[SEARCH_KEY_NAME] = " ".join(str(r[f]) for f in SEARCH_KEY_FIELDS if r.get(f))
        return rows

    genomes = []  # list of (summary_row, rows)

    for path in kept:
        sample = sample_name(path)
        fname = os.path.basename(path)
        print(f"[i] Parsing {sample}  ({fname}) ...", flush=True)

        type_counts = defaultdict(int)
        contigs = set()
        sample_rows = []

        for contig, features, seq in iter_records(path):
            contigs.add(contig)
            rows, gene_info, primary_tags = features_to_rows(
                sample, fname, contig, features, seq, args.sequence)
            sample_rows.extend(rows)

            # Orphan gene features (e.g. pseudogenes with no CDS partner).
            if args.genes:
                for feat in features:
                    if feat.type != "gene":
                        continue
                    lt = first(feat.quals, "locus_tag")
                    if lt and lt in primary_tags:
                        continue
                    # emit a row directly from the gene feature
                    start, end, strand, _ = parse_location(feat.location)
                    is_pseudo = bool(feat.quals.get("pseudo") or feat.quals.get("pseudogene"))
                    notes = list(feat.quals.get("note", []))
                    if feat.quals.get("pseudogene"):
                        notes.append("pseudogene=" + ";".join(feat.quals["pseudogene"]))
                    r = {c: "" for c in columns}
                    r.update({
                        "Sample": sample, "File": fname, "Contig": contig,
                        "Type": "gene", "Locus_Tag": lt,
                        "Gene": first(feat.quals, "gene"),
                        "Start": start if start is not None else "",
                        "End": end if end is not None else "",
                        "Strand": strand,
                        "Length_bp": (end - start + 1) if (start and end) else "",
                        "Pseudo": "yes" if is_pseudo else "",
                        "Note": " | ".join(notes),
                        "Nt_Sequence": extract_feature_seq(seq, feat.location) if args.sequence else "",
                    })
                    sample_rows.append(r)

        for r in sample_rows:
            type_counts[r["Type"]] += 1
        n_amr = sum(1 for r in sample_rows
                    if amr_classify(r.get("Gene", ""), r.get("Product", ""))[0])
        genomes.append(({
            "sample": sample, "file": fname, "contigs": len(contigs),
            "total": len(sample_rows), "types": dict(type_counts),
        }, sample_rows))
        print(f"    {len(sample_rows):,} features, {len(contigs)} contig(s), {n_amr} AMR gene(s).")

    # ---- COMBINED mode: everything in one workbook -----------------------
    if args.combined:
        all_rows, summary_rows = [], []
        for summ, rws in genomes:
            all_rows.extend(rws)
            summary_rows.append(summ)
        finalize(all_rows)
        output = args.output or default_output(args.inputs, kept)
        print(f"[i] Writing combined workbook ({len(all_rows):,} rows) -> {output} ...",
              flush=True)
        write_workbook(all_rows, columns, summary_rows, output, dupes)
        print(f"[+] Done. Wrote {output}")
        if args.open_after:
            open_path(output)
        return 0

    # ---- DEFAULT mode: one workbook per genome ---------------------------
    out_dir = per_genome_out_dir(args)
    os.makedirs(out_dir, exist_ok=True)
    written, used = [], set()
    for summ, rws in genomes:
        finalize(rws)
        stem = safe_name(summ["sample"])
        name, n = f"{stem}_annotations.xlsx", 2
        while name.lower() in used:           # guard against duplicate sample names
            name, n = f"{stem}_{n}_annotations.xlsx", n + 1
        used.add(name.lower())
        path = os.path.join(out_dir, name)
        print(f"[i] Writing {summ['sample']} ({len(rws):,} rows) ...", flush=True)
        write_workbook(rws, columns, [summ], path, [])
        written.append(path)

    print(f"[+] Done. Wrote {len(written)} workbook(s) into {os.path.abspath(out_dir)}")
    if args.open_after:
        # open the folder when there are several files; the file itself if just one
        open_path(written[0] if len(written) == 1 else out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
