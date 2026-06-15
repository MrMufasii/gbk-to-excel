# gbk-to-excel

Turn GenBank annotation files (`.gbk` / `.gbff` / `.gb` / `.genbank`) into a **searchable Excel workbook** — one row per gene, with a built-in search box, per-column filters, and an automatic antimicrobial-resistance (AMR) screen.

Built for **Bakta / Prokka**-style bacterial annotation output, but it works on any GenBank file.

---

## What you get

Each workbook contains these sheets:

| Sheet | What it's for |
|-------|---------------|
| **Lookup** | A single yellow search box — type a gene, product keyword, EC/KEGG/COG id, locus tag or sample and matching genes appear below. Works in Excel 2016+ **and** LibreOffice (no dynamic-array support needed). |
| **Annotations** | The master table — one row per feature (CDS, tRNA, rRNA, ncRNA, tmRNA, …), with a filter dropdown on every column. |
| **AMR_overview** | A *Sample × resistance-class* matrix: how many AMR genes each isolate carries. |
| **AMR_genes** | The detailed, filterable list of every flagged AMR gene, including **what triggered the flag** (gene symbol vs. product keyword) so you can sanity-check each hit. |
| **Summary** | Per-sample feature counts. |
| **README** | In-workbook help. |

Each gene row pulls out the gene name, product, coordinates, strand, length, every cross-reference of interest (EC, COG, KEGG, GO, UniRef, RefSeq, RFAM, SO), notes/inference, and — by default — the gene's nucleotide sequence (minus-strand genes reverse-complemented to coding-strand 5′→3′).

> ⚠️ **AMR caveat:** the AMR screen is a fast keyword / gene-symbol heuristic, **not** a validated resistance caller. For anything you'll publish, confirm hits with [AMRFinderPlus](https://github.com/ncbi/amr), [ResFinder](https://cge.food.dtu.dk/services/ResFinder/), or [CARD-RGI](https://card.mcmaster.ca/).

---

## Requirements

- **Python 3.7+**
- **[openpyxl](https://pypi.org/project/openpyxl/)** — the only third-party dependency:

  ```bash
  pip install openpyxl
  ```

  (The Windows installer below installs this for you automatically.)

---

## Usage

### Command line

```bash
# One workbook per genome, written next to the input
python gbk_to_excel.py path/to/genomes/

# A single file
python gbk_to_excel.py sample.gbff

# Merge every genome into ONE workbook
python gbk_to_excel.py path/to/genomes/ --combined -o all_genomes.xlsx

# Glob pattern (recursive)
python gbk_to_excel.py "results/**/*.gbff"
```

If you give no input at all, it searches the current directory and its parent for GenBank files.

Inputs may be **files, directories** (searched recursively), or **glob patterns**. Byte-identical duplicates (e.g. a root `.gbk` and a copy under `annotation/`) are detected and skipped automatically.

#### Options

| Option | Description |
|--------|-------------|
| `-o`, `--output PATH` | Output folder for per-genome files, or — with `--combined` — the single `.xlsx` path. Default: next to the input. |
| `--combined` | Merge all genomes into one workbook instead of one file per genome. |
| `--open` | Open the workbook (or folder) when finished. *(Windows)* |
| `--no-sequence` | Omit the per-gene nucleotide sequence column (smaller files). |
| `--translation` | Also include the full amino-acid translation column. |
| `--no-genes` | Drop orphan `gene` features (kept by default). |
| `-h`, `--help` | Show help. |

By default the tool writes **one workbook per genome**, named `<sample>_annotations.xlsx`. The sample label is derived from the folder layout (the directory above an `annotation/` folder, otherwise the immediate parent folder, otherwise the filename).

---

## Windows right-click integration (optional)

On Windows you can add an **"Extract annotations to Excel"** option to the right-click menu so you never touch a terminal.

**Install** (per-user, **no admin rights needed**):

```powershell
powershell -ExecutionPolicy Bypass -File .\Install-RightClick.ps1
```

…or just right-click `Install-RightClick.ps1` → **Run with PowerShell**.

This:
- detects Python on the machine and makes sure `openpyxl` is installed,
- registers a right-click verb for `.gbk` / `.gbff` / `.gb` / `.genbank` files and for folders,
- adds a **Send to → Extract GBK to Excel** shortcut (visible in the normal Windows 11 menu).

Then, on any GenBank file or a folder of them:

> Right-click → **Send to** → **Extract GBK to Excel**
> *(or: right-click → Show more options → Extract annotations to Excel)*

The workbook is written next to the input and opens automatically.

**Uninstall:**

```powershell
powershell -ExecutionPolicy Bypass -File .\Uninstall-RightClick.ps1
```

This removes the menu entries and the Send-to shortcut; your `.py` / `.cmd` files are left in place.

---

## Portability

The repo is self-contained and machine-independent. `gbk_to_excel.cmd` locates Python at runtime (via the `py` launcher, falling back to `python` on `PATH`), so you can copy / clone / zip this folder onto any Windows machine with Python 3 installed and it just works — no per-machine path editing. Recipients only need to run `Install-RightClick.ps1` once to wire up the right-click menu and install the dependency.

---

## Files

| File | Purpose |
|------|---------|
| `gbk_to_excel.py` | The converter (cross-platform). |
| `gbk_to_excel.cmd` | Portable Windows wrapper used by the right-click menu. |
| `Install-RightClick.ps1` | Adds the right-click / Send-to integration (per-user). |
| `Uninstall-RightClick.ps1` | Removes it. |

---

## License

Released under the [MIT License](LICENSE).
