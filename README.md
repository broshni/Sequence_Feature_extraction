# Sequence Feature Extraction Pipeline (SFVT)

This repository builds a standardized **Sequence Feature (SF) table** for a virus by combining
functional/structural feature annotations from **UniProt** and epitope annotations from
**IEDB**, mapping every feature onto a common reference (GenBank) sequence via pairwise
alignment, collapsing duplicate features reported by more than one source, and assigning
stable per-gene identifiers.

The pipeline is six scripts, run in order, each taking the previous script's output as input.
Every script is a standalone command-line tool (`--help` for full options) and writes its
output to a directory you specify — nothing is hardcoded to a particular virus or machine.

## Pipeline overview

```
SF_Table_Uniprot.py  ──┐
                        ├──> Edit_SFTable.py ──> mapping.py ──> merge_same_SF_coordinates.py ──> build_sf_ids_names.py
SF_Table_Epitope.py ──┘
```

1. **`SF_Table_Uniprot.py`** — given a UniProt taxonomic-search TSV export, filters to
   reviewed (Swiss-Prot) entries only, pulls per-entry feature annotations from the UniProt
   REST API, and builds a standardized SF table.
2. **`SF_Table_Epitope.py`** — given an IEDB epitope search export (`.xlsx`), reads the first
   sheet, cleans up IEDB's double-row header, drops "Discontinuous peptide on multi chain"
   entries and linear peptides with missing position data (tracked in a separate CSV), and
   builds a standardized SF table.
3. **`Edit_SFTable.py`** — merges the Uniprot and IEDB SF tables into one table. Optionally
   backfills `genbank_accession` (and `gene_ref`/`product_ref`) from a reference CSV.
4. **`mapping.py`** — for each row, fetches its source sequence (UniProt or GenBank) and its
   reference sequence (GenBank, via `genbank_accession`), aligns them with MAFFT, and maps the
   feature's coordinates onto the reference. Populates `start`/`end`/`segments` and adds
   `Source_Sequence`/`Ref_sequence`. Rows that can't be mapped (missing data, failed lookups,
   unmappable coordinates) are left blank rather than dropped, and the run continues.
5. **`merge_same_SF_coordinates.py`** — collapses rows that represent the same feature (same
   `gene` + `segments` + `sf_category`) reported by more than one source, merging the extra
   sources into the surviving row's `publication` column.
6. **`build_sf_ids_names.py`** — adds `sf_id` (e.g. `Measles_N_SF3`) and `sf_name` (e.g.
   `Measles_N_Helix_18(3)`) columns, producing the final SF table.

Steps 5 and 6 validate that their required columns are fully populated before running, and
will stop with a clear report (rather than silently producing bad output) if not — see
**Data quality checks** below.

## Requirements

- Python 3.9+
- Python packages: see `requirements.txt` (`pandas`, `requests`, `biopython`, `openpyxl`)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) — external command-line aligner, required
  by `mapping.py`. Must be available on your `PATH` (e.g. `brew install mafft` on macOS).

```
pip install -r requirements.txt
```

`mapping.py` also needs network access to the UniProt REST API and NCBI Entrez.

## Usage

Each script takes its input file(s) as positional arguments, an output directory (`-o`/
`--outdir`, created automatically if missing), and a name prefix (`-n`/`--name`) used to build
the output filename. Run any script with `--help` for the full list of options.

### 1. `SF_Table_Uniprot.py`

```
python SF_Table_Uniprot.py <uniprot_taxonomy_search>.tsv -o <outdir> -n <name>
```

- Input: a TSV export from a UniProt taxonomic search (must include a `Reviewed` column).
- Output: `<name>_Uniprot_SF.csv`, plus an intermediate `merged_features.csv`
  (override with `--merged-filename`).

### 2. `SF_Table_Epitope.py`

```
python SF_Table_Epitope.py <epitope_export>.xlsx -o <outdir> -n <name>
```

- Input: the raw `.xlsx` epitope table exported from IEDB (multi-sheet; the epitope table is
  read from the first sheet by default — override with `--sheet`).
- Output: `<name>_Epitope_SF.csv`, plus `<name>_Epitope_dropped_missing_position.csv` for any
  linear peptides that had to be excluded for missing position data.

### 3. `Edit_SFTable.py`

```
python Edit_SFTable.py <uniprot_sf>.csv <epitope_sf>.csv -o <outdir> -n <name> [--reference-file reference.csv]
```

- Inputs: the two SF tables from steps 1 and 2.
- `--reference-file` (optional): a CSV with `product`, `gene`, `protein_accession` columns,
  used to backfill `genbank_accession`/`gene_ref`/`product_ref` by matching product names.
  Skipped if omitted.
- Output: `<name>_SF_merged.csv`.

### 4. `mapping.py`

```
python mapping.py <merged_sf>.csv -o <outdir> -n <name> [--email you@example.com]
```

- Input: the merged SF table from step 3.
- Output: `<name>_SF_mapped.csv`. Scratch alignment files (per-row FASTA inputs/outputs used by
  MAFFT) are written to `<outdir>/_mapping_scratch/` rather than cluttering your working
  directory.
- Useful options: `--email` (required by NCBI Entrez), `--retries`/`--retry-delay`/`--timeout`
  (network resilience for UniProt/GenBank lookups), `--row-delay` (pause between rows to avoid
  rate limits), `--encoding` (input CSV encoding, default `latin-1`).

### 5. `merge_same_SF_coordinates.py`

```
python merge_same_SF_coordinates.py <mapped_sf>.csv -o <outdir> -n <name>
```

- Input: the mapped SF table from step 4 (**must** have `segments` populated — run `mapping.py`
  first).
- Output: `<name>_dedup.csv`.

### 6. `build_sf_ids_names.py`

```
python build_sf_ids_names.py <dedup_sf>.csv -o <outdir> -n <name>
```

- Input: the deduplicated SF table from step 5.
- Output: `<name>_SF_Final.csv` — the finished SF table.

## Data quality checks

`merge_same_SF_coordinates.py` and `build_sf_ids_names.py` both check, before doing any work,
that the columns they depend on (`gene`, `segments`, `sf_category`, `Source`, `source_id` for
the former; `virus_type`, `gene`, `sf_type`, `start`, `length` for the latter) have no empty
cells. If any are found, the script prints exactly which column(s) and which row indices are
affected and exits without writing an output file, instead of silently producing incorrect or
non-unique results (e.g. rows with a missing `gene` would otherwise all collide on the same
`sf_id`).

`mapping.py` takes the opposite approach for per-row failures: a row that can't be mapped
(missing `genbank_accession`, a failed network fetch after retries, unparseable coordinates,
etc.) is logged and left blank in the output rather than stopping the whole run, so one bad row
never costs you the rest of the table.

## Other files in this repository

The remaining top-level scripts (`extract_comments_uniprot.py`, `extract_features_Uniprot.py`,
`extract_references.py`, `get_IEBD_ID_epitope.py`, `get_source_accession_epitope.py`,
`Parse_Evidence_Uniprot.py`, `Parse_SF_category_Uniprot.py`,
`Parse_location_length_seq_epitope.py`, `Epitope_location and sequence_parse.py`) are early
prototypes whose logic has since been consolidated into the six pipeline scripts above.
`main.py` is an unused PyCharm project template.
