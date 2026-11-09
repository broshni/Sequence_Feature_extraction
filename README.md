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

Built to work for any virus: the virus is a search term and an accession, nothing in the code
is specific to one. Tested end to end on Measles, RSV (non-segmented) and influenza A
(segmented, 8 segments), and on name matching for Zika, dengue, HIV-1, hepatitis B and
SARS-CoV-2 (see **What works for any virus** below).

```
find_reference.py  (virus name -> reference genome accessions)
        |
build_reference.py  ->  suggest_aliases.py  ->  SF_Table_Pubmed.py --+
                                   |                                 |
SF_Table_Uniprot.py  --+           v                                 v
SF_Table_Epitope.py  --+--> Edit_SFTable.py --> filter_organism.py --> mapping.py
                                                                      |
                                       merge_same_SF_coordinates.py <-----+
                                                  |
                                       build_sf_ids_names.py
```

1. **`SF_Table_Uniprot.py`** — given a UniProt taxonomic-search TSV export, filters to
   reviewed (Swiss-Prot) entries only, pulls per-entry feature annotations from the UniProt
   REST API (retrying failures, and stopping loudly if any entry still cannot be retrieved),
   and builds a standardized SF table.
2. **`SF_Table_Epitope.py`** — given an IEDB epitope search export (`.xlsx`), reads the first
   sheet, cleans up IEDB's double-row header, drops "Discontinuous peptide on multi chain"
   entries and linear peptides with missing position data (tracked in a separate CSV), and
   builds a standardized SF table.
3. **`build_reference.py`** — reads the reference proteins (name, gene, accession) off one or
   more reference genome records. A segmented virus gives one accession per segment.
4. **`suggest_aliases.py`** — works out which product-name spellings in the source tables mean
   which reference protein, so no alias file is written by hand.
5. **`SF_Table_Pubmed.py`** — the third source. Searches PubMed for papers on the virus that
   report mutations, residues or epitopes, reads their abstracts, and extracts features. Every
   row is labelled `Source = PubMed miner`, cites its paper in `publication` as `PMID: n`, and
   keeps the sentence it came from in `evidence_sentence`. On by default in the UI.
6. **`Edit_SFTable.py`** — merges the source SF tables into one table and backfills
   `genbank_accession`, `gene_ref`, `product_ref` from the reference.
7. **`filter_organism.py`** — drops IEDB epitopes whose source organism is not the virus being
   built (related viruses, host proteins), judged by NCBI taxonomy lineage.
8. **`mapping.py`** — for each row, fetches its source sequence (UniProt or GenBank) and its
   reference sequence, aligns them with MAFFT, and maps the feature onto the reference in
   reference residue positions. Populates `start`/`end`/`segments` and adds
   `Source_Sequence`/`Ref_sequence`. Rows that can't be mapped are left blank rather than
   dropped, and the run continues.
9. **`merge_same_SF_coordinates.py`** — collapses rows that represent the same feature (same
   reference gene + `segments` + `sf_category`) reported by more than one source, merging the
   extra sources into the surviving row's `publication` column.
10. **`build_sf_ids_names.py`** — adds `sf_id` (e.g. `Measles_N_SF3`) and `sf_name` (e.g.
   `Measles_N_Helix_18(3)`) columns, producing the final SF table.

The last two steps validate that their required columns are fully populated before running,
and will stop with a clear report (rather than silently producing bad output) if not — see
**Data quality checks** below.

## Quick start: the UI

```
pip install -r requirements.txt
streamlit run app.py
```

Enter your email and the virus's name or taxonomy ID and press **Find reference genomes**: the
reference genome(s) are looked up at NCBI (a segmented virus is listed per strain, with every
segment's accession). Then upload the UniProt export and the IEDB export. The **PubMed miner** switch (on by default)
adds the third source; optionally upload a list of papers (a CSV with a `PMID` column, such as a
`pubmed_results.csv`) to mine as well as the automatic search. Press **Run pipeline**. Each step runs in order with its output wired into the
next; you can watch each one's log, and download every intermediate file at the end. A failed
step stops the run and shows its output rather than letting later steps run on bad input.

The look is in `ui_theme.py` (one stylesheet plus a few markup helpers), kept separate so
`app.py` stays about the pipeline. Theme defaults live in `.streamlit/config.toml`.

`run_pipeline.py` does the same thing from the command line:

```
python find_reference.py "measles virus" --email you@example.com      # lists the accessions
python run_pipeline.py <uniprot>.tsv <iedb>.xlsx NC_001498.1 -o <outdir> -n Measles \
    --pubmed-findings <miner_findings>.csv --email you@example.com
# segmented virus: one accession per segment, comma-separated
python run_pipeline.py <uniprot>.tsv <iedb>.xlsx NC_002016.1,NC_002017.1,... -o <outdir> -n Flu
```

The individual scripts below remain usable on their own, and are what both of these call.

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

### 2b. `build_reference.py` (optional, feeds step 3)

```
python build_reference.py <reference_genome_accession> -o <outdir> -n <name> [--aliases aliases.csv]
```

Builds the reference CSV that step 3 uses to backfill `genbank_accession`, by reading the CDS
features off an annotated reference genome instead of looking each protein accession up by
hand. For measles the reference genome is `NC_001498.1`, which yields 8 CDS (N, P, V, C, M, F,
H, L).

- `--aliases` (optional): a CSV with `alias` and `product` columns adding extra spellings for
  products already in the reference, e.g. `nucleoprotein` -> `nucleocapsid protein`. NCBI,
  UniProt and IEDB name the same protein differently. You normally do not need this: the
  pipeline works the aliases out for you (see `suggest_aliases.py` below).
- Output: `<name>_reference.csv`, ready to pass to step 3 as `--reference-file`.

### 2c. `suggest_aliases.py` (automatic)

```
python suggest_aliases.py <name>_reference.csv <uniprot_sf>.csv <epitope_sf>.csv -o <outdir> -n <name>
```

NCBI, UniProt and IEDB name one protein differently ("large polymerase protein" /
"RNA-directed RNA polymerase L" / "RNA dependent RNA polymerase"). This works out which names
mean which reference protein, straight from your source tables, so no alias file has to be
written by hand.

- Names are compared on the words that identify a protein, ignoring generic ones ("protein",
  "glycoprotein", "precursor"...). A bare letter such as `C` or `V` only counts when it is the
  whole name, so `Myelin basic protein` is never taken for the C protein.
- UniProt lists a protein's synonyms in brackets (`Nucleoprotein (Nucleocapsid protein)`);
  those link differently-spelled IEDB names to it.
- Anything ambiguous or unrecognised is left unmatched, never guessed. Those rows are skipped
  by mapping, and are usually host proteins or proteins the reference genome does not encode.
- On Measles it reaches 1707 of 1713 rows with no hand-written input, identical row-for-row to
  the hand-written aliases. On Zika (a polyprotein) it adds nothing and matches 97.5% through
  the single `polyprotein` reference.
- Output: `<name>_aliases.csv`, `<name>_aliases_unresolved.csv`, and
  `<name>_reference_aliased.csv` (the reference plus the aliases), which step 3 then uses as
  its `--reference-file`.

### 2d. `find_reference.py`

```
python find_reference.py "influenza A virus" --email you@example.com
```

Looks a virus up by name or NCBI taxonomy ID and lists its RefSeq reference genomes. Records
are grouped by taxon, so one group is one strain's full set of segments. Records of near-equal
length under one taxon (Zika has two) are alternative reference genomes rather than segments,
and are listed separately: choosing a reference is a scientific decision, so with several
candidates the UI makes you pick rather than defaulting.

### 3. `Edit_SFTable.py`

```
python Edit_SFTable.py <sf_table>.csv [<sf_table>.csv ...] -o <outdir> -n <name> [--reference-file reference.csv]
```

- Inputs: the per-source SF tables to merge — normally the two from steps 1 and 2, but any
  number may be given. The merge is a plain concat that never looks at which source a table
  came from, so a further source's table is just another argument; columns absent from one
  table are filled blank rather than dropped.
- `--reference-file` (optional): a CSV with `product`, `gene`, `protein_accession` columns,
  used to backfill `genbank_accession`/`gene_ref`/`product_ref` by matching product names
  (case-insensitive, on whole words). Generate it with `build_reference.py` rather than by
  hand. Skipped if omitted, which leaves `genbank_accession` empty and makes step 4 skip
  every row.
- Output: `<name>_SF_merged.csv`.

### 3b. `filter_organism.py`

```
python filter_organism.py <merged_sf>.csv --reference NC_001498.1 -o <outdir> -n <name>
```

An IEDB export for one virus also holds epitopes from related viruses and from host proteins
(a Zika export carries Dengue, West Nile and Yellow fever epitopes). Mapped onto this virus's
reference they would get coordinates that mean nothing there. Each IEDB row's source organism
is looked up in NCBI Taxonomy and kept only if it falls within the reference genome's species
(`--taxon` overrides). Only IEDB rows are judged; an organism NCBI cannot identify is kept, not
dropped. Dropped rows go to `<name>_dropped_other_organism.csv`.

### 4. `mapping.py`

```
python mapping.py <merged_sf>.csv -o <outdir> -n <name> [--email you@example.com]
```

- Input: the organism-filtered SF table (or any merged SF table).
- Output: `<name>_SF_mapped.csv`. Scratch alignment files are written to
  `<outdir>/_mapping_scratch/` rather than cluttering your working directory.
- `start`/`end`/`segments` are **reference residue positions**, not alignment columns. The two
  differ whenever a strain has an insertion relative to the reference: every column after the
  gap is shifted, so reading columns as positions misplaces every feature after it. (On
  Measles one F source with 3 extra N-terminal residues had all its features offset by 3 until
  this was fixed.) Residues inside an insertion have no reference position and are skipped.
- A row is skipped, with the reason logged, when its source sequence is less than
  `--min-identity` (default 0.5) identical to the reference protein: that is a different
  protein or a different subtype (an H4N6 neuraminidase against an N1 reference is 36%), and
  alignment coordinates between them would not mean anything.
- Each (source, reference) pair is fetched and aligned once and reused for every row that
  shares it, so a table of thousands of rows costs a few dozen alignments.
- Useful options: `--email` (required by NCBI Entrez), `--min-identity`,
  `--retries`/`--retry-delay`/`--timeout` (network resilience), `--row-delay` (pause between
  rows to avoid rate limits), `--encoding` (input CSV encoding, default `latin-1`).

### 5. `merge_same_SF_coordinates.py`

```
python merge_same_SF_coordinates.py <mapped_sf>.csv -o <outdir> -n <name> [--drop-unmapped]
```

- Input: the mapped SF table from step 4 (**must** have `segments` populated — run `mapping.py`
  first).
- `--drop-unmapped` (optional): drop rows `mapping.py` could not place on the reference
  (empty `segments`) instead of refusing to run. These are typically features on host or
  otherwise unrelated proteins with no counterpart on the reference. They are written to
  `<name>_dropped_unmapped.csv` for review.
- Rows are grouped by reference gene (`gene_ref`, falling back to `gene`), mapped segments and
  `sf_category` — the same gene resolution step 6 uses, so what counts as one feature here is
  what gets one `sf_id` there.
- Output: `<name>_dedup.csv`.

### 6. `build_sf_ids_names.py`

```
python build_sf_ids_names.py <dedup_sf>.csv -o <outdir> -n <name> [--virus-name Measles]
```

- Input: the deduplicated SF table from step 5.
- `--virus-name` (optional, defaults to `--name`): the short virus name used in `sf_id` and
  `sf_name`, e.g. `Measles_N_SF3` and `Measles_N_Helix_18(3)`. The raw `virus_type` column is
  deliberately not used — it holds a full strain description that varies per source entry
  (`Measles virus (strain Edmonston) (MeV) (...)`), which would give the same feature a
  different `sf_id` per strain.
- Genes come from `gene_ref` where available, so one reference gene gets one namespace rather
  than being split across the sources' own spellings (`N` vs `N NP`, `P/V` vs `P/V/C`).
- Output: `<name>_SF_Final.csv` — the finished SF table.

## 2c. `SF_Table_Pubmed.py` (the PubMed miner, third source)

```
python SF_Table_Pubmed.py <name>_reference.csv -o <outdir> -n <name> --virus "measles virus" \
    --accept-all [--pmid-file pubmed_results.csv] [--aliases aliases.csv] --email you@example.com
```

Mines sequence features out of PubMed abstracts. The virus's protein names come from the
reference CSV, so nothing is hardcoded per virus. The pipeline runs it with `--accept-all`.

**What lands in the table.** Every mined row has `Source = PubMed miner`, `source_id` = the
PMID, `virus_type` = the reference genome's virus name (as the IEDB rows have), and
`publication = PMID: <n>`, so it is always identifiable and traceable. The sentence it
was extracted from is kept in `evidence_sentence`, and `review_note` records the assumption
behind it (coordinates taken as relative to the reference protein, or a protein-attribution
conflict). When a feature is also reported by UniProt or IEDB it is kept once, `Source` lists every
source that reported it (`Uniprot, PubMed miner`), and the mined paper is added to that row's
`publication` as `PMID: <n>` (each paper listed once). Nothing the miner found goes missing:
on Measles all 79 mined features (66 distinct) and all 39 papers are in the final table.

- **The search query matters more than the parser.** A plain virus-name search returns reviews
  and epidemiology: on an 88-paper general Measles list, zero abstracts contained a single
  residue-level annotation. `--virus` builds a feature-targeted query instead (mutation,
  substitution, epitope, amino acid, residue, escape), which on Measles returns 642 papers and
  mines ~79 features from 39 of them; on RSV, 52 features from 31 papers. `--pmid-file` adds
  papers of your own (a CSV with a `PMID` column) to the search; a paper in both is mined once.
- **Every hit must be gene-attributed.** A mutation only counts if one of the virus's proteins
  is named in the same sentence. Without this, host genes mentioned in passing (BRAF `V600E` in
  an oncolytic-therapy paper, RNF213 `R4810K` in a vasculopathy study) become fake viral features.
- `sf_type` is `mutation`, `epitope` or `region`, mapping onto the existing `sf_category`
  vocabulary (`Sequence_Alterations`, `Epitope`, `Functional`). `Evidence_code` is `TAS`.
- **Ambiguous attributions are kept, flagged.** A paper that names several proteins can yield
  one position under each, and at most one is right (a paper on the matrix protein's F50P, in a
  sentence that also names the hemagglutinin). Every such row stays in the table with
  `review_note = protein attribution conflict: F, M`, so nothing the miner found is silently
  discarded and the doubtful rows are easy to filter. (Counting which protein a paper names
  most does not reliably pick the right one, so no guess is made.)
- **It is text mining, so it is not always right.** Mined rows are marked as such precisely so
  they can be filtered (`Source`) and checked (`evidence_sentence`). If PubMed cannot be
  reached, the run continues with an empty PubMed table rather than losing the other sources.
- **It needs a reference with one entry per protein.** A polyprotein virus (Zika) has a single
  "polyprotein" reference, so papers that say "NS1" or "envelope protein" are not recognised
  and very little is mined.

**Curated route.** Without `--accept-all`, candidates are written to
`<name>_Pubmed_SF_review.csv` with `needs_review=yes` instead. Check each against its
`evidence_sentence`, delete the wrong ones, set `needs_review` to `no` on the rest, then:

```
python SF_Table_Pubmed.py -o <outdir> -n <name> --reviewed <name>_Pubmed_SF_review.csv
```

and pass the resulting `<name>_Pubmed_SF.csv` to `run_pipeline.py --pubmed-sf` to use it in
place of mining.

## What works for any virus

Reference genomes are found by name or taxonomy ID, protein names are matched from your own
tables, and organisms are compared by NCBI lineage, so a new virus needs no code or
configuration. What has been checked:

| Virus | Kind | Checked |
|---|---|---|
| Measles | non-segmented | full pipeline |
| RSV | non-segmented | full pipeline; 1347 of 1347 mapped features' reference residues match `start..end` |
| Influenza A | 8 segments, 12 proteins | full pipeline on 30 real UniProt entries |
| Zika, dengue | polyprotein | name matching; organism filter on Zika (drops the 97 non-Zika epitopes) |
| HIV-1, hepatitis B, SARS-CoV-2 | various | name matching on the real UniProt exports (100%, 98.7%, 71%) |

Behaviour worth knowing:

- **Polyprotein viruses** have one reference protein (the polyprotein), so every row maps onto
  it. Rows that are not really that protein are rejected by the identity check at mapping time.
  Mature peptides (E, NS1...) are not used as references: a feature in a polyprotein source
  is numbered along the whole polyprotein.
- **Unmatched rows are left out, never guessed.** The UI lists them. They are usually host
  proteins, or accessory proteins the reference genome does not annotate (SARS-CoV-2 ORF3b/9b,
  RSV M2-2). If one is a protein the reference does have under another name, add the pair
  under *Extra product aliases* and run again.
- **The identity threshold is a judgement.** 0.5 rejects other subtypes and other proteins; a
  highly diverse virus (HIV-1 subtypes, HCV genotypes) may need `--min-identity` lowered.
- **A gene symbol can be a word pandas treats as missing** (`NA`, influenza neuraminidase). All
  tables are read through `sf_columns.read_csv`, where only an empty cell is missing.

## Coordinate semantics: ranges vs. residue pairs

`source_sf_location` is written `start-end` for almost every feature, but that notation means
two different things depending on the feature type:

- **Most types are contiguous.** `Helix 159-163`, `Chain 1-525`, `Transmembrane 38-58` and IEDB
  epitopes cover every residue in the span. The IEDB rows confirm this independently: each
  epitope's peptide sequence has exactly `end - start + 1` residues.
- **Disulfide bonds and cross-links are residue pairs.** UniProt annotates a disulfide bond as
  the two cysteines that are bonded, so `188-606` means Cys188 bonded to Cys606 — not a
  419-residue feature. These types are listed in `PAIR_SF_TYPES` in `mapping.py`.

For the pair types, `length` is 2 (or 1 if both positions are the same), `segments` is written
comma-separated (`188,606`) so nothing downstream reads it as a span, and the extracted
`Source_Sequence`/`Ref_sequence` hold just the two bonded residues (`CC`). Add a feature type
to `PAIR_SF_TYPES` if it should be read the same way.

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

## Tests

```
python test_build_reference.py
```

Asserts the logic that does not need the network: reference building, product matching (word
and hyphen boundaries, single-protein references, gene symbols such as `NA`), automatic alias
suggestion, the organism filter, reference-genome grouping, and coordinate mapping (including
a source with an N-terminal insertion, where alignment columns and reference positions
differ).

## Other files in this repository

`run_pipeline.py` runs the whole pipeline; `app.py` with `ui_theme.py` and `.streamlit/` is the
Streamlit front end. `sf_columns.py` holds the column conventions every step shares (how tables
are read and written), `taxonomy.py` the NCBI Taxonomy lookups, and `find_reference.py` finds a
virus's reference genome by name. `test_build_reference.py` is the self-check suite.
