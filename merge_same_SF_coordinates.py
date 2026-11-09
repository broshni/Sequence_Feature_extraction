#collapse SF-table rows that represent the same feature (same gene, mapped coordinates, and
#SF category) reported by more than one source, merging the extra sources' identifiers into
#the surviving row's 'publication' column. Run this after mapping.py has populated 'segments'.
import argparse
import os
import re
import sys
import pandas as pd

# Same gene resolution build_sf_ids_names.py uses, so a feature that gets one sf_id there is
# one group here. Grouping on the raw 'gene' instead splits a reference gene across the
# sources' own spellings ('P/V' vs 'P/V/C' vs 'P/V P') and leaves real duplicates unmerged.
from build_sf_ids_names import resolve_gene
from sf_columns import publication_label, read_csv, write_csv

# This script identifies duplicate features by grouping on 'gene' + 'segments' + 'sf_category',
# and pulls 'Source'/'source_id' into the merged row. If any of those cells are empty (e.g.
# 'segments' hasn't been populated by mapping.py yet), rows would silently get grouped
# together incorrectly - e.g. every row with a blank gene AND a blank segments would look
# like the "same" feature. So check upfront that these cells are populated, and refuse to
# proceed (rather than produce a wrong merge) if they're not.
REQUIRED_COLUMNS = ['gene', 'segments', 'sf_category', 'Source', 'source_id']


def drop_unmapped_rows(df):
    """Drop rows mapping.py could not place on the reference (empty 'segments'). They have no
    reference coordinates, so they can't be grouped or merged as features at all."""
    is_unmapped = df['segments'].isna() | (df['segments'].astype(str).str.strip() == '')
    return df[~is_unmapped].reset_index(drop=True), df[is_unmapped]


def find_incomplete_columns(df, required_columns):
    """Return a list of human-readable problem descriptions for any required column that is
    missing entirely or has empty/blank cells. Empty list means everything checks out."""
    problems = []
    for col in required_columns:
        if col not in df.columns:
            problems.append(f"  - '{col}': column not found in the input table at all")
            continue
        is_empty = df[col].isna() | (df[col].astype(str).str.strip() == '')
        empty_count = int(is_empty.sum())
        if empty_count:
            example_rows = df.index[is_empty].tolist()[:5]
            more = '' if empty_count <= 5 else f" (+{empty_count - 5} more)"
            problems.append(
                f"  - '{col}': {empty_count} of {len(df)} row(s) are empty - example row "
                f"index(es): {example_rows}{more}"
            )
    return problems


def already_cited(label, cited_text):
    """True when `label` ('PMID: 123' or 'IEDB:45') is already in a publication string. A PMID
    counts as cited wherever its number appears: UniProt writes 'PMID: 111, 222', so the second
    paper has no 'PMID:' of its own."""
    if label.startswith('PMID: '):
        return label[len('PMID: '):] in re.findall(r'\d+', cited_text)
    return label in cited_text


#function to merge rows that share the same gene/segments/sf_category, combining the
#merged-away rows' 'source:id' values into the surviving row's 'publication' column.
def merge_duplicate_sf_rows(df):
    df = df.copy()
    df['gene_segments_category'] = (
        resolve_gene(df) + '_' + df['segments'].astype(str) + '_' + df['sf_category'].astype(str)
    )
    df['source:id'] = [publication_label(s, i) for s, i in zip(df['Source'], df['source_id'])]

    # Ensure the 'publication' column exists so we always merge into the same column
    if 'publication' not in df.columns:
        df['publication'] = None

    merged_count = 0
    for gene_segment_category, group in df.groupby('gene_segments_category'):
        if len(group) > 1:
            # Keep the first row and merge the 'source:id' values of the rest into its
            # existing 'publication' value
            first_row_index = group.index[0]
            remaining_source_ids = group['source:id'].iloc[1:].tolist()

            existing_publication = df.at[first_row_index, 'publication']
            existing_publication = '' if pd.isna(existing_publication) else str(existing_publication)

            # Combine existing publication value (if any) with the merged-away rows' source:id
            # values, citing each paper once: UniProt and the PubMed miner often cite the same one.
            # Source and source_id stay the surviving row's own - one row names exactly one
            # source, and its id is that source's identifier. The sources that were merged away
            # are recorded in 'publication' instead.
            merged_values = [existing_publication] if existing_publication else []
            for label in remaining_source_ids:
                if not already_cited(label, ', '.join(merged_values)):
                    merged_values.append(label)
            df.at[first_row_index, 'publication'] = ', '.join(merged_values)

            # Drop the remaining rows
            df = df.drop(group.index[1:])
            merged_count += len(group) - 1

    df = df.drop(columns=['gene_segments_category', 'source:id'])
    df = df.reset_index(drop=True)
    return df, merged_count


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collapse SF table rows that represent the same feature (same gene, "
                    "mapped segments, and sf_category) reported by more than one source. "
                    "Requires 'segments' to already be populated (e.g. via mapping.py)."
    )
    parser.add_argument(
        "input_file",
        help="Path to the mapped SF table CSV to deduplicate (e.g. output of mapping.py).",
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Directory where the deduplicated CSV will be written (created if missing).",
    )
    parser.add_argument(
        "-n", "--name",
        default="output",
        help="Prefix used to name the output file, e.g. '<name>_dedup.csv' (default: 'output').",
    )
    parser.add_argument(
        "--drop-unmapped",
        action="store_true",
        help="Drop rows that mapping.py left unmapped (empty 'segments') instead of refusing "
             "to run. These are typically features on host or unrelated proteins that have no "
             "counterpart on the reference sequence. Dropped rows are written to "
             "'<name>_dropped_unmapped.csv' so they can be reviewed.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = read_csv(args.input_file)
    print(df.columns)

    if args.drop_unmapped and 'segments' in df.columns:
        df, dropped = drop_unmapped_rows(df)
        if len(dropped):
            dropped_path = os.path.join(args.outdir, f"{args.name}_dropped_unmapped.csv")
            write_csv(dropped, dropped_path)
            print(f"Dropped {len(dropped)} unmapped row(s) -> {dropped_path}")

    problems = find_incomplete_columns(df, REQUIRED_COLUMNS)
    if problems:
        print("Cannot proceed - the following important cell(s) need to be populated before merging duplicate SF coordinates:")
        print('\n'.join(problems))
        print("\nPopulate these columns (e.g. run mapping.py first to fill in 'segments') and re-run this script.")
        print("If the empty cells are rows mapping.py could not place on the reference, "
              "re-run with --drop-unmapped to exclude them.")
        sys.exit(1)

    df, merged_count = merge_duplicate_sf_rows(df)

    print(df)
    if merged_count:
        print(f"Merged away {merged_count} duplicate row(s); {len(df)} row(s) remain.")
    else:
        print(f"No duplicate rows found; all {len(df)} row(s) kept as-is.")

    output_file = os.path.join(args.outdir, f"{args.name}_dedup.csv")
    write_csv(df, output_file)
    print(f"Deduplicated SF table saved to {output_file}")


if __name__ == "__main__":
    main()
