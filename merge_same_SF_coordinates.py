#collapse SF-table rows that represent the same feature (same gene, mapped coordinates, and
#SF category) reported by more than one source, merging the extra sources' identifiers into
#the surviving row's 'publication' column. Run this after mapping.py has populated 'segments'.
import argparse
import os
import sys
import pandas as pd

# This script identifies duplicate features by grouping on 'gene' + 'segments' + 'sf_category',
# and pulls 'Source'/'source_id' into the merged row. If any of those cells are empty (e.g.
# 'segments' hasn't been populated by mapping.py yet), rows would silently get grouped
# together incorrectly - e.g. every row with a blank gene AND a blank segments would look
# like the "same" feature. So check upfront that these cells are populated, and refuse to
# proceed (rather than produce a wrong merge) if they're not.
REQUIRED_COLUMNS = ['gene', 'segments', 'sf_category', 'Source', 'source_id']


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


#function to merge rows that share the same gene/segments/sf_category, combining the
#merged-away rows' 'source:id' values into the surviving row's 'publication' column.
def merge_duplicate_sf_rows(df):
    df = df.copy()
    df['gene_segments_category'] = (
        df['gene'].astype(str) + '_' + df['segments'].astype(str) + '_' + df['sf_category'].astype(str)
    )
    df['source:id'] = df['Source'].astype(str) + ':' + df['source_id'].astype(str)

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

            # Combine existing publication value (if any) with the merged-away rows' source:id values
            merged_values = [v for v in [existing_publication] + remaining_source_ids if v]
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
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input_file)
    print(df.columns)

    problems = find_incomplete_columns(df, REQUIRED_COLUMNS)
    if problems:
        print("Cannot proceed - the following important cell(s) need to be populated before merging duplicate SF coordinates:")
        print('\n'.join(problems))
        print("\nPopulate these columns (e.g. run mapping.py first to fill in 'segments') and re-run this script.")
        sys.exit(1)

    df, merged_count = merge_duplicate_sf_rows(df)

    print(df)
    if merged_count:
        print(f"Merged away {merged_count} duplicate row(s); {len(df)} row(s) remain.")
    else:
        print(f"No duplicate rows found; all {len(df)} row(s) kept as-is.")

    output_file = os.path.join(args.outdir, f"{args.name}_dedup.csv")
    df.to_csv(output_file, index=False)
    print(f"Deduplicated SF table saved to {output_file}")


if __name__ == "__main__":
    main()
