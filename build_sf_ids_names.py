#read the SF table and add sf_id/sf_name columns, producing the final SF table.
import argparse
import os
import sys
import pandas as pd

# create_sf_ids() groups by 'virus_type'+'gene', and create_sf_names() reads 'virus_type',
# 'gene', 'sf_type', 'start' and 'length' for every row. A blank 'gene' is especially
# dangerous here: pandas' groupby drops NaN-key rows from cumcount() by default, so every
# row with a missing gene would collide on the same non-unique sf_id (e.g.
# "Measles_nan_SFnan" for all of them). So check upfront that these cells are populated,
# and refuse to proceed (rather than silently produce broken/duplicate IDs) if they're not.
REQUIRED_COLUMNS = ['virus_type', 'gene', 'sf_type', 'start', 'length']


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


#create sf_id in the format: virus_type_gene_SF followed by a unique number. The unique number should be grouped by the
#genes and virus_type. For example, if the virus_type is SARS-CoV-2 and the gene is S, then the SF_ids should be
#SARS-CoV-2_S_SF1, SARS-CoV-2_S_SF2, etc. for each SF in the gene S. Repeat for all genes and virus_types in the table.
def create_sf_ids(df):
    df['sf_id'] = df.groupby(['virus_type', 'gene']).cumcount() + 1
    df['sf_id'] = df.apply(lambda x: f"{x['virus_type']}_{x['gene']}_SF{x['sf_id']}", axis=1)
    df['sf_id'] = df['sf_id'].str.replace(r'\.0$', '', regex=True)
    return df


# create sf_names in the format: virus_type_gene_sf_type_start(length).
def create_sf_names(df):
    df['sf_name'] = df.apply(
        lambda x: f"{x['virus_type']}_{x['gene']}_{x['sf_type']}_{x['start']}({x['length']})",
        axis=1
    )
    return df


#function to remove .0 from the sf_name example change "Measles Virus_N_Compositional bias_433.0(20)" to "Measles Virus_N_Compositional bias_433(20)"
def remove_decimal_from_sf_name(df):
    df['sf_name'] = df['sf_name'].str.replace(r'\.0\(', '(', regex=True)
    return df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add sf_id and sf_name columns to a deduplicated SF table, producing the "
                    "final SF table. Requires virus_type/gene/sf_type/start/length to already "
                    "be fully populated (e.g. via mapping.py and merge_same_SF_coordinates.py)."
    )
    parser.add_argument(
        "input_file",
        help="Path to the deduplicated SF table CSV (e.g. output of merge_same_SF_coordinates.py).",
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Directory where the final CSV will be written (created if missing).",
    )
    parser.add_argument(
        "-n", "--name",
        default="output",
        help="Prefix used to name the output file, e.g. '<name>_SF_Final.csv' (default: 'output').",
    )
    parser.add_argument(
        "--encoding",
        default="latin-1",
        help="Text encoding to use when reading the input CSV (default: 'latin-1').",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sf_table = pd.read_csv(args.input_file, encoding=args.encoding)

    problems = find_incomplete_columns(sf_table, REQUIRED_COLUMNS)
    if problems:
        print("Cannot proceed - the following important cell(s) need to be populated before building sf_id/sf_name:")
        print('\n'.join(problems))
        print("\nPopulate these columns and re-run this script.")
        sys.exit(1)

    sf_table = create_sf_ids(sf_table)
    sf_table = create_sf_names(sf_table)
    sf_table = remove_decimal_from_sf_name(sf_table)
    print(sf_table)

    output_file = os.path.join(args.outdir, f"{args.name}_SF_Final.csv")
    sf_table.to_csv(output_file, index=False)
    print(f"Final SF table saved to {output_file}")


if __name__ == "__main__":
    main()
