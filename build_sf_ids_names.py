#read the SF table and add sf_id/sf_name columns, producing the final SF table.
import argparse
import os
import sys
import pandas as pd

from sf_columns import PUBMED_SOURCE, read_csv, write_csv

# create_sf_ids() groups by gene, and create_sf_names() reads 'gene', 'sf_type', 'start' and
# 'length' for every row. A blank 'gene' is especially dangerous here: pandas' groupby drops
# NaN-key rows from cumcount() by default, so every row with a missing gene would collide on
# the same non-unique sf_id (e.g. "Measles_nan_SFnan" for all of them). So check upfront that
# these cells are populated, and refuse to proceed (rather than silently produce
# broken/duplicate IDs) if they're not.
REQUIRED_COLUMNS = ['gene', 'sf_type', 'start', 'length']


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


def as_int(value):
    """'525.0' -> '525'. Coordinates and lengths are whole numbers; they only arrive looking
    like floats because pandas reads any column containing a blank cell as float64."""
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def resolve_gene(df):
    """The gene each sf_id is keyed on. Prefer 'gene_ref' (the reference genome's gene name):
    sf_ids describe a position on the reference, while the sources' own gene strings vary by
    entry for what is one reference gene ('N' vs 'N NP', 'P' vs 'P/V' vs 'P/V P'), which would
    otherwise split one gene's IDs across several namespaces."""
    if 'gene_ref' not in df.columns:
        return df['gene'].astype(str)
    gene_ref = df['gene_ref']
    has_ref = gene_ref.notna() & (gene_ref.astype(str).str.strip() != '')
    return gene_ref.where(has_ref, df['gene']).astype(str)


#create sf_id in the format: virus_gene_SF followed by a number unique within that gene. For
#example, with virus name SARS-CoV-2 and gene S the sf_ids are SARS-CoV-2_S_SF1,
#SARS-CoV-2_S_SF2, etc. The virus name is one constant for the whole table, so the same
#feature does not get a different id per source strain.
def create_sf_ids(df, virus_name):
    gene = resolve_gene(df)
    numbers = df.groupby(gene).cumcount() + 1
    df['sf_id'] = [f"{virus_name}_{g}_SF{n}" for g, n in zip(gene, numbers)]
    return df


# create sf_names in the format: virus_gene_sf_type_start(length).
def create_sf_names(df, virus_name):
    gene = resolve_gene(df)
    df['sf_name'] = [
        f"{virus_name}_{g}_{t}_{as_int(start)}({as_int(length)})"
        for g, t, start, length in zip(gene, df['sf_type'], df['start'], df['length'])
    ]
    return df


def name_virus_column(df, virus_name):
    """Give the table one clean virus_type - the same name the sf_ids use - and keep each row's
    original organism string as virus_type_source.

    The sources name the organism in their own way, at strain level and with their own
    punctuation ('Ebola virus (strain Mayinga-76)', 'Zaire ebolavirus'), so the column cannot
    be used to say which virus the table is for. Renaming it keeps that provenance while
    virus_type becomes a single, uniform value.
    """
    # on a re-run the table already has virus_type_source; renaming again would duplicate it
    if 'virus_type_source' not in df.columns:
        df = df.rename(columns={'virus_type': 'virus_type_source'})
    df = df.drop(columns=['virus_type'], errors='ignore')
    df.insert(0, 'virus_type', virus_name)
    return df


def add_pubmed_flag(df):
    """A plain yes/no column saying whether a row came from the PubMed miner, placed directly
    after source_id so the provenance reads together: Source, source_id, and this."""
    flag = 'sourced from pubmed miner'
    df = df.drop(columns=[flag], errors='ignore')
    values = (df['Source'].astype(str).str.strip() == PUBMED_SOURCE).map({True: 'yes', False: 'no'})
    df.insert(df.columns.get_loc('source_id') + 1, flag, values)
    return df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add sf_id and sf_name columns to a deduplicated SF table, producing the "
                    "final SF table. Requires gene/sf_type/start/length to already "
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
        "--virus-name",
        default=None,
        help="Short virus name used to build sf_id/sf_name, e.g. 'Measles' giving "
             "'Measles_N_SF3'. Defaults to --name. The raw 'virus_type' column is not used: "
             "it carries a full strain description that differs per source entry.",
    )
    parser.add_argument(
        "--encoding",
        default=None,
        help="Force a text encoding for the input CSV. By default it is read as UTF-8, "
             "falling back per byte to Windows-1252.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sf_table = read_csv(args.input_file, **({'encoding': args.encoding} if args.encoding else {}))

    problems = find_incomplete_columns(sf_table, REQUIRED_COLUMNS)
    if problems:
        print("Cannot proceed - the following important cell(s) need to be populated before building sf_id/sf_name:")
        print('\n'.join(problems))
        print("\nPopulate these columns and re-run this script.")
        sys.exit(1)

    virus_name = args.virus_name or args.name
    sf_table = name_virus_column(sf_table, virus_name)
    sf_table = add_pubmed_flag(sf_table)
    sf_table = create_sf_ids(sf_table, virus_name)
    sf_table = create_sf_names(sf_table, virus_name)
    print(sf_table)

    output_file = os.path.join(args.outdir, f"{args.name}_SF_Final.csv")
    write_csv(sf_table, output_file)
    print(f"Final SF table saved to {output_file}")


if __name__ == "__main__":
    main()
