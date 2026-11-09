#edit the SF table by adding the references for each protein
import argparse
import os
import re
import pandas as pd

from sf_columns import read_csv, write_csv

#stack any number of per-source SF tables (Uniprot, IEDB, and any further source) into one.
#The merge is a plain concat and never looks at which source a table came from, so a new
#source needs no change here - just pass its table as another argument.
def merge_sf_tables(sf_files):
    frames = []
    for path in sf_files:
        df = read_csv(path)
        print(f"  {os.path.basename(path)}: {len(df)} row(s)")
        frames.append(df)
    # columns absent from one table are filled with NaN rather than dropped, so a source that
    # carries an extra column does not quietly lose it
    return pd.concat(frames, ignore_index=True)

# write a function to check the reference file and if the gene column from reference file is in the
# gene column of the merged DataFrame, then add the value in the protein_accession value to genbank_accession column i
# in the merged DataFrame

def name_occurs(name, text):
    """True when `name` appears in `text` as a whole name.

    Word-boundary match, not a bare substring: a plain 'in' lets an unrelated product like
    'myelin basic protein' match the reference product 'C protein' (via 'basi-c protein').
    A hyphen counts as part of a word as well, so 'PB1' is not found inside 'PB1-F2', which is
    a different protein. suggest_aliases.py uses this same rule, so an alias it writes is
    guaranteed to match here.
    """
    return re.search(rf"(?<![\w-]){re.escape(str(name).lower())}(?![\w-])",
                     str(text).lower()) is not None


def find_reference_matches(product, df_reference):
    """Every reference row whose product name appears in `product` as a whole name."""
    return [ref_row for _, ref_row in df_reference.iterrows()
            if name_occurs(ref_row['product'], product)]


def check_and_update_genbank_accession(merged_df, reference_file):
    df_reference = read_csv(reference_file)
    single_protein = df_reference['protein_accession'].nunique() == 1
    fallback_rows = 0
    for col in ['genbank_accession', 'gene_ref', 'product_ref', 'gene']:
        if col not in merged_df.columns:
            merged_df[col] = None
        # An all-empty column read from CSV (or freshly created above) infers as float64,
        # which raises a TypeError the first time a string accession/gene is written into
        # it below. Force it to a string-friendly dtype up front.
        merged_df[col] = merged_df[col].astype('object')
    for index, row in merged_df.iterrows():
        matched_accession = None
        matched_gene = None
        matched_product = None
        matches = find_reference_matches(row['product'], df_reference)
        if not matches and single_protein:
            # One-protein reference (a polyprotein virus): every protein of the virus is part
            # of it, whatever name a source gives it. A row that is not really that protein is
            # rejected later, when mapping finds its sequence does not resemble the reference.
            matches = [df_reference.iloc[0]]
            fallback_rows += 1
        if matches:
            # several reference rows can match; the last one wins
            ref_row = matches[-1]
            matched_accession = ref_row['protein_accession']
            matched_gene = ref_row['gene']
            matched_product = ref_row['product']
        if matched_accession:
            merged_df.at[index, 'genbank_accession'] = matched_accession
            merged_df.at[index, 'gene_ref'] = matched_gene
            merged_df.at[index, 'product_ref'] = matched_product
            # IEDB rows arrive with no gene at all, and downstream steps group on 'gene'
            # and refuse to run while it's empty. The reference match just told us which
            # gene this row belongs to, so use it - without overwriting a gene the source
            # already supplied.
            existing_gene = merged_df.at[index, 'gene'] if 'gene' in merged_df.columns else None
            if pd.isna(existing_gene) or not str(existing_gene).strip():
                merged_df.at[index, 'gene'] = matched_gene
    if fallback_rows:
        print(f"  {fallback_rows} row(s) matched by name to nothing and were assigned to the "
              f"reference's only protein ({df_reference.iloc[0]['product']})")
    return merged_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge per-source SF tables into one combined table, optionally "
                    "backfilling genbank_accession from a reference file."
    )
    parser.add_argument(
        "sf_tables",
        nargs='+',
        help="Paths to the per-source SF table CSVs to merge, e.g. the outputs of "
             "SF_Table_Uniprot.py and SF_Table_Epitope.py. Any number may be given, so a "
             "further source's table is just another argument.",
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Directory where the merged CSV will be written (created if missing).",
    )
    parser.add_argument(
        "-n", "--name",
        default="output",
        help="Prefix used to name the merged output file, e.g. '<name>_SF_merged.csv' "
             "(default: 'output').",
    )
    parser.add_argument(
        "--reference-file",
        default=None,
        help="Optional reference CSV (with 'product', 'gene', 'protein_accession' columns) "
             "used to backfill genbank_accession/gene_ref/product_ref by matching product "
             "names. If omitted, this step is skipped.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Merging {len(args.sf_tables)} SF table(s):")
    merged_df_all = merge_sf_tables(args.sf_tables)

    if args.reference_file:
        merged_df_all = check_and_update_genbank_accession(merged_df_all, args.reference_file)
        print(merged_df_all[['gene', 'gene_ref', 'genbank_accession']])
    else:
        print(merged_df_all[['gene', 'genbank_accession']])

    # Save the updated merged DataFrame to a new CSV file
    merged_dataframe_file_path = os.path.join(args.outdir, f"{args.name}_SF_merged.csv")
    write_csv(merged_df_all, merged_dataframe_file_path)
    print(f"Merged SF table saved to {merged_dataframe_file_path}")


if __name__ == "__main__":
    main()
