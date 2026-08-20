#edit the SF table by adding the references for each protein
import argparse
import os
import pandas as pd

#read two SF tables (Uniprot and IEDB) and merge them
def merge_sf_tables(uniprot_file, iedb_file):
    # Read the Uniprot and IEDB SF tables
    df_uniprot = pd.read_csv(uniprot_file)
    df_iedb = pd.read_csv(iedb_file)

    # Merge the two DataFrames on 'SF_name'
    merged_df_All = pd.concat([df_uniprot, df_iedb], ignore_index=True)

    # Save the merged DataFrame to a new CSV file
    #merged_df.to_csv(output_file, index=False)
    #print(f"Merged SF table saved to {output_file}")
    return merged_df_All

# write a function to check the reference file and if the gene column from reference file is in the
# gene column of the merged DataFrame, then add the value in the protein_accession value to genbank_accession column i
# in the merged DataFrame

def check_and_update_genbank_accession(merged_df, reference_file):
    df_reference = pd.read_csv(reference_file)
    for col in ['genbank_accession', 'gene_ref', 'product_ref']:
        if col not in merged_df.columns:
            merged_df[col] = None
        # An all-empty column read from CSV (or freshly created above) infers as float64,
        # which raises a TypeError the first time a string accession/gene is written into
        # it below. Force it to a string-friendly dtype up front.
        merged_df[col] = merged_df[col].astype('object')
    for index, row in merged_df.iterrows():
        merged_product = str(row['product']).lower()
        matched_accession = None
        matched_gene = None
        matched_product = None
        for _, ref_row in df_reference.iterrows():
            ref_product = str(ref_row['product']).lower()
            if ref_product in merged_product:
                matched_accession = ref_row['protein_accession']
                matched_gene = ref_row['gene']
                matched_product = ref_row['product']
        if matched_accession:
            merged_df.at[index, 'genbank_accession'] = matched_accession
            merged_df.at[index, 'gene_ref'] = matched_gene
            merged_df.at[index, 'product_ref'] = matched_product
    return merged_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge a Uniprot SF table and an IEDB epitope SF table into one combined "
                    "table, optionally backfilling genbank_accession from a reference file."
    )
    parser.add_argument(
        "uniprot_file",
        help="Path to the Uniprot SF table CSV (e.g. output of SF_Table_Uniprot.py).",
    )
    parser.add_argument(
        "iedb_file",
        help="Path to the IEDB epitope SF table CSV (e.g. output of SF_Table_Epitope.py).",
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

    merged_df_all = merge_sf_tables(args.uniprot_file, args.iedb_file)

    if args.reference_file:
        merged_df_all = check_and_update_genbank_accession(merged_df_all, args.reference_file)
        print(merged_df_all[['gene', 'gene_ref', 'genbank_accession']])
    else:
        print(merged_df_all[['gene', 'genbank_accession']])

    # Save the updated merged DataFrame to a new CSV file
    merged_dataframe_file_path = os.path.join(args.outdir, f"{args.name}_SF_merged.csv")
    merged_df_all.to_csv(merged_dataframe_file_path, index=False)
    print(f"Merged SF table saved to {merged_dataframe_file_path}")


if __name__ == "__main__":
    main()
