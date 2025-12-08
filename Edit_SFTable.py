#edit the SF table by adding the references for each protein
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


merged_df_All = merge_sf_tables('/Users/rbhattac/Desktop/Measles_SFVT/Measles_Uniprot.csv', '/Users/rbhattac/Desktop/Measles_SFVT/Measles_epitope.csv')
#print(merged_df_All)

#change the virus_type column entries to our desired format
merged_df_All['virus_type'] = 'Measles Virus'
# Update the genbank_accession column based on the reference file
merged_df_All = check_and_update_genbank_accession(merged_df_All, '/Users/rbhattac/Desktop/Measles_SFVT/reference.csv')


#print the gene product and genbank_accession columns
print(merged_df_All[['gene', 'gene_ref', 'genbank_accession']])

# Save the updated merged DataFrame to a new CSV file
merged_dataframe_file_path = '/Users/rbhattac/Desktop/Measles_SFVT/Measles_SF_table_ALL.csv'
merged_df_All.to_csv(merged_dataframe_file_path, index=False)
print(f"Merged SF table saved to {merged_dataframe_file_path}")
