import pandas as pd
import requests
import time
import os
import ast

#main_folder: path to the main folder where the output file will be saved
main_folder = "/Users/rbhattac/Desktop/Chikungunya_SFVT"

#function to extract UniProt features from a TSV file obtained from taxonomic search of Uniprot
# and output a csv file with the all extracted features from multiple UniProt IDs.
#tsv_file_path: path to the TSV file containing UniProt IDs from taxonomic search

def extract_and_merge_uniprot_features(tsv_file_path):
    # Read TSV file for UniProt IDs
    uniprot_info = pd.read_csv(tsv_file_path, delimiter="\t")
    uniprot_ids = uniprot_info['Entry']

    base_url = "https://rest.uniprot.org/uniprotkb/{}.json"
    fields = (
        "ft_var_seq,ft_variant,ft_non_cons,ft_non_std,ft_non_ter,ft_conflict,"
        "ft_unsure,ft_act_site,ft_binding,ft_dna_bind,ft_site,ft_mutagen,"
        "ft_intramem,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,"
        "ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,ft_propep,"
        "ft_signal,ft_transit,ft_strand,ft_helix,ft_turn,ft_coiled,ft_compbias,"
        "ft_domain,ft_motif,ft_region,ft_repeat,ft_zn_fing"
    )

    global_df = pd.DataFrame()

    for uniprot_id in uniprot_ids:
        query_url = base_url.format(uniprot_id) + "?fields=" + fields
        try:
            response = requests.get(query_url)
            response.raise_for_status()
            data = response.json()
            if 'features' in data:
                dict_df = pd.json_normalize(data['features'])
                dict_df['uniprot_id'] = uniprot_id
                print(f"Data for {uniprot_id} retrieved successfully.")
                global_df = pd.concat([global_df, dict_df], ignore_index=True)
        except Exception as e:
            print(f"An error occurred for {uniprot_id}: {e}")
        time.sleep(1)

    merged_df = pd.merge(global_df, uniprot_info, left_on='uniprot_id', right_on='Entry')
    return merged_df

# function to extract publication information from the 'evidences' column of the merged DataFrame
def parse_evidence_column(df):
    results = []
    for idx, row in enumerate(df['evidences']):
        # Check for NaN only if row is a float (typical for NaN in pandas)
        if isinstance(row, float) and pd.isna(row):
            #print(f"NaN encountered at row index: {idx}")
            results.append({})
            continue
        # Only evaluate if row is a string
        if isinstance(row, str):
            row = ast.literal_eval(row)
        pubmed_ids = []
        for evidence in row:
            if 'source' in evidence and 'id' in evidence:
                if evidence['source'].lower() in ['pubmed', 'pmid']:
                    pubmed_ids.append(evidence['id'])
        result_row = {}
        if pubmed_ids:
            result_row['Publication'] = 'PMID: ' + ', '.join(pubmed_ids)
        results.append(result_row)
    results_df = pd.DataFrame(results)
    return results_df

# function to populate the SF_Category based on the 'type' column in the merged DataFrame
def parse_sf_category(merged_df):
    sf_categories = []
    for index, row in merged_df.iterrows():
        feature_type = row['type']
        if feature_type in ['Beta strand', 'Helix', 'Turn', 'Coiled coil']:
            sf_categories.append('Structural')
        elif feature_type == 'Mutagenesis':
            sf_categories.append('Sequence_Alterations')
        else:
            sf_categories.append('Functional')
    return pd.Series(sf_categories)

# function to extract comments from the merged DataFrame
def extract_comments(merged_df):
    comments = []
    for index, row in merged_df.iterrows():
        comment_parts = []
        comment_parts.append(str(row['type']))
        if pd.notna(row['description']):
            comment_parts.append(row['description'])
        if pd.notna(row['alternativeSequence.originalSequence']):
            comment_parts.append('original sequence:' + str(row['alternativeSequence.originalSequence']))
        alt_seq = row['alternativeSequence.alternativeSequences']
        if isinstance(alt_seq, (list, tuple)):
            if len(alt_seq) > 0:
                comment_parts.append('variant sequence:' + str(alt_seq))
        elif pd.notna(alt_seq):
            comment_parts.append('variant sequence:' + str(alt_seq))
        if pd.notna(row['ligand.name']):
            comment_parts.append('Binding site for:' + str(row['ligand.name']))

        #if row['type'] == 'Natural variant':
            #comment_parts.insert(0, 'Natural variation observed')

        comments.append(', '.join(comment_parts))

    return pd.Series(comments)

   

merged_df = extract_and_merge_uniprot_features('/Users/rbhattac/Desktop/Chikungunya_SFVT/Chikungunya_uniprot.tsv')
# Save the merged DataFrame to a CSV file
merged_dataframe_file_path = os.path.join(main_folder, "merged_features.csv")
merged_df.to_csv(merged_dataframe_file_path, index=False)
# print(f"All features saved to {merged_dataframe_file_path}.")
# print("Processing completed.")

#create an empty dataframe for the new SF table with the features extracted from UniProt
columns = ['virus_type', 'subtype', 'segment', 'genbank_accession', 'gene', 'product', 'sf_category',
           'sf_type', 'start', 'end','segments', 'Source', 'source_id',
           'source_strain', 'source_aa_sequence', 'source_sf_location','Evidence_code',
           'additional_metadata', 'comments']

SF_Table_Uniprot = pd.DataFrame(columns=columns)
#print(SF_Table_Uniprot)
#print(merged_df.columns)

# Ensure SF_Table_Uniprot has the same number of rows as merged_df
SF_Table_Uniprot = SF_Table_Uniprot.reindex(range(len(merged_df)))

# Assign values to the new DataFrame from the merged DataFrame
SF_Table_Uniprot['sf_type'] = merged_df['type'].values
SF_Table_Uniprot['Source'] = "Uniprot"
SF_Table_Uniprot['source_id'] = merged_df['uniprot_id'].values
SF_Table_Uniprot['gene'] = merged_df['Gene Names'].values
SF_Table_Uniprot['product'] = merged_df['Protein names'].values
SF_Table_Uniprot['virus_type'] = merged_df['Organism'].values
SF_Table_Uniprot['Evidence_code'] = "EXP"

#calculate the length of the feature
SF_Table_Uniprot['length'] = (merged_df['location.end.value'] - merged_df['location.start.value'])+1

#fill the source_sf_location column with start and end positions
SF_Table_Uniprot['source_sf_location'] = merged_df['location.start.value'].astype(str) + '-' + merged_df['location.end.value'].astype(str)

#call function to parse the evidence column and get the publication information
evidence_df = parse_evidence_column(merged_df)
SF_Table_Uniprot['publication'] = evidence_df['Publication'].values
#print(evidence_df)

#call function for SF_category
SF_Table_Uniprot['sf_category'] = parse_sf_category(merged_df)

#call function to extract comments
SF_Table_Uniprot['comments'] = extract_comments(merged_df)


#write SF_Table_Uniprot to a csv file
SF_Table_Uniprot_file_path = os.path.join(main_folder, "Chikungunya_Uniprot_SF.csv")
SF_Table_Uniprot.to_csv(SF_Table_Uniprot_file_path, index=False)
print(f"Measles_Uniprot.csv saved to {SF_Table_Uniprot_file_path}.")


