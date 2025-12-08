import pandas as pd
import json
import requests
import os
import time

main_folder = "/Users/rbhattac/Desktop/Chikungunya_SFVT"

# Read a TSV file as DataFrame for UniProt IDs
uniport_file_path = os.path.join(main_folder, "Chikungunya_uniprot.tsv")
uniprot_info = pd.read_csv(uniport_file_path, delimiter="\t")
uniprot_ids = uniprot_info['Entry']

print(uniprot_ids)

# Base URL for UniProt API
base_url = "https://rest.uniprot.org/uniprotkb/{}.json"

# Fields to retrieve
fields = (
    "ft_var_seq,ft_variant,ft_non_cons,ft_non_std,ft_non_ter,ft_conflict,"
    "ft_unsure,ft_act_site,ft_binding,ft_dna_bind,ft_site,ft_mutagen,"
    "ft_intramem,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,"
    "ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,ft_propep,"
    "ft_signal,ft_transit,ft_strand,ft_helix,ft_turn,ft_coiled,ft_compbias,"
    "ft_domain,ft_motif,ft_region,ft_repeat,ft_zn_fing"
)

# Initialize an empty global DataFrame
global_df = pd.DataFrame()

# Loop through each UniProt ID
for uniprot_id in uniprot_ids:
    print(f"Processing UniProt ID: {uniprot_id}")

    query_url = base_url.format(uniprot_id) + "?fields=" + fields

    try:
        response = requests.get(query_url)
        response.raise_for_status()  # Raise an error for HTTP errors

        # Read the JSON content
        data = response.json()

        # Normalize JSON data to a DataFrame
        if 'features' in data:
            dict_df = pd.json_normalize(data['features'])
            dict_df['uniprot_id'] = uniprot_id  # Add the UniProt ID column
            print(f"Data for {uniprot_id} retrieved successfully.")

            # Append the DataFrame to the global DataFrame
            global_df = pd.concat([global_df, dict_df], ignore_index=True)

        else:
            print(f"No 'features' found in response for {uniprot_id}.")

    except Exception as e:
        print(f"An error occurred for {uniprot_id}: {e}")

    # Sleep to avoid hitting API limits
    time.sleep(1)

# Merge global_df and uniprot_info on 'uniprot_id' and 'Entry'
merged_df = pd.merge(global_df, uniprot_info, left_on='uniprot_id', right_on='Entry')

# Save the merged DataFrame to a CSV file
merged_dataframe_file_path = os.path.join(main_folder, "merged_features.csv")
merged_df.to_csv(merged_dataframe_file_path, index=False)
print(f"All features saved to {merged_dataframe_file_path}.")
print("Processing completed.")