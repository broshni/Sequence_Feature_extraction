import argparse
import pandas as pd
import requests
import time
import os
import ast

# shared so step 1's 'length' and step 4's coordinate mapping agree on which types are pairs
from sf_columns import PAIR_SF_TYPES, read_csv, write_csv

#function to filter a UniProt TSV export down to only reviewed (Swiss-Prot) entries,
#based on the 'Reviewed' column populated by UniProt taxonomic-search exports.
def filter_reviewed_entries(uniprot_info):
    if 'Reviewed' not in uniprot_info.columns:
        raise ValueError(
            "Expected a 'Reviewed' column in the input TSV to filter on, but it was not found. "
            f"Available columns: {list(uniprot_info.columns)}"
        )

    total_count = len(uniprot_info)
    is_reviewed = uniprot_info['Reviewed'].astype(str).str.strip().str.lower() == 'reviewed'
    reviewed_df = uniprot_info[is_reviewed].reset_index(drop=True)

    print(
        f"Filtered to reviewed entries: {len(reviewed_df)} of {total_count} entries kept "
        f"({total_count - len(reviewed_df)} unreviewed entries excluded)."
    )
    return reviewed_df

def get_json_with_retries(url, retries=4, timeout=60):
    """GET a UniProt JSON document, retrying rate limits and server errors with a growing wait.
    Without this a single transient failure silently drops a protein's features."""
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=timeout)
            if response.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {response.status_code}")
            response.raise_for_status()
            return response.json()
        except Exception as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(2 ** attempt)
    raise last_error


#function to extract UniProt features from a TSV file obtained from taxonomic search of Uniprot
# and output a csv file with the all extracted features from multiple UniProt IDs.
#tsv_file_path: path to the TSV file containing UniProt IDs from taxonomic search

def extract_and_merge_uniprot_features(tsv_file_path):
    # Read TSV file for UniProt IDs
    uniprot_info = read_csv(tsv_file_path, delimiter="\t")

    # Keep only reviewed (Swiss-Prot) entries; unreviewed entries do not proceed downstream
    uniprot_info = filter_reviewed_entries(uniprot_info)
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
    failed = []

    for uniprot_id in uniprot_ids:
        query_url = base_url.format(uniprot_id) + "?fields=" + fields
        try:
            data = get_json_with_retries(query_url)
            if 'features' in data:
                dict_df = pd.json_normalize(data['features'])
                dict_df['uniprot_id'] = uniprot_id
                print(f"Data for {uniprot_id} retrieved successfully.")
                global_df = pd.concat([global_df, dict_df], ignore_index=True)
            else:
                print(f"{uniprot_id} has no annotated features.")
        except Exception as e:
            print(f"An error occurred for {uniprot_id}: {e}")
            failed.append(uniprot_id)
        time.sleep(1)

    if failed:
        # Carrying on would hand every later step a table silently missing these proteins.
        raise SystemExit(
            f"{len(failed)} of {len(uniprot_ids)} UniProt entries could not be retrieved even "
            f"after retries, so the feature table would be incomplete: {', '.join(failed)}. "
            "Re-run to try them again."
        )
    if global_df.empty:
        raise SystemExit(
            "None of the reviewed UniProt entries has annotated features, so there is "
            "nothing to build a table from."
        )

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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a standardized Sequence Feature (SF) table from a UniProt "
                    "taxonomic-search TSV export."
    )
    parser.add_argument(
        "tsv_file",
        help="Path to the input TSV file of UniProt IDs (from a UniProt taxonomic search).",
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Directory where output CSV files will be written (created if missing).",
    )
    parser.add_argument(
        "-n", "--name",
        default="output",
        help="Prefix used to name the output SF table, e.g. '<name>_Uniprot_SF.csv' "
             "(default: 'output').",
    )
    parser.add_argument(
        "--merged-filename",
        default="merged_features.csv",
        help="Filename for the intermediate merged UniProt features CSV "
             "(default: 'merged_features.csv').",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    merged_df = extract_and_merge_uniprot_features(args.tsv_file)

    # Save the merged DataFrame to a CSV file
    merged_dataframe_file_path = os.path.join(args.outdir, args.merged_filename)
    write_csv(merged_df, merged_dataframe_file_path)
    print(f"All features saved to {merged_dataframe_file_path}.")

    #create an empty dataframe for the new SF table with the features extracted from UniProt
    columns = ['virus_type', 'subtype', 'segment', 'genbank_accession', 'gene', 'product', 'sf_category',
               'sf_type', 'start', 'end', 'segments', 'Source', 'source_id',
               'source_strain', 'source_aa_sequence', 'source_sf_location', 'Evidence_code',
               'additional_metadata', 'comments']

    sf_table_uniprot = pd.DataFrame(columns=columns)

    # Ensure sf_table_uniprot has the same number of rows as merged_df
    sf_table_uniprot = sf_table_uniprot.reindex(range(len(merged_df)))

    # Assign values to the new DataFrame from the merged DataFrame
    sf_table_uniprot['sf_type'] = merged_df['type'].values
    sf_table_uniprot['Source'] = "Uniprot"
    sf_table_uniprot['source_id'] = merged_df['uniprot_id'].values
    sf_table_uniprot['gene'] = merged_df['Gene Names'].values
    sf_table_uniprot['product'] = merged_df['Protein names'].values
    sf_table_uniprot['virus_type'] = merged_df['Organism'].values
    sf_table_uniprot['Evidence_code'] = "EXP"

    #calculate the length of the feature. For a disulfide bond or cross-link, 'start-end' is
    #the pair of residues that are bonded rather than a span, so the feature is the 2 residues
    #themselves (1 if UniProt reported the same position twice) - not end-start+1, which would
    #call a bond between Cys188 and Cys606 a 419-residue feature.
    span_length = (merged_df['location.end.value'] - merged_df['location.start.value']) + 1
    is_pair = merged_df['type'].isin(PAIR_SF_TYPES)
    pair_length = (merged_df['location.end.value'] != merged_df['location.start.value']).map({True: 2, False: 1})
    sf_table_uniprot['length'] = span_length.where(~is_pair, pair_length)

    #fill the source_sf_location column with start and end positions
    sf_table_uniprot['source_sf_location'] = (
        merged_df['location.start.value'].astype(str) + '-' + merged_df['location.end.value'].astype(str)
    )

    #call function to parse the evidence column and get the publication information
    evidence_df = parse_evidence_column(merged_df)
    sf_table_uniprot['publication'] = evidence_df['Publication'].values

    #call function for SF_category
    sf_table_uniprot['sf_category'] = parse_sf_category(merged_df)

    #call function to extract comments
    sf_table_uniprot['comments'] = extract_comments(merged_df)

    #write sf_table_uniprot to a csv file
    sf_table_uniprot_file_path = os.path.join(args.outdir, f"{args.name}_Uniprot_SF.csv")
    write_csv(sf_table_uniprot, sf_table_uniprot_file_path)

    print(f"{args.name}_Uniprot_SF.csv saved to {sf_table_uniprot_file_path}.")


if __name__ == "__main__":
    main()

