# #parse the evidence column from merged_feature file and obtain the appropriate fields
# #input example for each row: [{'evidenceCode': 'ECO:0000255', 'source': 'PROSITE-ProRule', 'id': 'PRU01344'}, {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '32358203'}, {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '32526208'}, {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '33232691'}, {'evidenceCode': 'ECO:0007744', 'source': 'PDB', 'id': '7CYQ'}]
# #output : PubMed: 32358203, 32526208, 33232691 and if there is no pubmed then it will be empty

import pandas as pd
import ast

def parse_evidence_column(df):
    results = []
    for idx, row in enumerate(df['evidences']):
        if pd.isna(row):
            print(f"NaN encountered at row index: {idx}")
            results.append({})
            continue
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
    return results_df['Publication']

# Read the input DataFrame
df = pd.read_csv('/Users/rbhattac/Desktop/SFVT_pipeline/merged_features.csv')
# Parse the evidence column
parsed_df = parse_evidence_column(df)
print(parsed_df)
# Write the parsed DataFrame to a CSV file
#parsed_df.to_csv('/Users/rbhattac/Desktop/evidences_parsed.csv', index=False)