#extract source accession and epitope information from a CSV file containing epitope data

import pandas as pd

#read the CSV file containing epitope data
df_temp = pd.read_csv("/Users/rbhattac/Desktop/SFVT_pipeline/Sars_Epitope.csv")

# Remove columns containing "Related Object" in their header
df_epitope = df_temp.loc[:, ~df_temp.columns.str.contains("Related Object")]
#print(SF_Table_epitope.columns)

#remove the first row and make the second row the header
df_epitope.columns = df_epitope.iloc[0]
df_epitope = df_epitope[1:]
#reset the index
df_epitope.reset_index(drop=True, inplace=True)

# Function to extract source accession from "Source Molecule IRI" column. Split the string by '/' and
# and consider the last two parts. if the second last part is "uniprot" then return in format uniprot:last part value.
# if the second last part is "protein" then return in format NCBI:last part value.
#return a column with the source accession for each epitope.

def extract_source_accession(df):
    if 'Source Molecule IRI' in df.columns:
        return df['Source Molecule IRI'].apply(
            lambda x: f'uniprot:{x.split("/")[-1]}' if 'uniprot' in str(x)
            else f'NCBI:{x.split("/")[-1]}' if 'protein' in str(x)
            else None
        )
    else:
        raise ValueError("Column 'Source Molecule IRI' not found in the DataFrame.")

df_epitope['accession'] = extract_source_accession(df_epitope)
print(df_epitope['accession'])


