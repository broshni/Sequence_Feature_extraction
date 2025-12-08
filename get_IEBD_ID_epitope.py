#code to extract the IEBD ID for each epitope

import pandas as pd

df_epitope = pd.read_csv("/Users/rbhattac/Desktop/SFVT_pipeline/Sars_Epitope.csv")

df_epitope.columns = df_epitope.iloc[0]
df_epitope = df_epitope[1:]
#reset the index
df_epitope.reset_index(drop=True, inplace=True)

#read the "IEDB IRI" column from the DataFrame. split the string by '/' and take the last part as the IEBD ID.
# return a column with the IEBD ID for each epitope

def extract_iedb_id(df):
    if 'IEDB IRI' in df.columns:
        iedb_ids = df['IEDB IRI'].apply(lambda x: x.split('/')[-1] if isinstance(x, str) else None)
        return iedb_ids
    else:
        raise ValueError("Column 'IEDB IRI' not found in the DataFrame.")

iedb_ids = extract_iedb_id(df_epitope)
print(iedb_ids)




