#read the SF table and make a new table with SF_names and SF_ids
import pandas as pd

#read the SF table
SF_table = pd.read_csv('/Users/rbhattac/Desktop/Chikungunya_SFVT/CHIKV_SFTable_editduplicates.csv')
#remove sf_id column
#SF_table = SF_table.drop(columns=['sf_id'])
#print(SF_table)

#create sf_id in the format: virus_type_gene_SF followed by a unique number. The unique number should be grouped by the
#genes and virus_type. For example, if the virus_type is SARS-CoV-2 and the gene is S, then the SF_ids should be
#SARS-CoV-2_S_SF1, SARS-CoV-2_S_SF2, etc. for each SF in the gene S. Repeat for all genes and virus_types in the table.

def create_sf_ids(df):
    df['sf_id'] = df.groupby(['virus_type', 'gene']).cumcount() + 1
    df['sf_id'] = df.apply(lambda x: f"{x['virus_type']}_{x['gene']}_SF{x['sf_id']}", axis=1)
    df['sf_id'] = df['sf_id'].str.replace(r'\.0$', '', regex=True)
    return df

# create sf_names in the format: virus_type_gene_sf_type_start(length).

def create_sf_names(df):
    df['sf_name'] = df.apply(
        lambda x: f"{x['virus_type']}_{x['gene']}_{x['sf_type']}_{x['start']}({x['length']})",
        axis=1
    )
    return df

#function to remove .0 from the sf_name example change "Measles Virus_N_Compositional bias_433.0(20)" to "Measles Virus_N_Compositional bias_433(20)"
def remove_decimal_from_sf_name(df):
    df['sf_name'] = df['sf_name'].str.replace(r'\.0\(', '(', regex=True)
    return df

SF_table = create_sf_ids(SF_table)
# Create a new DataFrame with SF_names and SF_ids
#SF_names_ids = SF_table[['SF_name', 'SF_ids']].copy()
#print(SF_table)

SF_table = create_sf_names(SF_table)
SF_table = remove_decimal_from_sf_name(SF_table)
print(SF_table)

#Save the new DataFrame to a CSV file
SF_table.to_csv('/Users/rbhattac/Desktop/Chikungunya_SFVT/CHIKV_SFTable_Final_20Nov.csv', index=False)