# code to derive and populate the SF_Category from merged_df
import pandas as pd

#read the merged features file
df = pd.read_csv('/Users/rbhattac/Desktop/SFVT_pipeline/merged_features.csv')

# Function to parse the 'sf_category' column and extract relevant information
#input is merged_df
#the function will check the type column in merged_df and return a column with the SF_Category.
#if the type is 'Beta strand','Helix', 'turn', 'Coiled coil' then SF_category will be "structural"
#if the type is 'Mutagenesis' then SF_category will be "Sequence_Alterations"
#for all other types, the SF_category will be "Functional"

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

# Apply the function to create a new column 'sf_category'
df['sf_category'] = parse_sf_category(df)
print(df[['type', 'sf_category']])
