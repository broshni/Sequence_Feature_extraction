import pandas as pd

# # Read the CSV file into a DataFrame
# file_path = '/Users/rbhattac/Desktop/SFVT/SFVT_data/Dengue/Dengue1_SFTable.csv'
# df = pd.read_csv(file_path)
#
# # Create a new column that combines 'gene', 'segments', and 'sf_category' values
# df['gene_segments_category'] = df['gene'].astype(str) + '_' + df['segments'].astype(str) + '_' + df[
#     'sf_category'].astype(str)
# df['source:id'] = df['Source'].astype(str) + ':' + df['source_id'].astype(str)
#
# # Group by 'gene_segments_category' and process matching lines for each unique value
# result = []
# for gene_segment_category, group in df.groupby('gene_segments_category'):
#     if len(group) > 1:
#         print(f"Gene_Segment_category: {gene_segment_category}")
#         print("Matching line numbers:", group.index.tolist())
#
#         # Keep the first row in main dataframe "df" and copy the 'source:id' values of the rest to the first row
#         first_row = group.iloc[0].copy()
#         first_row['source:id'] = ','.join(group['source:id'].iloc[1:].tolist())
#
#         # Append the modified first row to the result list
#         result.append(first_row)
#
# # Create a new DataFrame from the result list
# result_df = pd.DataFrame(result)
#
# # Print the result DataFrame
# print(result_df)

import pandas as pd
#import chardet

# Read the CSV file into a DataFrame
#file_path = '/Users/rbhattac/Desktop/SARS-CoV2_SF_BV-BRC_edit.csv'
#df = pd.read_csv(file_path)
#read the csv file as a dataframe
#with open('/Users/rbhattac/Desktop/SARS-CoV2_SF_BV-BRC_edit.csv', 'rb') as file:
   # result = chardet.detect(file.read())

# file_encoding = result['encoding']
# print(f"Detected encoding: {file_encoding}")

# read the file
df = pd.read_csv('/Users/rbhattac/Desktop/Chikungunya_SFVT/Chikungunya_SFTable_ALL.csv')
print(df.columns)

# Create a new column that combines 'gene', 'segments', and 'sf_category' values
df['gene_segments_category'] = df['gene'].astype(str) + '_' + df['segments'].astype(str) + '_' + df[
    'sf_category'].astype(str)
df['source:id'] = df['Source'].astype(str) + ':' + df['source_id'].astype(str)

# Group by 'gene_segments_category' and process matching lines for each unique value
for gene_segment_category, group in df.groupby('gene_segments_category'):
    if len(group) > 1:
        # Keep the first row and copy the 'source:id' values of the rest to the first row
        first_row_index = group.index[0]
        #remove the source:id of the first row from group
        remaining_source_ids = group['source:id'].iloc[1:].tolist()
        #df.at[first_row_index, 'source:id'] = ','.join(group['source:id'].tolist())
        #df.at[first_row_index, 'Publication'] = ','.join(group['source:id'].tolist())
        df.at[first_row_index, 'Publication'] = ','.join(remaining_source_ids)

        # Drop the remaining rows
        df = df.drop(group.index[1:])

# Reset the index of the DataFrame
df = df.reset_index(drop=True)

# Print the updated DataFrame
print(df)

#writing to a csv file
df.to_csv('/Users/rbhattac/Desktop/Chikungunya_SFVT/CHIKV_SFTable_editduplicates.csv', index=False)