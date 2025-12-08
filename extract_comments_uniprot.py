#code to return comments column from the merged_df
import pandas as pd

#read the merged features file
df= pd.read_csv('/Users/rbhattac/Desktop/SFVT_pipeline/merged_features.csv')

# Function to extract comments from the merged DataFrame
# input is merged_df
# the function will check the 'description', 'alternativeSequence.originalSequence', 'alternativeSequence.alternativeSequences
# and 'ligand.name' columns in merged_df and combine them in comma separated format to form the
# comments column.
#the comments column will be a string with all the values combined. Format will be: value in the "description" column, value in the "alternativeSequence.originalSequence" column,
# value in the "alternativeSequence.alternativeSequences" column, value in the "ligand.name" column.
# if any of these columns are NaN, it will be skipped.
#if all the columns are NaN, then the comments column will the value in the 'type' column.

# def extract_comments(merged_df):
#     comments = []
#     for index, row in merged_df.iterrows():
#         comment_parts = []
#         if pd.notna(row['description']):
#             comment_parts.append(row['description'])
#         if pd.notna(row['alternativeSequence.originalSequence']):
#             comment_parts.append('original sequence:' + str(row['alternativeSequence.originalSequence']))
#         if pd.notna(row['alternativeSequence.alternativeSequences']):
#             comment_parts.append('variant sequence:' + str(row['alternativeSequence.alternativeSequences']))
#         if pd.notna(row['ligand.name']):
#             comment_parts.append('Binding site for:' + str(row['ligand.name']))
#
#         if not comment_parts:
#             comment_parts.append(row['type'])
#
#         if row['type'] == 'Natural variant':
#             comment_parts.insert(0, 'Natural variation observed')
#
#         comments.append(', '.join(comment_parts))
#
#     return pd.Series(comments)

def extract_comments(merged_df):
    comments = []
    for index, row in merged_df.iterrows():
        comment_parts = []
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

        if row['type'] == 'Natural variant':
            comment_parts.insert(0, 'Natural variation observed')

        if not comment_parts:
            comments.append(str(row['type']))
        else:
            comments.append(', '.join(comment_parts))

    return pd.Series(comments)

# Apply the function to create a new column 'comments'
df['comments'] = extract_comments(df)
print(df[['type', 'comments']])

# Save the updated DataFrame with comments to a new CSV file
df.to_csv('/Users/rbhattac/Desktop/SFVT_pipeline/merged_features_with_comments.csv', index=False)

