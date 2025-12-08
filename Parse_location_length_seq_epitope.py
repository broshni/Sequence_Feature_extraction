# extract source_sf_location, length and source_aa_sequence from the epitope file

import pandas as pd
import re

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

# function to extract source_sf_location, length and source_aa_sequence for the linear peptide epitope.
# the function takes a DataFrame as input and returns 3 columns: source_sf_location, length and source_aa_sequence
# if the "Object Type" is "Linear peptide" then combine the "Starting Position" and "Ending Position" columns
# and return the string in the format "Starting Position-Ending Position" for source_sf_location.
# for length, calculate the difference between "Ending Position" and "Starting Position" and add 1.
# for source_aa_sequence, extract the amino acid sequence from the "Name" column.

def extract_linear_info(df):
    if 'Object Type' in df.columns and 'Starting Position' in df.columns and 'Ending Position' in df.columns and 'Name' in df.columns:
        source_sf_location = df.apply(
            lambda x: f"{x['Starting Position']}-{x['Ending Position']}" if x['Object Type'] == 'Linear peptide' else None,
            axis=1
        )
        length = df.apply(
            lambda x: int(x['Ending Position']) - int(x['Starting Position']) + 1 if x['Object Type'] == 'Linear peptide' else None,
            axis=1
        )
        source_aa_sequence = df.apply(
            lambda x: x['Name'] if x['Object Type'] == 'Linear peptide' else None,
            axis=1
        )
        return source_sf_location, length, source_aa_sequence
    else:
        raise ValueError("Required columns not found in the DataFrame.")

# # Call the function to extract the information
# source_sf_location, length, source_aa_sequence = extract_linear_info(df_epitope)
# print(source_sf_location, length, source_aa_sequence)

# Function to extract source_sf_location, length and source_aa_sequence for the Discontinuous peptide epitope.
# if the "Object Type" is "Discontinuous peptide" the read the "Name" column and parse it to separate
# the amino acid sequence and the positions.for example, if the "Name" column is "G104, E126, W231" then
# the source_sf_location should be "104, 126, 231", length should be 3 and source_aa_sequence should be "G, E, W".

def extract_discontinuous_info(df):
    if 'Object Type' in df.columns and 'Name' in df.columns:
        def parse_location(name):
            locations = []
            for part in str(name).split(','):
                part = part.strip()
                nums = re.findall(r'\d+', part)
                locations.extend(nums)
            return ', '.join(locations) if locations else None

        def parse_length(name):
            return len([p for p in str(name).split(',') if p.strip()])

        def parse_sequence(name):
            aas = []
            for part in str(name).split(','):
                part = part.strip()
                match = re.match(r'^([A-Za-z])', part)
                if match:
                    aas.append(match.group(1))
            return ', '.join(aas) if aas else None

        source_sf_location = df.apply(
            lambda x: parse_location(x['Name']) if x['Object Type'] == 'Discontinuous peptide' else None,
            axis=1
        )
        length = df.apply(
            lambda x: parse_length(x['Name']) if x['Object Type'] == 'Discontinuous peptide' else None,
            axis=1
        )
        source_aa_sequence = df.apply(
            lambda x: parse_sequence(x['Name']) if x['Object Type'] == 'Discontinuous peptide' else None,
            axis=1
        )
        return source_sf_location, length, source_aa_sequence
    else:
        raise ValueError("Required columns not found in the DataFrame.")

# # Call the function to extract the information for discontinuous peptides
# source_sf_location_discontinuous, length_discontinuous, source_aa_sequence_discontinuous = extract_discontinuous_info(df_epitope)
# print(source_sf_location_discontinuous, length_discontinuous, source_aa_sequence_discontinuous)

#now write a function to parse a df which has both linear and discontinuous peptides and give a DataFrame with the source_sf_location, length and source_aa_sequence for both types of peptides.
def parse_epitope_info(df):
    if 'Object Type' in df.columns:
        linear_sf_location, linear_length, linear_sequence = extract_linear_info(df)
        discontinuous_sf_location, discontinuous_length, discontinuous_sequence = extract_discontinuous_info(df)

        # Combine results into a new DataFrame
        parsed_df = pd.DataFrame({
            'source_sf_location': linear_sf_location.combine_first(discontinuous_sf_location),
            'length': linear_length.combine_first(discontinuous_length),
            'source_aa_sequence': linear_sequence.combine_first(discontinuous_sequence)
        })
        return parsed_df
    else:
        raise ValueError("Column 'Object Type' not found in the DataFrame.")

# Call the function to parse the epitope information
parsed_epitope_info = parse_epitope_info(df_epitope)
# Print the parsed DataFrame
print(parsed_epitope_info)

