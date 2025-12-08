# read a csv file obtained by doing a search on IEDB for a virus
import pandas as pd
import re

#read the "IEDB IRI" column from the DataFrame. split the string by '/' and take the last part as the IEBD ID.
# return a column with the IEBD ID for each epitope

def extract_iedb_id(df):
    if 'IEDB IRI' in df.columns:
        iedb_ids = df['IEDB IRI'].apply(lambda x: x.split('/')[-1] if isinstance(x, str) else None)
        return iedb_ids
    else:
        raise ValueError("Column 'IEDB IRI' not found in the DataFrame.")

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

#function to extract source_sf_location, length and source_aa_sequence for the linear peptide epitope.
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


# Function to extract source_sf_location, length and source_aa_sequence for the Discontinuous peptide epitope.

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


#function to parse a df which has both linear and discontinuous peptides and give a DataFrame with the source_sf_location, length and source_aa_sequence for both types of peptides.
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


#read the CSV file containing epitope data
df_temp = pd.read_csv("/Users/rbhattac/Desktop/Chikungunya_SFVT/Chikungunya_IEDB.csv")

# Remove columns containing "Related Object" in their header
df_epitope = df_temp.loc[:, ~df_temp.columns.str.contains("Related Object")]
#print(SF_Table_epitope.columns)

#remove the first row and make the second row the header
df_epitope.columns = df_epitope.iloc[0]
df_epitope = df_epitope[1:]
#reset the index
df_epitope.reset_index(drop=True, inplace=True)

print(df_epitope.columns)

#create an empty dataframe for the new SF table with the features extracted from UniProt
columns = ['virus_type', 'subtype', 'segment', 'genbank_accession', 'gene', 'product', 'sf_category',
           'sf_type', 'start', 'end','segments', 'Source', 'source_id',
           'source_strain', 'source_aa_sequence', 'source_sf_location','Evidence_code',
           'additional_metadata', 'comments']

SF_Table_epitope = pd.DataFrame(columns=columns)

# Ensure SF_Table_epitope has the same number of rows as df_epitope
SF_Table_epitope = SF_Table_epitope.reindex(range(len(df_epitope)))

# Assign values to the new DataFrame from the df_epitope DataFrame
SF_Table_epitope['sf_type'] = "epitope"
SF_Table_epitope['virus_type'] = df_epitope['Species'].values
SF_Table_epitope['product'] = df_epitope['Source Molecule'].values
SF_Table_epitope['sf_category'] = 'Epitope'
SF_Table_epitope['Source'] = "IEDB"
SF_Table_epitope['evidence_code'] = "EXP"
SF_Table_epitope['comments'] = "Epitope, " + df_epitope['Object Type'].fillna('').astype(str)

#call function to obtain the IEDB ID for each epitope
SF_Table_epitope['source_id'] = extract_iedb_id(df_epitope)

# call function to obtain the source accession for each epitope
SF_Table_epitope['additional_metadata'] = extract_source_accession(df_epitope)

#call function to parse the epitope information and obtain source_sf_location, length and source_aa_sequence
parsed_epitope_info = parse_epitope_info(df_epitope)
print(parsed_epitope_info)
SF_Table_epitope['source_sf_location'] = parsed_epitope_info['source_sf_location']
SF_Table_epitope['length'] = parsed_epitope_info['length']
SF_Table_epitope['source_aa_sequence'] = parsed_epitope_info['source_aa_sequence']


#write SF_Table_epitope to a csv file
SF_Table_epitope.to_csv("/Users/rbhattac/Desktop/Chikungunya_SFVT/Chikungunya_epitope_SF.csv", index=False)



