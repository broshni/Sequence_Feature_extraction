# read an Excel export obtained by doing a search on IEDB for a virus
import argparse
import os
import pandas as pd

from sf_columns import write_csv
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
#Linear peptide rows with a missing/non-numeric Starting or Ending Position are left blank
#(instead of crashing) and reported so the underlying data gap can be investigated separately.
def extract_linear_info(df):
    required_columns = {'Object Type', 'Starting Position', 'Ending Position', 'Name'}
    if not required_columns.issubset(df.columns):
        raise ValueError("Required columns not found in the DataFrame.")

    missing_position_rows = []

    def _row_info(x):
        if x['Object Type'] != 'Linear peptide':
            return pd.Series({'source_sf_location': None, 'length': None, 'source_aa_sequence': None})

        start = pd.to_numeric(x['Starting Position'], errors='coerce')
        end = pd.to_numeric(x['Ending Position'], errors='coerce')

        if pd.isna(start) or pd.isna(end):
            missing_position_rows.append(str(x.get('Name', '')))
            return pd.Series({
                'source_sf_location': None,
                'length': None,
                'source_aa_sequence': x['Name'],
            })

        start_i, end_i = int(start), int(end)
        return pd.Series({
            'source_sf_location': f"{start_i}-{end_i}",
            'length': end_i - start_i + 1,
            'source_aa_sequence': x['Name'],
        })

    result = df.apply(_row_info, axis=1)

    if missing_position_rows:
        sample = ', '.join(missing_position_rows[:5])
        more = '' if len(missing_position_rows) <= 5 else f" (+{len(missing_position_rows) - 5} more)"
        print(
            f"Warning: {len(missing_position_rows)} 'Linear peptide' row(s) missing/non-numeric "
            f"Starting or Ending Position - source_sf_location and length left blank for these. "
            f"Examples: {sample}{more}"
        )

    return result['source_sf_location'], result['length'], result['source_aa_sequence']


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
        # Use a nullable integer dtype so whole-number lengths render as "15" rather than
        # "15.0" (plain int dtype can't hold NaN for any row that ends up with no length).
        parsed_df['length'] = parsed_df['length'].astype('Int64')
        return parsed_df
    else:
        raise ValueError("Column 'Object Type' not found in the DataFrame.")


#function to drop epitopes whose Object Type is "Discontinuous peptide on multi chain" -
#these don't have a single-chain location/sequence that this pipeline can represent.
#Matching is fuzzy (case-insensitive, ignores spaces/hyphens/underscores) so variants like
#"Discontinuous peptide on multichain", "Multi-Chain", extra whitespace, etc. are all caught.
def filter_multichain_discontinuous(df):
    if 'Object Type' not in df.columns:
        raise ValueError("Column 'Object Type' not found in the DataFrame.")

    def _normalize(value):
        return re.sub(r'[^a-z0-9]', '', str(value).lower())

    target = _normalize('Discontinuous peptide on multi chain')

    total_count = len(df)
    is_multichain = df['Object Type'].apply(lambda x: _normalize(x) == target)
    filtered_df = df[~is_multichain].reset_index(drop=True)

    print(
        f"Removed {int(is_multichain.sum())} 'Discontinuous peptide on multi chain' entries "
        f"out of {total_count} total rows; {len(filtered_df)} rows remain."
    )
    return filtered_df


#function to drop 'Linear peptide' epitopes with a missing/non-numeric Starting or Ending
#Position - these can't be assigned a real source_sf_location/length, so rather than leaving
#blank values in the final table they're excluded and returned separately for review.
def filter_missing_linear_positions(df):
    required_columns = {'Object Type', 'Starting Position', 'Ending Position'}
    if not required_columns.issubset(df.columns):
        raise ValueError("Required columns not found in the DataFrame.")

    start = pd.to_numeric(df['Starting Position'], errors='coerce')
    end = pd.to_numeric(df['Ending Position'], errors='coerce')
    is_linear = df['Object Type'] == 'Linear peptide'
    is_missing_position = is_linear & (start.isna() | end.isna())

    total_count = len(df)
    dropped_df = df.loc[is_missing_position].reset_index(drop=True)
    filtered_df = df.loc[~is_missing_position].reset_index(drop=True)

    if len(dropped_df):
        names = dropped_df['Name'].astype(str).tolist() if 'Name' in df.columns else []
        sample = ', '.join(names[:5])
        more = '' if len(names) <= 5 else f" (+{len(names) - 5} more)"
        print(
            f"Removed {len(dropped_df)} 'Linear peptide' row(s) missing/non-numeric Starting or "
            f"Ending Position out of {total_count} total rows; {len(filtered_df)} rows remain. "
            f"Examples: {sample}{more}"
        )

    return filtered_df, dropped_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a standardized Sequence Feature (SF) table from an IEDB epitope "
                    "search export (.xlsx)."
    )
    parser.add_argument(
        "xlsx_file",
        help="Path to the IEDB epitope export Excel file (.xlsx).",
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Directory where the output CSV will be written (created if missing).",
    )
    parser.add_argument(
        "-n", "--name",
        default="output",
        help="Prefix used to name the output SF table, e.g. '<name>_Epitope_SF.csv' "
             "(default: 'output').",
    )
    parser.add_argument(
        "--sheet",
        default=0,
        help="Sheet name or 0-based sheet index to read from the Excel file "
             "(default: 0, the first sheet).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Sheet can be given as a name (string) or a 0-based index (int) - try to interpret
    # numeric-looking input as an index, otherwise treat it as a sheet name.
    sheet = args.sheet
    if isinstance(sheet, str) and sheet.isdigit():
        sheet = int(sheet)

    #read the Excel file containing epitope data (IEDB exports come as .xlsx with multiple
    #sheets; the epitope table itself is on the first sheet by default).
    df_temp = pd.read_excel(args.xlsx_file, sheet_name=sheet)

    # Remove columns containing "Related Object" in their header
    df_epitope = df_temp.loc[:, ~df_temp.columns.str.contains("Related Object")]

    #remove the first row and make the second row the header
    df_epitope.columns = df_epitope.iloc[0]
    df_epitope = df_epitope[1:]
    #reset the index
    df_epitope.reset_index(drop=True, inplace=True)

    print(df_epitope.columns)

    #drop "Discontinuous peptide on multi chain" entries before they proceed downstream
    df_epitope = filter_multichain_discontinuous(df_epitope)

    #drop 'Linear peptide' entries with a missing/non-numeric Starting or Ending Position;
    #keep the dropped rows so they can be reviewed instead of silently disappearing
    df_epitope, dropped_missing_position = filter_missing_linear_positions(df_epitope)
    if len(dropped_missing_position):
        dropped_file_path = os.path.join(args.outdir, f"{args.name}_Epitope_dropped_missing_position.csv")
        write_csv(dropped_missing_position, dropped_file_path)
        print(f"Rows dropped for missing Starting/Ending Position saved to {dropped_file_path}.")

    #create an empty dataframe for the new SF table with the features extracted from UniProt
    columns = ['virus_type', 'subtype', 'segment', 'genbank_accession', 'gene', 'product', 'sf_category',
               'sf_type', 'start', 'end', 'segments', 'Source', 'source_id',
               'source_strain', 'source_aa_sequence', 'source_sf_location', 'Evidence_code',
               'additional_metadata', 'comments']

    sf_table_epitope = pd.DataFrame(columns=columns)

    # Ensure sf_table_epitope has the same number of rows as df_epitope
    sf_table_epitope = sf_table_epitope.reindex(range(len(df_epitope)))

    # Assign values to the new DataFrame from the df_epitope DataFrame
    sf_table_epitope['sf_type'] = "epitope"
    sf_table_epitope['virus_type'] = df_epitope['Species'].values
    sf_table_epitope['product'] = df_epitope['Source Molecule'].values
    sf_table_epitope['sf_category'] = 'Epitope'
    sf_table_epitope['Source'] = "IEDB"
    sf_table_epitope['Evidence_code'] = "EXP"
    sf_table_epitope['comments'] = "Epitope, " + df_epitope['Object Type'].fillna('').astype(str)

    #call function to obtain the IEDB ID for each epitope
    sf_table_epitope['source_id'] = extract_iedb_id(df_epitope)

    # call function to obtain the source accession for each epitope
    sf_table_epitope['additional_metadata'] = extract_source_accession(df_epitope)

    #call function to parse the epitope information and obtain source_sf_location, length and source_aa_sequence
    parsed_epitope_info = parse_epitope_info(df_epitope)
    print(parsed_epitope_info)
    sf_table_epitope['source_sf_location'] = parsed_epitope_info['source_sf_location']
    sf_table_epitope['length'] = parsed_epitope_info['length']
    sf_table_epitope['source_aa_sequence'] = parsed_epitope_info['source_aa_sequence']

    #write sf_table_epitope to a csv file
    sf_table_epitope_file_path = os.path.join(args.outdir, f"{args.name}_Epitope_SF.csv")
    write_csv(sf_table_epitope, sf_table_epitope_file_path)

    print(f"{args.name}_Epitope_SF.csv saved to {sf_table_epitope_file_path}.")


if __name__ == "__main__":
    main()

