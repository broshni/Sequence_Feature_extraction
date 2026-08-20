#parse the SFVT merged SF table and map each row's source-sequence coordinates onto its
#reference (GenBank) sequence via a pairwise alignment. Populates the start/end/segments
#columns already present in the merged table and adds Source_Sequence/Ref_sequence columns.
#Every original column and row is preserved in the output - rows that can't be mapped for
#any reason are left blank in the mapping columns rather than dropped, and processing moves
#on to the next row.
import argparse
import os
import socket
import time
from urllib.error import HTTPError, URLError
from Bio import AlignIO
import pandas as pd
from Bio import Entrez, SeqIO
import subprocess
import requests as r
from io import StringIO

# Transient network errors worth retrying (read timeouts, dropped/reset connections, etc.)
# rather than giving up on a row after a single hiccup.
TRANSIENT_NETWORK_ERRORS = (
    r.exceptions.Timeout,
    r.exceptions.ConnectionError,
    URLError,
    HTTPError,
    socket.timeout,
    ConnectionError,
)


def fetch_uniprot(uni_id, retries=3, retry_delay=5, timeout=30):
    """Fetch a protein sequence from UniProt by accession. Retries on transient network
    errors (timeouts, connection resets) with a backing-off delay before raising."""
    base_url = "http://www.uniprot.org/uniprot/"
    current_url = f"{base_url}{uni_id}.fasta"

    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            response = r.get(current_url, timeout=timeout)
            response.raise_for_status()
            records = list(SeqIO.parse(StringIO(response.text), "fasta"))
            if not records:
                raise ValueError(f"No FASTA record returned for UniProt ID '{uni_id}'")
            return records[0].seq
        except TRANSIENT_NETWORK_ERRORS as e:
            last_exc = e
            if attempt < retries:
                wait = retry_delay * attempt
                print(f"    UniProt fetch for '{uni_id}' failed ({e}); retrying in {wait}s (attempt {attempt}/{retries})...")
                time.sleep(wait)
    raise last_exc


def fetch_genbank_record(accession_number, email="rbhattac@jcvi.org", retries=3, retry_delay=5):
    """Fetch a GenBank protein record by accession. Retries on transient network errors
    (timeouts, connection resets) with a backing-off delay before raising."""
    Entrez.email = email

    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            handle = Entrez.efetch(db="protein", id=accession_number, rettype="gb", retmode="text")
            try:
                return SeqIO.read(handle, "genbank")
            finally:
                handle.close()
        except TRANSIENT_NETWORK_ERRORS as e:
            last_exc = e
            if attempt < retries:
                wait = retry_delay * attempt
                print(f"    GenBank fetch for '{accession_number}' failed ({e}); retrying in {wait}s (attempt {attempt}/{retries})...")
                time.sleep(wait)
    raise last_exc


def align_sequences(ref_seq, otherseq, output_file):
    #mafft_command = f" mafft --6merpair --anysymbol --keeplength --addfragments {otherseq} {ref_seq} > {output_file}"
    mafft_command = f"mafft --6merpair --anysymbol --addfragments {otherseq} {ref_seq} > {output_file}"
    result = subprocess.run(mafft_command, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"mafft alignment failed (exit code {result.returncode})")


def parse_source_coordinates(coordinates):
    if "-" in coordinates:
        coordinate_list = coordinates.split("-")
        coordinate_list = [int(i) for i in coordinate_list]
        coordinate_list.sort()
    elif "," in coordinates:
        coordinate_list = coordinates.split(",")
        coordinate_list = [int(i) for i in coordinate_list]
        coordinate_list.sort()
    else:
        coordinate_list = [int(coordinates)]
    return coordinate_list #return a list of coordinates


def get_ref_pos(alignment, num, coordinates):
    # access the source through index
    source_seq = alignment[num].seq

    pos = 0  # set the starting pos
    ref_position = []
    for index, amino in enumerate(source_seq):
        if amino != "-":
            pos += 1
            if pos in coordinates:
                ref_position.append(index)

    return ref_position


#function to slice the alignment down to just the feature's mapped columns and return the
#(source_seq, ref_seq) pair for those columns. Handles contiguous ('-'), discontinuous (','),
#and single-position coordinate formats.
def extract_sf_alignment(alignment, source_coordinates, new_list):
    if not new_list:
        raise ValueError("No aligned reference positions found for this feature")

    if "-" in source_coordinates:
        sf_alignment = alignment[:, (new_list[0] - 1):new_list[-1]]
    elif "," in source_coordinates:
        # extract each mapped column individually and concatenate them together
        sf_alignment = None
        for i in new_list:
            col = alignment[:, (i - 1):i]
            sf_alignment = col if sf_alignment is None else sf_alignment + col
    else:
        # single coordinate
        sf_alignment = alignment[:, (new_list[0] - 1):new_list[0]]

    return sf_alignment[1].seq, sf_alignment[0].seq


#function to map one SF table row's source_sf_location onto its reference sequence.
#Returns (start, end, segments, source_seq, ref_seq); raises on any failure so the caller
#can decide how to handle it (see main(), which leaves the row blank and moves on).
def map_sf_row(row, ref_fasta_path, source_fasta_path, alignment_output_path,
               email="rbhattac@jcvi.org", retries=3, retry_delay=5, timeout=30):
    source_coordinates = row['source_sf_location']
    if pd.isna(source_coordinates) or not str(source_coordinates).strip():
        raise ValueError("Missing source_sf_location")
    source_coordinates = str(source_coordinates)

    ref_accession = row['genbank_accession']
    if pd.isna(ref_accession) or not str(ref_accession).strip():
        raise ValueError("Missing genbank_accession")

    ref_record = fetch_genbank_record(ref_accession, email=email, retries=retries, retry_delay=retry_delay)
    with open(ref_fasta_path, "w") as f:
        f.write(f">{ref_accession}\n{ref_record.seq}\n")

    if str(row['Source']).lower() == 'uniprot':
        source_id = row['source_id']
        source_seq = fetch_uniprot(source_id, retries=retries, retry_delay=retry_delay, timeout=timeout)
    else:
        metadata = str(row['additional_metadata'])
        parts = metadata.split(":")
        if len(parts) < 2:
            raise ValueError(f"Could not parse additional_metadata '{metadata}'")
        source_type, source_id = parts[0], parts[1]
        if source_type.upper() == 'NCBI':
            source_record = fetch_genbank_record(source_id, email=email, retries=retries, retry_delay=retry_delay)
            source_seq = source_record.seq
        elif source_type.lower() == 'uniprot':
            source_seq = fetch_uniprot(source_id, retries=retries, retry_delay=retry_delay, timeout=timeout)
        else:
            raise ValueError(f"Unrecognized source type '{source_type}' in additional_metadata")

    with open(source_fasta_path, "w") as f:
        f.write(f">{source_id}\n{source_seq}\n")

    align_sequences(ref_fasta_path, source_fasta_path, alignment_output_path)
    alignment = AlignIO.read(alignment_output_path, "fasta")

    coordinates = parse_source_coordinates(source_coordinates)
    ref_pos = get_ref_pos(alignment, 1, coordinates)
    new_list = [x + 1 for x in ref_pos]

    if not new_list:
        raise ValueError("Could not map any source coordinates onto the reference alignment")

    start = new_list[0]
    end = new_list[-1]
    segments = ','.join(map(str, new_list)) if len(new_list) > 2 else f"{start}-{end}"

    source_seq_sf, ref_seq_sf = extract_sf_alignment(alignment, source_coordinates, new_list)

    return start, end, segments, str(source_seq_sf), str(ref_seq_sf)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map each row of a merged SF table onto its reference (GenBank) sequence "
                    "via a pairwise alignment, populating start/end/segments and adding "
                    "Source_Sequence/Ref_sequence columns. Rows that can't be mapped are left "
                    "blank rather than dropped."
    )
    parser.add_argument(
        "input_file",
        help="Path to the merged SF table CSV to map (e.g. output of Edit_SFTable.py).",
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Directory where the mapped CSV (and scratch alignment files) will be written "
             "(created if missing).",
    )
    parser.add_argument(
        "-n", "--name",
        default="output",
        help="Prefix used to name the mapped output file, e.g. '<name>_SF_mapped.csv' "
             "(default: 'output').",
    )
    parser.add_argument(
        "--encoding",
        default="latin-1",
        help="Text encoding to use when reading the input CSV (default: 'latin-1').",
    )
    parser.add_argument(
        "--email",
        default="rbhattac@jcvi.org",
        help="Email address to send with NCBI Entrez requests, as required by NCBI "
             "(default: 'rbhattac@jcvi.org').",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of attempts per network fetch (UniProt/GenBank) before giving up on "
             "that row - helps ride out transient read timeouts/connection resets "
             "(default: 3).",
    )
    parser.add_argument(
        "--retry-delay",
        type=float,
        default=5,
        help="Base delay in seconds between retries of a failed fetch; the wait grows with "
             "each attempt (default: 5).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30,
        help="Timeout in seconds for each UniProt HTTP request (default: 30).",
    )
    parser.add_argument(
        "--row-delay",
        type=float,
        default=1,
        help="Delay in seconds after each row, successful or not, to go easier on the "
             "UniProt/NCBI APIs and reduce the chance of rate-limit-triggered errors "
             "(default: 1).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Scratch files MAFFT/Entrez read and write per row - kept in their own subfolder of
    # outdir so repeated runs don't collide with unrelated files on the Desktop.
    scratch_dir = os.path.join(args.outdir, "_mapping_scratch")
    os.makedirs(scratch_dir, exist_ok=True)
    ref_fasta_path = os.path.join(scratch_dir, "ref.fasta")
    source_fasta_path = os.path.join(scratch_dir, "sequence1.fasta")
    alignment_output_path = os.path.join(scratch_dir, "alignment_output.fasta")

    output_file = os.path.join(args.outdir, f"{args.name}_SF_mapped.csv")

    #read the merged SF table as a dataframe - the output keeps every column from this table
    sf_table = pd.read_csv(args.input_file, encoding=args.encoding)

    # 'start', 'end' and 'segments' already exist (empty) in the merged SF table schema;
    # Source_Sequence/Ref_sequence are new. Force all five to an object dtype so an
    # all-empty column (which pandas would otherwise infer as float64) can safely hold the
    # ints/strings we write into it below.
    for col in ['start', 'end', 'segments', 'Source_Sequence', 'Ref_sequence']:
        if col not in sf_table.columns:
            sf_table[col] = None
        sf_table[col] = sf_table[col].astype('object')

    failed_rows = []
    for index, row in sf_table.iterrows():
        print(f"Processing row {index} (source_id={row.get('source_id', '')})...")
        try:
            start, end, segments, source_seq, ref_seq = map_sf_row(
                row, ref_fasta_path, source_fasta_path, alignment_output_path,
                email=args.email, retries=args.retries, retry_delay=args.retry_delay,
                timeout=args.timeout,
            )
        except Exception as e:
            print(
                f"Row {index}: mapping failed - {e}. Leaving start/end/segments/"
                f"Source_Sequence/Ref_sequence blank for this row and moving on."
            )
            failed_rows.append(index)
            if args.row_delay:
                time.sleep(args.row_delay)
            continue

        sf_table.at[index, 'start'] = start
        sf_table.at[index, 'end'] = end
        sf_table.at[index, 'segments'] = segments
        sf_table.at[index, 'Source_Sequence'] = source_seq
        sf_table.at[index, 'Ref_sequence'] = ref_seq

        if args.row_delay:
            time.sleep(args.row_delay)

    sf_table.to_csv(output_file, index=False)
    print(f"Mapping results saved to {output_file}")
    if failed_rows:
        print(
            f"{len(failed_rows)} of {len(sf_table)} row(s) could not be mapped "
            f"(see log above for reasons); their mapping columns were left blank: {failed_rows}"
        )


if __name__ == "__main__":
    main()
