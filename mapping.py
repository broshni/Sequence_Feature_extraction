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

# PAIR_SF_TYPES: feature types whose 'a-b' location is a residue pair, not a span.
from sf_columns import PAIR_SF_TYPES, read_csv, write_csv
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


def is_pair_feature(sf_type):
    """True if this feature type's 'a-b' location means two discrete residues, not a range."""
    return str(sf_type).strip() in PAIR_SF_TYPES


def wanted_source_residues(coordinates, is_pair=False):
    """The source residues a location refers to, and whether they form a contiguous span.

    'a-b' is every residue from a to b - except for residue-pair types (disulfide bonds), where
    it is just residues a and b. 'a,b,c' is exactly those residues. 'a' is that one residue.
    """
    text = str(coordinates).strip()
    if '-' in text:
        ends = sorted(int(part) for part in text.split('-'))
        if len(ends) != 2:
            raise ValueError(f"Could not parse location '{coordinates}'")
        if ends[0] == ends[1]:
            return [ends[0]], True          # one residue, e.g. a disulfide reported as '139-139'
        return (ends, False) if is_pair else (list(range(ends[0], ends[1] + 1)), True)
    if ',' in text:
        return sorted({int(part) for part in text.split(',')}), False
    return [int(text)], True


def residue_map(alignment):
    """For every source residue (1-based): (alignment column, reference residue position).

    The reference position is None where the source residue sits in an insertion, a stretch
    the reference has no residues for. Alignment columns are NOT reference positions: the
    reference row has gaps wherever the source has an insertion, and every column after a gap
    is shifted by it, so reading columns as positions misplaces every feature after one.
    """
    ref_row, source_row = str(alignment[0].seq), str(alignment[1].seq)
    placed, ref_position, source_position = {}, 0, 0
    for column, (ref_char, source_char) in enumerate(zip(ref_row, source_row), start=1):
        if ref_char != '-':
            ref_position += 1
        if source_char != '-':
            source_position += 1
            placed[source_position] = (column, ref_position if ref_char != '-' else None)
    return placed


def percent_identity(alignment):
    """Fraction of residues identical where source and reference are both present."""
    both = [(r, s) for r, s in zip(str(alignment[0].seq), str(alignment[1].seq))
            if r != '-' and s != '-']
    if not both:
        return 0.0
    return sum(1 for r, s in both if r.upper() == s.upper()) / len(both)


def map_feature(alignment, source_coordinates, is_pair=False):
    """Place a source feature on the reference.

    Returns (start, end, segments, source_seq, ref_seq), with start/end/segments in REFERENCE
    residue positions and the two sequences cut from the alignment for just this feature.
    """
    residues, contiguous = wanted_source_residues(source_coordinates, is_pair)
    placed = residue_map(alignment)
    hits = [placed[r] for r in residues if r in placed and placed[r][1] is not None]
    if not hits:
        raise ValueError("None of the feature's residues could be placed on the reference")

    ref_row, source_row = str(alignment[0].seq), str(alignment[1].seq)
    columns = [column for column, _ in hits]
    positions = [position for _, position in hits]
    if contiguous:
        first, last = min(columns), max(columns)
        source_seq, ref_seq = source_row[first - 1:last], ref_row[first - 1:last]
        start, end = min(positions), max(positions)
        segments = f"{start}-{end}"
    else:
        # discrete residues: just those columns, comma-separated so nothing reads them as a span
        source_seq = ''.join(source_row[c - 1] for c in columns)
        ref_seq = ''.join(ref_row[c - 1] for c in columns)
        start, end = positions[0], positions[-1]
        segments = ','.join(map(str, positions))
    return start, end, segments, source_seq, ref_seq


def new_cache():
    return {'sequences': {}, 'alignments': {}, 'failures': {}}


def load_sequence(kind, identifier, cache, email, retries, retry_delay, timeout):
    """A sequence from GenBank ('genbank') or UniProt ('uniprot'), fetched once per run.
    A failed fetch is remembered too, so the dozens of rows sharing one bad accession do not
    each spend a full round of retries rediscovering it."""
    key = (kind, identifier)
    if key in cache['failures']:
        raise ValueError(cache['failures'][key])
    if key not in cache['sequences']:
        try:
            if kind == 'genbank':
                record = fetch_genbank_record(identifier, email=email, retries=retries,
                                              retry_delay=retry_delay)
                cache['sequences'][key] = str(record.seq)
            else:
                cache['sequences'][key] = str(fetch_uniprot(
                    identifier, retries=retries, retry_delay=retry_delay, timeout=timeout))
        except Exception as exc:
            cache['failures'][key] = f"could not fetch {kind} sequence {identifier}: {exc}"
            raise ValueError(cache['failures'][key]) from exc
    return cache['sequences'][key]


#function to map one SF table row's source_sf_location onto its reference sequence.
#Returns (start, end, segments, source_seq, ref_seq); raises on any failure so the caller
#can decide how to handle it (see main(), which leaves the row blank and moves on).
def map_sf_row(row, ref_fasta_path, source_fasta_path, alignment_output_path,
               email="rbhattac@jcvi.org", retries=3, retry_delay=5, timeout=30,
               min_identity=0.5, cache=None):
    cache = cache if cache is not None else new_cache()

    source_coordinates = row['source_sf_location']
    if pd.isna(source_coordinates) or not str(source_coordinates).strip():
        raise ValueError("Missing source_sf_location")
    source_coordinates = str(source_coordinates)

    ref_accession = row['genbank_accession']
    if pd.isna(ref_accession) or not str(ref_accession).strip():
        raise ValueError("Missing genbank_accession")

    if str(row['Source']).lower() == 'uniprot':
        source_kind, source_id = 'uniprot', row['source_id']
    else:
        metadata = str(row['additional_metadata'])
        parts = metadata.split(":")
        if len(parts) < 2:
            if metadata.strip() and metadata.strip().lower() != 'nan':
                raise ValueError(f"Could not parse additional_metadata '{metadata}'")
            # No source sequence named: the feature's coordinates are already numbered along
            # the reference protein, as a paper's residue numbers are. Aligning the reference
            # to itself maps them through unchanged and still records the residues involved.
            parts = ['NCBI', str(ref_accession)]
        source_type, source_id = parts[0], parts[1]
        if source_type.upper() == 'NCBI':
            source_kind = 'genbank'
        elif source_type.lower() == 'uniprot':
            source_kind = 'uniprot'
        else:
            raise ValueError(f"Unrecognized source type '{source_type}' in additional_metadata")

    # one alignment per (reference, source) pair - a table has thousands of rows but only a
    # few dozen distinct source proteins, so aligning per row repeats the same work
    pair = (ref_accession, source_kind, source_id)
    if pair not in cache['alignments']:
        ref_seq = load_sequence('genbank', ref_accession, cache, email, retries, retry_delay, timeout)
        source_seq = load_sequence(source_kind, source_id, cache, email, retries, retry_delay, timeout)
        with open(ref_fasta_path, "w") as f:
            f.write(f">{ref_accession}\n{ref_seq}\n")
        with open(source_fasta_path, "w") as f:
            f.write(f">{source_id}\n{source_seq}\n")
        align_sequences(ref_fasta_path, source_fasta_path, alignment_output_path)
        cache['alignments'][pair] = AlignIO.read(alignment_output_path, "fasta")
    alignment = cache['alignments'][pair]

    # A source that barely resembles the reference is a different protein (or another
    # organism's), and the "mapping" would be coordinates on a meaningless alignment.
    identity = percent_identity(alignment)
    if identity < min_identity:
        raise ValueError(f"source {source_id} is only {identity:.0%} identical to reference "
                         f"{ref_accession} (minimum {min_identity:.0%}) - not the same protein")

    return map_feature(alignment, source_coordinates, is_pair_feature(row.get('sf_type')))


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
        default=None,
        help="Force a text encoding for the input CSV. By default it is read as UTF-8, "
             "falling back per byte to Windows-1252.",
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
        "--min-identity",
        type=float,
        default=0.5,
        help="Skip a row when its source sequence is less than this fraction identical to the "
             "reference protein it is being mapped onto (default: 0.5). Catches another "
             "organism's protein, or the wrong protein, rather than mapping onto it.",
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
    sf_table = read_csv(args.input_file, **({'encoding': args.encoding} if args.encoding else {}))

    # 'start', 'end' and 'segments' already exist (empty) in the merged SF table schema;
    # Source_Sequence/Ref_sequence are new. Force all five to an object dtype so an
    # all-empty column (which pandas would otherwise infer as float64) can safely hold the
    # ints/strings we write into it below.
    for col in ['start', 'end', 'segments', 'Source_Sequence', 'Ref_sequence']:
        if col not in sf_table.columns:
            sf_table[col] = None
        sf_table[col] = sf_table[col].astype('object')

    failed_rows = []
    cache = new_cache()
    for index, row in sf_table.iterrows():
        print(f"Processing row {index} (source_id={row.get('source_id', '')})...")
        try:
            start, end, segments, source_seq, ref_seq = map_sf_row(
                row, ref_fasta_path, source_fasta_path, alignment_output_path,
                email=args.email, retries=args.retries, retry_delay=args.retry_delay,
                timeout=args.timeout, min_identity=args.min_identity, cache=cache,
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

    write_csv(sf_table, output_file)
    print(f"Mapping results saved to {output_file}")
    if failed_rows:
        print(
            f"{len(failed_rows)} of {len(sf_table)} row(s) could not be mapped "
            f"(see log above for reasons); their mapping columns were left blank: {failed_rows}"
        )


if __name__ == "__main__":
    main()
