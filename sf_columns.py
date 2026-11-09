#shared column conventions for the SF pipeline. Kept in its own module so every step can use
#them without importing a step that drags in Bio/MAFFT.
import codecs
import csv
import re

import pandas as pd

# Positions and lengths are whole numbers. pandas reads any column containing a blank cell as
# float64, so without this they round-trip through the pipeline as '525.0' / '1.0'.
INT_COLUMNS = ('start', 'end', 'length')

# The canonical SF table column set every source table is built to. Sources may carry extra
# columns beyond these (the merge keeps them), but these are the ones the pipeline reads.
SF_TABLE_COLUMNS = [
    'virus_type', 'subtype', 'segment', 'genbank_accession', 'gene', 'product',
    'sf_category', 'sf_type', 'start', 'end', 'segments', 'Source', 'source_id',
    'source_strain', 'source_aa_sequence', 'source_sf_location', 'Evidence_code',
    'additional_metadata', 'comments',
]

# Feature types whose 'start-end' location is the two residues that are bonded rather than a
# span: UniProt annotates a disulfide bond as the pair of cysteines ('188-606' = Cys188 bonded
# to Cys606), and a cross-link the same way. Every other range type ('Helix 159-163',
# 'Chain 1-525', an IEDB epitope) really is contiguous.
PAIR_SF_TYPES = {'Disulfide bond', 'Cross-link'}


def to_int_columns(df):
    """Return df with INT_COLUMNS as nullable Int64, so they write as '525' not '525.0'.

    A column is left untouched if any non-blank cell isn't a number, or if any value isn't a
    whole number - better to write the odd-looking value than to silently coerce or round it
    away.
    """
    df = df.copy()
    for col in INT_COLUMNS:
        if col not in df.columns:
            continue
        numeric = pd.to_numeric(df[col], errors='coerce')
        was_filled = df[col].notna() & (df[col].astype(str).str.strip() != '')
        lost_values = (was_filled & numeric.isna()).any()
        has_fractions = not (numeric.dropna() % 1 == 0).all()
        if lost_values or has_fractions:
            continue
        df[col] = numeric.astype('Int64')
    return df


# The Source value of every row that came from the PubMed miner, and the one place it is
# spelled. A row's Source names exactly one source and its source_id is that source's own
# identifier - a UniProt accession, an IEDB id, or the paper's numeric PMID.
PUBMED_SOURCE = 'PMID'


def publication_label(source, source_id):
    """How a source record is cited in the 'publication' column: PubMed rows by their PMID
    ('PMID: 40493685'), every other source as 'Source:id' ('IEDB:2322', 'Uniprot:P04851')."""
    if str(source).strip().lower() in ('pmid', 'pubmed', 'pubmed miner'):
        return f"PMID: {source_id}"
    return f"{source}:{source_id}"


# Values like '1-36' or '70-77' are what Excel silently turns into dates ('1-Mar', '70-Jul').
# Quoting does not stop it - Excel strips the quotes and converts anyway - so those cells are
# written as ="1-36", which Excel shows as the text it is. read_csv unwraps them again, so the
# pipeline and any other reader see the plain value.
DATE_LIKE = re.compile(r'^\s*\d{1,4}\s*-\s*\d{1,4}\s*$')
EXCEL_GUARD = re.compile(r'^="(.*)"$')


def is_text(series):
    """True for a column that holds text. Checked by dtype kind rather than `== object`:
    pandas 3 gives string columns a 'str' dtype, which is not object."""
    return pd.api.types.is_string_dtype(series) or series.dtype == object


def excel_guard(value):
    return f'="{value}"' if isinstance(value, str) and DATE_LIKE.match(value) else value


def strip_excel_guard(value):
    if isinstance(value, str):
        match = EXCEL_GUARD.match(value)
        if match:
            return match.group(1)
    return value


def _windows_fallback(error):
    """Decode a byte UTF-8 rejected as Windows-1252 rather than discarding it."""
    return error.object[error.start:error.end].decode('cp1252', errors='replace'), error.end


codecs.register_error('windows_fallback', _windows_fallback)


def read_csv(path, **kwargs):
    """pd.read_csv for SF tables: only an empty cell counts as missing, and text decodes the
    same way everywhere.

    pandas also treats the text 'NA', 'N/A', 'NULL', 'None' and 'nan' as missing - and 'NA' is
    the gene symbol of influenza neuraminidase, so every NA protein would silently lose its
    gene (and with it the reference match, the sf_id and the dedup grouping) at each step that
    reads a table. Reading the symbol as the text it is avoids that.
    """
    kwargs.setdefault('keep_default_na', False)
    kwargs.setdefault('na_values', [''])
    # Read as UTF-8, falling back per byte to Windows-1252. Tables carry quoted sentences from
    # papers, which hold en-dashes and Greek letters; reading those as latin-1 (as the steps
    # used to) turns every one into mojibake, and plain UTF-8 stops dead on a stray byte.
    if not kwargs.get('encoding'):
        kwargs['encoding'] = 'utf-8'
        kwargs.setdefault('encoding_errors', 'windows_fallback')
    df = pd.read_csv(path, **kwargs)
    for column in df.columns:
        if is_text(df[column]):
            df[column] = df[column].map(strip_excel_guard)
    return df


def write_csv(df, path, excel_safe=True, **kwargs):
    """to_csv for SF tables: integer columns as integers, every field quoted, no index column.

    With excel_safe (the default), coordinate ranges are written as ="1-36" so Excel does not
    read them as dates. Pass excel_safe=False for a plain file.
    """
    out = to_int_columns(df)
    if excel_safe:
        out = out.copy()
        for column in out.columns:
            if is_text(out[column]):
                out[column] = out[column].map(excel_guard)
    kwargs.setdefault('quoting', csv.QUOTE_ALL)
    out.to_csv(path, index=False, **kwargs)
