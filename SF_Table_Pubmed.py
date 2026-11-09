#build a standardized SF table from a PubMed-miner findings export, so mined features can be
#merged alongside the UniProt and IEDB tables as a third source. The mining itself is done
#elsewhere; this reads what it produced and puts it into the SF table's shape.
import argparse
import os
import re

import pandas as pd

from sf_columns import PUBMED_SOURCE, SF_TABLE_COLUMNS, read_csv, write_csv

# columns the findings export must carry; anything else in it is optional
REQUIRED_COLUMNS = ['pmid', 'virus', 'protein', 'sequence_feature_name', 'position', 'quote']

# sf_category values the UniProt and IEDB tables already use, so mined rows share the vocabulary
STRUCTURAL_RE = re.compile(r'\b(?:helix|helices|strand|sheet|turn|coil|coiled)\b', re.I)
EPITOPE_RE = re.compile(r'\bepitopes?\b', re.I)


def read_findings(path, encoding=None):
    """Read a findings export. Miner exports are typically UTF-8 but carry occasional stray
    Windows-1252 bytes, which the shared reader handles; `encoding` forces one instead."""
    return read_csv(path, **({'encoding': encoding} if encoding else {}))


def sf_category(target_type, feature_name):
    """The miner's own classification, mapped onto the SF table's category vocabulary."""
    if str(target_type).strip().lower() == 'mutation':
        return 'Sequence_Alterations'
    name = str(feature_name)
    if EPITOPE_RE.search(name):
        return 'Epitope'
    if STRUCTURAL_RE.search(name):
        return 'Structural'
    return 'Functional'


def clean_position(position):
    """The miner's 'position' written the way the SF table writes locations: '20-48' for a span,
    '21,22,24' for discrete residues.

    'position' is the authoritative field: start_position/end_position hold only its first and
    last value, and are blank on the rows that list discrete residues.
    """
    text = re.sub(r'\s+', '', str(position))
    if ',' not in text:
        return text
    residues = set()
    for part in (p for p in text.split(',') if p):
        if '-' in part:
            # a span inside a list: keep every residue it covers
            first, last = sorted(int(end) for end in part.split('-'))
            residues.update(range(first, last + 1))
        else:
            residues.add(int(part))
    return ','.join(str(residue) for residue in sorted(residues))


def split_features(feature_name, position):
    """The (name, position) pairs a findings row describes.

    A row can hold several features at once, as parallel comma-separated lists:
    'Long loop, L209, E211, G215' with '146-340, 209, 211, 215' is four features, not one
    feature at 198 residues. They are split only when the two lists are the same length - a
    single name against several positions ('interface' with '21, 22, 24, 26, 37') is one
    feature that spans those residues, and stays one row.
    """
    names = [part.strip() for part in str(feature_name).split(',')]
    positions = [part.strip() for part in str(position).split(',')]
    if len(positions) > 1 and len(names) == len(positions):
        return list(zip(names, positions))
    return [(str(feature_name), str(position))]


def feature_length(position, given):
    """How many residues the feature covers. The export's own length is used when it has one;
    otherwise it is counted from the position, which is right for a span and for a residue list."""
    if given is not None and pd.notna(given) and str(given).strip():
        return int(float(given))
    text = clean_position(position)
    if ',' in text:
        return len(text.split(','))
    ends = [int(part) for part in text.split('-') if part]
    return ends[-1] - ends[0] + 1 if len(ends) == 2 else 1


def convert(findings):
    """A findings export (one row per mined feature) as a standardized SF table."""
    missing = [column for column in REQUIRED_COLUMNS if column not in findings.columns]
    if missing:
        raise SystemExit(
            f"The findings file has no {', '.join(repr(c) for c in missing)} column(s). "
            f"It needs at least: {', '.join(REQUIRED_COLUMNS)}."
        )

    rows = []
    for finding in findings.itertuples():
      for feature_name, position in split_features(finding.sequence_feature_name,
                                                   finding.position):
        rows.append({
            'virus_type': finding.virus,
            'subtype': '',
            'segment': '',
            'genbank_accession': '',
            # the miner names the protein rather than the gene; the reference match fills gene_ref
            'gene': finding.protein,
            'product': finding.protein,
            'sf_category': sf_category(getattr(finding, 'target_type', ''), feature_name),
            'sf_type': feature_name,
            'start': '',
            'end': '',
            'segments': '',
            'Source': PUBMED_SOURCE,
            'source_id': str(finding.pmid).strip(),
            'source_strain': '',
            'source_aa_sequence': '',
            'source_sf_location': clean_position(position),
            'Evidence_code': 'EXP',
            'additional_metadata': '',
            # the sentence the feature was read from is its evidence, so it is its comment
            'comments': finding.quote,
            'length': feature_length(position, getattr(finding, 'length', None)),
            'publication': f'PMID: {str(finding.pmid).strip()}',
            # kept alongside 'comments' so the evidence stays identifiable as the miner's own
            # quotation even if a comment is later edited
            'evidence_sentence': finding.quote,
            'review_note': 'coordinates assumed relative to the reference protein',
        })
    return pd.DataFrame(rows, columns=SF_TABLE_COLUMNS + [
        'length', 'publication', 'evidence_sentence', 'review_note'])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a standardized SF table from a PubMed-miner findings export, so its "
                    "features can be merged alongside the UniProt and IEDB tables."
    )
    parser.add_argument(
        "findings_file",
        help="The miner's findings CSV, with pmid, virus, protein, sequence_feature_name, "
             "position and quote columns.",
    )
    parser.add_argument("-o", "--outdir", required=True,
                        help="Directory where the output CSV will be written (created if missing).")
    parser.add_argument("-n", "--name", default="output",
                        help="Prefix used to name the output file (default: 'output').")
    parser.add_argument("--encoding", default=None,
                        help="Force a text encoding for the findings file. By default it is "
                             "read as UTF-8, falling back per byte to Windows-1252, which "
                             "handles the mixed encoding these exports usually have.")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    findings = read_findings(args.findings_file, args.encoding)
    sf_table = convert(findings)

    out_path = os.path.join(args.outdir, f"{args.name}_Pubmed_SF.csv")
    write_csv(sf_table, out_path)

    print(f"{len(findings)} finding(s) from {findings['pmid'].nunique()} paper(s) "
          f"-> {len(sf_table)} SF row(s)")
    if len(sf_table):
        print("  by category:", sf_table['sf_category'].value_counts().to_dict())
        print("  by protein :", sf_table['gene'].value_counts().to_dict())
    print(f"PubMed SF table saved to {out_path}")


if __name__ == "__main__":
    main()
