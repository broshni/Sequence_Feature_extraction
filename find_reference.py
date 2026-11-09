#find the reference genome(s) for a virus from its name or NCBI taxonomy ID, so nobody has to
#know an accession number. Segmented viruses have one RefSeq record per segment, so records are
#grouped by taxon: one group is one strain's full set of segments.
import argparse
import sys
from collections import OrderedDict

from Bio import Entrez

import taxonomy


def group_records(summaries, names):
    """Group RefSeq nucleotide summaries into one entry per taxon.

    summaries: esummary records (need AccessionVersion, TaxId, Length).
    names: {taxid: scientific name}.
    Only NC_ accessions are kept - those are RefSeq genomes and genome segments; other prefixes
    are transcripts and gene records that are not references for mapping.
    Returns [{'taxid', 'name', 'accessions', 'records', 'bp'}], largest group first.
    """
    groups = OrderedDict()
    for record in summaries:
        accession = str(record['AccessionVersion'])
        if not accession.startswith('NC_'):
            continue
        taxid = str(int(record['TaxId']))  # Entrez hands back a special int type
        group = groups.setdefault(taxid, {
            'taxid': taxid, 'name': names.get(taxid, taxid), 'entries': []})
        group['entries'].append((accession, int(record['Length'])))

    ordered = []
    for group in groups.values():
        entries = sorted(group['entries'])
        lengths = [bp for _, bp in entries]
        # Segments of one genome differ a great deal in length (influenza: 890-2341 bp). Records
        # of near-equal length under one taxon are instead alternative reference genomes for the
        # same virus (Zika has two), and using them together would give two competing rows for
        # every protein - so each is its own candidate.
        alternatives = len(entries) > 1 and min(lengths) >= 0.95 * max(lengths)
        sets = [[e] for e in entries] if alternatives else [entries]
        for entry_set in sets:
            ordered.append({
                'taxid': group['taxid'], 'name': group['name'],
                'accessions': [a for a, _ in entry_set],
                'records': len(entry_set), 'bp': sum(bp for _, bp in entry_set)})
    return sorted(ordered, key=lambda g: (-g['records'], g['name'], g['accessions']))


def find_reference_genomes(query, email, retmax=200):
    """(taxon, groups) for a virus name or taxonomy ID. `taxon` is what the query resolved to,
    so the caller can show it; `groups` are the candidate reference sets."""
    taxonomy.configure(email)
    first = taxonomy.lookup(query)
    if first is None:
        raise ValueError(f"NCBI Taxonomy has no entry for '{query}'")

    # A generic name ('respiratory syncytial virus') can resolve to a taxon that has no RefSeq
    # genome of its own while a closely related one does, so try the other matches in turn.
    candidates = [first['taxid']]
    if not str(query).strip().isdigit():
        candidates += [t for t in taxonomy.search(query) if t != first['taxid']]

    for taxid in candidates:
        taxon = first if taxid == first['taxid'] else taxonomy.lookup(taxid)
        term = f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP]"
        handle = Entrez.esearch(db='nuccore', term=term, retmax=retmax)
        ids = Entrez.read(handle)['IdList']
        handle.close()
        if not ids:
            continue
        handle = Entrez.esummary(db='nuccore', id=','.join(ids))
        summaries = Entrez.read(handle)
        handle.close()
        names = taxonomy.names_for({str(int(s['TaxId'])) for s in summaries})
        groups = group_records(summaries, names)
        if groups:
            return taxon, groups

    # Last resort: search the genome records' own names. Catches a generic name whose RefSeq
    # genomes are filed under a differently-named taxon ('respiratory syncytial virus').
    words = ' AND '.join(f'{w}[Title]' for w in str(query).split() if w.isalnum())
    if words:
        term = f'({words}) AND srcdb_refseq[PROP]'
        handle = Entrez.esearch(db='nuccore', term=term, retmax=retmax)
        ids = Entrez.read(handle)['IdList']
        handle.close()
        if ids:
            handle = Entrez.esummary(db='nuccore', id=','.join(ids))
            summaries = Entrez.read(handle)
            handle.close()
            names = taxonomy.names_for({str(int(s['TaxId'])) for s in summaries})
            groups = group_records(summaries, names)
            if groups:
                return first, groups
    return first, []


def parse_args():
    parser = argparse.ArgumentParser(
        description="List the RefSeq reference genomes for a virus, from its name or NCBI "
                    "taxonomy ID. A segmented virus appears as one group per strain, listing "
                    "every segment's accession."
    )
    parser.add_argument("query", help="Virus name or NCBI taxonomy ID, e.g. 'measles virus' or 11234.")
    parser.add_argument("--email", default="rbhattac@jcvi.org",
                        help="Email address sent with NCBI Entrez requests.")
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        taxon, groups = find_reference_genomes(args.query, args.email)
    except ValueError as exc:
        sys.exit(str(exc))
    print(f"'{args.query}' -> {taxon['name']} (taxid {taxon['taxid']}, {taxon['rank'] or 'no rank'})")
    if not groups:
        sys.exit("No RefSeq reference genome found for this taxon. Try a more specific name, "
                 "or give a genome accession directly.")
    for number, group in enumerate(groups, 1):
        print(f"\n{number}. {group['name']}  [taxid {group['taxid']}]")
        print(f"   {group['records']} record(s), {group['bp']:,} bp: {', '.join(group['accessions'])}")


if __name__ == "__main__":
    main()
