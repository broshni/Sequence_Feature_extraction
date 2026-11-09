#drop IEDB epitopes whose source organism is not the virus being built. An IEDB export for one
#virus also holds epitopes from related viruses and from host proteins (a Zika export carries
#Dengue, West Nile and Yellow fever epitopes), and mapping those onto this virus's reference
#would place them at coordinates that mean nothing there. Organisms are compared by NCBI
#taxonomy lineage, so no virus names are written into the code.
import argparse
import os
import sys

import re

import pandas as pd
from Bio import Entrez

import taxonomy
from sf_columns import read_csv, write_csv


def target_taxon(reference_accession, taxon_override=None):
    """The taxon the table is being built for: the species of the reference genome, or the
    `taxon_override` (name or ID) if one is given. Returns the lookup dict, or None."""
    if taxon_override:
        return taxonomy.lookup(taxon_override)
    handle = Entrez.esummary(db='nuccore', id=reference_accession)
    taxid = str(int(Entrez.read(handle)[0]['TaxId']))
    handle.close()
    return taxonomy.lookup(taxid)


def target_species_taxid(taxon):
    """Species taxid of `taxon`; the taxon itself when NCBI gives it no species rank."""
    return taxon['species_taxid'] or taxon['taxid']


def organism_agreement(organism_counts, species_taxid, lookup=taxonomy.lookup, limit=8):
    """Fraction of rows whose organism falls within `species_taxid`, judged on the most common
    organisms (up to `limit` lookups). organism_counts: {organism name: rows}. Organisms NCBI
    cannot identify are left out of the fraction rather than counted against it."""
    judged = agree = 0
    for name, count in sorted(organism_counts.items(), key=lambda kv: -kv[1])[:limit]:
        taxon = lookup(name)
        if taxon is None:
            continue
        judged += count
        agree += count if species_taxid in taxon['lineage'] else 0
    return agree / judged if judged else None


# words that appear in any virus name and so identify nothing on their own
GENERIC_NAME_WORDS = {'virus', 'viruses', 'strain', 'isolate', 'type', 'human', 'species',
                      'genus', 'variant', 'subtype'}


def name_words(name):
    """The words of an organism name that actually identify it."""
    return {word for word in re.findall(r'[a-z]+', str(name).lower())
            if len(word) > 3 and word not in GENERIC_NAME_WORDS}


def names_agree(name, target_words):
    """True when an organism name shares an identifying word with the virus being built.

    Used only for names NCBI cannot resolve. Matching on prefixes lets 'Ebola' agree with
    'ebolavirus', while 'Marburg' agrees with neither - NCBI has no searchable 'Marburg virus'
    (the species is now Orthomarburgvirus marburgense), so without this the row would be kept.
    """
    for word in name_words(name):
        for target in target_words:
            if word.startswith(target[:5]) or target.startswith(word[:5]):
                return True
    return False


def filter_rows(df, species_taxid, lookup=taxonomy.lookup, target_lineage=(), target_words=()):
    """Split rows by whether their source organism is the virus being built.

    IEDB and PubMed rows are judged, because both name an organism that may not be this virus:
    an IEDB export carries epitopes from related viruses, and a paper about Ebola also discusses
    Marburg and Reston. UniProt rows are not judged - they came from a taxonomy search for this
    virus already.

    A row is kept when its organism is the target species, something below it (a strain), or
    something above it: a paper that says only 'Ebolavirus' means the genus this virus sits in,
    not a different virus. A sibling species (Reston, Bundibugyo) is on neither path and is
    dropped. An organism NCBI cannot identify at all is kept unless its name shares no word
    with the virus being built - the benefit of the doubt, but not to a name that plainly
    belongs to something else.

    Returns (kept, dropped, unknown) where unknown lists the organism names that could not be
    checked.
    """
    judged = df['Source'].astype(str).str.strip().str.lower().isin(('iedb', 'pmid'))
    target_lineage = set(target_lineage) | {str(species_taxid)}
    verdict, unknown = {}, []
    for name in df.loc[judged, 'virus_type'].dropna().astype(str).unique():
        taxon = lookup(name)
        if taxon is None:
            unknown.append(name)
            verdict[name] = names_agree(name, target_words) if target_words else True
        else:
            verdict[name] = (species_taxid in taxon['lineage']
                             or taxon['taxid'] in target_lineage)

    def belongs(row):
        if not judged[row.name]:
            return True
        return verdict.get(str(row['virus_type']), True)

    keep = df.apply(belongs, axis=1)
    return df[keep].reset_index(drop=True), df[~keep].reset_index(drop=True), unknown


def parse_args():
    parser = argparse.ArgumentParser(
        description="Drop IEDB epitopes whose source organism is not the virus being built "
                    "(related viruses, host proteins), judged by NCBI taxonomy lineage."
    )
    parser.add_argument("input_file", help="Merged SF table CSV (output of Edit_SFTable.py).")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Directory for the outputs (created if missing).")
    parser.add_argument("-n", "--name", default="output",
                        help="Prefix for the output files (default: 'output').")
    parser.add_argument("--reference", default=None,
                        help="Reference genome accession; its species is the virus being built. "
                             "With several (a segmented virus) the first is used.")
    parser.add_argument("--taxon", default=None,
                        help="Virus name or taxonomy ID to keep, instead of deriving it from "
                             "--reference.")
    parser.add_argument("--email", default="rbhattac@jcvi.org",
                        help="Email address sent with NCBI Entrez requests.")
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.reference and not args.taxon:
        sys.exit("Give --reference (a genome accession) or --taxon (a virus name or ID).")
    os.makedirs(args.outdir, exist_ok=True)
    taxonomy.configure(args.email)

    df = read_csv(args.input_file)
    for column in ('Source', 'virus_type'):
        if column not in df.columns:
            sys.exit(f"{args.input_file} has no '{column}' column - is it a merged SF table?")

    reference = args.reference.replace(',', ' ').split()[0] if args.reference else None
    taxon = target_taxon(reference, args.taxon)
    out = lambda suffix: os.path.join(args.outdir, f"{args.name}_{suffix}")

    if taxon is None:
        print("Could not identify the virus's taxon, so no epitopes were filtered.")
        write_csv(df, out('SF_filtered.csv'))
        write_csv(df.iloc[0:0], out('dropped_other_organism.csv'))
        return

    species = target_species_taxid(taxon)
    kept, dropped, unknown = filter_rows(df, species, target_lineage=taxon['lineage'],
                                         target_words=name_words(taxon['name']))

    write_csv(kept, out('SF_filtered.csv'))
    write_csv(dropped, out('dropped_other_organism.csv'))

    print(f"Keeping organisms within {taxon['name']} (species taxid {species})")
    print(f"{len(df)} row(s) in, {len(kept)} kept, {len(dropped)} dropped as another organism")
    if len(dropped):
        print("Dropped by source organism:")
        for organism, count in dropped['virus_type'].value_counts().items():
            print(f"  {organism}: {count}")
    if unknown:
        print(f"Could not identify {len(unknown)} organism name(s) at NCBI: "
              f"{', '.join(unknown[:8])} (kept unless the name shares no word with "
              f"{taxon['name']})")
    print(f"\nFiltered table saved to {out('SF_filtered.csv')}")


if __name__ == "__main__":
    main()
