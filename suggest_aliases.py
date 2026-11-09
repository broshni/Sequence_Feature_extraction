#work out which product-name spellings in the source SF tables mean which reference protein, so
#nobody has to write an alias file by hand. NCBI, UniProt and IEDB name one protein differently
#('large polymerase protein' / 'RNA-directed RNA polymerase L' / 'RNA dependent RNA polymerase');
#this finds the links from the data itself and writes them as aliases for Edit_SFTable.py.
import argparse
import os
import re
import sys

import pandas as pd

from Edit_SFTable import find_reference_matches, name_occurs
from sf_columns import read_csv, write_csv

# Words that say what kind of thing a protein is, not which protein it is. They are dropped
# before two names are compared, so 'Hemagglutinin glycoprotein' and 'hemagglutinin protein'
# both reduce to {hemagglutinin}.
GENERIC_WORDS = {
    'protein', 'glycoprotein', 'precursor', 'partial', 'putative', 'uncharacterized',
    'large', 'small', 'structural', 'non-structural', 'nonstructural', 'subunit', 'chain',
    'peptide', 'virus', 'viral', 'major', 'minor',
}

MIN_ALIAS_LENGTH = 3

# When a product could be two reference proteins, the better one has to beat the other by this
# factor or the product is left unmatched rather than guessed.
AMBIGUITY_MARGIN = 1.5


def name_phrases(product):
    """The name phrases in a product string: the main name plus any parenthetical synonyms.

    Bracketed sections ('[Cleaved into: ...]', '[Measles morbillivirus]') name chains or
    organisms rather than this protein, so they are dropped.
    """
    text = re.sub(r'\[[^\]]*\]', ' ', str(product))
    parts = [re.sub(r'\([^()]*\)', ' ', text)] + re.findall(r'\(([^()]*)\)', text)
    phrases = []
    for part in parts:
        part = re.sub(r'\s+', ' ', part).strip(' ,;:.')
        if len(part) >= MIN_ALIAS_LENGTH and part.lower() not in {p.lower() for p in phrases}:
            phrases.append(part)
    return phrases


def distinct_tokens(text):
    """The words in a name that identify the protein. Falls back to every word when a name is
    nothing but generic ones, so 'polyprotein' still has something to match on."""
    tokens = set(re.findall(r"[a-z0-9]+(?:-[a-z0-9]+)*", str(text).lower()))
    return (tokens - GENERIC_WORDS) or tokens


def same_protein(name_tokens, phrase_tokens):
    # A bare letter ('C', 'V', 'M') only names a protein when it is the whole name; as a
    # subset it would match any product that happens to contain that letter.
    if all(len(token) == 1 for token in name_tokens):
        return name_tokens == phrase_tokens
    return name_tokens <= phrase_tokens


def gene_tokens(value):
    """The gene symbol(s) in a gene field: only the first word, which is the gene's primary
    name - the rest are synonyms and the names of other products of the same gene ('PA PA-X'
    is the PA gene, which also makes PA-X). 'P/V/C' is three symbols in one gene."""
    words = str(value).split() if value is not None else []
    if not words or words[0].lower() == 'nan':
        return set()
    return {t.upper() for t in re.split(r'/', words[0]) if t}


def gene_index(reference):
    """{gene symbol: reference rows carrying it}. 'P/V/C' puts one row under each of P, V, C."""
    index = {}
    for row in reference.to_dict('records'):
        for token in gene_tokens(row['gene']):
            index.setdefault(token, []).append(row)
    return index


def resolve_by_gene(product, genes, index):
    """(main-name phrase, reference_row) when the product's gene symbols single out exactly one
    reference protein. Gene symbols agree across NCBI and UniProt even where the protein names
    are spelled differently ('haemagglutinin' / 'Hemagglutinin'), so this catches what the name
    comparison cannot. A symbol shared by several proteins (P/V/C) singles out none."""
    accessions = {}
    for token in genes:
        for row in index.get(token, []):
            accessions[row['protein_accession']] = row
    if len(accessions) != 1:
        return None
    phrases = name_phrases(product)
    if not phrases:
        return None
    return phrases[0], next(iter(accessions.values()))


def known_names(reference):
    """(tokens, reference_row) for every name the reference already gives a protein: its
    product name, and its gene symbol where only one protein has that gene ('M protein')."""
    rows = reference.to_dict('records')
    gene_counts = reference['gene'].astype(str).value_counts()
    names = []
    for row in rows:
        # RefSeq names carry synonyms too: 'non-structural protein 1 (1C)', 'Nucleoprotein (N)'
        for phrase in name_phrases(row['product']) or [str(row['product'])]:
            names.append((distinct_tokens(phrase), row))
        gene = str(row['gene']).strip().lower()
        if gene and gene != 'nan' and gene_counts[str(row['gene'])] == 1:
            names.append(({gene}, row))
    return names


def resolve(product, names):
    """(phrase, reference_row) for the reference protein this product names, or None when it
    names none, or two that cannot be told apart."""
    best = {}
    for phrase in name_phrases(product):
        tokens = distinct_tokens(phrase)
        for name_tokens, ref in names:
            if same_protein(name_tokens, tokens):
                score = len(name_tokens) / len(tokens)
                key = ref['protein_accession']
                if key not in best or score > best[key][0]:
                    best[key] = (score, phrase, ref)
    if not best:
        return None
    ranked = sorted(best.values(), key=lambda hit: -hit[0])
    if len(ranked) > 1 and ranked[0][0] < ranked[1][0] * AMBIGUITY_MARGIN:
        return None
    return ranked[0][1], ranked[0][2]


def suggest(reference, product_rows, product_genes=None):
    """Aliases for every product the reference does not already match.

    product_rows maps each source product string to the number of rows that carry it, and
    product_genes (optional) to the gene symbols those rows carry. Returns (aliases,
    unresolved): aliases is [{'alias', 'product', 'rows'}], unresolved the products that could
    not be placed.
    """
    names = known_names(reference)
    index = gene_index(reference)
    product_genes = product_genes or {}

    unresolved = []
    for product in product_rows:
        matches = find_reference_matches(product, reference)
        accessions = {m['protein_accession'] for m in matches}
        if len(accessions) == 1:
            # Already matched. UniProt lists a protein's synonyms in brackets, and those
            # synonyms are how a differently-spelled IEDB product gets linked to it.
            for phrase in name_phrases(product):
                names.append((distinct_tokens(phrase), matches[-1]))
        elif not accessions:
            unresolved.append(product)

    aliases = {}
    progress = True
    while progress:
        progress = False
        for product in list(unresolved):
            hit = resolve(product, names)
            if hit is None:
                hit = resolve_by_gene(product, product_genes.get(product, ()), index)
            if hit is None:
                continue
            phrase, ref = hit
            if not name_occurs(phrase, product):
                continue
            entry = aliases.setdefault(phrase.lower(), {
                'alias': phrase, 'product': ref['product'], 'rows': 0})
            entry['rows'] += product_rows[product]
            unresolved.remove(product)
            # a newly placed product can link the next one, so go round again
            for synonym in name_phrases(product):
                names.append((distinct_tokens(synonym), ref))
            progress = True
    return list(aliases.values()), unresolved


def expand_reference(reference, aliases):
    """The reference with each alias added as a row of its own, as build_reference does."""
    by_product = {str(r['product']).lower(): r for r in reference.to_dict('records')}
    rows = reference.to_dict('records')
    for entry in aliases:
        target = by_product[str(entry['product']).lower()]
        rows.append({'product': entry['alias'], 'gene': target['gene'],
                     'protein_accession': target['protein_accession']})
    return pd.DataFrame(rows, columns=['product', 'gene', 'protein_accession'])


def coverage(reference, product_rows):
    """(rows matched, rows ambiguous) when step 4 matches products against `reference`."""
    if reference['protein_accession'].nunique() == 1:
        # Edit_SFTable assigns every row to a one-protein reference, whatever its name
        return sum(product_rows.values()), 0
    matched = ambiguous = 0
    for product, count in product_rows.items():
        accessions = {m['protein_accession'] for m in find_reference_matches(product, reference)}
        if len(accessions) == 1:
            matched += count
        elif len(accessions) > 1:
            ambiguous += count
    return matched, ambiguous


def read_product_rows(sf_tables):
    """({product: row count}, total rows, {product: gene symbols seen with it})."""
    counts, genes = {}, {}
    total = 0
    for path in sf_tables:
        df = read_csv(path)
        if 'product' not in df.columns:
            raise SystemExit(f"{path} has no 'product' column - is it a source SF table?")
        total += len(df)
        for product, count in df['product'].dropna().astype(str).value_counts().items():
            counts[product] = counts.get(product, 0) + int(count)
        if 'gene' in df.columns:
            for product, gene in df[['product', 'gene']].dropna().astype(str).drop_duplicates().values:
                genes.setdefault(product, set()).update(gene_tokens(gene))
    return counts, total, genes


def parse_args():
    parser = argparse.ArgumentParser(
        description="Suggest product aliases automatically: work out which product-name "
                    "spellings in the source SF tables mean which reference protein, so no "
                    "alias file has to be written by hand."
    )
    parser.add_argument("reference_file",
                        help="Reference CSV (product, gene, protein_accession) from "
                             "build_reference.py.")
    parser.add_argument("sf_tables", nargs='+',
                        help="Source SF table CSVs whose product names should be matched, "
                             "e.g. the outputs of SF_Table_Uniprot.py and SF_Table_Epitope.py.")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Directory for the outputs (created if missing).")
    parser.add_argument("-n", "--name", default="output",
                        help="Prefix for the output files (default: 'output').")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    reference = read_csv(args.reference_file)
    for column in ('product', 'gene', 'protein_accession'):
        if column not in reference.columns:
            raise SystemExit(f"{args.reference_file} is missing a '{column}' column")

    product_rows, total, product_genes = read_product_rows(args.sf_tables)
    before, _ = coverage(reference, product_rows)

    aliases, unresolved = suggest(reference, product_rows, product_genes)
    expanded = expand_reference(reference, aliases)
    after, ambiguous = coverage(expanded, product_rows)

    out = lambda suffix: os.path.join(args.outdir, f"{args.name}_{suffix}")
    write_csv(pd.DataFrame(aliases, columns=['alias', 'product', 'rows']), out('aliases.csv'))
    write_csv(expanded, out('reference_aliased.csv'))
    write_csv(pd.DataFrame(
        [{'product': p, 'rows': product_rows[p]} for p in unresolved],
        columns=['product', 'rows']).sort_values('rows', ascending=False),
        out('aliases_unresolved.csv'))

    print(f"{total} source row(s), {len(product_rows)} distinct product name(s)")
    print(f"Matched to a reference protein: {before} row(s) before aliases, {after} after "
          f"({after / total:.1%})" if total else "No source rows")
    print(f"\nAliases added ({len(aliases)}):")
    for entry in sorted(aliases, key=lambda e: -e['rows']):
        print(f"  {entry['alias']!r} -> {entry['product']}  ({entry['rows']} row(s))")

    unresolved_rows = sum(product_rows[p] for p in unresolved)
    blank_rows = total - sum(product_rows.values())
    unresolved_rows += blank_rows
    print(f"\nLeft unmatched: {len(unresolved)} product name(s)"
          f"{f' plus {blank_rows} row(s) with no product name' if blank_rows else ''}, "
          f"{unresolved_rows} row(s) in all. "
          "Mapping skips these rows; they are usually host proteins or proteins the "
          "reference genome does not encode.")
    for product in sorted(unresolved, key=lambda p: -product_rows[p])[:8]:
        print(f"  {product[:90]!r}  ({product_rows[product]} row(s))")
    if len(unresolved) > 8:
        print(f"  ... and {len(unresolved) - 8} more in {out('aliases_unresolved.csv')}")

    if ambiguous:
        print(f"\nWARNING: {ambiguous} row(s) have a product name that fits more than one "
              "reference protein; step 4 will pick the last one listed in the reference.")

    print(f"\nAliases saved to {out('aliases.csv')}")
    print(f"Reference with aliases saved to {out('reference_aliased.csv')}")


if __name__ == "__main__":
    main()
