#build the reference CSV (product, gene, protein_accession) that Edit_SFTable.py consumes,
#by reading the CDS features off a reference genome record instead of copying them by hand
import argparse
import os
import time
import pandas as pd
from sf_columns import read_csv
from Bio import Entrez, SeqIO


def cds_rows(record):
    """One row per CDS in the record: its product name, gene, and protein accession."""
    rows = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        q = feature.qualifiers
        protein_accession = q.get("protein_id", [None])[0]
        if not protein_accession:
            continue
        # /product is what we match UniProt/IEDB product strings against; fall back to the
        # gene name so a CDS is never dropped just for lacking a product qualifier.
        gene = q.get("gene", [None])[0] or q.get("locus_tag", [None])[0]
        product = q.get("product", [None])[0] or gene
        if not product:
            continue
        rows.append({"product": product, "gene": gene, "protein_accession": protein_accession})
    return rows


def add_aliases(rows, aliases_file):
    """Add extra product spellings, so UniProt/IEDB names that don't match NCBI's wording
    still resolve. The alias file maps 'alias' -> an existing reference 'product'; the
    aliased row reuses that product's gene and accession."""
    aliases = read_csv(aliases_file)
    missing = {"alias", "product"} - set(aliases.columns)
    if missing:
        raise SystemExit(f"Alias file {aliases_file} is missing column(s): {sorted(missing)}")

    by_product = {r["product"].lower(): r for r in rows}
    extra = []
    for _, a in aliases.iterrows():
        target = by_product.get(str(a["product"]).lower())
        if target is None:
            raise SystemExit(
                f"Alias '{a['alias']}' points at product '{a['product']}', which is not one of "
                f"the reference products: {sorted(by_product)}"
            )
        extra.append({"product": a["alias"], "gene": target["gene"],
                      "protein_accession": target["protein_accession"]})
    return rows + extra


def fetch_genome(accession, email, retries=3, retry_delay=5):
    Entrez.email = email
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            return record
        except Exception as exc:
            last_error = exc
            print(f"  fetch of {accession} failed (attempt {attempt}/{retries}): {exc}")
            if attempt < retries:
                time.sleep(retry_delay)
    raise RuntimeError(f"Could not fetch {accession} after {retries} attempts: {last_error}")


def check_reference_is_this_virus(accession, sf_table, email):
    """Stop when the reference genome belongs to a different virus than the source table.

    Everything downstream depends on it: with another virus's reference, product names match
    the wrong proteins, the organism filter drops the right epitopes, and mapping rejects every
    row as too dissimilar - leaving a near-empty table with no error at all."""
    import taxonomy
    from filter_organism import organism_agreement, target_species_taxid, target_taxon
    from sf_columns import read_csv

    taxonomy.configure(email)
    taxon = target_taxon(accession)
    organisms = read_csv(sf_table)['virus_type'].dropna().astype(str).value_counts().to_dict()
    agreement = organism_agreement(organisms, target_species_taxid(taxon))
    if agreement is None:
        print(f"Could not identify the organisms in {sf_table}; reference not checked.")
        return
    if agreement < 0.5:
        top = ', '.join(list(organisms)[:3])
        raise SystemExit(
            f"The reference genome {accession} is {taxon['name']}, but {1 - agreement:.0%} of the "
            f"rows in {os.path.basename(sf_table)} are from another virus ({top}). Choose this "
            "virus's reference genome (find_reference.py lists them) and run again.")
    print(f"Reference {accession} ({taxon['name']}) matches {agreement:.0%} of the source rows")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build the reference CSV (product, gene, protein_accession) used by "
                    "Edit_SFTable.py --reference-file, from the CDS features of one or more "
                    "reference genome accessions (e.g. NC_001498.1 for measles)."
    )
    parser.add_argument(
        "accession",
        nargs='+',
        help="Reference genome nucleotide accession(s), e.g. 'NC_001498.1'. A segmented virus "
             "has one accession per segment: give them all, space- or comma-separated "
             "(find_reference.py lists them).",
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Directory where the reference CSV will be written (created if missing).",
    )
    parser.add_argument(
        "-n", "--name",
        default="output",
        help="Prefix used to name the output file, i.e. '<name>_reference.csv' "
             "(default: 'output').",
    )
    parser.add_argument(
        "--aliases",
        default=None,
        help="Optional CSV with 'alias' and 'product' columns, adding extra product spellings "
             "for products already in the reference (e.g. alias 'nucleoprotein' -> product "
             "'nucleocapsid protein'), so UniProt/IEDB naming still matches.",
    )
    parser.add_argument(
        "--check-against",
        default=None,
        help="A source SF table (e.g. <name>_Uniprot_SF.csv). Stops with an error when most of "
             "its rows are from a different virus than the reference genome - a reference left "
             "over from another virus would otherwise silently map nothing.",
    )
    parser.add_argument(
        "--email",
        default="rbhattac@jcvi.org",
        help="Email address to send with NCBI Entrez requests, as required by NCBI.",
    )
    parser.add_argument(
        "--retries", type=int, default=3,
        help="How many times to retry the NCBI fetch before giving up (default: 3).",
    )
    parser.add_argument(
        "--retry-delay", type=float, default=5,
        help="Seconds to wait between fetch retries (default: 5).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    accessions = [a for chunk in args.accession for a in chunk.replace(',', ' ').split()]
    rows = []
    for accession in accessions:
        record = fetch_genome(accession, args.email, args.retries, args.retry_delay)
        found = cds_rows(record)
        print(f"Fetched {record.id}: {record.description} ({len(found)} protein(s))")
        rows.extend(found)
    if not rows:
        raise SystemExit(
            f"No CDS features with a protein accession found in {', '.join(accessions)}. "
            "Check that these are annotated genome records."
        )
    if args.check_against:
        check_reference_is_this_virus(accessions[0], args.check_against, args.email)

    # the same protein can be annotated twice (e.g. overlapping records); keep one row each
    seen, unique = set(), []
    for row in rows:
        key = (row['protein_accession'], str(row['product']).lower())
        if key not in seen:
            seen.add(key)
            unique.append(row)
    rows = unique

    if args.aliases:
        rows = add_aliases(rows, args.aliases)

    df = pd.DataFrame(rows, columns=["product", "gene", "protein_accession"])
    out_path = os.path.join(args.outdir, f"{args.name}_reference.csv")
    df.to_csv(out_path, index=False)
    print(df.to_string(index=False))
    print(f"\nReference table ({len(df)} rows) saved to {out_path}")


if __name__ == "__main__":
    main()
