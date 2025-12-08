from Bio import Entrez, SeqIO
import csv
import pandas as pd

#fetch the genbank record
def fetch_genbank_record(accession_number):
    Entrez.email = "bvbrcacapria@gmail.com"     #put your email
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching GenBank record: {e}")
        return None


def extract_cds_info(record):
    cds_list = []
    for feature in record.features:
        if feature.type == "CDS":
            gene = feature.qualifiers.get("gene", [""])[0]
            product = feature.qualifiers.get("product", [""])[0]
            protein_id = feature.qualifiers.get("protein_id", [""])[0]
            cds_list.append({
                "gene": gene,
                "product": product,
                "protein_accession": protein_id
            })
    return cds_list

# Example usage:
record = fetch_genbank_record("NC_001498.1")
cds_info = extract_cds_info(record)
for cds in cds_info:
     print(cds)

#write the CDS information to a CSV file
def write_cds_to_csv(cds_info, output_file):
    with open(output_file, "w", newline='') as csvfile:
        fieldnames = ["gene", "product", "protein_accession"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for cds in cds_info:
            writer.writerow(cds)
write_cds_to_csv(cds_info, "/Users/rbhattac/Desktop/Measles_SFVT/reference.csv")