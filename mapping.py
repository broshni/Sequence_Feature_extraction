#parse the SFVT csv file and extract the source and reference information
from Bio import AlignIO
import pandas as pd
from Bio import Entrez, SeqIO
import csv
import subprocess
import requests as r
from io import StringIO
import math

def fetch_uniprot(uni_id):
    cID = uni_id
    try:
        baseUrl = "http://www.uniprot.org/uniprot/"
        currentUrl = baseUrl + cID + ".fasta"
        response = r.post(currentUrl)
        cData = ''.join(response.text)

        Seq = StringIO(cData)
        pSeq = list(SeqIO.parse(Seq, 'fasta'))
        #print(pSeq[0].seq)
        return pSeq[0].seq
    except Exception as e:
        print(f"Error fetching uniprot record: {e}")
        return None

def fetch_genbank_record(accession_number):
    Entrez.email = "rbhattac@jcvi.org"     #put your email
    try:
        handle = Entrez.efetch(db="protein", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching GenBank record: {e}")
        return None

def align_sequences(ref_seq, otherseq, output_file):
    #mafft_command = f" mafft --6merpair --anysymbol --keeplength --addfragments {otherseq} {ref_seq} > {output_file}"
    mafft_command = f" mafft --6merpair --anysymbol --addfragments {otherseq} {ref_seq} > {output_file}"
    subprocess.call(mafft_command, shell=True)

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
    print(alignment[num].id)
    source_seq = alignment[num].seq

    pos = 0  # set the starting pos
    ref_position = []
    for index, amino in enumerate(source_seq):
        # print(index, amino)

        if amino != "-":
            # print(amino)

            pos += 1
            if pos in coordinates:
                print(pos)
                print(index, amino)
                ref_position.append(index)

    print(ref_position)
    return ref_position


def main():

    #read the csv file as a dataframe
    SF_Table = pd.read_csv("/Users/rbhattac/Desktop/Chikungunya_SFVT/map_left.csv")
    #SF_Table = pd.read_csv("/Users/rbhattac/Desktop/results_code_don_1.csv")
    #print column names of the dataframe
    print(SF_Table['additional_metadata'])

    #open the csv file to write
    with open("/Users/rbhattac/Desktop/Chikungunya_SFVT/map_left_mapping.csv", "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header
        csv_writer.writerow(
            ["Source Accession", "Source Coordinates", "start", "end", "segments", "Source_Sequence", "Ref_Sequence"])

        #access rows of the dataframe one by one
        for index, row in SF_Table.iterrows():
            print(index)
            source=""
            ref_seq=""
            source_seq=""
            source_coordinates=row["source_sf_location"]

            #assign the ref sequence
            #ref_seq = "/Users/rbhattac/Desktop/Dengue1_ref.fasta"
            ref=row['genbank_accession']
            try:
                record = fetch_genbank_record(ref)
                ref_seq_temp = record.seq
            except Exception as e:
                print(f"Error fetching GenBank record: {e}")
                continue
            # write ref_seq_temp to the fasta file
            print(str(ref_seq_temp))
            ref_seq= "/Users/rbhattac/Desktop/ref.fasta"
            file = open(ref_seq, "w+")
            file.write(">"+str(ref)+"\n")
            file.write(str(ref_seq_temp))
            file.close()


            #ref = row['genbank_accession']
            if row['Source'].lower() =='uniprot':
                source = row['source_id']
                source_seq = fetch_uniprot(source)
            else:
                temp_source = row['additional_metadata']
                #temp_source = str(temp_source)
                temp_source = temp_source.split(":")
                if temp_source[0] == 'NCBI':
                    source = temp_source[1]
                    record= fetch_genbank_record(source)
                    source_seq =record.seq
                    print(source_seq)
                elif temp_source[0] == 'Uniprot' or temp_source[0] == 'uniprot':
                    source = temp_source[1]
                    source_seq = fetch_uniprot(source)
                    print(source_seq)
                #print(source_seq)

            # write source_seq to the fasta file
            other_seq = "/Users/rbhattac/Desktop/sequence1.fasta"
            file1 = open(other_seq, "w+")
            file1.write(">"+str(source)+"\n")
            file1.write(str(source_seq))
            file1.close()

            #align the source sequence with the reference sequence
            output_file = "/Users/rbhattac/Desktop/output_test.fasta"
            print("starting to calculate the alignment...")
            align_sequences(ref_seq, other_seq, output_file)

            #read the alignment file
            alignment = AlignIO.read(output_file, "fasta")
            #
            # #set the source seq coordinates in the alignment
            # seq_source = alignment[1].seq

            #extract and parse the source coordinates
            try:
                coordinates = parse_source_coordinates(source_coordinates)
            except Exception as e:
                print(f"Error parsing source coordinates: {e}")
                continue

            #get the reference position
            ref_pos = get_ref_pos(alignment, 1, coordinates)

            #from the ref_pos extract start and end coordinates
            segment = ""
            start = ""
            end = ""
            SF_seq = ""
            new_list = [x + 1 for x in ref_pos]
            if len(new_list) == 2:
                start = str(new_list[0])
                end = str(new_list[-1])
                segment = f"{start}-{end}"
            elif len(new_list) > 2:
                start = str(new_list[0])
                end = str(new_list[-1])
                segment_temp = ','.join(map(str, new_list))
                segment = '"' + segment_temp + '"'
            elif len(new_list) == 1:
                start = str(new_list[0])
                end = str(new_list[0])
                segment = f"{start}-{end}"

            #extract the SF sequence from the alignment
            SF_seq=""
            ref_SF=""

            #check if length of new_list  is greater than 0
            if len(new_list)> 0:

                if "-" in source_coordinates:
                    SF1_seq = alignment[:, (new_list[0] - 1):new_list[-1]]
                    print(SF1_seq)
                    SF_seq=SF1_seq[1].seq

                elif "," in coordinates:
                    # extract each coordinate from alignment and concatenate them together
                    for i in new_list:
                        SF1_temp = alignment[:, (i - 1):i]
                        if i == new_list[0]:
                            SF1_seq = SF1_temp
                        else:
                            SF1_seq += SF1_temp
                        print(SF1_seq)
                        #ref_SF = SF1_seq[0].seq
                        SF_seq = SF1_seq[1].seq





            #segment = ','.join(map(str, new_list))
            print(f"Start: {start}, End: {end}")
            print(f"Segment: {segment}")
            print(f"SF_seq: {SF_seq}")
            print(f"Reference SF: {SF1_seq[0].seq}")

            #extract the sequence

            #write to the csv file
            csv_writer.writerow([source, source_coordinates, start, end, segment, SF_seq, SF1_seq[0].seq, SF1_seq])
            #csv_writer.writerow([source, source_coordinates, start, end, segment, SF_seq])
    #         row['start'] = start
    #         row['end'] = end
    #         row['segments'] = segment
    #         #write dataframe to a csv file
    # SF_Table.to_csv('/Users/rbhattac/Desktop/out.csv', index=False)






if __name__ == "__main__":
    main()
