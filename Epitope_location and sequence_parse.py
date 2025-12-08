#write a function which reads a column of csv file and outputs another csv file with columns after running the function
#input: an csv file with a column of epitope coordinates
#output: an csv file with three columns, one with the residue names and the other with the residue numbers and another with length
#import openpyxl

#pip install openpyxl
import pandas as pd
import csv
def parse_epitope_coordinates_file(input_file_path):


    # Read the input file
    epitope_df = pd.read_csv(input_file_path)
    epitopes = epitope_df['Epitope'].tolist()
    #print(epitope_df)
    print(epitopes)
    return epitopes




#write a function to parse a discontinuos epitope coordinate file
#input: a text line from the file
#output: a comma separated string of the start and end coordinates separated by tab comma separated sequences
#input example: "G104, E126, W231"
#output example: "G,E,W" "104,126,231"

def parse_epitope_coordinates(line):
    # Split the input line by commas
    residues = line.split(', ')
    #print(residues)
    # Extract the residue names and numbers
    residue_names = [residue[0] for residue in residues]
    residue_numbers = [residue[1:] for residue in residues]
    print(residue_names)
    print(residue_numbers)

    # Join the residue names and numbers into comma-separated strings
    #residue_names_str = ','.join(residue_names)
    #residue_numbers_str = ','.join(str(residue_numbers))
    #print(residue_names_str)
    #print(residue_numbers_str)
    #return residue_names_str, residue_numbers_str
    return residue_names, residue_numbers

# Test the function with an example input
#example_line = "G304, K305, F306, K307, V308, V309, K310, I312, R323, Q325, P364, K388, L389, N390"
#print(parse_epitope_coordinates(example_line))  # Expected output: ("G,E,W", "104,126,231")

def main():
    # Parse the epitope coordinates from the input file
    epitopes = parse_epitope_coordinates_file("/Users/rbhattac/Desktop/Sars_epitope.csv")
    #print(epitopes)

    #open the output file to write
    with open("/Users/rbhattac/Desktop/Sars_parsed_epitope.csv", "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write header
        csv_writer.writerow(["Sequence", "Positions", "Length"])
        for line in epitopes:
            residue_names, residue_numbers = parse_epitope_coordinates(line)

            # Join the residue names and numbers into comma-separated strings
            residue_names_str = ", ".join(residue_names)
            residue_numbers_str = ", ".join(residue_numbers)
            #calculate how many residues are there
            length_epitope = len(residue_names)

            # Write data to CSV
            csv_writer.writerow([residue_names_str, residue_numbers_str, length_epitope])

    # Write the parsed epitope coordinates to the output file
    #parse_epitope_coordinates_file(epitopes, "/Users/rbhattac/Desktop/parsed_epitope.txt")
#Test the function with an example input
#input_file_path = "/Users/rbhattac/Desktop/epitope.csv"
#output_file_path = "/Users/rbhattac/Desktop/parsed_epitope.txt"
#parse_epitope_coordinates_file(input_file_path, output_file_path)

if __name__ == "__main__":
    main()