

# This script will parse through CSV files to capture how many polymorphisms,
# per individual, were found in each peptide. This will be presented as both 
# an aboslute count as well as a ratio. 
#
# Run the script using a command like this:
# python3 count_polymorphisms.py -peptide_map peptide_map.csv -sample_map sample_map.csv
#
# Author: James Matsumura

import pandas,argparse

def main():

    parser = argparse.ArgumentParser(description='Script to count polymorphisms given a peptide map and sample data.')
    parser.add_argument('-peptide_map', type=str, required=True, help='Path to peptide file.')
    parser.add_argument('-sample_maps', type=str, required=True, help='Comma separated list of sample data.')
    args = parser.parse_args()

    possible_positions = set() # store each valid possible AA polymorphic position

    samples = args.sample_maps.split(',')

    # First do a quick pull from the column names of the sample data files and
    # find the valid PolyAA positions.
    for sample in samples:
        sample_df = pandas.read_csv(sample)
        for col in list(sample_df): # extract the column names
            if col.startswith("AA"):
                possible_positions.add(col[2:])

    peptide_df = pandas.read_csv(args.peptide_map) # open up the peptide file

    # Now build a dictionary which contains a dictionary. This is so we can map
    # a peptide ID to a list of its possible polymorphisms.
    peptide_dict = {} # will contain the peptide name, a list of all its polymorphic 
                      # AAs+positions, and the length of the peptide.
    for index,row in peptide_df.iterrows():

        # Only keep info on those peptides with polymorphic AAs
        if row["Polymorphic AA only"] != "":

            start_end = "{0}:{1}".format(row["StartAA"],row["EndAA"])

            peptide_dict[start_end] = {
                'name': row["Name"],
                'length': row["#polyAA"],
                'poly_aas': extract_polymorphic_aas(row["StartAA"],row["EndAA"],row["Polly AA Start"],row["Polly AA End"],row["Sequence"],row["Polymorphic AA only"])
            }

# Function which takes in a range for the peptide sequence, range for the 
# polymorphic AAs, the peptide sequence, the polypeptides, and a set for
# possible positions for these poly AAs to reside.
def extract_polymorphic_aas(pep_start,pep_stop,aa_start,aa_stop,aas,poly_aas):


if __name__ == '__main__':
    main()