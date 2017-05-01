

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

    # row names are peptides, col names are sample IDs by their Subj+Timepoint
    rownames,colnames = ([] for i in range(2))

    peptide_df = pandas.read_csv(args.peptide_map) # open up the peptide file
    peptide_dict = {}

    for index,row in peptide_df.iterrows():

        peptide_dict[row["Name"]] = {
            'name': row["Name"],
            'seq': row["Sequence"],
            'seq_s': row["StartAA"],
            'seq_e': row["EndAA"]
        }

        rownames.append(row["Name"]) # maintain dict order via this

    samples = args.sample_maps.split(',')
    sample_dict = {}
    for sample in samples:
        sample_df = pandas.read_csv(sample)
        # want these ordered so that we can build the output in order
        poly_aa_pos = [] 

        for col in list(sample_df):
            if col.startswith("AA"):
                # this tells us which polymorphic AA positions we will check 
                # against each peptide fragment
                poly_aa_pos.append(col.split('AA')[1])

        # each sample will have a dif set of fragments we want to check, note which
        relevant_peptides = set() 
        for fragment in rownames:
            start = int(peptide_dict[fragment]['seq_s'])
            end = int(peptide_dict[fragment]['seq_e'])

            for pos in poly_aa_pos:
                if start <= int(pos) <= end:
                    relevant_peptides.add(fragment)
                    break

        # Build a list where each element will be a list that represents the 
        # sample v peptideS alignment variables. 
        final_matrix = [] 

        for index,row in sample_df.iterrows():

            curr_sample = []

            unique_id = "{0}.{1}".format(row["Subj"],row["Timepoint"])
            colnames.append(unique_id)

            for fragment in rownames:
                # If this fragment has no poly AAs, then give it a perfect score
                if fragment not in relevant_peptides: 
                    curr_sample.append("1.000")

                # Time to inspect whether or not the polymorphisms are present
                else:
                    relevant_positions = []
                    for pos in poly_aa_pos:
                        if start <= int(pos) <= end:
                            relevant_positions.append(pos)
                        elif end < int(pos):
                            break



if __name__ == '__main__':
    main()