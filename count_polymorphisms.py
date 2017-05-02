

# This script will parse through CSV files to capture how many polymorphisms,
# per individual, were found in each peptide. This will be presented as either
# an absolute count of how many poly AAs are present or a ratio for how many
# AAs matched to the reference. 
#
# Run the script using a command like this:
# python3 count_polymorphisms.py -peptide_map peptide_map.csv -sample_map sample_map.csv -rel_or_abs rel -outfile outfile.csv
#
# Author: James Matsumura

import pandas,argparse

def main():

    parser = argparse.ArgumentParser(description='Script to count polymorphisms given a peptide map and sample data.')
    parser.add_argument('-peptide_map', type=str, required=True, help='Path to peptide file.')
    parser.add_argument('-sample_map', type=str, required=True, help='Path to sample file.')
    parser.add_argument('-rel_or_abs', type=str, required=True, help='Whether to output a matrix of relative match % or absolute counts for how many polymorphisms in that peptide.')
    parser.add_argument('-outfile', type=str, required=True, help='Name of the output file.')
    args = parser.parse_args()

    # row names are peptides, col names are sample IDs by their Subj+Timepoint
    rownames,colnames = ([] for i in range(2))

    peptide_df = pandas.read_csv(args.peptide_map) # open up the peptide file
    peptide_dict,sample_dict = ({} for i in range(2))

    for index,row in peptide_df.iterrows():

        peptide_dict[row["Name"]] = {
            'name': row["Name"],
            'seq': row["Sequence"],
            'seq_s': row["StartAA"],
            'seq_e': row["EndAA"]
        }

        rownames.append(row["Name"]) # maintain dict order via this

    sample_df = pandas.read_csv(args.sample_map)
    # want these ordered so that we can build the output in order
    poly_aa_pos = [] 

    for col in list(sample_df):
        if col.startswith("AA"):
            # this tells us which polymorphic AA positions we will check 
            # against each peptide fragment
            poly_aa_pos.append(int(col.split('AA')[1]))

    # each sample will have a dif set of fragments we want to check, note which
    relevant_peptides = set() 
    for fragment in rownames:
        start = int(peptide_dict[fragment]['seq_s'])
        end = int(peptide_dict[fragment]['seq_e'])

        for pos in poly_aa_pos:
            if start <= pos <= end:
                relevant_peptides.add(fragment)
                break

    # Build a list where each element will be a list that represents the 
    # sample v peptideS alignment variables. 
    final_matrix = [] 

    for index,row in sample_df.iterrows():

        unique_id = "{0}.{1}".format(row["Subj"],row["Timepoint"])
        colnames.append(unique_id)

        curr_sample = [] # a list of all the values found for this sample

        for fragment in rownames:
            # If this fragment has no poly AAs, then give it a perfect score
            if fragment not in relevant_peptides: 
                if args.rel_or_abs == 'rel':
                    curr_sample.append("1.00000")
                elif args.rel_or_abs == 'abs':
                    curr_sample.append("0")

            # Time to inspect whether or not the polymorphisms are present
            else:
                relevant_positions = []
                num_of_poly = 0
                seq = peptide_dict[fragment]['seq']
                start = int(peptide_dict[fragment]['seq_s'])
                end = int(peptide_dict[fragment]['seq_e'])

                for pos in poly_aa_pos:
                    if start <= pos <= end:
                        relevant_positions.append(pos)
                    elif end < int(pos): # got em all, do a check
                        break

                if len(relevant_positions) > 0:
                    for pos in relevant_positions:

                        pos_with_aa = "AA{0}".format(pos)

                        # need to offset based on where the fragment starts
                        adj_pos = pos - start

                        if seq[adj_pos] != row[pos_with_aa]:
                            num_of_poly += 1

                if args.rel_or_abs == 'rel':
                    seq_len = int(len(seq))
                    num_match = seq_len - num_of_poly
                    curr_sample.append("{0:.5f}".format(num_match/seq_len))

                elif args.rel_or_abs == 'abs':
                    curr_sample.append("{0}".format(num_of_poly))

        # now we have all the data for this sample column
        final_matrix.append(curr_sample) 

    rows = pandas.Index(rownames,name="rows")
    cols = pandas.Index(colnames,name="columns")
    df = pandas.DataFrame(data=final_matrix,index=cols,columns=rows)
    df = df.transpose()
    df.to_csv(args.outfile,sep='\t')


if __name__ == '__main__':
    main()