

# This script will parse through the output of a Clustal Omega nt alignment 
# (http://www.ebi.ac.uk/Tools/msa/clustalo/) and yields a matrix for 
# base frequency at each position. The rows are A-T-C-G while the columns
# are each base in the sequence. 
#
# This script is meant for a single MSA alignment for analyzing primers so it
# assumes it can fit the file into memory. 
#
# Run the script using a command like this:
# python3 analyze_clustal_omega.py -input_file clustal_omega.out -seq_length 15 -output_file out.matrix
#
# Author: James Matsumura

import argparse

def main():

    parser = argparse.ArgumentParser(description='Script to analyze multiple outputs from Primer Search so that different primer sets can be compared.')
    parser.add_argument('-input_file', type=str, required=True, help='Path to the output from primersearch_to_tsv.py.')
    parser.add_argument('-output_file', type=str, required=True, help='Name/location of the output file.')
    args = parser.parse_args()

    nt_order = ["A","T","C","G","-"] # guarantee this for the output
    final_list = []
    num_of_entries = 0 # need to keep track of total base possibilities per position
 
    # Parse the input to build a string to represent each base position
    with open(args.input_file,'r') as co_in:

        alignments = co_in.readlines()[3:-1]
        num_of_entries = len(alignments)
        # Initialize the final list with one element per base being aligned
        final_list = [""] * len(alignments[0].split()[1].strip())
        
        for entry in alignments: 
            entry = entry.strip()
            seq = entry.split()[1]

            for idx in range(0,len(seq)):
                final_list[idx] += seq[idx]

    # Now do counts / ratios for each base at each position and print
    counts,ratios = ([] for i in range(2))

    with open(args.output_file,'w') as out:

        for x in range(0,len(final_list)):

            out.write("\t{0}".format(x+1)) # don't use python's 0-based indexing for the output
            pos = final_list[x]

            for base in nt_order:
                counts.append(pos.count(base))
                ratios.append("{0:.2f}".format(pos.count(base)/num_of_entries*100))

        out.write("\n")

        # Reverse the last loop to print by base X pos
        for j in range(0,len(nt_order)):

            base = nt_order[j]

            out.write("{0}".format(base))
            
            # The current data structure has each base count per position every n
            # steps away depending on the length of nt_order. 
            for x in range(j,len(counts),len(nt_order)):

                out.write("\t{0} ({1})".format(ratios[x],counts[x]))

            out.write("\n")


if __name__ == '__main__':
    main()