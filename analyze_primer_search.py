

# This script will parse through numerous outputs from EMBOSS's Primer Search. 
# These files must have been formatted by primersearch_to_tsv.py in order to 
# work. The input for this script is a list file where each line locates an
# individual *.tsv file from primersearch_to_tsv.py as well as the FASTA 
# reference file used to align to by Primer Search (the -seqall option).
#
# The output will be a matrix 
#
# Run the script using a command like this:
# python3 analyze_primer_search.py -list /path/to/list_of_tsvs -reference /path/to/ref.fsa -mismatch 3 -min_amp 50 -max_amp 3000 -output_file /path/to/out_matrix.tsv
#
# Author: James Matsumura

import argparse
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to analyze multiple outputs from Primer Search so that different primer sets can be compared.')
    parser.add_argument('-list', type=str, required=True, help='Path to a list file where each line denotes a path to an output file from primersearch_to_tsv.py.')
    parser.add_argument('-reference', type=str, required=True, help='Reference FASTA sequences provided for the -seqall option of Primer Search.')
    parser.add_argument('-mismatch', type=int, required=True, help='Maximum number of mismatches for one of the primers.')
    parser.add_argument('-min_amp', type=int, required=True, help='Minimum valid amplimer length.')
    parser.add_argument('-max_amp', type=int, required=True, help='Maximum valid amplimer length.')
    parser.add_argument('-output_file', type=str, required=True, help='Location for the output file.')
    args = parser.parse_args()
 
    valid_pair_dict,valid_amp_dict = (defaultdict(list) for i in range(2)) # new list with excessive mismatches removed
    ps_out,order = ([] for i in range(2))

    # Iterate over the list file to isolate each individual primer search output file
    with open(args.list,'r') as ps_outputs:
        for file in ps_outputs:
            file = file.rstrip()
            ps_out.append(file)

    # Establish the maximum mismatch, by default the worst value will be 2*mismatch+1
    # since we should not consider those alignments with mismatch values larger than 
    # specified. Thus, a "no_match" value will be set dependent on the input mismatch 
    # amount. 
    no_match = args.mismatch*2+1
    print("A value of {0} in the output matrix is considered to not align well enough and should not be considered a match.".format(no_match))

    # Now grab all the IDs possible for the primers to align to 
    with open(args.reference,'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                id = line.split(' ')[0][1:]
                valid_pair_dict[id] = [no_match] * len(ps_out)
                valid_amp_dict[id] = [0] * len(ps_out)
                order.append(id) # establish ordering using the FASTA file

    # Now identify the best possible hits for each set of primers
    for file in ps_out:

        index = ps_out.index(file)

        with open(file,'r') as results:
            for line in results:
                elements = line.split('\t')
                id = elements[1]
                amp = int(elements[2].split(' ')[0])
                f_mm = int(elements[3])
                r_mm = int(elements[4])

                if f_mm > args.mismatch or r_mm > args.mismatch:
                    continue # skip any with too many mismatches
                elif amp < args.min_amp or amp > args.max_amp:
                    continue # skip any outside the desired amplimer range

                # If we've made it here, there are not too many mismatches.
                # Now we want to establish the "best" case which will maintain
                # the lowest total number of mismatches for a given primer search
                # hit (e.g. F has 1 mismatch and R has 1 mismatch will be 
                # considered a better result than if just F or R have 3 
                # mismatches with the other having 0).

                previous_best = valid_pair_dict[id][index]
                curr_val = f_mm+r_mm

                if curr_val < previous_best:
                    valid_pair_dict[id][index] = curr_val
                    valid_amp_dict[id][index] = amp

    # Reassign all the names of the input files to something more legible by
    # eliminating the bulk of the path
    for x in ps_out:
        filename = x.split('/')[-1]
        ps_out[ps_out.index(x)] = filename

    # Have all the data we need, build a final output matrix
    with open(args.output_file,'w') as out:
        # Write out the column header which denotes the individual file name
        for primer_set in ps_out:
            out.write("\t{0}".format("{0} total mismatch".format(primer_set)))
        for primer_set in ps_out:
            out.write("\t{0}".format("{0} amplimer".format(primer_set)))
        out.write("\n")
        # Each row will be an ID from the reference
        for id in order:
            out.write("{0}\t{1}\t".format(id,("\t").join(map(str,valid_pair_dict[id]))))
            out.write("{0}\n".format(("\t").join(map(str,valid_amp_dict[id]))))
                

if __name__ == '__main__':
    main()