

# This script will parse through the output of primersearch_to_tsv.py
# and extract sequences from the original FASTA entries that represent
# where these primers have aligned. 
#
# The output will be a FASTA file with F/R primer alignment sections. 
#
# Run the script using a command like this:
# python3 extract_primer_regions.py -input /path/to/primersearch_to_tsv.out.tsv -buffer 5 -min_match 0 -max_match 3 -min_amp 50 -max_amp 3000 -reference /path/to/ref.fsa -output_fsa /path/to/out_fsa
#
# Author: James Matsumura

import argparse

def main():

    parser = argparse.ArgumentParser(description='Script to analyze multiple outputs from Primer Search so that different primer sets can be compared.')
    parser.add_argument('-input', type=str, required=True, help='Path to the output from primersearch_to_tsv.py.')
    parser.add_argument('-reference', type=str, required=True, help='Reference FASTA sequences provided for the -seqall option of Primer Search.')
    parser.add_argument('-min_match', type=int, required=True, help='Minimum number of total mismatches.')
    parser.add_argument('-max_match', type=int, required=True, help='Maximum number of total mismatches.')
    parser.add_argument('-min_amp', type=int, required=True, help='Minimum valid amplimer length.')
    parser.add_argument('-max_amp', type=int, required=True, help='Maximum valid amplimer length.')
    parser.add_argument('-f_len', type=int, required=True, help='Length of forward primer.')
    parser.add_argument('-r_len', type=int, required=True, help='Length of reverse primer.')
    parser.add_argument('-buffer', type=int, required=True, help='How many bases to add to each end of the primer aligned region.')
    parser.add_argument('-output_dir', type=str, required=True, help='Location for the output F/R FASTA files.')
    args = parser.parse_args()
 
    valid_pair_dict,valid_pos_dict,fasta_dict = ({} for i in range(3))
    seq,id = ("" for i in range(2))
    order = []

    # Now grab all the IDs possible for the primers to align to 
    with open(args.reference,'r') as fasta:
        for line in fasta:

            line = line.strip()

            if line.startswith('>'):

                if len(seq) > 0: # have a sequence add to dict
                    fasta_dict[id] = seq
                    seq = "" # reinit for next seq

                id = line.split(' ')[0][1:]
                valid_pair_dict[id] = 9999 # arbitrarily high max
                valid_pos_dict[id] = ""
                order.append(id) # establish ordering using the FASTA file  

            else:
                seq += line

        fasta_dict[id] = seq

    with open(args.input,'r') as results:
        for line in results:
            line = line.strip()
            elements = line.split('\t')
            id = elements[1]
            amp = int(elements[2].split(' ')[0])
            f_mm = int(elements[3])
            r_mm = int(elements[4])
            tot_mm = f_mm+r_mm

            if not (args.min_match <= tot_mm <= args.max_match):
                continue # skip any not within specified range
            elif not (args.min_amp <= amp <= args.max_amp):
                continue # skip any outside the desired amplimer range

            # Now we want to establish the "best" case which will maintain
            # the lowest total number of mismatches for a given primer search
            # hit (e.g. F has 1 mismatch and R has 1 mismatch will be 
            # considered a better result than if just F or R have 3 
            # mismatches with the other having 0).

            previous_best = valid_pair_dict[id]
            curr_val = tot_mm

            if curr_val < previous_best: # found a new best
                valid_pair_dict[id] = curr_val
                valid_pos_dict[id] = "{0}:{1}".format(elements[5],elements[6][1:-1])

    forward = "{0}/forward.fsa".format(args.output_dir)
    reverse = "{0}/reverse.fsa".format(args.output_dir)

    with open(forward,'w') as f_out:
        with open(reverse,'w') as r_out:
            for id in order:
                if valid_pos_dict[id] != "": # only want those that found valid positions
                    seq_len = len(fasta_dict[id])

                    # Grab the ranges for primer extraction
                    f_start_pos,f_end_pos,r_start_pos,r_end_pos = (0 for i in range(4))

                    f_pos = int(valid_pos_dict[id].split(':')[0]) - 1 # correct for python 0-based indexing
                    r_pos = int(valid_pos_dict[id].split(':')[1]) - 1

                    if (f_pos - args.buffer) < 0:
                        print("Warning, buffer has reached the end of sequence {0}".format(id))
                        f_start_pos = 0
                    else:
                        f_start_pos = (f_pos - args.buffer)

                    f_end_pos = (f_pos + args.buffer + args.f_len)

                    if (r_pos + args.buffer) > seq_len:
                        print("Warning, buffer has reached the end of sequence {0}".format(id))
                        r_pos = seq_len
                    else:
                        r_end_pos = (seq_len - r_pos + args.buffer)

                    r_start_pos = (seq_len - r_pos - args.buffer - args.r_len)

                    f_out.write(">{0}\n{1}\n".format(id,fasta_dict[id][f_start_pos:f_end_pos]))
                    r_out.write(">{0}\n{1}\n".format(id,fasta_dict[id][r_start_pos:r_end_pos]))


if __name__ == '__main__':
    main()