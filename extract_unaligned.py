

# This script will parse through the output of primersearch_to_tsv.py
# and extract those sequences which do not meet the set parameters of
# mismatch/amp length. This can output either DNA or Protein sequences 
# depending on the option passed.
#
# The output will be a FASTA file. 
#
# Run the script using a command like this:
# python3 extract_unaligned.py -input /path/to/primersearch_to_tsv.out.tsv -min_match 0 -max_match 3 -min_amp 50 -max_amp 3000 -nt_or_aa aa -reference /path/to/ref.fsa -output_fsa /path/to/out_fsa
#
# Author: James Matsumura

import argparse
from Bio import SeqIO

def main():

    parser = argparse.ArgumentParser(description='Script to extract those sequences which are not being hit by the primer set.')
    parser.add_argument('-input', type=str, required=True, help='Path to the output from primersearch_to_tsv.py.')
    parser.add_argument('-reference', type=str, required=True, help='Reference FASTA sequences provided for the -seqall option of Primer Search.')
    parser.add_argument('-min_match', type=int, required=True, help='Minimum number of mismatches for either F or R alignment.')
    parser.add_argument('-max_match', type=int, required=True, help='Maximum number of mismatches for either F or R alignment.')
    parser.add_argument('-min_amp', type=int, required=True, help='Minimum valid amplimer length.')
    parser.add_argument('-max_amp', type=int, required=True, help='Maximum valid amplimer length.')
    parser.add_argument('-nt_or_aa', type=str, required=True, help='Nucleotides or Amino Acids for the output.')
    parser.add_argument('-output', type=str, required=True, help='Location for the output F/R FASTA files.')
    args = parser.parse_args()
 
    aligned = set()

    # Store the FASTA file in memory
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.reference,"fasta"))

    # Find all those sequences that DO meet the criteria
    with open(args.input,'r') as results:
        for line in results:
            line = line.strip()
            elements = line.split('\t')
            id = elements[1]
            amp = int(elements[2].split(' ')[0])
            f_mm = int(elements[3])
            r_mm = int(elements[4])

            # Pass if either the F or R alignment has mismatches outside the desired range
            if not (args.min_match <= f_mm <= args.max_match):
                continue
            elif not (args.min_match <= r_mm <= args.max_match):
                continue
            elif not (args.min_amp <= amp <= args.max_amp):
                continue # skip any outside the desired amplimer range

            # We only care whether or not this sequence had a valid target region
            aligned.add(id)

    # Now only write out those that didn't make it
    with open(args.output,'w') as out:
        for entry in seq_dict:
            if seq_dict[entry].id not in aligned:
                if args.n_or_aa == "aa":
                    seq_dict[entry].seq = seq_dict[entry].seq.translate()
                    SeqIO.write(seq_dict[entry],out,"fasta")
                else:
                    SeqIO.write(seq_dict[entry],out,"fasta")


if __name__ == '__main__':
    main()