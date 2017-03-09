

# This script will parse through a SAM file resulting from a paired-end 
# Bowtie2 alignment for primers and result in a TSV matrix for references (Y) 
# and primers (X). The user can specify the maximum amount of mismatches 
# allowed.
#
# Run the script using a command like this:
# python3 analyze_sam.py -sam /path/to/bowtie_out.sam -mismatch 3 -output_file /path/to/out.tsv
#
# Author: James Matsumura

import argparse,pysam
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to isolate all reads and where they aligned to given a BAM file.')
    parser.add_argument('-sam', type=str, required=True, help='Path to a SAM file derived from Bowtie2.')
    parser.add_argument('-mismatch', type=int, required=False, help='Maximum number of mismatches allowed to be considered a valid alignment. Defaults to 0.')
    parser.add_argument('-output_file', type=str, required=True, help='Location for the output file.')
    args = parser.parse_args()
 
    ref_dict = {} # capture all reference sequences

    input = pysam.AlignmentFile(args.sam,'r')

    # First grab all the headers. Probably some way to grab this via pysam
    # but it's not readily clear.
    with open(args.sam,'r') as header:
        for line in header:
            if line.startswith('@PG'): # leave if done getting headers
                break

            elif line.startswith('@SQ'):
                elements = line.split('\t')
                header = elements[1].split(':')[1]
                ref_dict[header] = []

    max_mm = 0  # maximum number of mismatches allowed for an alignment
    if args.mismatch: # user specificying threshold
        max_mm = args.mismatch

    for alignment in input.fetch(until_eof=True): # iterate over all reads in the SAM file.

        current_mm = alignment.get_tag('XM') # grab number of mismatches
        primer = alignment.query_name # grab primer/query name
        ref = alignment.reference_name # grab reference sequence name

        ref_dict[ref].append("{0}::{1}".format(primer,current_mm))

    input.close()

    # At this point we have a reference dictionary where each reference 
    # either has an empty list or a list that is comprised of pairs of 
    # elements (one element for each read). 

    # First iterate over the pairs to get rid of any that exceed the 
    # specified mismatch value. 
    valid_pair_dict = defaultdict(list) # new list with excessive mismatches removed
    valid_pair_keys = defaultdict(set) # use this to isolate best laignment per reference
    for reference in ref_dict:
        for f,r in read_pairs(ref_dict[reference]):

            if int(f.split('::')[1]) > max_mm:
                continue
            elif int(r.split('::')[1]) > max_mm:
                continue
            else:
                valid_pair_dict[reference].append("{0}\t{1}".format(f,r))
                valid_pair_keys[reference].add(f[:-1]) # grab the primer


def read_pairs(iterable):
    a = iter(iterable)
    return zip(a,a)


if __name__ == '__main__':
    main()