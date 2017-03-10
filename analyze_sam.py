

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
    ref_order = [] # establish an order for the final list based on reference order

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
                ref_order.append(header)

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
    valid_pair_keys = defaultdict(set) # use this to isolate best alignment per reference
    primer_order = set()
    for reference in ref_dict:
        for f,r in read_pairs(ref_dict[reference]):

            if int(f.split('::')[1]) > max_mm:
                continue
            elif int(r.split('::')[1]) > max_mm:
                continue
            else:
                
                # If this is the first of this primer for this reference, arbitrarily
                # consider it to be the best for now. 
                if f.split('::')[0] not in valid_pair_keys:
                    valid_pair_keys[reference].add(f.split('::')[0]) # grab the primer
                    valid_pair_dict[reference].append("{0}\t{1}".format(f,r))
                    primer_order.add(f.split('::')[0])
                    primer_order.add(r.split('::')[0])

                # If a primer pair already exists, see if this set has a lower
                # total mismatch amount and perhaps is a better representation
                # of the alignment. 
                else:
                    new_mm = (int(f.split('::')[1])+int(r.split('::')[1]))

                    copy_vp_dict = valid_pair_dict # need a copy so we can modify original

                    for pair in copy_vp_dict[reference]:
                        if pair.startswith(f[:-1]): # found the pair to compare to
                            fwd = pair.split('\t')[0]
                            rev = pair.split('\t')[1] 
                            established_mm = (int(fwd.split('::')[1])+int(rev.split('::')[1]))

                            # If we've found a pair with less total mismatches, make
                            # these the representative.
                            if new_mm < established_mm:
                                idx = copy_vp_dict[reference].index(pair)
                                valid_pair_dict[reference].pop(idx)
                                valid_pair_dict[reference].append("{0}\t{1}".format(f,r))

    primer_idxs = sorted(primer_order)
    final_dict = defaultdict(list)
    for reference in ref_order: # we want all references even those that did not align
        # initialize however many slots there are for the number of different primers
        final_dict[reference] = ['-1']*(len(primer_idxs)) 
        
        # Grab all relevant pair mismatch values for this reference
        for pairs in valid_pair_dict[reference]:

            split_pairs = pairs.split('\t')
            f_ref = split_pairs[0].split('::')[0]
            f_mm = split_pairs[0].split('::')[1]
            r_ref = split_pairs[1].split('::')[0]
            r_mm = split_pairs[1].split('::')[1]
            f_idx = primer_idxs.index(f_ref)
            r_idx = primer_idxs.index(r_ref)

            final_dict[reference][f_idx] = str(f_mm)
            final_dict[reference][r_idx] = str(r_mm)

    with open(args.output_file,'w') as out:
        out.write("\t{0}\n".format(("\t").join(primer_idxs))) # header is primer IDs

        # subsequent lines will be the reference and all the MM values for those primers
        for reference in ref_order:
            out.write("{0}\t{1}\n".format(reference,("\t".join(final_dict[reference]))))


def read_pairs(iterable):
    a = iter(iterable)
    return zip(a,a)


if __name__ == '__main__':
    main()