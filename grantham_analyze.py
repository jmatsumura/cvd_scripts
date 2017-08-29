#! /usr/bin/env python3

# This script will calculate the GV score of all peptides present in a file
# slicing at the specified start/stop positions (entered as 1-based index).
# The peptides are assumed to all be of the same length and already aligned.
#
# Run the script using a command like this:
# grantham_analyze.py -p peptides.csv -b 30 -e 50 -baa
#
# Author: James Matsumura

import argparse
from collections import defaultdict, OrderedDict

def main():

    parser = argparse.ArgumentParser(description='Script to calculate the GV score for peptides.')
    
    parser.add_argument('-p','--peptides', type=str, required=True, help='Path to a file containing a list of peptides.')
    parser.add_argument('-b','--beginning', type=int, default=1, required=False, help='Beginning position of sequence (first AA is considered position 1), defaults to beginning of the peptide.')
    parser.add_argument('-e','--end', type=int, default=1000000, required=False, help='Ending position of sequence (first AA is considered position 1), defaults to end of the peptide.')
    parser.add_argument('-o','--outfile', type=str, required=True, help='Name/location of the file to write out to.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-ap','--all_pairwise', action='store_true', help='Compute GD score for the entire sequence of all possible pairs, all GD values are added to a total sum then divided by total length.')
    group.add_argument('-sp','--single_pairwise', action='store_true', help='Obtain the most extreme GD score found at a single position for the entire sequence of all possible pairs.')
    group.add_argument('-baa','--by_amino_acid', action='store_true', help='Obtain the most extreme GV score found at a single position across ALL sequences, this is not pairwise.')
    args = parser.parse_args()

    # Establish the values from Table 1 of Grantham's work: https://www.ncbi.nlm.nih.gov/pubmed/4843792
    gv_vals = {
        'A': AminoAcid(0.0,8.1,31.0),
        'I': AminoAcid(0.0,5.2,111.0),
        'L': AminoAcid(0.0,4.9,111.0),
        'V': AminoAcid(0.0,5.9,84.0),
        'F': AminoAcid(0.0,5.2,132.0),
        'W': AminoAcid(0.13,5.4,170.0),
        'Y': AminoAcid(0.20,6.2,136.0),
        'N': AminoAcid(1.33,11.6,56.0),
        'C': AminoAcid(2.75,5.5,55.0),
        'Q': AminoAcid(0.89,10.5,85.0),
        'M': AminoAcid(0.0,5.7,105.0),
        'S': AminoAcid(1.42,9.2,32.0),
        'T': AminoAcid(0.71,8.6,61.0),
        'R': AminoAcid(0.65,10.5,124.0),
        'H': AminoAcid(0.58,10.4,96.0),
        'K': AminoAcid(0.33,11.3,119.0),
        'D': AminoAcid(1.38,13.0,54.0),
        'E': AminoAcid(0.92,12.3,83.0),
        'G': AminoAcid(0.74,9.0,3.0),
        'P': AminoAcid(0.39,8.0,32.5)
    }

    if args.beginning == 1 and args.end == 1000000:
        print("Calculating GV scores for each pair using the entire length of the peptide sequence.")
    else:
        print("Calculating GV scores for each pair between AAs {} and {} of the peptide sequence.".format(args.beginning,args.end))
 
    # grab all the peptides, use the ID in the input as the ID in the hash 
    peptides = OrderedDict()
    with open(args.peptides,'r') as infile:
        for line in infile:
            peptides[line.split(',')[0]] = line.split(',')[1].strip()[(args.beginning-1):(args.end)]
    
    key_lookup = list(peptides.keys()) # based on index, get the key for this dict

    if args.all_pairwise: # pairwise comparisons

        # store the results of all the comparisons in a giant list, order is 
        # guaranteed to be 1v2,1v3,1v4....2v3,2v4,2v5....(n-1)v(n)
        results = [] 

        cur_base_sample_pos = 0 # start at 0 and work up

        for pep1 in range(0,len(peptides)):

            results.append('BREAKPOINT') # break on this for outputting

            pep1_seq = peptides[key_lookup[pep1]] 

            for pep2 in range(cur_base_sample_pos+1,len(peptides)): # never repeat a comparison
                
                pep2_seq = peptides[key_lookup[pep2]]

                results.append(pairwise_gv(gv_vals,pep1_seq,pep2_seq))

            cur_base_sample_pos += 1

    elif args.by_amino_acid: # find worst case at a position and note it

        unique_sets = defaultdict(set)

        results = []

        for j in range(0,len(peptides[key_lookup[0]])):
            unique_sets[j] = set([ele[j] for ele in peptides.values()])

        for j in range(0,len(peptides[key_lookup[0]])):
            if len(unique_sets[j]) == 1: # only one AA found in this position
                results.append('0.0')
            else:
                cmin,pmin,vmin = (1000.0 for i in range(3))
                cmax,pmax,vmax = (0.0 for i in range(3))

                for x in unique_sets[j]:

                    if gv_vals[x].c < cmin:
                        cmin = gv_vals[x].c
                    if gv_vals[x].c > cmax:
                        cmax = gv_vals[x].c

                    if gv_vals[x].p < pmin:
                        pmin = gv_vals[x].p
                    if gv_vals[x].p > pmax:
                        pmax = gv_vals[x].p

                    if gv_vals[x].v < vmin:
                        vmin = gv_vals[x].v
                    if gv_vals[x].v > vmax:
                        vmax = gv_vals[x].v

                results.append(single_gv(cmax,cmin,pmax,pmin,vmax,vmin))

        with open(args.outfile,'w') as out:
            out.write("\t".join(str(x) for x in range(args.beginning,args.end+1)))
            out.write("\n{}".format("\t".join(results)))
            

class AminoAcid(object):
    """
    A particular amino acid which has three main properties relevant to GV calculation:
        c: composition
        p: polarity
        v: molecular volume
    """

    def __init__(self, c, p, v):
        self.c = c
        self.p = p
        self.v = v

# Takes in the GV data information in addition to the two seqs to compare, returns
# the overall GV score across the sequence
def pairwise_gv(gv_vals,seq1,seq2):

    pass 

# Takes in (c|p|v)-(max|min) and calculates the Grantham difference
def single_gv(cmax,cmin,pmax,pmin,vmax,vmin):

    c = (1.833*(cmax-cmin))**2
    p = (0.1018*(pmax-pmin))**2
    v = (0.000399*(vmax-vmin))**2

    return "{0:.3f}".format(50.723*((c+p+v)**(.5)))


if __name__ == '__main__':
    main()