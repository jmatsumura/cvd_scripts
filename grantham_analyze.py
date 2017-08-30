#! /usr/bin/env python3

# This script performs different Grantham-based analyses on a set of aligned
# peptides of the same length. These sequences can be subset into particular
# start/stop positions (entered as a 1-based index). For instance, 1 and 10
# set as the -b and -e arguments respectively will analyze AAs 1-10 in each
# sequence. 
#
# Grantham Difference (GD) is calculated when doing pairwise comparisons. One 
# can either isolate the most egregious difference (-sp option) or calculate a
# a total difference that is "normalized" by the length of the sequence (-ap option). 
# The division by length helps differentiate high concentrations of close proximity
# polymorphisms from a more dispersed set of polymorphisms. 
# Further reading: https://www.ncbi.nlm.nih.gov/pubmed/4843792 
#
# For formatting GD/pairwise outputs, one can use forceSymmetric in R if 
# all values need to be present
# Further reading: https://cran.r-project.org/web/packages/Matrix/Matrix.pdf
#
# Grantham Variation (GV) is calculated when looking at each position across 
# all sequences and isolating the mins/maxs for c/p/v values.
# Futher reading: https://www.ncbi.nlm.nih.gov/pubmed/16014699
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
    parser.add_argument('-s','--subset', type=str, required=False, help='Path to a file that includes IDs, one per line, for those peptides to include from the peptide file.')
    parser.add_argument('-o','--outfile', type=str, required=True, help='Name/location of the file to write out to.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-ap','--all_pairwise', action='store_true', help='Compute GD score for the entire sequence of all possible pairs, all GD values are added to a total sum then divided by total length.')
    group.add_argument('-sp','--single_pairwise', action='store_true', help='Obtain the most extreme GD score found at a single position for the entire sequence between all pair combinations.')
    group.add_argument('-baa','--by_amino_acid', action='store_true', help='Obtain the most extreme GV score found at a single position across ALL sequences, this is not pairwise.')
    args = parser.parse_args()

    # Establish the values from Table 1 of Grantham's work: https://www.ncbi.nlm.nih.gov/pubmed/4843792
    gd_vals = {
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

    subset_us = set()

    if args.subset:
        with open(args.subset,'r') as subset_list:
            for line in subset_list:
                subset_us.add(line.strip())

    # grab all the peptides, use the ID in the input as the ID in the hash 
    peptides = OrderedDict()
    with open(args.peptides,'r') as infile:
        for line in infile:
            if args.subset:
                if line.split(',')[0] not in subset_us:
                    continue

            peptides[line.split(',')[0]] = line.split(',')[1].strip()[(args.beginning-1):(args.end)]
    
    key_lookup = list(peptides.keys()) # based on index, get the key for this dict

    if args.by_amino_acid: # calculate GV for each position (worst case across all sequences)

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

                # Find the most extreme value for every Grantham metric
                for x in unique_sets[j]:

                    if gd_vals[x].c < cmin:
                        cmin = gd_vals[x].c
                    if gd_vals[x].c > cmax:
                        cmax = gd_vals[x].c

                    if gd_vals[x].p < pmin:
                        pmin = gd_vals[x].p
                    if gd_vals[x].p > pmax:
                        pmax = gd_vals[x].p

                    if gd_vals[x].v < vmin:
                        vmin = gd_vals[x].v
                    if gd_vals[x].v > vmax:
                        vmax = gd_vals[x].v

                results.append(single_gd(cmax,cmin,pmax,pmin,vmax,vmin))

        with open(args.outfile,'w') as out:
            out.write("\t".join(str(x) for x in range(args.beginning,args.beginning+len(results))))
            out.write("\n{}".format("\t".join(results)))

    else: # pairwise

        # store the results of all the comparisons in a giant list, order is 
        # guaranteed to be 1v2,1v3,1v4....2v3,2v4,2v5....(n-1)v(n)
        results = [] 

        cur_base_sample_pos = 0 # start at 0 and work up

        for pep1 in range(0,len(peptides)):

            subresult = []
            subresult.append("0.000") # sample compared to itself will always be 0

            pep1_seq = peptides[key_lookup[pep1]] 

            for pep2 in range(cur_base_sample_pos+1,len(peptides)): # never repeat a comparison
                
                pep2_seq = peptides[key_lookup[pep2]]

                if args.single_pairwise: # only 
                    subresult.append(pairwise_gd(gd_vals,pep1_seq,pep2_seq,True))
                else:
                    subresult.append(pairwise_gd(gd_vals,pep1_seq,pep2_seq,False))

            results.append(subresult)

            cur_base_sample_pos += 1

        with open(args.outfile,'w') as out:
            out.write("\t{}".format("\t".join(key_lookup)))

            for i,res in enumerate(results):
                filler = "\t"
                if i != 0:
                    filler = "\t{}\t".format("\t".join(["NA"]*(len(results[0])-len(results[i]))))
                out.write("\n{}{}{}".format(key_lookup[i],filler,"\t".join(res)))
            

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

# Takes in the Grantham data in addition to the two seqs to compare, returns
# either the largest GD found throughout if largest_only == True or all GD
# scores summed and divided by the length of the sequence.
def pairwise_gd(g_vals,seq1,seq2,largest_only):

    overall_gd = 0.0 # store all GDs found
    max_gd = 0.0 # since we have to iterate through already, just keep track manually

    for i,aa in enumerate(seq1):
        if aa == seq2[i]:
            pass
        else:
            cur_gd = float(single_gd(gd_vals[aa].c,g_vals[seq2[i]].c,gd_vals[aa].p,g_vals[seq2[i]].p,gd_vals[aa].v,gd_vals[seq2[i]].v))
            if cur_gd > max_gd:
                max_gd = cur_gd
            
            overall_gd += cur_gd

    if largest_only:
        return "{0:.3f}".format(max_gd)
    else:
        return "{0:.3f}".format(overall_gd/len(seq1))

# Takes in two sets of (c|p|v) values and calculates the Grantham difference
def single_gd(c1,c2,p1,p2,v1,v2):

    c = (1.833*(c1-c2))**2
    p = (0.1018*(p1-p2))**2
    v = (0.000399*(v1-v2))**2

    return "{0:.3f}".format(50.723*((c+p+v)**(.5)))


if __name__ == '__main__':
    main()