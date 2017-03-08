
# Takes the output from EMBOSS' primer search program and produces a TSV file
# with the following columns:
#
# (1) Primer name
# (2) Sequence ID
# (3) Amplimer length
# (4) Forward # of mismatches
# (5) Reverse # of mismatches
# (6) Forward position
# (7) Reverse position
#
# This script takes two arguments: 
# -input_file: the output generated from EMBOSS' primer search
# -output_file: name of the TSV file here to write output
#
# Run the script using a command like this:
# python3 primersearch_to_tsv.py -input_file primersearch.out -output_file primersearch_out.tsv
#
# Author: James Matsumura
# Contact: jmatsumura@som.umaryland.edu

import argparse,os,re

def main():

    parser = argparse.ArgumentParser(description="Script to convert the output of EMBOSS' primer search program into a TSV file.")
    parser.add_argument('-input_file', type=str, required=True, help='Location of the output from primer search.')
    parser.add_argument('-output_file', type=str, required=True, help='Location/name of the TSV converted file.')
    args = parser.parse_args()

    regex_for_fr_s = r'strand\sat\s(\[?\d+\]?)\swith\s(\d+)' # capture position and mismatches

    # Assign variables for the columns noted in the description above
    pn,sid,al = ("" for i in range(3))
    fm,rm,fp,rp = (0 for i in range(4))
    # booleans for tracking the forward/reverse regions
    forward_section,reverse_section = (False for i in range(2))

    with open(args.output_file,'w') as outfile:
        with open(args.input_file,'r') as infile:
            for line in infile:
                line = line.strip()

                if line.startswith('Primer name'):
                    pn = line.split(' ',2)[2]
                
                elif line.startswith('Sequence:'):
                    sid = line.split(' ')[1]
                    forward_section = True

                elif forward_section == True and line.startswith('G'):
                    fm = re.search(regex_for_fr_s,line).group(2)
                    fp = re.search(regex_for_fr_s,line).group(1)
                    forward_section = False
                    reverse_section = True

                elif reverse_section == True:
                    rm = re.search(regex_for_fr_s,line).group(2)
                    rp = re.search(regex_for_fr_s,line).group(1)
                    reverse_section = False

                elif line.startswith('Amplimer length:'):
                    al = line.split(' ',2)[2]
                    
                    # Once we've made it here we need to output the current entry
                    outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(pn,sid,al,fm,rm,fp,rp))
            

if __name__ == '__main__':
    main()