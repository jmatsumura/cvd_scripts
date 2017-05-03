

# 
#
# Run the script using a command like this:
# 
#
# Author: James Matsumura

import pandas,argparse

def main():

    parser = argparse.ArgumentParser(description='Script to convert the AMA1 data into a "long" format.')
    parser.add_argument('-peptide_map', type=str, required=True, help='Path to peptide file.')
    parser.add_argument('-sample_map', type=str, required=True, help='Path to sample file.')
    parser.add_argument('-raw_file', type=str, required=True, help='Path to microarray results.')
    parser.add_argument('-abs_matrices', type=str, required=True, help='Name of absolute value matrices.')
    parser.add_argument('-rel_matrices', type=str, required=True, help='Name of relative ratio value matrices.')
    parser.add_argument('-outfile', type=str, required=True, help='Name of the output file.')
    args = parser.parse_args()

    peptide_order,sample_order = ([] for i in range(2))
    peptide_dict,sample_dict = ({} for i in range(2))
    pep_map_back = {}

    peptide_df = pandas.read_csv(args.peptide_map,sep="\t") 
    sample_df = pandas.read_csv(args.sample_map,sep="\t") 
    raw_df = pandas.read_csv(args.raw_file,sep="\t") 

    abs = args.abs_matrices.split(',')
    rel = args.rel_matrices.split(',')
    sample_names = set()

    # Only need to iterate over one of either abs/rel as the names are the same for both
    for r in rel:
        r_df = pandas.read_csv(r,sep="\t")
        for x in list(r_df):
            sample_names.add(x)

    rel_values = {}
    # now grab the relative/abs counts for each peptide+sample ID pair
    for r in rel:
        r_df = pandas.read_csv(r,sep="\t")
        cols = list(r_df)
        for index,row in r_df.iterrows():
            for j in range(1,len(row)):
                key = "{0}:{1}".format(row[0],cols[j])
                rel_values[key] = row[j]

    abs_values = {}
    for a in abs:
        a_df = pandas.read_csv(a,sep="\t")
        cols = list(a_df)
        for index,row in a_df.iterrows():
            for j in range(1,len(row)):
                key = "{0}:{1}".format(row[0],cols[j])
                abs_values[key] = row[j]

    for index,row in sample_df.iterrows():
        if row["sample_id"] in sample_names:
            sample_order.append(row["sample_id"])
            sample_dict[row["sample_id"]] = {'sample_idx':row["sample_index"],'sample_pad':row["slide_pad"],'study_idx':row["study_index"]}

    peptides = set()
    for index,row in peptide_df.iterrows():
        peptide_order.append(row["peptide_id"])
        peptide_dict[row["peptide_id"]] = {'peptide_idx':row["peptide_index"],'row':[]}

        if 'peptide' in row["peptide_id"]:
            peptides.add(row["peptide_id"])
            pep_map_back[row["peptide_id"]] = row["peptide_id"]
        else:
            oh_one = "{0}_0.01".format(row["peptide_id"])
            oh_oh_three = "{0}_0.003".format(row["peptide_id"])
            peptides.add(oh_one)
            pep_map_back[oh_one] = row["peptide_id"]
            peptides.add(oh_oh_three)
            pep_map_back[oh_oh_three] = row["peptide_id"]

    slides = set() # differentiate between those with G25 and those without
    for col in list(raw_df):
        if col.startswith("Slide") and "G25" not in col:
            slides.add(col)

    for index,row in raw_df.iterrows():
        if "peptide_0" in row["ID"]:
            row["ID"] = row["ID"].replace('0','')
        if row["ID"] in peptides: # a relevant peptide row
            original_peptide_id = pep_map_back[row["ID"]]
            peptide_dict[original_peptide_id]['row'].append(row)

    with open (args.outfile,'w') as out:
        out.write("peptide_idx\tsample_idx\tstudy_idx\tdilution\tseroreactivity\trel_match\tpoly_aa_count\n")

        for pep in peptide_order:
            for sam in sample_order:

                pep_idx = peptide_dict[pep]['peptide_idx']
                sam_idx = sample_dict[sam]['sample_idx']
                stu_idx = sample_dict[sam]['study_idx']
                abs = abs_values["{0}:{1}".format(pep,sam)]
                rel = rel_values["{0}:{1}".format(pep,sam)]

                for x in peptide_dict[pep]['row']:
                    dilution = x['Dilution']
                    slide = sample_dict[sam]['sample_pad']
                    seroreactivity = ""
                    if slide in slides:
                        seroreactivity = x[slide]
                    
                    out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(pep_idx,sam_idx,stu_idx,dilution,seroreactivity,rel,abs))
    


if __name__ == '__main__':
    main()