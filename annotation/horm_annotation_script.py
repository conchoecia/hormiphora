#!/usr/bin/env python
import pandas as pd
import numpy as np
import gzip
import collections

workingdir="/Users/darrin/git/hormiphora/annotation/"
filepath = "{}/raw_files/Hcal_annotation_v5.csv".format(workingdir)

df = pd.read_csv(filepath, header=0, sep=',')
df["checked"] = "none"

#First make sure that someone has checked the gene
for i, row in df.iterrows():
    C1 = False
    C2 = False
    if type(row["DTS_checked"]) == str:
        C1 = row["DTS_checked"].strip().lower() in ['y', 'yes']
    if type(row["WRF_checked"]) == str:
        C2 = row["WRF_checked"].strip().lower() in ['y', 'yes']
    df.at[i,'checked'] = C1 or C2

# First make sure that each row has one of DTS_checked or WRF_checked
t1 = df.loc[df['checked'] == False, ]
assert len(t1) == 0

# Now make sure that each row has something in stringtie_id, isoseq_hq_id, pinfish_id, or isoseq_singleton_id
#  don't count rows that also have isoseq reads containing ID m64069
df["one_row_one_gene"] = "none"
for i, row in df.iterrows():
    C1_ST = False
    C2_IS = False
    C3_PF = False
    C4_SI = False
    C5_CO = False
    if type(row["stringtie_id"]) == str:
        C1_ST = row["stringtie_id"].strip().lower() != ""
    if type(row["isoseq_hq_id"]) == str:
        C2_IS = row["isoseq_hq_id"].strip().lower() != ""
    if type(row["pinfish_id"]) == str:
        C3_PF = row["pinfish_id"].strip().lower() != ""
    if type(row["isoseq_singleton_id"]) == str:
        C4_SI = row["isoseq_singleton_id"].strip().lower() != ""
    if type(row["comment"]) == str:
        if 'm64069' in row["comment"].strip().lower():
            C5_CO = True
    df.at[i,'one_row_one_gene'] = C1_ST or C2_IS or C3_PF or C4_SI or C5_CO

t1 = df.loc[df['one_row_one_gene'] == False, ]
assert len(t1) == 0

# Now get a list of genes left over that still need annotation from reads
#  m64069
df["only_minimap"] = "none"
for i, row in df.iterrows():
    C1_ST = False
    C2_IS = False
    C3_PF = False
    C4_SI = False
    C5_CO = False
    if type(row["stringtie_id"]) == str:
        C1_ST = row["stringtie_id"].strip().lower() != ""
    if type(row["isoseq_hq_id"]) == str:
        C2_IS = row["isoseq_hq_id"].strip().lower() != ""
    if type(row["pinfish_id"]) == str:
        C3_PF = row["pinfish_id"].strip().lower() != ""
    if type(row["isoseq_singleton_id"]) == str:
        C4_SI = row["isoseq_singleton_id"].strip().lower() != ""
    df.at[i,'only_minimap'] = C1_ST or C2_IS or C3_PF or C4_SI

t1 = df.loc[df['only_minimap'] == False, ]
print(t1)
# just skip this for now since it's in-progress
#assert len(t1) == 0

# Now we make sure that all of the pinfish files have unique IDs
#  if they are a hash, all the IDs will be unique
P1="/Users/darrin/Downloads/HC_spliced/UCSC_Hcal_v1_B1_LR.pinfish_clusters.gff.gz"
P2="/Users/darrin/Downloads/HC_spliced/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c7p10.gff.gz"
P3="/Users/darrin/Downloads/HC_spliced/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c2p20.gff.gz"

ids = {}
Files = [P1, P2, P3]
for thisfile in Files:
    with gzip.open(thisfile, "rb") as this_f:
        for line in this_f:
            line = line.decode("utf-8")
            if "pinfish	transcript" in line:
                print(line.split('\t'))
                sys.exit()
                ID = line.split('\t')[8].split(';')[0].replace("ID=","")
                if ID in ids:
                    ids[ID] += 1
                else:
                    ids[ID] = 1
assert len(np.unique(list(ids.values()))) == 1
#yes, all the hashes are the same

def main():
    pass

if __name__== "__main__":
    main()
