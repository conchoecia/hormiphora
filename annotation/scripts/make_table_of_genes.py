#!/usr/bin/env python

seqs_orig="h1flnc_to_Jan24_annoth1.pilon/h1flnc_to_Jan24_annoth1.pilon.fasta"
h1_fasta="prottrans_predict/h1.prottrans.fasta"
h2_fasta="prottrans_predict/h2.prottrans.fasta"

seqlist = set()
from Bio import SeqIO
import pandas as pd
import sys

gene_dict = {}
for record in SeqIO.parse(seqs_orig, "fasta"):
    seq_record = ".".join(record.id.split('.')[0:5])
    seqlist.add(seq_record)
    this_chr = seq_record.split('.')[2]
    this_gene= int(seq_record.split('.')[3].replace('g',''))
    if this_chr not in gene_dict:
        gene_dict[this_chr] = set()
    gene_dict[this_chr].add(this_gene)

# get the number of genes on each chr
for key in gene_dict:
    gene_dict[key] = max(gene_dict[key])


# get the genes in the h1 file
h1_dict = {x: -1 for x in seqlist}
for record in SeqIO.parse(h1_fasta, "fasta"):
    seqid = ".".join(record.id.split(".")[0:5])
    h1_dict[seqid] = len(record.seq)

# get the genes in the h2 file
h2_dict = {x: -1 for x in seqlist}
for record in SeqIO.parse(h2_fasta, "fasta"):
    seqid = ".".join(record.id.split(".")[0:5])
    h2_dict[seqid] = len(record.seq)

# make the dataset
df = pd.DataFrame([h1_dict, h2_dict])
df = df.T
df.columns = ["h1", "h2"]
convert_dict = {'h1': int, 'h2': int}
df = df.fillna(-2)
df = df.astype(convert_dict)


df["keep"] = "NO"
for i, row in df.iterrows():
    # both genes are the same length
    if row['h1'] == row['h2']:
        # if they're both 1, pass
        if row['h1'] == -1 and row["h2"] == -1:
            pass
        # something else is larger than -1
        else:
            df.at[i, "keep"] = "h1"
    # They're not the same length
    else:
        if row['h1'] > row['h2']:
            df.at[i, "keep"] = "h1"
        else:
            df.at[i, "keep"] = "h2"

print(df)
df.to_csv("prottrans_predict/diff_lengths.txt", sep='\t')

# we need this later to sort the output
df["isoform_index"] = -1
for i, row in df.iterrows():
    df.at[i, "isoform_index"] = int(i.split('.')[4].replace('i', ''))

df["fasta"] = ""

# add the sequence to the pandas df
col_to_file = {"h1": h1_fasta, "h2": h2_fasta}
for key in col_to_file:
    for record in SeqIO.parse(col_to_file[key], "fasta"):
        seqid = ".".join(record.id.split(".")[0:5])
        #this_chr = seqid.split('.')[2]
        #this_gene= int(seqid.split('.')[3].replace('g',''))
        one_row = df.loc[seqid]
        assert len(one_row) == 5
        if one_row["keep"] == key:
            fasta_string = ">{}\n{}".format(seqid, record.seq)
            #print("fasta string ", fasta_string)
            df.at[seqid, "fasta"] = fasta_string

# print the fasta file out in order
#first get the num Cs
num_c = max([int(key.replace('c','')) for key in gene_dict if key[0] == 'c'])
num_s = max([int(key.replace('sca','')) for key in gene_dict if 'sca' in key])
sort_list = []
for i in range(1, num_c+1):
    sort_list.append("c{}".format(i))
for i in range(1, num_s+1):
    sort_list.append("sca{}".format(i))

new_fasta = open("prottrans_predict/Hcv1.1.pep", "w")
for this_c in sort_list:
    if this_c in gene_dict:
        for this_gene in range(1, gene_dict[this_c]+1):
            this_id = "Hcv1.1.{}.g{}.".format(this_c, this_gene)
            print_these_df = df.filter(like=this_id, axis=0).sort_values(by ="isoform_index")
            if len(print_these_df) > 0:
                for i, row in print_these_df.iterrows():
                    #print the isoforms to the fasta
                    if row["keep"] != "NO":
                        print(row["fasta"], file=new_fasta)
new_fasta.close()
