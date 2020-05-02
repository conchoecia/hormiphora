#!/usr/bin/env python
"""
This program makes a table of protein sizes between two haplotypes,
  then makes a single file with all the proteins
"""

import sys

# the length must be 6
# 0 is the program name
# 1 is one of the original nucleotide fasta files
# 2 is the h1 pep file
# 3 is the h2 pep file      (output)
# 4 is the diff lengths csv (output)
if len(sys.argv) != 6:
    print("""# the length must be 6
# 0 is the program name
# 1 is one of the original nucleotide fasta files
# 2 is the h1 pep file
# 3 is the h2 pep file
# 4 is the diff lengths csv   (output)
# 5 is the model protein file (output)""")
    sys.exit()

seqs_orig        = sys.argv[1]
h1_fasta         = sys.argv[2]
h2_fasta         = sys.argv[3]
diff_lengths_csv = sys.argv[4]
model_pep_file   = sys.argv[5]


seqlist = set()
from Bio import SeqIO
import pandas as pd
import sys

print(" - Reading in the sequences from the transcripts from the gff", file=sys.stderr)
gene_dict = {}
for record in SeqIO.parse(seqs_orig, "fasta"):
    seq_record = ".".join(record.id.split('.')[0:5])
    seqlist.add(seq_record)
    this_chr = seq_record.split('.')[2]
    this_gene= int(seq_record.split('.')[3].replace('g',''))
    if this_chr not in gene_dict:
        gene_dict[this_chr] = set()
    gene_dict[this_chr].add(this_gene)

print(" - Getting the number of genes on each chromosome.", file=sys.stderr)
# get the number of genes on each chr
for key in gene_dict:
    gene_dict[key] = max(gene_dict[key])

print(" - Getting the h1 genes.", file=sys.stderr)
# get the genes in the h1 file
h1_dict = {x: -1 for x in seqlist}
for record in SeqIO.parse(h1_fasta, "fasta"):
    seqid = ".".join(record.id.split(".")[0:5])
    h1_dict[seqid] = len(record.seq)

print(" - Getting the h2 genes.", file=sys.stderr)
# get the genes in the h2 file
h2_dict = {x: -1 for x in seqlist}
for record in SeqIO.parse(h2_fasta, "fasta"):
    seqid = ".".join(record.id.split(".")[0:5])
    h2_dict[seqid] = len(record.seq)

# make the dataset
print(" - Making a table of which isoform to keep.", file=sys.stderr)
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
df.to_csv(diff_lengths_csv, sep='\t', index_label="gene")

print(" - Getting the isoform indices.", file=sys.stderr)
# we need this later to sort the output
df["isoform_index"] = -1
for i, row in df.iterrows():
    df.at[i, "isoform_index"] = int(i.split('.')[4].replace('i', ''))

df["fasta"] = ""

print(" - Adding the protein sequences to the pandas dataframe.", file=sys.stderr)
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

print(" - Getting number of chromosomes.", file=sys.stderr)
# print the fasta file out in order
#first get the num Cs
num_c = max([int(key.replace('c','')) for key in gene_dict if key[0] == 'c'])
num_s = max([int(key.replace('sca','')) for key in gene_dict if 'sca' in key])
sort_list = []
for i in range(1, num_c+1):
    sort_list.append("c{}".format(i))
for i in range(1, num_s+1):
    sort_list.append("sca{}".format(i))

print(" - Printing out the fasta file.", file=sys.stderr)
new_fasta = open(model_pep_file, "w")
for this_c in sort_list:
    if this_c in gene_dict:
        for this_gene in range(1, gene_dict[this_c]+1):
            this_id = ".{}.g{}.".format(this_c, this_gene)
            print_these_df = df.filter(like=this_id, axis=0).sort_values(by ="isoform_index")
            if len(print_these_df) > 0:
                for i, row in print_these_df.iterrows():
                    #print the isoforms to the fasta
                    if row["keep"] != "NO":
                        print(row["fasta"], file=new_fasta)
new_fasta.close()
