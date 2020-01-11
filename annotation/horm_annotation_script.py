#!/usr/bin/env python3
import ast
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter
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
assert len(t1) == 0

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
                ID = line.split('\t')[8].split(';')[0].replace("ID=","")
                if ID in ids:
                    ids[ID] += 1
                else:
                    ids[ID] = 1
assert len(np.unique(list(ids.values()))) == 1
#yes, all the hashes are the same

class Gene:
    """
    just a container for holding information about genes.
    mostly useful to hold information to add to the flags,
    like "interesting", "operon", et cetera
    """
    
class Transcript:
    """
    Instance attributes:
     - filename
        - The filename will be used to look up the transcript information.
     - scaffold
        - The scaffold is the name of the scaffold on which the transcript is
          located. Like 'c1', 'c2', 'M', et cetera
     - num_on_scaffold
        - This is the index of the transcript on that scaffold. For instance 
          if this gene is the 117th gene that occurs when looking at the chromosome
          from 5' to 3', then this number is 117
     - ID
        - The unique identifier of this gene, or specific isoform in the file
          provided in filename. Depending on the filetype this will change how things
          are parsed
     - ftype
        - The program that was used to generate the transcripts. Can be one of 
          the following types.
          - "stringtie"
          - "manual"
          - "isoseq"
          - "pinfish"
        - We will use grep to find these IDs. This is useful because it allows
          us to control if we want specific isoforms of a gene, or every single isoform.
          This is only possible for stringtie and IsoSeq IDs.
    """
    def __init__(self, filename, scaffold, num_on_scaffold, ID, ftype):
        self.filename = filename
        self.scaffold = scaffold
        self.num_on_scaffold = num_on_scaffold
        self.ID = ID
        self.ftype = ftype

#There is a unique way to parse the IsoSeq files or Stringtie files.
# The ID can be like this, "PB.2026", with no period-delimited additional index.
# The ID can also be like this, "PB.2026.1", with an additional index, ".1".
# If the ID provided has no additional index, this means to turn every
#  isoform of this gene into a transcript for the new Hormiphora genome.
# If there is an additional index, that means to accept only that specific isoform
#  for the transcript.
# Because of this complicating factor, we first have to parse the IsoSeq
#  and stringtie files to build an index of which genes have which transcripts.

class gffFile:
    """
    This is just a way to keep track of what transcripts are in what file.
    It will contain all of the necessary information to pull out whole genes
     or individual transcripts from a gff file.
    
     - The program that was used to generate the transcripts. Can be one of 
        the following types.
        - "stringtie"
        - "manual"
        - "isoseq"
        - "pinfish"
    """
    def __init__(self, filename, filetype):
        # the init store the filename, and also opens up the file and parses it for all of
        #  the pertinent information
        
        # if we're parsing pinfish genes, it is better to put everything into
        # a single gffFile class to make it easier to parse later on.
        self.transcript_id_to_string = {}
        self.filename = filename
        if type(filename) == list and filetype == pinfish:
            for thisfile in self.filename:
                if thisfile.endswith(".gff.gz"):
                    f = gzip.open(thisfile, "rb")
                elif thisfile.endswith(".gff"):
                    f = open(thisfile, "r")
                else:
                    raise Error("dunno what kind of file this is. must be .gff or .gff.gz")

        
