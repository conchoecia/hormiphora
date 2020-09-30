#!/usr/bin/env python

"""
This makes plots to try to find the centromeres in all of the sequences
"""

#!/usr/bin/env python3
import ast
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns; sns.set()
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import struct

# set seaborn stuff
#sns.set(rc={'text.usetex' : True})
sns.set_style("ticks", {'font.family': ['sans-serif'],
                            'font.sans-serif': ['Helvetica'],
                            'grid.color': '.95'})

# Preserve the vertical order of embedded images:
matplotlib.rcParams['image.composite_image'] = False
# text as font in pdf
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import sys
from Bio import SeqIO

fasta_file = sys.argv[1]

seq_to_length = {}

for record in SeqIO.parse(fasta_file, "fasta"):
    seq_to_length[str(record.id)] = len(record.seq)

for entry in seq_to_length:
    print("{}\t{}".format(entry, seq_to_length[entry]), file=sys.stderr)

seq_to_size_class = {}
seq_to_pos_width = {}
current_sequence = ""
for line in sys.stdin:
    line = line.strip()
    if line:
        splitd = line.split(" ")
        isint = False
        try:
            testint = int(splitd[0])
            isint = True
        except:
            pass
        if splitd[0] == "Sequence:":
            # we have found a new sequence
            if splitd[1] != current_sequence:
                current_sequence = splitd[1]
                seq_to_pos_width[current_sequence] = [[],[],[]] #pos, percent_of_genome, width
        elif isint:
            seq_to_pos_width[current_sequence][0].append( int(splitd[0]) )
            seq_to_pos_width[current_sequence][1].append( 100 * ((int(splitd[1]) - int(splitd[0]))/seq_to_length[current_sequence]) )
            seq_to_pos_width[current_sequence][2].append( int(splitd[1]) - int(splitd[0]) )
            print("{}\t{}\t{}".format(current_sequence, int(splitd[1])-int(splitd[0]), "\t".join(splitd)))

for seqid in seq_to_pos_width:
    print(seqid, file=sys.stderr)
    x     = seq_to_pos_width[seqid][0]
    y     = seq_to_pos_width[seqid][1]
    width = seq_to_pos_width[seqid][2]
    plt.figure()
    plt.bar(x, y, width=width, bottom=None, align='center', color = "black", lw=0)
    plt.savefig("{}_repeats.pdf".format(seqid))
    plt.cla()
