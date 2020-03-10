#!/usr/bin/env python

#prints out the transcripts that are unique to this file.
import sys

#structure is key: stem     Hcv1.1.c9.g786
#           value: absolute Hcv1.1.c9.g786.i2
all_genes = {}

for line in sys.stdin:
    line = line.strip()
    stem = ".".join(line.split('.')[0:4])
    if stem not in all_genes:
        all_genes[stem] = []
    all_genes[stem].append(line)

for key in all_genes:
    if len(all_genes[key]) == 1:
        print(all_genes[key][0])
