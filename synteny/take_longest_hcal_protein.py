#!/usr/bin/env python

# take_longest_hcal_protein.py

'''
./take_longest_hcal_protein.py Hcv1a1d20200309_model_proteins.pep > Hcv1a1d20200309_model_proteins.pep.longest_only.fa
'''

import sys
from Bio import SeqIO

if len(sys.argv) < 2:
	sys.exit( __doc__ )
else:
	lastgene = ""
	secrec_dict = {} # key is length, value is seqrecord object
	for seqrec in SeqIO.parse(sys.argv[1],"fasta"):
		seqid = seqrec.id
		geneid = seqid.rsplit(".",1)[0]
	#	print >> sys.stderr, seqid, geneid, secrec_dict
		if geneid==lastgene:
			seqlength = len(seqrec.seq)
			secrec_dict[seqlength] = seqrec
		else: # moved on to a different gene
			if secrec_dict: # if previous gene has entries, print the longest
				longest_isoform_len = max(secrec_dict.keys())
				sys.stdout.write( secrec_dict[longest_isoform_len].format("fasta") )
				secrec_dict = {} # print old entry, reset dict, add new entry
				seqlength = len(seqrec.seq)
				secrec_dict[seqlength] = seqrec
			else: # first gene
				seqlength = len(seqrec.seq)
				secrec_dict[seqlength] = seqrec
		lastgene = geneid
	else: # for last gene
		longest_isoform_len = max(secrec_dict.keys())
		sys.stdout.write( secrec_dict[longest_isoform_len].format("fasta") )
		secrec_dict = {}
