#!/usr/bin/env python
import pandas as pd
import numpy as np
import gzip
import collections
import transcript_classes as tc

def main():
    annotation_spreadsheet = "raw_files/Hcal_annotation_v5.csv"
    pinfish_files = ["raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters.gff.gz",
                     "raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c7p10.gff.gz",
                     "raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c2p20.gff.gz"]
    # Now we make sure that all of the pinfish files have unique IDs
    #  if they are a hash, all the IDs will be unique

    # read in the spreadsheet
    df = pd.read_csv(annotation_spreadsheet, header=0, sep=',')
    #First make sure that someone has checked the gene
    df = tc.sumone_has_checked(df)
    # now make sure that each row has something (a gene/transcript)
    df = tc.each_row_has_something(df)
    # now make sure that there are no more genes that still need a transcript,
    #  but that have a minimap ID
    tc.still_needs_transcript(df)
    # make a dict to store the GFF files
    GFFs = {}
    GFFs["pinfish"] = tc.gffFile(pinfish_files, "pinfish")
    print(GFFs["pinfish"].IDTS)

if __name__== "__main__":
    main()
