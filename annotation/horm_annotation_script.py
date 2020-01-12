#!/usr/bin/env python
import pandas as pd
import gzip
import sys
import transcript_classes as tc

def main():
    annotation_spreadsheet = "raw_files/Hcal_annotation_v5.csv"
    pinfish_files = ["raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters.gff.gz",
                     "raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c7p10.gff.gz",
                     "raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c2p20.gff.gz"]
    stringtie = ["raw_files/UCSC_Hcal_v1_B1_LR.stringtie_f01.gff.gz"]
    isoseq_hq = ["raw_files/GLO64_isoseq.collapsed.filtered.gff.gz"]
    isoseq_singletons = ["raw_files/GLO64_singletons.collapsed.gff.gz"]
    # Now we make sure that all of the pinfish files have unique IDs
    #  if they are a hash, all the IDs will be unique

    # read in the spreadsheet
    df = pd.read_csv(annotation_spreadsheet, header=0, sep=',')
    #First make sure that someone has checked the gene
    df = tc.sumone_has_checked(df)
    # now make sure that each row has something (a gene/transcript)
    df = tc.each_row_has_something(df)
    print(df.columns)
    # now make sure that there are no more genes that still need a transcript,
    #  but that have a minimap ID
    tc.still_needs_transcript(df)
    # make a dict to store the GFF files
    GFFs = {}
    # now get all of the transcripts
    GFFs["pinfish"] = tc.gffFile(pinfish_files, "pinfish")
    GFFs["stringtie"] = tc.gffFile(stringtie, "stringtie")
    GFFs["isoseq_hq"] = tc.gffFile(isoseq_hq, "isoseq_hq")
    GFFs["isoseq_singletons"] = tc.gffFile(isoseq_singletons, "isoseq_singletons")

    # now parse the spreadsheet and print out new transcripts using the GFFs dict
#ndex(['chromosome', 'stringtie_id', 'DTS_checked', 'spliced_in_intron',
#             'WRF_checked', 'isoseq_hq_id', 'pinfish_id', 'isoseq_singleton_id',
#             'remove_st', 'interesting', 'comment', 'time', 'Unnamed: 12', 'checked',
#             'one_row_one_gene'],
    column_name_to_GFF_map = {"stringtie_id": "stringtie",
                              "isoseq_hq_id": "isoseq_hq",
                              "isoseq_singleton_id": "isoseq_singletons",
                              "pinfish_id": "pinfish"
                              }
    parse_spreadsheet(df, GFFs, colum_name_to_GFF_map)

if __name__== "__main__":
    main()
