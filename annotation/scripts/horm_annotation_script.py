#!/usr/bin/env python
import pandas as pd
import gzip
import sys
import transcript_classes as tc

def parse_file_keyword(thisfile, keyword):
    # now get all of the transcripts
    print("parsing {}".format(keyword), file=sys.stderr)
    return_me = tc.gffFile(thisfile, keyword)
    print("  - found {} genes encompassing {} isoforms".format(
        len(return_me.GTT), sum([len(return_me.GTT[key]) for key in return_me.GTT])),
          file=sys.stderr)
    if keyword == "manual":
        for gene in return_me.GTT:
            for tx in return_me.GTT[gene]:
                print("    - {}".format(tx), file=sys.stderr)
    return(return_me)

def main():
    annotation_spreadsheet = sys.argv[1]
    pinfish_files = ["raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters.gff.gz",
                     "raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c7p10.gff.gz",
                     "raw_files/UCSC_Hcal_v1_B1_LR.pinfish_clusters_c2p20.gff.gz"]
    stringtie = ["raw_files/UCSC_Hcal_v1_B1_LR.stringtie_f01.gff.gz"]
    augustus  = ["raw_files/horm_augustus.gff"]
    isoseq_hq = ["raw_files/GLO64_isoseq.collapsed.filtered.gff.gz"]
    isoseq_singletons = ["raw_files/GLO64_singletons.collapsed.gff.gz"]
    manual  = ["raw_files/horm_manual.gff"]
    # Now we make sure that all of the pinfish files have unique IDs
    #  if they are a hash, all the IDs will be unique
    # now parse the spreadsheet and print out new transcripts using the GFFs dict

    column_name_to_GFF_map = {"stringtie_id": "stringtie",
                              "isoseq_hq_id": "isoseq_hq",
                              "isoseq_singleton_id": "isoseq_singletons",
                              "pinfish_id": "pinfish",
                              "augustus": "augustus",
                              "manual": "manual"
                              }

    # read in the spreadsheet
    df = pd.read_csv(annotation_spreadsheet, header=0, sep=',', comment = '#')
    #First make sure that someone has checked the gene
    df = tc.sumone_has_checked(df)
    # now make sure that each transcript has a sensible chromosome
    chr_list = ["c{}".format(i) for i in range(1,14)] + ["sca{}".format(i) for i in range(1,33)] + ["M"]
    indices = tc.sensible_chromosomes(df, chr_list)
    if len(indices) != 0:
        print("some of the rows didn't have a chromosomes", file=sys.stderr)
        print(indices, file=sys.stderr)
        raise Exception("missing chromosomes in rows")
    # now make sure that each row has something (a gene/transcript)
    df = tc.each_row_has_something(df, column_name_to_GFF_map)
    #print(df.columns, file=sys.stderr)
    # now make sure that there are no rows that have no stringtie but have a delete
    indices = tc.delete_but_no_stringtie(df)
    if len(indices) != 0:
        print("some of the rows had no stringtie, but are marked for deletion", file=sys.stderr)
        print(indices, file=sys.stderr)
        raise Exception("marked for deletion but missing stringtie")


    # PASSED CHECKS. NOW ANNOTATE.
    # make a dict to store the GFF files
    GFFs = {}
    # now get all of the transcripts
    parse_these = {"pinfish": pinfish_files,
                   "stringtie": stringtie,
                   "isoseq_hq": isoseq_hq,
                   "isoseq_singletons": isoseq_singletons,
                   "augustus": augustus,
                   "manual": manual
                  }

    for key in parse_these:
        keyword = key
        filelist = parse_these[key]
        GFFs[keyword] = parse_file_keyword(filelist, keyword)
        #if key == "augustus":
        #    print(GFFs["augustus"].IDTS, file=sys.stderr)
        #    print(GFFs["augustus"].GTT, file=sys.stderr)

    # Now make sure that there are no gene IDs shared between any of the GFFs
    all_ids = dict()
    #print(GFFs, file = sys.stderr)
    for key in GFFs:
        for thisgene in GFFs[key].GTT:
            #print(thisgene, file = sys.stderr)
            if thisgene not in all_ids:
                all_ids[thisgene] = 1
            else:
                all_ids[thisgene] += 1
    #check stringtime_manual
    print_this = False
    print_message = "The following genes were duplicates across multiple GFF files.\n"
    for thisgene in all_ids:
        if all_ids[thisgene] > 1:
            print_this = True
            print_message = print_message + "  - {} - {}\n".format(thisgene, all_ids[thisgene])
    if print_this:
        print(print_message, file = sys.stderr)
        raise IOError("See above message.")


    tc.parse_spreadsheet(df, GFFs, column_name_to_GFF_map)

if __name__== "__main__":
    main()
