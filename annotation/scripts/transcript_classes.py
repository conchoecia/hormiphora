#!/usr/bin/env python
import gzip
import os
import pandas as pd
import numpy as np
import sys
pd.options.display.width = 0

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

        #IDTS = transcript ID To String
        self.IDTS = {}
        # GTT genes to transcripts
        self.GTT  = {}
        #not sure what this is
        self.used_isoforms = dict()
        self.filename = filename
        self.filetype = filetype
        if type(filename) == str:
            print(filename, file=sys.stderr)
            raise Exception("Files should be passed in a list.")
        elif type(filename) == list:
            for thisfile in self.filename:
                if not os.path.exists(thisfile):
                    raise IOError("file {} does not exist".format(thisfile))
                isgz = False
                if thisfile.endswith(".gff.gz"):
                    f = gzip.open(thisfile, "rb")
                    isgz = True
                elif thisfile.endswith(".gff"):
                    f = open(thisfile, "r")
                else:
                    print("bad filename -> ", filename, file=sys.stderr)
                    raise Exception("""dunno what kind of file this is.
                                must be .gff or .gff.gz""")
                for line in f:
                    nl = ""
                    if isgz:
                        nl = line.decode("utf-8")
                    else:
                        nl=line
                    if nl and (nl.strip() != ""):
                        splitd = nl.split('\t')
                        #print("printing nl: ", nl, file=sys.stderr)
                        #print("printing splitd: ", splitd, file=sys.stderr)
                        if str(splitd[2]).strip() in ["transcript", "mRNA"]:
                            # we have just found a new transcript.
                            # make sure that the input is legal
                            if splitd[1] in  ["pinfish", "StringTie", "custom", "AUGUSTUS"]:
                                if not splitd[8].split(';')[0].startswith("ID="):
                                    # should start with ID=
                                    print(line, file=sys.stderr)
                                    raise Exception("""There is some input error. We found
                                    a line that doesn't have field 9 starting with ID=.
                                    all pinfish/StringTie transcripts start with this""")
                            elif splitd[1] in ["PacBio"]:
                                if not splitd[8].split(';')[0].startswith("gene_id"):
                                    # should start with gene_id
                                    print(line, file=sys.stderr)
                                    raise Exception("""There is some input error. We found
                                    a line that doesn't have field 9 starting with gene_id.
                                    all PacBio transcripts start with this""")
                            else:
                                print(nl, file=sys.stderr)
                                raise IOError("Encountered some unknown while parsing gene type.")

                            # now that we made sure the input is legal, let's parse
                            #  the gene id and the transcript ID.
                            if splitd[1] == "pinfish":
                                tID = splitd[8].split(';')[0].replace("ID=","").strip()
                                gID = tID
                            elif splitd[1] in ["StringTie", "AUGUSTUS", "custom"]:
                                tID = [x.replace("ID=","").strip() for x in splitd[8].split(';') if "ID=" in x][0]
                                gID = ".".join(tID.split('.')[0:-1])
                            elif splitd[1] == "PacBio":
                                tID = splitd[8].split(';')[1].replace("transcript_id","").strip().replace("\"", "")
                                gID = splitd[8].split(';')[0].replace("gene_id","").strip().replace("\"", "")
                            else:
                                raise IOError("Encountered some unknown while parsing gene id.")

                            if tID in self.IDTS:
                                raise Exception("""This transcript ({}) is already in the map.
                                            we shouldn't see it here yet.""".format(
                                            tID))
                            # Now that we have the transcript ID and the geneID
                            #  store them in the instance
                            if splitd[1] == "pinfish":
                                # every pinfish transcript is also its own gene.
                                #  each transcript just has its own hash.
                                assert tID not in self.IDTS
                                assert gID not in self.GTT
                                self.GTT[gID] = [tID]
                            elif splitd[1] in ["StringTie", "PacBio", "AUGUSTUS", "custom"]:
                                assert tID not in self.IDTS
                                if gID not in self.GTT:
                                    self.GTT[gID] = [tID]
                                else:
                                    self.GTT[gID].append(tID)
                            self.IDTS[tID] = nl
                        elif splitd[2] == "exon":
                            # we have found an exon for this transcript. We should
                            #  have already found the transcript itself.

                            #make sure that the input is legal
                            if splitd[1] in  ["pinfish"]:
                                if not splitd[8].split(';')[0].startswith("Parent="):
                                    print(line, file=sys.stderr)
                                    raise Exception("""There is some input error. We found
                                  a line that doesn't have field 9 starting with Parent=.
                                    all pinfish exons start with this""")
                            elif splitd[1] in  ["StringTie", "AUGUSTUS"]:
                                if not splitd[8].split(';')[0].startswith("ID="):
                                    print(line, file=sys.stderr)
                                    raise Exception("""There is some input error. We found
                                  a line that doesn't have field 9 starting with ID=.
                                    all StringTie and AUGUSTUS exons start with this""")
                            elif splitd[1] in  ["PacBio"]:
                                if not splitd[8].split(';')[0].startswith("gene_id \""):
                                    print(line, file=sys.stderr)
                                    raise Exception("""There is some input error. We found
                                  a line that doesn't have field 9 starting with gene_id.
                                  all PacBio exons start with this""")
                            elif splitd[1] in ["custom"]:
                                # the input format is variable, but should have parent
                                if "Parent=" not in splitd[8]:
                                    print(line, file=sys.stderr)
                                    raise Exception("""There is some input error. We found
                                    a line for a custom gene that doesn't have field 9 
                                    containing Parent=
                                    all custom exons contain this""")
                            else:
                                raise IOError("Encountered some unknown while parsing exons")

                            # now that we made sure the input is legal, let's parse
                            #  the gene id and the transcript ID.
                            # this block is messy and needs to be reworked and refactored. Redundant code.
                            if splitd[1] == "pinfish":
                                tID = splitd[8].split(';')[0].replace("Parent=","").strip()
                            elif str(splitd[1]).strip() in ["StringTie", "AUGUSTUS"]:
                                tID = splitd[8].split(';')[1].replace("Parent=","").strip()
                            elif splitd[1] == "PacBio":
                                tID = splitd[8].split(';')[1].replace("transcript_id","").strip().replace("\"", "")
                            elif splitd[1] == "custom":
                                temp = splitd[8].split(';')
                                parent_index = 0
                                for i in range(len(temp)):
                                    if "Parent=" in temp[i]:
                                        parent_index=i
                                tID = temp[parent_index].replace("Parent=","").strip()
                            else:
                                raise IOError("Encountered some unknown while parsing transcript IDs")


                            # Now that we have the transcript ID
                            #  store them in the instance
                            if tID not in self.IDTS:
                                # we should have already seen the transcript
                                #  if the gff file is sorted properly
                                print("offending ID: ", tID, file=sys.stderr)
                                print("offending file:", thisfile, file = sys.stderr)
                                raise Exception("""For some reason we found an exon for a
                                transcript before we found the transcript itself.
                                The GFF file should have all of the transcripts first.""")
                            self.IDTS[tID] += nl

def DoL_empty(DoL):
    """
    checks if a dictionary of lists (DoL) is empty.
    For example
     { key: [], key2: []} is empty
     { key: [1,2], key2: []} is not

    Returns True if empty, false if there is something in the DoL
    """
    for key in DoL:
        if len(DoL[key]) > 0:
            return False
    return True

def parse_spreadsheet(df, GFFs, CTGm):
    """Go through the spreadsheet,
      one row at a time, and construct transcripts

    GFFs are the GFF file objects as a dict with a lookup key
      and a value of the GFF object

    CTGm is the column_name_to_GFF_map
    """
    the_source = "Hcv1"
    this_chromosome = ""
    gene_counter = ""
    isoform_counter = 0
    # this is used to catch the error mode in which the isoform counter is not incrementing for each transcript
    isoforms_printed = set()
    for i, row in df.iterrows():
        #first, determine which chromosome we are on, and what transcript
        row_chr = str(row["chromosome"])
        if row_chr != this_chromosome:
            #we have either just started, or transitioned to a new chromosome
            this_chromosome = row_chr
            gene_counter = 0

        #make sure that the gene_name string doesn't have any illegal characters
        illegal_chars = [",", ";", "="]
        if type(row["gene_name"]) == str:
            for thischar in illegal_chars:
                #print("going to print row",file=sys.stderr)
                #print(row, file=sys.stderr)
                if thischar in row["gene_name"]:
                    raise IOError("{}\nthis row's gene_name field has an illegal character: {}".format(row, thischar))
        #NOW WE PARSE THE ROW - EACH ROW IS A GENE
        # first we check if this row has a delete flag or not.

        # we have flagged it for deletion
        # if there is nothing else in all the columns aside from stringtie,
        #  then just skip this. Darrin probably made this row.
        #  The annotation pattern for the first four chromosomes was to
        #  mark the row for deletion if the stringtie model was bad, then
        #  insert a new row below and use the correct model. WRF put the
        # correct model on the same line, then DTS started doing that on
        # chrs c5 and c6. Hence the special parsing.

        #there is either a stringtie ID, or there is not one.
        isoforms_in_this_gene = {CTGm[key]: [] for key in CTGm}
        # look in every single possible place we may have drawn a transcript from
        for key in CTGm:
            skip = False
            if key == "stringtie_id":
                if str(row["remove_st"]).strip().lower() in ["y", "yes"]:
                    skip = True
            if not skip and not pd.isnull(row[key]):
                # we have found a gene for this GFF file
                splitd = str(row[key]).split(",")
                for tx in splitd:
                    isoforms_in_this_gene[CTGm[key]].append(tx.strip())
            if len(isoforms_in_this_gene[CTGm[key]]) ==  0:
                isoforms_in_this_gene.pop(CTGm[key])
        # we now have all the isoforms_in_this_gene. now print them out
        if not DoL_empty(isoforms_in_this_gene):
            #there are some transcripts here.
            gene_counter += 1
            isoform_counter = 1
            this_gene = "Hcv1.1.{}.g{}".format(this_chromosome, gene_counter)

            # look for "FORWARD STRAND" and "REVERSE STRAND" to override
            #  the strand for everything
            CATCH_EM_ALL = ["FORWARD STRAND", "REVERSE STRAND"]
            strand = ""
            # this block sets
            for this_one in CATCH_EM_ALL:
                if this_one in str(row["comment"]).strip():
                    matching_cases = 0
                    if "FORWARD" in str(row["comment"]).strip():
                        strand = "+"
                        matching_cases += 1
                    if "REVERSE" in str(row["comment"]).strip():
                        strand = "-"
                        matching_cases += 1
                    #check if something weird happened
                    if matching_cases == 0:
                        raise Exception("""We shouldn't have found a 0 here.
                        This means the genes comment says this should be both
                        forward strand and reverse strand.
                        Consult your local programmer to debug.""")
                    if matching_cases > 1:
                        raise Exception("""Matched to multiple cases for strand.
                        This shouldn'ta happened. You comment should either
                        contain FORWARD STRAND or contain REVERSE STRAND.""")
            #add them to the gene and print
            print_buffer = ""
            gene_coords = [-1, -1]
            this_chr = ""
            for key in isoforms_in_this_gene:
                for this_isoform_ID in isoforms_in_this_gene[key]:
                    # the transcript ID might be for a single isoform, or it might
                    #   be for a bunch of isoforms. To determine which, we will
                    #   look in the GFF file's GTT instance attribute
                    #   (Gene-to-transcript)
                    look_these_up = []
                    if this_isoform_ID in GFFs[key].GTT:
                        #This is a gene and there are other transcripts to look up
                        look_these_up = GFFs[key].GTT[this_isoform_ID]
                    elif this_isoform_ID in GFFs[key].IDTS:
                        #this is a single transcript that we're adding
                        look_these_up = [this_isoform_ID]
                    else:
                        #We should have encountered something. Raise an error if
                        #  we didn't find it.
                        print(isoforms_in_this_gene, file=sys.stderr)
                        print("We couldn't find: ", this_isoform_ID, " in ", key, file=sys.stderr)
                        raise Exception("Couldn't find the transcript in the GFF file object")
                    # make sure that this isoform isn't a duplicate
                    for lookup in look_these_up:
                        if lookup not in GFFs[key].used_isoforms:
                            GFFs[key].used_isoforms[lookup] = 1
                        else:
                            GFFs[key].used_isoforms[lookup] += 1
                    #now that we have a list of isoforms to add to this gene,
                    #  start printing things out
                    for lookup in look_these_up:
                        #HCv1.1.c1.0000
                        this_isoform = "{}.i{}".format(this_gene,isoform_counter)
                        isoform_counter += 1
                        iso_string = GFFs[key].IDTS[lookup]
                        lines_split = iso_string.split('\n')
                        # for now just print everything out
                        exon_counter=1
                        for line in lines_split:
                            if line.strip():
                                gff_split = line.split('\t')
                                if len(gff_split) != 9:
                                    print("erroneous line: ", file = sys.stderr)
                                    print(gff_split, file = sys.stderr)
                                    raise IOError("This line was too long.")
                                #get the info if the gene is spliced in an intron
                                SII=""
                                if not pd.isnull(row["spliced_in_intron"]):
                                    SII="y"
                                else:
                                    SII="n"

                                #col0 sequence
                                if this_chr == "":
                                    this_chr = gff_split[0].strip()
                                #col1 source
                                gff_split[1] = the_source
                                #col2 feature
                                #col3 start
                                if (int(gff_split[3]) < gene_coords[0]) or (gene_coords[0] == -1):
                                    gene_coords[0] = int(gff_split[3])
                                #col4 end
                                if (int(gff_split[4]) > gene_coords[1]) or (gene_coords[1] == -1):
                                    gene_coords[1] = int(gff_split[4])
                                #col5 score
                                gff_split[5] = "."
                                #col6 strand
                                if strand == "":
                                    # take the first strand value if we haven't
                                    #  seen it yet
                                    strand = gff_split[6].strip()
                                # seems redundant but important for overriding strand
                                #  when necessary
                                gff_split[6] = strand
                                #col7 phase
                                #col9 attribute

                                #is the gene interesting
                                INT=""
                                if not pd.isnull(row["interesting"]):
                                    INT="y"
                                else:
                                    INT="n"
                                #print("GFFSPLIT")
                                #print(gff_split)
                                if str(gff_split[2]).strip() in ["transcript", "mRNA"]:
                                    gff_split[2] = "transcript"
                                    if this_isoform in isoforms_printed:
                                        print(isoforms_in_this_gene,
                                              file=sys.stderr)
                                        raise Exception("We already printed the isoform {}".format(this_isoform))
                                    else:
                                        isoforms_printed.add(this_isoform)
                                    comment="ID={0};Parent={1};source_program={2};source_ID={3}".format(this_isoform,this_gene,key,this_isoform_ID)
                                elif gff_split[2] == "exon":
                                    comment="Parent={0}".format(this_isoform,
                                             exon_counter)
                                    exon_counter += 1
                                else:
                                    print(gff_split, file=sys.stderr)
                                    raise Exception("Encountered some type of GFF entry we don't konw")
                                if SII=="y":
                                    comment += ";SII=y"
                                if INT=="y":
                                    comment += ";INT=y"
                                gff_split[8] = comment
                                print_buffer += "{}\n".format("\t".join(gff_split))
            #now that we have looked at every isoform construct a gene line
            gene = [this_chr, the_source, "gene", str(gene_coords[0]),
                    str(gene_coords[1]), ".",
                    strand, ".",
                    "ID={0};Name={0}".format(this_gene)]
            if type(row["gene_name"]) == str:
                gene[-1] += ";Description={}".format(row["gene_name"].strip())
            print("\t".join(gene))
            print(print_buffer, end="")
    # now that everything has been parsed, we can check to see if any GFFs
    #  had an isoform used more than once
    print_message = """ - We found that the following isoforms were used more
       than one time. Please make sure that within a column, an isoform is only
       used once.\n"""
    print_yes = False
    for key in GFFs:
        print_message = print_message + "    - {}\n".format(key)
        for entry in GFFs[key].used_isoforms:
            if GFFs[key].used_isoforms[entry] > 1:
                print_yes = True
                print_message = print_message + "      - {} - {}\n".format(entry, GFFs[key].used_isoforms[entry])
    if print_yes:
        print(print_message, file = sys.stderr)
        raise IOError("Input error. See above message")

def sumone_has_checked(df):
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
    t1 = df.loc[df['checked'] == False, ]
    # program will crash if this fails
    if len(t1) != 0:
        print(t1)
        raise Exception("There are some genes that haven't been checked.")

    return(df)

def delete_but_no_stringtie(df):
    """
    make sure that there are no rows that are slated to delete while
     lacking stringtie IDs. That is probably an error.
    """
    error_rows = []
    for i, row in df.iterrows():
        if pd.isnull(row["stringtie_id"]) and not pd.isnull(row["remove_st"]):
            error_rows.append(i)
    return(error_rows)

def sensible_chromosomes(df, chr_list):
    """
    This function goes through the dataframe and makes sure that every
     row has a chromosome name that makes sense
    """
    indices = []
    for i, row in df.iterrows():
        if str(row["chromosome"]) not in chr_list:
            indices.append(i)
    return indices

def each_row_has_something(df, it_with_columns):
    # Now make sure that each row has something in stringtie_id,
    #  isoseq_hq_id, pinfish_id, or isoseq_singleton_id
    # don't count rows that also have isoseq reads containing ID m64069
    df["one_row_one_gene"] = "none"
    for i, row in df.iterrows():
        hasone = False
        for colname in it_with_columns:
            if type(row[colname]) == str:
                if row[colname].strip().lower() != "":
                    hasone = True
        # what is this doing? I don't know exactly
        if type(row["comment"]) == str:
            for this_thing in ["m64069", "manual", "augustus"]:
                if this_thing in row["comment"].strip().lower():
                    hasone = True
        df.at[i,'one_row_one_gene'] = True

    t1 = df.loc[df['one_row_one_gene'] == False, ]
    if len(t1) != 0:
        print(t1, file=sys.stderr)
        raise IOError("the rows above don't have any annotations")
    return(df)

# now make sure that each pinfish file has a unique hash
def check_pinfish_ids(pinfish_list):
    ids = {}
    for thisfile in pinfish_list:
        if thisfile.endswith(".gff.gz"):
            f = gzip.open(thisfile, "rb")
            isgz = True
        elif thisfile.endswith(".gff"):
            f = open(thisfile, "r")
        else:
            raise Exception("""dunno what kind of file this is.
                        must be .gff or .gff.gz""")
        for line in f:
            nl = ""
            if isgz:
                nl = line.decode("utf-8")
            else:
                nl=line
            if "transcript" in nl and "pinfish" in nl:
                ID = line.split('\t')[8].split(';')[0].replace("ID=","")
                if ID in ids:
                    ids[ID] += 1
                else:
                    ids[ID] = 1
    assert len(np.unique(list(ids.values()))) == 1
    #yes, all the hashes are the same
