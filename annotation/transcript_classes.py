#!/usr/bin/env python
import gzip
import os
import pandas as pd
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
                    splitd = nl.split('\t')
                    if str(splitd[2]).strip() in ["transcript", "mRNA"]:
                        # we have just found a new transcript.
                        # make sure that the input is legal
                        if splitd[1] in  ["pinfish", "StringTie"]:
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

                        # now that we made sure the input is legal, let's parse
                        #  the gene id and the transcript ID.
                        if splitd[1] == "pinfish":
                            tID = splitd[8].split(';')[0].replace("ID=","").strip()
                            gID = tID
                        elif splitd[1] in ["StringTie", "AUGUSTUS"]:
                            tID = splitd[8].split(';')[0].replace("ID=","").strip()
                            gID = ".".join(tID.split('.')[0:-1])
                        elif splitd[1] == "PacBio":
                            tID = splitd[8].split(';')[1].replace("transcript_id","").strip().replace("\"", "")
                            gID = splitd[8].split(';')[0].replace("gene_id","").strip().replace("\"", "")
                        if tID in self.IDTS:
                            raise Exception("""This transcript is already in the map.
                                        we shouldn't see it here yet.""")
                        # Now that we have the transcript ID and the geneID
                        #  store them in the instance
                        if splitd[1] == "pinfish":
                            # every pinfish transcript is also its own gene.
                            #  each transcript just has its own hash.
                            assert tID not in self.IDTS
                            assert gID not in self.GTT
                            self.GTT[gID] = [tID]
                        elif splitd[1] in ["StringTie", "PacBio", "AUGUSTUS"]:
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
                        elif splitd[1] in  ["StringTie"]:
                            if not splitd[8].split(';')[0].startswith("ID="):
                                print(line, file=sys.stderr)
                                raise Exception("""There is some input error. We found
                              a line that doesn't have field 9 starting with ID=.
                                all pinfish exons start with this""")
                        elif splitd[1] in  ["PacBio"]:
                            if not splitd[8].split(';')[0].startswith("gene_id \""):
                                print(line, file=sys.stderr)
                                raise Exception("""There is some input error. We found
                              a line that doesn't have field 9 starting with ID=.
                                all pinfish exons start with this""")
                        elif splitd[1] in ["AUGUSTUS"]:
                            # I don't want to implement for a single gene
                            pass

                        # now that we made sure the input is legal, let's parse
                        #  the gene id and the transcript ID.
                        if splitd[1] == "pinfish":
                            tID = splitd[8].split(';')[0].replace("Parent=","").strip()
                        elif str(splitd[1]).strip() in ["StringTie", "AUGUSTUS"]:
                            tID = splitd[8].split(';')[1].replace("Parent=","").strip()
                        elif splitd[1] == "PacBio":
                            tID = splitd[8].split(';')[1].replace("transcript_id","").strip().replace("\"", "")

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
    """
    this_chromosome = ""
    transcript_counter = ""
    isoform_counter = 0
    for i, row in df.iterrows():
        #first, determine which chromosome we are on, and what transcript
        row_chr = str(row["chromosome"])
        if row_chr != this_chromosome:
            #we have either just started, or transitioned to a new chromosome
            this_chromosome = row_chr
            transcript_counter = 0

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
        transcripts_in_this_gene = {CTGm[key]: [] for key in CTGm}
        for key in CTGm:
            skip = False
            if key == "stringtie_id":
                if str(row["remove_st"]).lower() in ["y", "yes"]:
                    skip = True
                    if "B1_LR" not in str(row["stringtie_id"]):
                        print("rejecting: ", row["stringtie_id"], file=sys.stderr)
            if not skip and not pd.isnull(row[key]):
                # we have found a gene for this GFF file
                splitd = str(row[key]).split(",")
                for tx in splitd:
                    transcripts_in_this_gene[CTGm[key]].append(tx.strip())
            if len(transcripts_in_this_gene[CTGm[key]]) ==  0:
                transcripts_in_this_gene.pop(CTGm[key])
        # we now have all the transcripts_in_this_gene. now print them out
        if not DoL_empty(transcripts_in_this_gene):
            #there are some transcripts here.
            transcript_counter += 1
            this_transcript = "Hcv1.1.{}.t{}".format(this_chromosome, transcript_counter)
            #add them to the gene and print
            isoform_counter = 1
            for key in transcripts_in_this_gene:
                for this_transcript_ID in transcripts_in_this_gene[key]:
                    #the transcript ID might be for a single isoform, or it might
                    #  be for a bunch of isoforms. To determine which, we will
                    #  look in the GFF file's GTT instance attribute (Gene-to-transcript)
                    look_these_up = []
                    if this_transcript_ID in GFFs[key].GTT:
                        #This is a gene and there are other transcripts to look up
                        look_these_up = GFFs[key].GTT[this_transcript_ID]
                    elif this_transcript_ID in GFFs[key].IDTS:
                        #this is a single transcript that we're adding
                        look_these_up = [this_transcript_ID]
                    else:
                        #We should have encountered something. Raise an error if
                        #  we didn't find it.
                        print(transcripts_in_this_gene)
                        print("We couldn't find: ", this_transcript_ID, " in ", key, file=sys.stderr)
                        raise Exception("Couldn't find the transcript in the GFF file object")
                    #now that we have a list of isoforms to add to this gene,
                    #  start printing things out
                    for lookup in look_these_up:
                        #HCv1.1.c1.0000
                        this_isoform = "{}.i{}".format(this_transcript,isoform_counter) 
                        iso_string = GFFs[key].IDTS[lookup]
                        lines_split = iso_string.split('\n')
                        # for now just print everything out
                        exon_counter=1
                        for line in lines_split:
                            if line.strip():
                                gff_split = line.split('\t')
                                #get the info if the gene is spliced in an intron
                                SII=""
                                if not pd.isnull(row["spliced_in_intron"]):
                                    SII="y"
                                else:
                                    SII="n"

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
                                    comment="ID={0};Name={0}".format(this_isoform)
                                elif gff_split[2] == "exon":
                                    comment="ID={0}.e{1};Parent={0}".format(this_isoform,
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
                                print("\t".join(gff_split))
                        isoform_counter += 1

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
    assert len(t1) == 0

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

def each_row_has_something(df):
    # Now make sure that each row has something in stringtie_id,
    #  isoseq_hq_id, pinfish_id, or isoseq_singleton_id
    # don't count rows that also have isoseq reads containing ID m64069
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
            for this_thing in ["m64069", "manual", "augustus"]:
                if this_thing in row["comment"].strip().lower():
                    C5_CO = True

        df.at[i,'one_row_one_gene'] = C1_ST or C2_IS or C3_PF or C4_SI or C5_CO

    t1 = df.loc[df['one_row_one_gene'] == False, ]
    print(t1, file=sys.stderr)
    assert len(t1) == 0
    return(df)

# print out list of genes that still need a transcript
def still_needs_transcript(df):
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
    print(t1, file=sys.stderr )
    # just skip this for now since it's in-progress
    #assert len(t1) == 0
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
