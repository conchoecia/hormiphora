#!/usr/bin/env python

import gzip


# I don't think I will actually use this
#class Transcript:
#    """
#    Instance attributes:
#     - filename
#        - The filename will be used to look up the transcript information.
#     - scaffold
#        - The scaffold is the name of the scaffold on which the transcript is
#          located. Like 'c1', 'c2', 'M', et cetera
#     - num_on_scaffold
#        - This is the index of the transcript on that scaffold. For instance
#          if this gene is the 117th gene that occurs when looking at the chromosome
#          from 5' to 3', then this number is 117
#     - ID
#        - The unique identifier of this gene, or specific isoform in the file
#          provided in filename. Depending on the filetype this will change how things
#          are parsed
#     - ftype
#        - The program that was used to generate the transcripts. Can be one of
#          the following types.
#          - "stringtie"
#          - "manual"
#          - "isoseq"
#          - "pinfish"
#        - We will use grep to find these IDs. This is useful because it allows
#          us to control if we want specific isoforms of a gene, or every single isoform.
#          This is only possible for stringtie and IsoSeq IDs.
#    """
#    def __init__(self, filename, scaffold, num_on_scaffold, ID, ftype):
#        self.filename = filename
#        self.scaffold = scaffold
#        self.num_on_scaffold = num_on_scaffold
#        self.ID = ID
#        self.ftype = ftype

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
            print(filename)
            raise Exception("Files should be passed in a list.")
        elif type(filename) == list:
            for thisfile in self.filename:
                isgz = False
                if thisfile.endswith(".gff.gz"):
                    f = gzip.open(thisfile, "rb")
                    isgz = True
                elif thisfile.endswith(".gff"):
                    f = open(thisfile, "r")
                else:
                    print("bad filename -> ", filename)
                    raise Exception("""dunno what kind of file this is.
                                must be .gff or .gff.gz""")
                for line in f:
                    nl = ""
                    if isgz:
                        nl = line.decode("utf-8")
                    else:
                        nl=line
                    splitd = nl.split('\t')
                    if splitd[2] == "transcript":
                        # we have just found a new transcript.
                        # make sure that the input is legal
                        if splitd[1] in  ["pinfish", "StringTie"]:
                            if not splitd[8].split(';')[0].startswith("ID="):
                                # should start with ID=
                                print(line)
                                raise Exception("""There is some input error. We found
                                a line that doesn't have field 9 starting with ID=.
                                all pinfish/StringTie transcripts start with this""")
                        elif splitd[1] in ["PacBio"]:
                            if not splitd[8].split(';')[0].startswith("gene_id"):
                                # should start with gene_id
                                print(line)
                                raise Exception("""There is some input error. We found
                                a line that doesn't have field 9 starting with gene_id.
                                all PacBio transcripts start with this""")

                        # now that we made sure the input is legal, let's parse
                        #  the gene id and the transcript ID.
                        if splitd[1] == "pinfish":
                            tID = splitd[8].split(';')[0].replace("ID=","").strip()
                            gID = tID
                        elif splitd[1] == "StringTie":
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
                        elif splitd[1] in ["StringTie", "PacBio"]:
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
                                print(line)
                                raise Exception("""There is some input error. We found
                              a line that doesn't have field 9 starting with Parent=.
                                all pinfish exons start with this""")
                        elif splitd[1] in  ["StringTie"]:
                            if not splitd[8].split(';')[0].startswith("ID="):
                                print(line)
                                raise Exception("""There is some input error. We found
                              a line that doesn't have field 9 starting with ID=.
                                all pinfish exons start with this""")
                        elif splitd[1] in  ["PacBio"]:
                            if not splitd[8].split(';')[0].startswith("gene_id \""):
                                print(line)
                                raise Exception("""There is some input error. We found
                              a line that doesn't have field 9 starting with ID=.
                                all pinfish exons start with this""")


                        # now that we made sure the input is legal, let's parse
                        #  the gene id and the transcript ID.
                        if splitd[1] == "pinfish":
                            tID = splitd[8].split(';')[0].replace("Parent=","").strip()
                        elif splitd[1] == "StringTie":
                            tID = splitd[8].split(';')[1].replace("Parent=","").strip()
                        elif splitd[1] == "PacBio":
                            tID = splitd[8].split(';')[1].replace("transcript_id","").strip().replace("\"", "")

                        # Now that we have the transcript ID
                        #  store them in the instance
                        if tID not in self.IDTS:
                            # we should have already seen the transcript
                            #  if the gff file is sorted properly
                            print("offending ID: ", tID)
                            print("offending file:", thisfile)
                            raise Exception("""For some reason we found an exon for a
                            transcript before we found the transcript itself.
                            The GFF file should have all of the transcripts first.""")
                        self.IDTS[tID] += nl

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
            if 'm64069' in row["comment"].strip().lower():
                C5_CO = True
        df.at[i,'one_row_one_gene'] = C1_ST or C2_IS or C3_PF or C4_SI or C5_CO

    t1 = df.loc[df['one_row_one_gene'] == False, ]
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
    print(t1)
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
                print(line.split('\t'))
                sys.exit()
                ID = line.split('\t')[8].split(';')[0].replace("ID=","")
                if ID in ids:
                    ids[ID] += 1
                else:
                    ids[ID] = 1
    assert len(np.unique(list(ids.values()))) == 1
    #yes, all the hashes are the same
