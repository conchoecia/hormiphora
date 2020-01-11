#!/usr/bin/env python

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

        #IDTS = transcript ID To String
        self.IDTS = {}
        self.filename = filename
        if type(filename) == list and filetype == pinfish:
            for thisfile in self.filename:
                isgz = False
                if thisfile.endswith(".gff.gz"):
                    f = gzip.open(thisfile, "rb")
                    isgz = True
                elif thisfile.endswith(".gff"):
                    f = open(thisfile, "r")
                else:
                    raise Error("""dunno what kind of file this is.
                                must be .gff or .gff.gz""")
                for line in f:
                    nl = ""
                    if isgz:
                        nl = line.decode("utf-8")
                    else:
                        nl=line
                    splitd = line.split('\t')
                    if splitd[1] == "pinfish" and splitd[2] == "transcript":
                        # we have just found a new transcript in pinfish
                        if not splitd[8].split(';')[0].startswith("ID="):
                            #make sure the input is legal
                            # should start with ID=
                            print(line)
                            raise Error("""There is some input error. We found
                            a line that doesn't have field 9 starting with ID=.
                            all pinfish transcripts start with this""")
                        ID = splitd[8].split(';')[0].replace("ID=","")
                        if ID in self.IDTS:
                            raise Error("""This transcript is already in the map.
                                        we shouldn't see it here yet.""")
                        #just add the line to the dict. we'll parse it later
                        self.IDTS[ID] = nl
                    elif splitd[1] == "pinfish" and splitd[2] == "exon":
                        # we have found an exon for this transcript. We should
                        #  have already found the transcript itself.
                        if not splitd[8].split(';')[0].startswith("Parent="):
                            #make sure the input is legal
                            # should start with Parent=
                            print(line)
                            raise Error("""There is some input error. We found
                          a line that doesn't have field 9 starting with Parent=.
                            all pinfish exons start with this""")
                        ID = splitd[8].split(';')[0].replace("Parent=","")
                        if ID not in self.IDTS:
                            # we should have already seen the transcript
                            #  if the gff file is sorted properly
                            raise Error("For some reason we found an exon for a
                            transcript before we found the transcript itself.
                            The GFF file should have all of the transcripts first.")
                        self.IDTS[ID] = nl


