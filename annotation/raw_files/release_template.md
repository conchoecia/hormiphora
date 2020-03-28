The genome file is not included in this release `.tar.gz`. Download the genome file here: [UCSC_Hcal_v1.fa.gz](https://github.com/conchoecia/hormiphora/blob/master/annotation/raw_files/UCSC_Hcal_v1.fa.gz)

This release contains annotation and protein files for the Hcalv1 genome. Most likely you will use files:

- `RELEASEPREFIX_release/RELEASEPREFIX_model_proteins.pep.gz`
  - The model proteins for each transcript. NB - not all transcripts had CDS.
- `RELEASEPREFIX_release/RELEASEPREFIX_transcripts.fasta.gz`
  - Transcript files generated directly from the genome. May contain prematurely truncated CDS.
- `RELEASEPREFIX_release/RELEASEPREFIX.gff.gz`
  - Genome annotation of transcripts.
- `RELEASEPREFIX_release/protein_size_table_RELEASEPREFIX.csv`
  - A table showing the protein size differences in the within-transcript-phased transcript haplotypes, as well as which was selected for the model proteins.
- `RELEASEPREFIX_release/partly_phased/`
  - `RELEASEPREFIX_release/partly_phased/h1_pilon_RELEASEPREFIX.fasta.gz`
    - Pseudohapltype h1 of within-transcript-phased transcripts. Each transcript is derived from a single haplotype, but it is not phased with respect to all other transcripts in the genome.
  - `RELEASEPREFIX_release/partly_phased/h1_RELEASEPREFIX.pep.gz`
    - Putative proteins from the above fasta file.
  - `RELEASEPREFIX_release/partly_phased/h2_pilon_RELEASEPREFIX.fasta.gz`
    - Pseudohaplotype h2 of the within-transcript-phased transcripts
  - `RELEASEPREFIX_release/partly_phased/h2_RELEASEPREFIX.pep.gz`
    - Putative proteins from the above fasta file.
- `RELEASEPREFIX_release/phased/`
  - `RELEASEPREFIX_release/phased/RELEASEPREFIX_h1_phased_nucl.fasta.gz`
    - Transcripts that are from h1. Matches the whole-genome phased vcf file.
  - `RELEASEPREFIX_release/phased/RELEASEPREFIX_h1_phased_protein.pep.gz`
    - Proteins from the above file.
  - `RELEASEPREFIX_release/phased/RELEASEPREFIX_h2_phased_nucl.fasta.gz`
    - Transcripts that are from h2. Matches the whole-genome phased vcf file.
  - `RELEASEPREFIX_release/phased/RELEASEPREFIX_h2_phased_protein.pep.gz`
    - Proteins from the above file.
  - `RELEASEPREFIX_release/phased/transcripts_unique_to_h1.RELEASEPREFIX.list`
    - Transcripts that were able to be assigned to haplotype 1 (h1) of the whole-genome phasing. You probably won't need this file.
  - `RELEASEPREFIX_release/phased/transcripts_unique_to_h2.RELEASEPREFIX.list`
    - Same as above, but to h2. You probably won't need this file.
  - `RELEASEPREFIX_release/phased/transcripts_shared_by_both_should_be_empty.RELEASEPREFIX.list`
    - Intermediate check that no transcripts are shared by both haplotypes. Should be empty.
  - `RELEASEPREFIX_release/phased/second_list_of_transcripts_shared_by_both_should_be_empty.RELEASEPREFIX.list`
    - Final check that no transcripts are shared by both haplotypes. Should be empty.
