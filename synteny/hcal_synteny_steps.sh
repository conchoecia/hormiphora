#!/bin/bash

# from protein set

./take_longest_hcal_protein.py Hcv1a1d20200309_model_proteins.pep > Hcv1a1d20200309_model_proteins.pep.longest_only.fa

# synteny with ML2

~/diamond-latest/diamond blastp -q Hcv1a1d20200309_model_proteins.pep.longest_only.fa -d ~/genomes/mnemiopsis_leidyi/ML2.2.aa -o Hcv1a1d20200309_model_proteins.pep.vs_ml2.tab

~/diamond-latest/diamond blastp -q Hcv1a1d20200309_model_proteins.pep.longest_only.fa -d ~/genomes/mnemiopsis_leidyi/MLRB2.2.aa -o Hcv1a1d20200309_model_proteins.pep.vs_ml2unfiltered.tab

~/git/genomeGTFtools/scaffold_synteny.py -f ~/genomes/hormiphora_californensis/UCSC_Hcal_v1.fa -F ~/genomes/mnemiopsis_leidyi/MlScaffold09.nt.gz -q Hcv1a1d20200309.gff -d ~/genomes/mnemiopsis_leidyi/ML2.2.gene_only.gff -b Hcv1a1d20200309_model_proteins.pep.vs_ml2.tab -G 20 -l 110 -L 60 > hcalv1_v_ml2_2d_synteny_points.tab

Rscript ~/git/genomeGTFtools/synteny_2d_plot.R hcalv1_v_ml2_2d_synteny_points.tab Mnemiopsis-leidyi Hormiphora-californensis

# locally sorted version

~/git/genomeGTFtools/scaffold_synteny.py -f ~/genomes/hormiphora_californensis/UCSC_Hcal_v1.fa -F ~/genomes/mnemiopsis_leidyi/MlScaffold09.nt.gz -q Hcv1a1d20200309.gff -d ~/genomes/mnemiopsis_leidyi/ML2.2.gene_only.gff -b Hcv1a1d20200309_model_proteins.pep.vs_ml2.tab -G 20 -l 110 -L 60 --local-positions > hcalv1_v_ml2_2d_synteny_points_local.tab

Rscript synteny_2d_plot_w_2nd_genome_sorted.R

# microsynteny to ML2

~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200309.gff -d ~/genomes/mnemiopsis_leidyi/ML2.2.gene_only.gff -b Hcv1a1d20200309_model_proteins.pep.vs_ml2.tab -G > hcalv1_v_ml2_microsynteny.gff

# Parsing Hcv1a1d20200309.gff  Sat Mar 14 13:26:03 2020
# Found 29925 genes  Sat Mar 14 13:26:03 2020
# Parsing /home/wrf/genomes/mnemiopsis_leidyi/ML2.2.gene_only.gff  Sat Mar 14 13:26:03 2020
# Found 16548 genes  Sat Mar 14 13:26:03 2020
# Parsing tabular blast output Hcv1a1d20200309_model_proteins.pep.vs_ml2.tab  Sat Mar 14 13:26:03 2020
# Found blast hits for 9724 query sequences  Sat Mar 14 13:26:03 2020
# Removed 797 hits by evalue, kept 51733 hits
# Names parsed as Hcv1.1.sca25.g1.i1 from Hcv1.1.sca25.g1.i1, and ML09422a from ML09422a
# make GFF output: True
# searching for colinear blocks of at least 3 genes, with up to 5 intervening genes
# Found 283 possible split genes  Sat Mar 14 13:26:03 2020
# Most genes on a query scaffold was 3823  Sat Mar 14 13:26:03 2020
# Found 564 total putative synteny blocks for 2229 genes  Sat Mar 14 13:26:03 2020
# Average block is 3.95, longest block was 12 genes  Sat Mar 14 13:26:03 2020
# Total block span was 13105231 bases  Sat Mar 14 13:26:03 2020
#3 303
#4 128
#5 57
#6 46
#7 13
#8 5
#9 7
#10 3
#11 1
#12 1


# microsynteny to Pbachei

~/diamond-latest/diamond blastp -q Hcv1a1d20200309_model_proteins.pep.longest_only.fa -d ~/genomes/pleurobrachia_bachei/pbachei_03_filtered_gene_models_transcripts_collapsed_only.prot.fa -o Hcv1a1d20200309_model_proteins.pep.vs_pbachei.tab

 ~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200309.gff -d ~/genomes/pleurobrachia_bachei/pbachei_03_filtered_gene_models_transcripts_collapsed_only.gff -b Hcv1a1d20200309_model_proteins.pep.vs_pbachei.tab -G > hcalv1_v_pbachei_microsynteny.gff

# Parsing Hcv1a1d20200309.gff  Sat Mar 14 13:57:56 2020
# Found 29925 genes  Sat Mar 14 13:57:56 2020
# Parsing /home/wrf/genomes/pleurobrachia_bachei/pbachei_03_filtered_gene_models_transcripts_collapsed_only.gff  Sat Mar 14 13:57:56 2020
# Found 7014 genes  Sat Mar 14 13:57:56 2020
# Parsing tabular blast output Hcv1a1d20200309_model_proteins.pep.vs_pbachei.tab  Sat Mar 14 13:57:56 2020
# Found blast hits for 6151 query sequences  Sat Mar 14 13:57:56 2020
# Removed 412 hits by evalue, kept 27663 hits
# Names parsed as Hcv1.1.sca23.g1.i1 from Hcv1.1.sca23.g1.i1, and scaffold907_1_size33252_gene_8440Barcelona from scaffold907_1_size33252_gene_8440Barcelona
# make GFF output: True
# searching for colinear blocks of at least 3 genes, with up to 5 intervening genes
# Found 402 possible split genes  Sat Mar 14 13:57:56 2020
# Most genes on a query scaffold was 3823  Sat Mar 14 13:57:56 2020
# Found 297 total putative synteny blocks for 1269 genes  Sat Mar 14 13:57:56 2020
# Average block is 4.27, longest block was 15 genes  Sat Mar 14 13:57:56 2020
# Total block span was 9562252 bases  Sat Mar 14 13:57:56 2020
#3 128
#4 78
#5 45
#6 20
#7 8
#8 5
#9 7
#10 2
#11 2
#12 1
#15 1


# vs triad v2

~/diamond-latest/diamond blastp -q Hcv1a1d20200309_model_proteins.pep.longest_only.fa -d ~/genomes/trichoplax_adhaerens/triad_augustus_t1_only.prot.fasta -o Hcv1a1d20200309_model_proteins.pep.vs_triad.tab

~/git/genomeGTFtools/scaffold_synteny.py -f ~/genomes/hormiphora_californensis/UCSC_Hcal_v1.fa -F ~/genomes/trichoplax_adhaerens/Triad1_genomic_scaffolds.fasta -q Hcv1a1d20200309.gff -d ~/genomes/trichoplax_adhaerens/Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff -b Hcv1a1d20200309_model_proteins.pep.vs_triad.tab -G 20 -l 110 -L 100 --blast-db-delimiter "__" > hcalv1_v_triad_2d_synteny_points.tab

Rscript ~/git/genomeGTFtools/synteny_2d_plot.R hcalv1_v_triad_2d_synteny_points.tab Trichoplax-adhaerens Hormiphora-californensis

#~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200309.gff -d ~/genomes/trichoplax_adhaerens/Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff -b Hcv1a1d20200309_model_proteins.pep.vs_triad.tab -G --blast-db-delimiter "__" > hcalv1_v_triad_microsynteny.gff

~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200309.gff -d ~/genomes/trichoplax_adhaerens/Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff -b Hcv1a1d20200309_model_proteins.pep.vs_triad.tab -G --blast-db-delimiter "__" -m 2 > hcalv1_v_triad_microsynteny_m2.gff








