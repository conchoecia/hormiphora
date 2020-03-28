#!/bin/bash

cd ~/genomes/hormiphora_californensis/hormiphora/annotation/Hcv1a1d20200325_release
gzip -dc Hcv1a1d20200325.gff.gz > Hcv1a1d20200325.gff
gzip -dc Hcv1a1d20200325_model_proteins.pep.gz > Hcv1a1d20200325_model_proteins.pep

# from protein set
cd ~/genomes/hormiphora_californensis/hormiphora/synteny
ln -s ../annotation/Hcv1a1d20200325_release/Hcv1a1d20200325_model_proteins.pep
ln -s ../annotation/Hcv1a1d20200325_release/Hcv1a1d20200325.gff
./take_longest_hcal_protein.py Hcv1a1d20200325_model_proteins.pep > Hcv1a1d20200325_model_proteins.pep.longest_only.fa

# synteny with ML2

~/diamond-latest/diamond blastp -q Hcv1a1d20200325_model_proteins.pep.longest_only.fa -d ~/genomes/mnemiopsis_leidyi/ML2.2.aa -o Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab

~/diamond-latest/diamond blastp -q Hcv1a1d20200325_model_proteins.pep.longest_only.fa -d ~/genomes/mnemiopsis_leidyi/MLRB2.2.aa -o Hcv1a1d20200325_model_proteins.pep.vs_ml2unfiltered.tab

~/git/genomeGTFtools/scaffold_synteny.py -f ~/genomes/hormiphora_californensis/UCSC_Hcal_v1.fa -F ~/genomes/mnemiopsis_leidyi/MlScaffold09.nt.gz -q Hcv1a1d20200325.gff -d ~/genomes/mnemiopsis_leidyi/ML2.2.gene_only.gff -b Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab -G 20 -l 110 -L 60 > hcalv1_v_ml2_2d_synteny_points.tab

Rscript ~/git/genomeGTFtools/synteny_2d_plot.R hcalv1_v_ml2_2d_synteny_points.tab Mnemiopsis-leidyi Hormiphora-californensis

# locally sorted version

~/git/genomeGTFtools/scaffold_synteny.py -f ~/genomes/hormiphora_californensis/UCSC_Hcal_v1.fa -F ~/genomes/mnemiopsis_leidyi/MlScaffold09.nt.gz -q Hcv1a1d20200325.gff -d ~/genomes/mnemiopsis_leidyi/ML2.2.gene_only.gff -b Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab -G 20 -l 110 -L 109 --local-positions > hcalv1_v_ml2_2d_synteny_points_local_ml109.tab

Rscript synteny_2d_plot_w_2nd_genome_sorted.R

# microsynteny to ML2

~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200325.gff -d ~/genomes/mnemiopsis_leidyi/ML2.2.gene_only.gff -b Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab -G > hcalv1_v_ml2_microsynteny.gff 2> Hcv1a1d20200325_hcalv1_v_ml2_microsynteny.log

# microsynteny to Pbachei

~/diamond-latest/diamond blastp -q Hcv1a1d20200325_model_proteins.pep.longest_only.fa -d ~/genomes/pleurobrachia_bachei/pbachei_03_filtered_gene_models_transcripts_collapsed_only.prot.fa -o Hcv1a1d20200325_model_proteins.pep.vs_pbachei.tab

 ~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200325.gff -d ~/genomes/pleurobrachia_bachei/pbachei_03_filtered_gene_models_transcripts_collapsed_only.gff -b Hcv1a1d20200325_model_proteins.pep.vs_pbachei.tab -G > hcalv1_v_pbachei_microsynteny.gff 2> Hcv1a1d20200325_hcalv1_v_pbachei_microsynteny.log

# vs triad v2

~/diamond-latest/diamond blastp -q Hcv1a1d20200325_model_proteins.pep.longest_only.fa -d ~/genomes/trichoplax_adhaerens/triad_augustus_t1_only.prot.fasta -o Hcv1a1d20200325_model_proteins.pep.vs_triad.tab

~/git/genomeGTFtools/scaffold_synteny.py -f ~/genomes/hormiphora_californensis/UCSC_Hcal_v1.fa -F ~/genomes/trichoplax_adhaerens/Triad1_genomic_scaffolds.fasta -q Hcv1a1d20200325.gff -d ~/genomes/trichoplax_adhaerens/Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff -b Hcv1a1d20200325_model_proteins.pep.vs_triad.tab -G 20 -l 110 -L 100 --blast-db-delimiter "__" > hcalv1_v_triad_2d_synteny_points.tab

Rscript ~/git/genomeGTFtools/synteny_2d_plot.R hcalv1_v_triad_2d_synteny_points.tab Trichoplax-adhaerens Hormiphora-californensis

#~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200325.gff -d ~/genomes/trichoplax_adhaerens/Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff -b Hcv1a1d20200325_model_proteins.pep.vs_triad.tab -G --blast-db-delimiter "__" > hcalv1_v_triad_microsynteny.gff

~/git/genomeGTFtools/microsynteny.py -q Hcv1a1d20200325.gff -d ~/genomes/trichoplax_adhaerens/Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff -b Hcv1a1d20200325_model_proteins.pep.vs_triad.tab -G --blast-db-delimiter "__" -m 2 > hcalv1_v_triad_microsynteny_m2.gff 2> Hcv1a1d20200325_hcalv1_v_triad_microsynteny.log








