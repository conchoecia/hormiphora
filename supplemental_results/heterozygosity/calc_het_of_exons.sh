#!/bin/bash
GFF=../../annotation/Hcv1a1d20200309_release/Hcv1a1d20200309.gff.gz
# Map the shotgun reads to the genome with bwa mem to generate this bam file
BAM=/path/to/shot_to_assem.sorted.bam
REF=../annotation/raw_files/UCSC_Hcal_v1.fa.gz
# download chep here https://github.com/conchoecia/chep
PTA=/path/to/chep_pileup_to_array
#CP=/home/dschultz/git/chep/bin/chep_plot

cp ${REF} ./
gunzip UCSC_Hcal_v1.fa.gz
REF=UCSC_Hcal_v1.fa

#get the exons
rm exons.bed
zcat ${GFF} | grep 'exon' | awk '{printf("%s\t%d\t%d\n", $1, $4 -1, $5 -1)}' > exons.bed
rm exons.sorted.bed
sort -k1,1 -k2n,2 exons.bed > exons.sorted.bed
rm exons.final.sorted.bed
bedtools merge -i exons.sorted.bed > exons.final.sorted.bed

RAW=mpileup_depth_v_ref.txt
samtools mpileup -f ${REF} -l exons.final.sorted.bed ${BAM} | \
  ${PTA} > ${RAW}
chep_plot -f mpileup_depth_v_ref.txt -x 0 -X 250
