#!/bin/bash
THREADS=90
REF=/bigdata/user/darrin/hormiphora/assembly_versions/UCSC_Hcal_v1.fa
GFF=../Hcalv1_annotation1_20190124.gff
FLNC=/bigdata/raw/rna_reads/PacBio/PB405_1plex_IsoSeq/final_bam/GLO64_flnc.fastq.gz
PILONJAR=/usr/local/bin/pilon-1.22.jar
PREFIX=Jan24_annot
RNAF=/bigdata/make_txomes/reads/CTE_Hormiphora_californensis_D914_T1_f.fastq.gz
RNAR=/bigdata/make_txomes/reads/CTE_Hormiphora_californensis_D914_T1_r.fastq.gz



#### first generate a sequence for each isoform
#gffread -w ${PREFIX}.transcripts.fasta -g ${REF} ${GFF}
#bwa index ${PREFIX}.transcripts.fasta
#
### now map illumina to each transcript
#bwa mem -t ${THREADS} ${PREFIX}.transcripts.fasta ${RNAF} ${RNAR} | \
#	samtools view -hb -@ ${THREADS} - | \
#        samtools sort -@ ${THREADS} > ILL_to_${PREFIX}.sorted.bam
#samtools index ILL_to_${PREFIX}.sorted.bam
#
## now map flnc to each transcript
#minimap2 -ax asm20 -t ${THREADS} ${PREFIX}.transcripts.fasta ${FLNC} | \
#	samtools view -hb -@ ${THREADS} - | \
#	samtools sort -@ ${THREADS} > flnc_to_${PREFIX}.sorted.bam
#samtools index flnc_to_${PREFIX}.sorted.bam
#
#
## merge bams
#samtools merge ILL_and_flnc.merged.bam ILL_to_${PREFIX}.sorted.bam flnc_to_${PREFIX}.sorted.bam
#samtools index ILL_and_flnc.merged.bam
#
# call a vcf
#samtools faidx ${PREFIX}.transcripts.fasta
#freebayes-parallel <(freebayes_fasta_generate_regions.py ${PREFIX}.transcripts.fasta.fai 100000) ${THREADS} \
#	    -f ${PREFIX}.transcripts.fasta ILL_and_flnc.merged.bam > calls.vcf
#bcftools mpileup -Ou -f ${PREFIX}.transcripts.fasta ILL_and_flnc.merged.bam | bcftools call -mv -Ov -o calls.vcf

#whatshap phase --indels --ignore-read-groups \
#	--reference ${PREFIX}.transcripts.fasta \
#	-o whatshap_phased.vcf \
#	calls.vcf ILL_to_${PREFIX}.sorted.bam flnc_to_${PREFIX}.sorted.bam
#bgzip -c whatshap_phased.vcf > whatshap_phased.vcf.gz
#tabix -p vcf whatshap_phased.vcf.gz
#whatshap haplotag --ignore-read-groups -o flnc_haplotagged.bam \
#	--reference ${PREFIX}.transcripts.fasta whatshap_phased.vcf.gz flnc_to_${PREFIX}.sorted.bam

## generate two references 
#for INDEX in 1 2
#do
#    bcftools consensus -H ${INDEX} -f ${PREFIX}.transcripts.fasta whatshap_phased.vcf.gz > ${PREFIX}.haplotype${INDEX}.fasta
#done

## now get one fastq file per haplotig
#for INDEX in 1 2
#do
#    samtools view -F 4 flnc_haplotagged.bam | \
#	 grep "HP:i:${INDEX}" | cut -f1 | sort | uniq > h${INDEX}.reads.txt
#    seqtk subseq ${FLNC} h${INDEX}.reads.txt >  h${INDEX}.flnc.fastq
#    gzip -f h${INDEX}.flnc.fastq
#done

# now map the reads to the haplotypes and correct with pilon
#for INDEX in 1 2
#do
#    minimap2 -ax asm20 -t ${THREADS} \
#	${PREFIX}.haplotype${INDEX}.fasta h${INDEX}.flnc.fastq.gz | \
#	samtools view -hb -@ ${THREADS} - | \
#	samtools sort -@ ${THREADS} - > h${INDEX}flnc_to_${PREFIX}h${INDEX}.sorted.bam 
#    samtools index h${INDEX}flnc_to_${PREFIX}h${INDEX}.sorted.bam
#    java -Xmx100G -jar ${PILONJAR} --genome ${PREFIX}.haplotype${INDEX}.fasta \
#        --unpaired h${INDEX}flnc_to_${PREFIX}h${INDEX}.sorted.bam \
#        --outdir h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon \
#        --output h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon \
#        --threads ${THREADS}
#done

# just rename the sequences
#for INDEX in 1 2
#do
#    sed -i "s/_pilon/.h${INDEX}/g" h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon.fasta
#done

# now map the FLNC reads back for QC
#for INDEX in 1 2
#do
#    minimap2 -ax asm20 -t ${THREADS} \
#        h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon.fasta \
#	h${INDEX}.flnc.fastq.gz | \
#	samtools view -hb -@ ${THREADS} - | \
#	samtools sort -@ ${THREADS} - > h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_h1_pilon.sorted.bam 
#    samtools index h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_h1_pilon.sorted.bam
#done

# now run Transdecoder.LongOrfs on the top strands
for INDEX in 1 2
do
    TransDecoder.LongOrfs \
        -t h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon.fasta \
	-m 50 -S
done


