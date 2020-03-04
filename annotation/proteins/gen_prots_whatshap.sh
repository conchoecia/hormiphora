#!/bin/bash
THREADS=45
REF=/bigdata/user/darrin/hormiphora/assembly_versions/UCSC_Hcal_v1.fa
GFF=../Hcalv1_annotation1_20190226.gff
FLNC=/bigdata/raw/rna_reads/PacBio/PB405_1plex_IsoSeq/final_bam/GLO64_flnc.fastq.gz
PILONJAR=/usr/local/bin/pilon-1.22.jar
PREFIX=Feb_24_annot
RNAF=/bigdata/make_txomes/reads/CTE_Hormiphora_californensis_D914_T1_f.fastq.gz
RNAR=/bigdata/make_txomes/reads/CTE_Hormiphora_californensis_D914_T1_r.fastq.gz
GVCF=/bigdata/user/darrin/hormiphora/phasing/attempt2/final_output/largest_blocks.hap.phased.VCF


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
## merge bams
#samtools merge ILL_and_flnc.merged.bam ILL_to_${PREFIX}.sorted.bam flnc_to_${PREFIX}.sorted.bam
#samtools index ILL_and_flnc.merged.bam

## call a vcf
#samtools faidx ${PREFIX}.transcripts.fasta
#freebayes-parallel <(freebayes_fasta_generate_regions.py ${PREFIX}.transcripts.fasta.fai 100000) ${THREADS} \
#	    -f ${PREFIX}.transcripts.fasta ILL_and_flnc.merged.bam > calls.vcf
#
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
#
## now get one fastq file per haplotig
#for INDEX in 1 2
#do
#    samtools view -F 4 flnc_haplotagged.bam | \
#	 grep "HP:i:${INDEX}" | cut -f1 | sort | uniq > h${INDEX}.reads.txt
#    seqtk subseq ${FLNC} h${INDEX}.reads.txt >  h${INDEX}.flnc.fastq
#    gzip -f h${INDEX}.flnc.fastq
#done
#
## now map the reads to the haplotypes and correct with pilon
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
#
## just rename the sequences
#for INDEX in 1 2
#do
#    sed -i "s/_pilon/.h${INDEX}/g" h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon.fasta
#done
#
## now map the FLNC reads back for QC
#for INDEX in 1 2
#do
#    minimap2 -ax asm20 -t ${THREADS} \
#        h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon.fasta \
#	h${INDEX}.flnc.fastq.gz | \
#	samtools view -hb -@ ${THREADS} - | \
#	samtools sort -@ ${THREADS} - > h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_h1_pilon.sorted.bam 
#    samtools index h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_h1_pilon.sorted.bam
#done
#
#rm -rf prottrans_predict
#mkdir prottrans_predict
### now run prottrans.py
#for INDEX in 1 2
#do
#    python2 prottrans.py -r -a 50 h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon.fasta | cut -d'_' -f1 > prottrans_predict/h${INDEX}.prottrans.fasta
#done
#
#python make_table_of_genes.py
#gzip prottrans_predict/Hcv1.1.pep

## this bit is all for phasing the transcripts
## map haplotype1 and haplotype2 back to the genome
#mkdir phased_to_ref
#for INDEX in 1 2
#do
#    minimap2 -ax splice -t ${THREADS} ${REF} \
#        h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon/h${INDEX}flnc_to_${PREFIX}h${INDEX}.pilon.fasta | \
#        samtools view -hb -@ ${THREADS} - | \
#        samtools sort -@ ${THREADS} - > phased_to_ref/pilon_tx_${INDEX}_to_ref.sorted.bam
#    samtools index phased_to_ref/pilon_tx_${INDEX}_to_ref.sorted.bam
#done

# merge the bams
samtools merge -f phased_to_ref/phased_transcripts_to_ref.bam phased_to_ref/pilon_tx_1_to_ref.sorted.bam phased_to_ref/pilon_tx_2_to_ref.sorted.bam
samtools index phased_to_ref/phased_transcripts_to_ref.bam

cp ${GVCF} temp_whole_genome.vcf
bgzip -c temp_whole_genome.vcf > temp_whole_genome.vcf.gz
tabix -p vcf temp_whole_genome.vcf.gz

whatshap haplotag --ignore-read-groups -o phased_to_ref/all_transcripts_haplotag.bam \
    --reference ${REF} temp_whole_genome.vcf.gz phased_to_ref/phased_transcripts_to_ref.bam

# now make a list of things that were phased
samtools view phased_to_ref/all_transcripts_haplotag.bam | grep 'HP:i:1' | cut -f1 | sort | uniq | sort > phased_to_ref/haplotype0_transcript.list
samtools view phased_to_ref/all_transcripts_haplotag.bam | grep 'HP:i:2' | cut -f1 | sort | uniq | sort > phased_to_ref/haplotype1_transcript.list

# a list of transcripts that are shared by both haplotypes. Should be nothing
comm -12 phased_to_ref/haplotype0_transcript.list phased_to_ref/haplotype1_transcript.list > phased_to_ref/transcripts_shared_by_both.list

# need to make a python script to parse this
# now get the transcripts that were able to be assigned uniquely to one haplotype or the other
for INDEX in 0 1
do
    cat phased_to_ref/haplotype${INDEX}_transcript.list | \
        python transcript_list_unique.py > phased_to_ref/tx_phased_to_h${INDEX}.list
done

# double check a list of transcripts that are shared by both haplotypes. Should be nothing
comm -12 phased_to_ref/haplotype0_transcript.list phased_to_ref/haplotype1_transcript.list > phased_to_ref/transcripts_shared_by_both_2.list

mkdir final_phased_transcripts
# generate transcript-specific haplotypes
for INDEX in 0 1
do
    OIND=1
    seqtk subseq h${OIND}flnc_to_${PREFIX}h${OIND}.pilon/h${OIND}flnc_to_${PREFIX}h${OIND}.pilon.fasta phased_to_ref/tx_phased_to_h${INDEX}.list > final_phased_transcripts/tx_phased_to_h${INDEX}.fasta
    seqtk subseq prottrans_predict/h${OIND}.prottrans.fasta phased_to_ref/tx_phased_to_h${INDEX}.list > final_phased_transcripts/tx_phased_to_h${INDEX}.pep
    OIND=2
    seqtk subseq h${OIND}flnc_to_${PREFIX}h${OIND}.pilon/h${OIND}flnc_to_${PREFIX}h${OIND}.pilon.fasta phased_to_ref/tx_phased_to_h${INDEX}.list >> final_phased_transcripts/tx_phased_to_h${INDEX}.fasta
    seqtk subseq prottrans_predict/h${OIND}.prottrans.fasta phased_to_ref/tx_phased_to_h${INDEX}.list >> final_phased_transcripts/tx_phased_to_h${INDEX}.pep
done

#now fix the headers
mkdir final_phased_transcripts
# generate transcript-specific haplotypes
for INDEX in 0 1
do
    TDIR=final_phased_transcripts
    cut -d'.' -f1,2,3,4,5 ${TDIR}/tx_phased_to_h${INDEX}.fasta | bioawk -v var=${INDEX} -cfastx '{printf(">%s.h%s\n%s\n", $name, var, $seq)}' > ${TDIR}/tx_phased_to_h${INDEX}.correct_headers.fasta
    cut -d'.' -f1,2,3,4,5 ${TDIR}/tx_phased_to_h${INDEX}.pep   | bioawk -v var=${INDEX} -cfastx '{printf(">%s.h%s\n%s\n", $name, var, $seq)}' > ${TDIR}/tx_phased_to_h${INDEX}.correct_headers.pep
    cut -d'.' -f1,2,3,4,5 ${TDIR}/tx_phased_to_h${INDEX}.fasta | bioawk -v var=${INDEX} -cfastx '{printf(">%s.h%s\n%s\n", $name, var, $seq)}' > ${TDIR}/tx_phased_to_h${INDEX}.correct_headers.fasta
    cut -d'.' -f1,2,3,4,5 ${TDIR}/tx_phased_to_h${INDEX}.pep   | bioawk -v var=${INDEX} -cfastx '{printf(">%s.h%s\n%s\n", $name, var, $seq)}' > ${TDIR}/tx_phased_to_h${INDEX}.correct_headers.pep
done
