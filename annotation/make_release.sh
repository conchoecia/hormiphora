#!/bin/bash
RELEASE=Hcv1a1d20200324
mkdir ${RELEASE}_release
find ./final_output/ -type f -name "*" -exec cp {} ${RELEASE}_release \;
rm ${RELEASE}_release/*vcf*
rm ${RELEASE}_release/*.amb
rm ${RELEASE}_release/*.ann
rm ${RELEASE}_release/*.bwt
rm ${RELEASE}_release/*.fai
rm ${RELEASE}_release/*.pac
rm ${RELEASE}_release/*.sa

# gzip fasta
for file in ./${RELEASE}_release/*.fasta
do
    echo "gzipping ${file}"
    gzip ${file}
done

# gzip pep
for file in ./${RELEASE}_release/*.pep
do
    echo "gzipping ${file}"
    gzip ${file}
done

# gzip gff
for file in ./${RELEASE}_release/*.gff
do
    gzip ${file}
done

#move phase files to extra
mkdir ${RELEASE}_release/phased
mv ${RELEASE}_release/*_phased_* ${RELEASE}_release/phased/
mv ${RELEASE}_release/*empty* ${RELEASE}_release/phased/
mv ${RELEASE}_release/*transcripts_unique_to_h* ${RELEASE}_release/phased/

mkdir ${RELEASE}_release/partly_phased
mv ${RELEASE}_release/h1_*.pep.gz ${RELEASE}_release/partly_phased/
mv ${RELEASE}_release/h2_*.pep.gz ${RELEASE}_release/partly_phased/
mv ${RELEASE}_release/h1_*.fasta.gz ${RELEASE}_release/partly_phased/
mv ${RELEASE}_release/h2_*.fasta.gz ${RELEASE}_release/partly_phased/

# add the relase notes
cat raw_files/release_template.md | sed 's/RELEASEPREFIX/${RELEASE}' > ${RELEASE}_release/README_${RELEASE}.md

tar -czvf ${RELEASE}_release.tar.gz ${RELEASE}_release
