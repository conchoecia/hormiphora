# hormiphora

This repo contains the code necessary to generate the annotation for the _Hormiphora californensis_ genome, the annotation as releases, and the code and results from the _H. californensis_ genome assembly paper.

# Directory

## The assembly - Hcv1

[Download the gzipped genome assembly fasta file here.](https://github.com/conchoecia/hormiphora/blob/master/annotation/raw_files/UCSC_Hcal_v1.fa.gz)

## Annotation

To download the latest annotation, navigate to the [releases](https://github.com/conchoecia/hormiphora/releases) page and download the `Hc[version]_release.tar.gz` file. The releases are editioned like so, and dot character, `.`, -delimited :

```
Hcv1.av93
 Hcv1 = Hormiphora californensis genome assembly version 1
 av93 = annotation version 93
```

- The releases contain specific documentation, but briefly, each release contains the three most important files:
  - Model proteins (use these for protein analyses - do not translate proteins from CDS sequences generated from the GFF file/assembly file.
  - The transcripts, may contain prematurely truncated CDS sequences. See above for getting the final model proteins.
  - GFF of the transcripts.

## [`supplemental_results`](https://github.com/conchoecia/hormiphora/tree/master/supplemental_results)

- Contains directories with supplementary files for the following analyses:
  - `heterozygosity` contains files that were generated in the process of calculating the heterozygosity of the _H. californensis_ genome.
  - `intergenic_antisense` contains a Snakefile and config file used to investigate nested intronic genes. Also includes the data output for Hormiphora.
