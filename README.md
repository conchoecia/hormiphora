# hormiphora

![hormiphora image](supplemental_results/pictures/Misc/Horm_0720_small.JPG)

This repo contains the code necessary to generate the annotation for the _Hormiphora californensis_ genome, the annotation as releases, and the code and results from the _H. californensis_ genome assembly paper.

This repo is maintained by [Darrin T. Schultz](https://github.com/conchoecia) and [Warren R. Francis](https://github.com/wrf).

# Directory

## The assembly - Hcv1

[Download the gzipped genome assembly fasta file here.](https://github.com/conchoecia/hormiphora/blob/master/annotation/raw_files/UCSC_Hcal_v1.fa.gz)

## [Annotation](https://github.com/conchoecia/hormiphora/tree/master/annotation)

To download the latest annotation, [navigate to the releases page and download the `Hc[version]_release.tar.gz` file](https://github.com/conchoecia/hormiphora/releases). The releases are editioned like so, and dot character, `.`, -delimited :

```
Hcv1.av93
 Hcv1 = Hormiphora californensis genome assembly version 1
 av93 = annotation version 93
```

- The releases contain specific documentation, but briefly, each release contains the three most important files:
  - Model proteins (use these for protein analyses - do not translate proteins from CDS sequences generated from the GFF file/assembly file.
  - The transcripts, may contain prematurely truncated CDS sequences. See above for getting the final model proteins.
  - GFF of the transcripts.

The actual annotation directory in the repo contains the files necessary to generate the current annotation version. The annotation can be reconstructed by running snakemake in that directory.

## [`supplemental_results`](https://github.com/conchoecia/hormiphora/tree/master/supplemental_results)

- Contains directories with supplementary files for the following analyses:
  - [`heterozygosity`](https://github.com/conchoecia/hormiphora/tree/master/supplemental_results/heterozygosity) contains files that were generated in the process of calculating the heterozygosity of the _H. californensis_ genome.
  - [`intergenic_antisense`](https://github.com/conchoecia/hormiphora/tree/master/supplemental_results/intergenic_antisense) contains a Snakefile and config file used to investigate nested intronic genes. Also includes the data output for Hormiphora.
  - [`centromere_plots`](https://github.com/conchoecia/hormiphora/tree/master/supplemental_results/centrometere_plots) contains an annotation of the repeats present in the genome, as well as a python file used to plot this in repeat frequency vs coordinate to look for repeat-rich regions, as well as the plots from this analysis.
  - `pictures`
    - [`Hc2`](https://github.com/conchoecia/hormiphora/tree/master/supplemental_results/pictures/Hc2) contains pictures of the Hc2 _H. californensis_ individual. It was collected with the [MBARI ROV Doc Ricketts](https://www.mbari.org/at-sea/vehicles/remotely-operated-vehicles/rov-doc-ricketts/)

## [TADs](https://github.com/conchoecia/hormiphora/tree/master/TADs)

This directory contains the TADs for _H. californensis_.

## [Phased genome assembly](https://github.com/conchoecia/hormiphora/tree/master/phased_genome_assembly)

Contains the whole genome assembly, converted into haplotype-specific fasta files using the phased VCF files.

