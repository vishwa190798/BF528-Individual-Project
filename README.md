# Project Description

The objectives of this study were to determine if myocytes revert the transcriptional phenotype to a less differentiated state during regeneration and to systematically interrogate the transcriptional data to identify and validate potential regulators of this process. A global gene expression pattern is profiled over the course of mouse cardiac myocyte differentiation both in vitro and in vivo to compare this transcriptional signature of differentiation to a cardiac myocyte explant model whereby the cardiac myocytes loose fully differentiated phenotype to identify genes and gene network that changed dynamically during these processes. The RNA seq datasets are interrogated to predict and validate the upstream and downstream regulators along with its associated pathways that can modulate cell cycle state of cardiac myocytes. The aim of this project was to replicate the findings of O’Meara et al using similar tools and methods.

# Respository Contents

## Programmer Files

### run_tophat(1).qsub

Input: P0_1_1.fastq, P0_1_2.fastq, mm9 reference genome

Dependencies: python2, Bowtie2, topHat, samtools

Execution: `qsub run_tophat(1).qsub`

Output: accepted_hits.bam file

### run_RSeQC.qsub

Input: accepted_hits.bam, accepted_hits.bam.bai, mm9.bed

Dependencies: python3, rseqc, samtools, R

Execution: `qsub run_RSeQC.qsub`

Output: geneBody_coverage plot, inner_distance plot, one quality control metrice

### run_cufflinks.qsub

Input: accepted_hits.bam file and mm9.gtf and mm9.fa

Dependencies: python3, cufflinks, R

Execution: `qsub run_cufflinks.qsub`

Output: gene.FPKM_tracing

### run_cuffdiff.qsub

Input: accepted_hits.bam files for P0_1, P0_2, AD_1, AD_2

Dependencies: cufflinks

Execution: `qsub run_cuffdiff.qsub`

Output: gene_exp.diff

### FPKMhistogram.R

Input: genes.FPKM_tracing

Dependencies: tidyverse package

Execution: It is recommended to run this code on R studio

Alternatively to run on command line:

      module load R/4.0.2

      Rscript FPKMhistogram.R

Output: A histogram for FPKM values of all genes

## Biologist files

### Biologist.R

Input: genes.FPKM_tracing files for P0_1, P0_2, AD_1, AD_2, P4_1, P4_2, P7_1, P7_2

Dependencies: pheatmap, tidyverse, gplots, gridExtra packages

Execution: It is recommended to run this code on R studio

Alternatively to run on command line:

      module load R/4.0.2

      Rscript Biologist.R

Output: A line graph with FPKM values of representative Sarcomere, Mitochondria and cell cycle genes which were significant and differentially expressed and a clustered heatmap of top 1000 DEG
