#!/bin/bash

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494940/SRR8494940_1.fastq.gz | gunzip > data/all_real_reads_ont.fastq

curl "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?report=fasta&id=26111730" > data/ref_e_coli_cft073.fasta
bwa index data/ref_e_coli_cft073.fasta

seqtk seq -a data/all_real_reads_ont.fastq > data/all_real_reads_ont.fasta
seqtk sample -s 42 data/all_real_reads_ont.fasta 0.1618 > data/real_reads_ont.fasta

snakemake --snakefile pipeline/uncorrected.snakefile quast/real_reads_ont.g500.c4.yacrd.miniasm/report.txt quast/real_reads_ont.g500.c4.yacrd.wtdbg2/report.txt mapping/real_reads_ont.g500.c4.yacrd.bam mapping/real_reads_ont.g500.c4.yacrd.paf

snakemake --snakefile pipeline/fpa.snakefile fpa/quast/real_reads_ont/report.txt fpa/quast/fpa_real_reads_ont/report.txt

snakemake --snakefile pipeline/combo.snakefile combo/real_reads_ont_mm.fasta combo/real_reads_ont_pymfm.fasta combo/quast/real_reads_ont_mm/report.txt combo/quast/real_reads_ont_pymfm/report.txt
