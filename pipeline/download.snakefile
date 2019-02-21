rule all:
    input:
        "data/ref_e_coli_cft073.fasta",
        "data/real_reads_pb.fasta",
        "data/real_reads_pb.fastq",
        "data/real_reads_ont.fasta",
        "data/real_reads_ont.fastq",

rule dl_reference:
    output:
        "data/ref_e_coli_cft073.fasta",
        
    shell:
        "curl ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_14_collection/escherichia_coli_cft073/dna/Escherichia_coli_cft073.ASM744v1.dna.chromosome.Chromosome.fa.gz | gunzip > {output}"

rule dl_reads_pb:
    output:     
        "data/real_reads_pb.fastq",
        
    shell:
        "curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/001/SRR8494911/SRR8494911_subreads.fastq.gz | gunzip > data/real_reads_pb.fastq"

rule dl_reads_ont:
    output:     
        "data/real_reads_ont.fastq",
        
    shell:
        "curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494940/SRR8494940.fastq.gz | gunzip > data/real_reads_ont.fastq"
        
rule fastq_to_fasta:
    input:
        "{name}.fastq",
        
    output:
        "{name}.fasta",
        
    shell:
        "sed -n '1~4s/^@/>/p;2~4p' {input} > {output}"
