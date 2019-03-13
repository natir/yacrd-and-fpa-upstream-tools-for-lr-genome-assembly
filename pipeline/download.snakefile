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
        " && ".join([
            "curl ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_14_collection/escherichia_coli_cft073/dna/Escherichia_coli_cft073.ASM744v1.dna.chromosome.Chromosome.fa.gz > {output}.gz",
            "gzip -dk {output}.gz > {output}"
            ])

rule dl_reference:
    output:
        "data/ref_d_melanogaster.fasta",

    shell:
        " && ".join([
            "curl ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz > {output}.gz",
            "gzip -dk {output}.gz > {output}"
        ])
        
rule dl_reads_pb:
    output:     
        "data/all_real_reads_pb.fastq.gz",
        
    shell:
        "curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/001/SRR8494911/SRR8494911_subreads.fastq.gz > {output}"
   
rule dl_reads_ont:
    output:     
        "data/all_real_reads_ont.fastq.gz",
        
    shell:
        "curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494940/SRR8494940.fastq.gz > {output}"

rule dl_reads_ont:
    output:     
        "data/d_melanogaster_reads_ont.fastq.gz",
        
    shell:
        "curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR670/003/SRR6702603/SRR6702603.fastq.gz > {output}"
        
rule subsampling_pb:
    input:
        "data/all_real_reads_pb.fastq.gz"
    output:
        "data/real_reads_pb.fastq"
    shell:
        "seqtk sample -s 42 {input} 0.18 > {output}"

rule subsampling_ont:
    input:
        "data/all_real_reads_ont.fastq.gz"
    output:
        "data/real_reads_ont.fastq"
    shell:
        "seqtk sample -s 42 {input} 0.16 > {output}"
        
        
rule gz_to_fastq:
    input:
        "{name}.fastq.gz",
        
    output:
        "{name}.fastq",
        
    shell:
        "gzip -d -c {input} > {output}"

        
rule fastq_to_fasta:
    input:
        "{name}.fastq",
        
    output:
        "{name}.fasta",
        
    shell:
        "sed -n '1~4s/^@/>/p;2~4p' {input} > {output}"
