genome_size = {"real_reads_pb": "5.2M", "real_reads_ont": "5.2M", "h_sapiens_chr1_reads_ont": "248M", "d_melanogaster_reads_ont": "143M"}

rule all:
    input:
        "scrubbing/real_reads_pb.raw.fasta",
        "scrubbing/real_reads_pb.yacrd.fasta",
        "scrubbing/real_reads_pb.yacrd2.fasta",
        "scrubbing/real_reads_pb.dascrubber.fasta",
        "scrubbing/real_reads_pb.miniscrub.fasta",
        "scrubbing/real_reads_ont.raw.fasta",
        "scrubbing/real_reads_ont.yacrd.fasta",
        "scrubbing/real_reads_ont.yacrd2.fasta",
        "scrubbing/real_reads_ont.dascrubber.fasta",
        "scrubbing/real_reads_ont.miniscrub.fasta",

rule yacrd_test:
    input:
        "scrubbing/real_reads_pb.yacrd.1.8.fasta",
        "scrubbing/real_reads_pb.yacrd.1.7.fasta",
        "scrubbing/real_reads_pb.yacrd.1.6.fasta",
        "scrubbing/real_reads_pb.yacrd.1.5.fasta",
        "scrubbing/real_reads_pb.yacrd.1.4.fasta",
        "scrubbing/real_reads_pb.yacrd.2.8.fasta",
        "scrubbing/real_reads_pb.yacrd.2.7.fasta",
        "scrubbing/real_reads_pb.yacrd.2.6.fasta",
        "scrubbing/real_reads_pb.yacrd.2.5.fasta",
        "scrubbing/real_reads_pb.yacrd.2.4.fasta",
        "scrubbing/real_reads_pb.yacrd.3.8.fasta",
        "scrubbing/real_reads_pb.yacrd.3.7.fasta",
        "scrubbing/real_reads_pb.yacrd.3.6.fasta",
        "scrubbing/real_reads_pb.yacrd.3.5.fasta",
        "scrubbing/real_reads_pb.yacrd.3.4.fasta",
        "scrubbing/real_reads_pb.yacrd.4.8.fasta",
        "scrubbing/real_reads_pb.yacrd.4.7.fasta",
        "scrubbing/real_reads_pb.yacrd.4.6.fasta",
        "scrubbing/real_reads_pb.yacrd.4.5.fasta",
        "scrubbing/real_reads_pb.yacrd.4.4.fasta",
        "scrubbing/real_reads_ont.yacrd.1.8.fasta",
        "scrubbing/real_reads_ont.yacrd.1.7.fasta",
        "scrubbing/real_reads_ont.yacrd.1.6.fasta",
        "scrubbing/real_reads_ont.yacrd.1.5.fasta",
        "scrubbing/real_reads_ont.yacrd.1.4.fasta",
        "scrubbing/real_reads_ont.yacrd.2.8.fasta",
        "scrubbing/real_reads_ont.yacrd.2.7.fasta",
        "scrubbing/real_reads_ont.yacrd.2.6.fasta",
        "scrubbing/real_reads_ont.yacrd.2.5.fasta",
        "scrubbing/real_reads_ont.yacrd.2.4.fasta",
        "scrubbing/real_reads_ont.yacrd.3.8.fasta",
        "scrubbing/real_reads_ont.yacrd.3.7.fasta",
        "scrubbing/real_reads_ont.yacrd.3.6.fasta",
        "scrubbing/real_reads_ont.yacrd.3.5.fasta",
        "scrubbing/real_reads_ont.yacrd.3.4.fasta",
        "scrubbing/real_reads_ont.yacrd.4.8.fasta",
        "scrubbing/real_reads_ont.yacrd.4.7.fasta",
        "scrubbing/real_reads_ont.yacrd.4.6.fasta",
        "scrubbing/real_reads_ont.yacrd.4.5.fasta",
        "scrubbing/real_reads_ont.yacrd.4.4.fasta",
        
rule d_melano:
    input:
        "scrubbing/d_melanogaster_reads_ont.raw.fasta",
        "scrubbing/d_melanogaster_reads_ont.yacrd.fasta",
        "scrubbing/d_melanogaster_reads_ont.yacrd2.fasta",
        "scrubbing/d_melanogaster_reads_ont.dascrubber.fasta",
        "scrubbing/d_melanogaster_reads_ont.miniscrub.fasta",
        
rule raw:
    input:
        "data/{file}.fasta"

    output:
        "scrubbing/{file}.raw.fasta"

    shell:
        "ln -s $(readlink -f {input}) {output}"
        
rule yacrd_pb:
    input:
        "data/{prefix}_pb.fasta",
        
    output:
        "scrubbing/{prefix}_pb.yacrd.fasta",

    benchmark:
        "benchmarks/{prefix}_pb.yacrd.txt",
        
    shell:
        " && ".join([
            "minimap -t 8 -x ava-pb {input} {input} | fpa -l 1000 -s -i > scrubbing/{wildcards.prefix}_pb.paf",
            "yacrd scrubbing -c 1 -m scrubbing/{wildcards.prefix}_pb.paf -r scrubbing/{wildcards.prefix}_pb.yacrd -s {input} -S {output}"
            ])

rule yacrd_ont:
    input:
        "data/{prefix}_ont.fasta",
        
    output:
        "scrubbing/{prefix}_ont.yacrd.fasta",
        
    benchmark:
        "benchmarks/{prefix}_ont.yacrd.txt",

    shell:
        " && ".join([
            "minimap -t 8 -x ava-ont {input} {input} | fpa -l 1000 -s -i > scrubbing/{wildcards.prefix}_ont.paf",
            "yacrd scrubbing -c 1 -m scrubbing/{wildcards.prefix}_ont.paf -r scrubbing/{wildcards.prefix}_ont.yacrd -s {input} -S {output}"
            ])

rule minimap_pb:
    input:
        reads="data/{prefix}_pb.fasta",

    output:
        overlaps="scrubbing/{prefix}_pb.paf"

    shell:
        "minimap -t 8 -x ava-pb {input} {input} | fpa -l 1000 -s -i > {output}"

rule minimap_ont:
    input:
        reads="data/{prefix}_ont.fasta",

    output:
        overlaps="scrubbing/{prefix}_ont.paf",

    shell:
        "minimap -t 8 -x ava-ont {input} {input} | fpa -l 1000 -s -i > {output}"
        
rule yacrd_combo:
    input:
        reads="data/{prefix}.fasta",
        overlaps="scrubbing/{prefix}.paf",
    output:
        "scrubbing/{prefix}.yacrd.{coverage}.{discard}.fasta",
        
    benchmark:
        "benchmarks/{prefix}.yacrd.{coverage}.{discard}.txt",

    shell:
        "yacrd scrubbing -c {wildcards.coverage} -n 0.{wildcards.discard} -m {input.overlaps} -r scrubbing/{wildcards.prefix}.{wildcards.coverage}.{wildcards.discard}.yacrd -s {input.reads} -S {output}"

            
rule dascrubber:
    input:
        "data/{prefix}.fasta",
        
    output:
        "scrubbing/{prefix}.dascrubber.fasta",

    params:
        genome_size=lambda wildcards, output: genome_size[wildcards.prefix]
        
    benchmark:
        "benchmarks/{prefix}.dascrubber.txt",        

    shell:
        "dascrubber_wrapper.py --dbsplit_options='-x1000' --daligner_options='-T8' --datander_options='-T8' -i {input} -g {params.genome_size} > {output}"

rule miniscrub:
    input:
        "data/{prefix}.fastq",
        
    output:
        "scrubbing/{prefix}.miniscrub.fasta",
        
    benchmark:
        "benchmarks/{prefix}.miniscrub.txt",
        
    shell:
        " && ".join([
            "run_miniscrub.sh --processes 8 --output scrubbing/{wildcards.prefix}.miniscrub.fastq {input}",
            "sed -n '1~4s/^@/>/p;2~4p' scrubbing/{wildcards.prefix}.miniscrub.fastq > {output}"
            ])

