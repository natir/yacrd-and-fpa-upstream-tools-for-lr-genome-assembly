rule all:
    input:
        "scrubbing/real_reads_pb.raw.fasta",
        "scrubbing/real_reads_pb.yacrd.fasta",
        "scrubbing/real_reads_pb.dascrubber.fasta",
        "scrubbing/real_reads_pb.miniscrub.fasta",
        "scrubbing/real_reads_ont.raw.fasta",
        "scrubbing/real_reads_ont.yacrd.fasta",
        "scrubbing/real_reads_ont.dascrubber.fasta",
        "scrubbing/real_reads_ont.miniscrub.fasta",

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
        "minimap -t 8 -x ava-pb {input} {input} | fpa -l 1000 -s -i > scrubbing/{wildcards.prefix}_pb.paf | yacrd -c 1 -o scrubbing/{wildcards.prefix}_pb.yacrd -s {input} --splited-suffix .yacrd && mv data/{wildcards.prefix}_pb.yacrd.fasta scrubbing/"

rule yacrd_ont:
    input:
        "data/{prefix}_pb.fasta",
        
    output:
        "scrubbing/{prefix}_ont.yacrd.fasta",
        
    benchmark:
        "benchmarks/{prefix}_ont.yacrd.txt",

    shell:
        "minimap -t 8 -x ava-ont {input} {input} | fpa -l 1000 -s -i > scrubbing/{wildcards.prefix}_ont.paf | yacrd -c 1 -o scrubbing/{wildcards.prefix}_ont.yacrd -s {input} --splited-suffix .yacrd && mv data/{wildcards.prefix}_ont.fasta scrubbing/"
        
rule dascrubber:
    input:
        "data/{prefix}.fasta",
        
    output:
        "scrubbing/{prefix}.dascrubber.fasta",

    benchmark:
        "benchmarks/{prefix}.dascrubber.txt",        

    shell:
        "dascrubber_wrapper.py --daligner_options='-T 8' --datander_options='-T 8' -i {input} -g 4.6M > {output}"

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

