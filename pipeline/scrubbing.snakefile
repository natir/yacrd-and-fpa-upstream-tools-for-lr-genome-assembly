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
        "data/{file}.{ext}"

    output:
        "scrubbing/{file}.raw.{ext}"

    shell:
        "ln -s {input} {output}"
        
rule yacrd_pb:
    input:
        "data/{prefix}_pb.fasta",
        
    output:
        "scrubbing/{prefix}_pb.yacrd.fasta",

    benchmark:
        "benchmarks/{prefix}_pb.yacrd.txt",
        
    shell:
        "minimap -x ava-pb {input} {input} | fpa -l 500 -i > scrubbing/{prefix}_pb.paf | yacrd -c 1 -i -o scrubbing/{prefix}_pb.yacrd -s {input} --splited-suffix .yacrd && mv data/{prefix}_pb.fasta scrubbing/"

rule yacrd_ont:
    input:
        "data/{prefix}_pb.fasta",
        
    output:
        "scrubbing/{prefix}_ont.yacrd.fasta",
        
    benchmark:
        "benchmarks/{prefix}_ont.yacrd.txt",

    shell:
        "minimap -x ava-ont {input} {input} | fpa -l 500 -i > scrubbing/{prefix}_pb.paf | yacrd -c 1 -i -o scrubbing/{prefix}_pb.yacrd -s {input} --splited-suffix .yacrd && mv data/{prefix}_pb.fasta scrubbing/"
        
rule dascrubber:
    input:
        "data/{prefix}_pb.fasta",
        
    output:
        "scrubbing/{prefix}.dascrubber.fasta",

    benchmark:
        "benchmarks/{prefix}.dascrubber.txt",        

    shell:
        "dascrubber_wrapper.py -i {input} -g 4.6M > {output}"

rule miniscrub:
    input:
        "data/{prefix}.fasta",
        
    output:
        "scrubbing/{prefix}.miniscrub.fastq",
        
    benchmark:
        "benchmarks/{prefix}.miniscrub.txt",
        
    shell:
        "miniscrub.py --output {output} {input}"
