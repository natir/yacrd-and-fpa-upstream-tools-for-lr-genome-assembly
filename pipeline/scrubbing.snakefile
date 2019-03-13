genome_size = {"real_reads_pb": "5.2M", "real_reads_ont": "5.2M", "h_sapiens_chr1_reads_ont": "248M", "d_melanogaster_reads_ont": "143M"}

read_type = {"ont": "ava-ont", "pb": "ava-pb"}

rule raw:
    input:
        "data/{file}.fasta"

    output:
        "scrubbing/{file}.raw.fasta"

    shell:
        "ln -s $(readlink -f {input}) {output}"


rule yacrd:
    input:
        reads="data/{prefix}.fasta"
        
    output:
        "scrubbing/{prefix}_{techno}.{coverage}.{discard}.yacrd.fasta"

    benchmark:
        "benchmarks/scrubbing/{prefix}_{techno}.{coverage}.{discard}.yacrd.txt"
        
    shell:
        " && ".join([
            "minimap -x ava-{wildcards.techno} {input.reads} {input.reads} | fpa -i -l 2000 > scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.paf",
            "yacrd -m scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.paf -s {input.reads} -r scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.yacrd -S {output}",
        ])

        
rule yacrd_precision:
    input:
        reads="data/{prefix}.fasta"
        
    output:
        "scrubbing/{prefix}_{techno}.{coverage}.{discard}.precision.yacrd.fasta"

    benchmark:
        "benchmarks/scrubbing/{prefix}_{techno}.{coverage}.{discard}.precision.yacrd.txt"
        
    shell:
        " && ".join([
            "minimap -x ava-{wildcards.techno} -g 1000 -n 3 {input.reads} {input.reads} > scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.precision.paf",
            "yacrd -m scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.precision.paf -s {input.reads} -r scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.yacrd -S {output}",
        ])
    
            
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

