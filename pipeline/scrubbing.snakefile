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
        reads="data/{prefix}_{techno}.fasta"
        
    output:
        "scrubbing/{prefix}_{techno}.{coverage}.{discard}.yacrd.fasta"

    benchmark:
        "benchmarks/scrubbing/{prefix}_{techno}.{coverage}.{discard}.yacrd.txt"
        
    shell:
        " && ".join([
            "minimap -t 16 -x ava-{wildcards.techno} {input.reads} {input.reads} | fpa -i -l 2000 > scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.paf",
            "yacrd -m scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.paf -s {input.reads} -r scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.yacrd -S {output}",
        ])

        
rule yacrd_precision:
    input:
        reads="data/{prefix}_{techno}.fasta"
        
    output:
        "scrubbing/{prefix}_{techno}.{coverage}.{discard}.precision.yacrd.fasta"

    benchmark:
        "benchmarks/scrubbing/{prefix}_{techno}.{coverage}.{discard}.precision.yacrd.txt"
        
    shell:
        " && ".join([
            "minimap -t16 -x ava-{wildcards.techno} -g 1000 -n 3 {input.reads} {input.reads} > scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.precision.paf",
            "yacrd -m scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.precision.paf -s {input.reads} -r scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.yacrd -S {output}",
        ])

        
coverage = {"real_reads_pb": "49", "real_reads_ont": "49", "d_melanogaster_reads_ont": "63"} 
rule dascrubber:
    input:
        "data/{prefix}.fasta",
        
    output:
        "scrubbing/{prefix}.dascrubber.fasta",

    params:
        coverage=lambda wildcards, output: coverage[wildcards.prefix]
        
    benchmark:
        "benchmarks/{prefix}.dascrubber.txt",        

    shell:
        " && ".join([
            "mkdir -p dascrubber/{wildcards.prefix}/"
            "cd dascrubber/{wildcards.prefix}/",

            "./script/rename_with_fake_pacbio.py {input} renamed_reads.fasta"
            "fasta2DB reads.db renamed_reads.fasta",
            "DBsplit -s200 -x100 reads",

            "mkdir align_temp",
            "HPC.daligner -v -M16 -Palign_temp -T16 reads | csh",
            "rm -r align_temp",

            "HPC.REPmask -v -c{params.coverage} reads reads.reads.las | csh",

            "mkdir align_temp",
            "datander -v -Palign_temp -T16 reads",
            "rm -r align_temp",

            "TANmask -v reads TAN.reads",

            "mkdir align_temp",
            "HPC.daligner -v -Palign_temp -mrep -mtan -T16 reads | csh",
            "rm -r align_temp",

            "HPC.DAScover -v reads reads.reads.las | csh",

            "DASqv -v -c{params.coverage} reads reads.reads.las",

            "DAStrim -v reads reads.reads.las",

            "DASpatch -v reads reads.reads.las",

            "DASedit '-v' reads patched_reads",
        ])
        
rule miniscrub:
    input:
        "data/{prefix}.fastq",
        
    output:
        "scrubbing/{prefix}.miniscrub.fasta",
        
    benchmark:
        "benchmarks/{prefix}.miniscrub.txt",
        
    shell:
        " && ".join([
            "module load tensorflow/1.12.0/anaconda3",
            "python3 /home/pierre.marijon/tools/jgi-miniscrub/miniscrub.py --processes 16 --output scrubbing/{wildcards.prefix}.miniscrub.fastq {input}",
            "sed -n '1~4s/^@/>/p;2~4p' scrubbing/{wildcards.prefix}.miniscrub.fastq > {output}"
            ])

