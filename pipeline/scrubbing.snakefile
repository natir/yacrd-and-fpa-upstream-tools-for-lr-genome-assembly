read_type = {"ont": "ava-ont", "pb": "ava-pb"}

rule raw:
    input:
        "data/{file}_{techno}.fasta"

    output:
        "scrubbing/{file}_{techno}.raw.fasta"

    wildcard_constraints:
        techno="[^\.]+"
        
    shell:
        "ln -s $(readlink -f {input}) {output}"


rule yacrd:
    input:
        reads="data/{prefix}_{techno}.fasta"
        
    output:
        "scrubbing/{prefix}_{techno}.{coverage,\d+}.{discard,\d+}.yacrd.fasta"

    benchmark:
        "benchmarks/{prefix}_{techno}.{coverage}.{discard}.yacrd.txt"
        
    shell:
        " && ".join([
            "minimap2 -t 16 -x ava-{wildcards.techno} {input.reads} {input.reads} | fpa drop -i -l 2000 > scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.paf",
            "yacrd scrubbing -m scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.paf -s {input.reads} -r scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.yacrd -S {output} -c {wildcards.coverage} -n 0.{wildcards.discard}",
        ])

        
rule yacrd_precision:
    input:
        reads="data/{prefix}_{techno}.fasta"
        
    output:
        "scrubbing/{prefix}_{techno}.{coverage}.{discard}.precision.yacrd.fasta"

    benchmark:
        "benchmarks/{prefix}_{techno}.{coverage}.{discard}.precision.yacrd.txt"
        
    shell:
        " && ".join([
            "minimap2 -t16 -x ava-{wildcards.techno} -g 1000 -n 3 {input.reads} {input.reads} > scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.precision.paf",
            "yacrd scrubbing -m scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.precision.paf -s {input.reads} -r scrubbing/{wildcards.prefix}_{wildcards.techno}.{wildcards.coverage}.{wildcards.discard}.yacrd -S {output} -c {wildcards.coverage} -n 0.{wildcards.discard}",
        ])

rule minimap_yacrd_gc:
    input:
        reads="data/{prefix}_{techno}.fasta"
        
    output:
        "scrubbing/{prefix}_{techno}.g{dist}.paf"

    benchmark:
        "benchmarks/{prefix}_{techno}.g{dist}.minimap.yacrd.txt"
        
    shell:
        "minimap2 -t16 -x ava-{wildcards.techno} -g {wildcards.dist} {input.reads} {input.reads} > {output}"
        
rule yacrd_gc:
    input:
        reads="data/{prefix}_{techno}.fasta",
        overlap="scrubbing/{prefix}_{techno}.g{dist}.paf"
    output:
        "scrubbing/{prefix}_{techno}.g{dist}.c{coverage}.yacrd.fasta"

    benchmark:
        "benchmarks/{prefix}_{techno}.g{dist}.c{coverage}.yacrd.txt"
        
    shell:
        "yacrd scrubbing -m {input.overlap} -s {input.reads} -r scrubbing/{wildcards.prefix}_{wildcards.techno}.g{wildcards.dist}.c{wildcards.coverage}.yacrd -S {output} -c {wildcards.coverage} -n 0.4"
        
        
rule dascrubber:
    input:
        "data/{prefix}.fasta",
        
    output:
        "scrubbing/{prefix}.dascrubber.fasta",

    params:
        coverage=lambda wildcards, output: config["coverage"][wildcards.prefix]
        
    benchmark:
        "benchmarks/{prefix}.dascrubber.txt",        

    shell:
        " && ".join([
            "rm -rf dascrubber/{wildcards.prefix}/",
            "mkdir -p dascrubber/{wildcards.prefix}/",
            "cd dascrubber/{wildcards.prefix}/",

            "./../../script/rename_with_fake_pacbio.py ../../{input} renamed_reads.fasta",

            "echo '---------- fasta2DB ----------'",
            "fasta2DB reads.db renamed_reads.fasta",

            "echo '---------- DBsplit ----------'",
            "DBsplit -s200 -x100 reads",

            "NB_BLOCK=$(cat reads.db | grep 'blocks' | cut -d' ' -f 3- | sed 's/\s//g')",

            "echo '---------- HPC.daligner ----------'",
            "mkdir align_temp",
            "HPC.daligner -v -M16 -Palign_temp -T16 reads | csh",
            "rm -r align_temp",

            "echo '---------- Maybe LAmerge ----------'",
            "if [ ${{NB_BLOCK}} != 1 ]; then LAmerge reads.las reads.*.las; fi",

            "echo '---------- HPC.REPmask ----------'",
            "if [ ${{NB_BLOCK}} = 1 ]; then BLOCK=1; else BLOCK=2; fi",
            "mkdir align_temp",
            "HPC.REPmask -v -Palign_temp -g${{BLOCK}} -T16 -c{params.coverage} reads 1-${{NB_BLOCK}} | csh",
            "rm -r align_temp",

            "echo '---------- HPC.TANmask ----------'",
            "mkdir align_temp",
            "HPC.TANmask -v reads -T16 -Palign_temp | csh",
            "rm -r align_temp",

            "echo '---------- HPC.daligner ----------'",
            "mkdir align_temp",
            "rm reads.*las",
            "HPC.daligner -v -Palign_temp -T16 -mrep -mtan -T16 reads | csh",
            "rm -r align_temp",

            "echo '---------- Maybe LAmerge ----------'",
            "if [ ${{NB_BLOCK}} != 1 ]; then LAmerge reads.las reads.*.las; fi",

            "echo '---------- DAScover ----------'",
            "DAScover -v reads reads.las",

            "echo '---------- DASqv ----------'",
            "DASqv -v -c{params.coverage} reads reads.las",

            "echo '---------- DAStrim ----------'",
            "DAStrim -v reads reads.las",

            "echo '---------- DASpatch ----------'",
            "DASpatch -v reads reads.las",

            "echo '---------- DASedit ----------'",
            "DASedit -v reads patched_reads",

            "echo '---------- DB2fasta ----------'",
            "mv renamed_reads.fasta temp.fasta",
            "DB2fasta -vU patched_reads",
            "seqtk seq renamed_reads.fasta > ../../{output}"
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
            "mkdir -p miniscrub/{wildcards.prefix}/",
            "cd miniscrub/{wildcards.prefix}/"
            "python3 /home/pierre.marijon/tools/jgi-miniscrub/miniscrub.py --processes 16 --output ../../scrubbing/{wildcards.prefix}.miniscrub.fastq ../../{input}",
            "sed -n '1~4s/^@/>/p;2~4p' ../../scrubbing/{wildcards.prefix}.miniscrub.fastq > ../../{output}"
            ])

rule miniscrub_cpu:
    input:
        "data/{prefix}_{techno}.fastq",
        
    output:
        "scrubbing/{prefix}_{techno}.miniscrub.cpu.fasta",
        
    benchmark:
        "benchmarks/{prefix}_{techno}.miniscrub.cpu.txt",
        
    shell:
        " && ".join([
            "mkdir -p miniscrub/{wildcards.prefix}_{wildcards.techno}_cpu/",
            "cd miniscrub/{wildcards.prefix}_{wildcards.techno}_cpu/",
            "/home/pierre.marijon/tools/jgi-miniscrub/venv_cpu/bin/python3 /home/pierre.marijon/tools/jgi-miniscrub/miniscrub.py --processes 16 --output ../../scrubbing/{wildcards.prefix}_{wildcards.techno}.miniscrub.cpu.fastq ../../{input}",
            "sed -n '1~4s/^@/>/p;2~4p' ../../scrubbing/{wildcards.prefix}_{wildcards.techno}.miniscrub.cpu.fastq > ../../{output}"
            ])

