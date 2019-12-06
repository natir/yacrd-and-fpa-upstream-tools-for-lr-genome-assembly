
rule miniasm:
    input:
        "scrubbing/{prefix}_{tech}.{scrubber}.fasta"

    output:
        ovl="assembly/{prefix}_{tech}.{scrubber}.miniasm.paf",
        asm="assembly/{prefix}_{tech}.{scrubber}.miniasm.fasta",
        graph="assembly/{prefix}_{tech}.{scrubber}.miniasm.gfa",

    benchmark:
        "benchmarks/{prefix}_{tech}.{scrubber}.miniasm.txt",

    params:
        tech=lambda wildcards, output: tech2tech(wildcards, output)
        
    shell:
        " && ".join([
            "minimap2 -t 16 -x ava-{params.tech} {input} {input} > {output.ovl}",
            "miniasm -f {input} {output.ovl} > {output.graph}",
            "./script/gfaminiasm2fasta.py {output.graph} {output.asm}"
        ])
        
rule wdbtg2:
    input:
        "scrubbing/{prefix}_{tech,[^\.]+}.{scrubber}.fasta"

    output:
        asm="assembly/{prefix}_{tech,[^\.]+}.{scrubber}.wtdbg2.fasta",
        layout="assembly/{prefix}_{tech,[^\.]+}.{scrubber}.wtdbg2.ctg.lay.gz"

    benchmark:
        "benchmarks/{prefix}_{tech}.{scrubber}.wdbtg2.txt",

    params:
        genome_size=lambda wildcards, output: config["genome_size"][wildcards.prefix],
        tech=lambda wildcards, output: config["prefix_tech2tech_wtdbg2"][wildcards.prefix + "_" + wildcards.tech]
        
    shell:
        " && ".join([
            "wtdbg2 -t 16 -g {params.genome_size} -x {params.tech} -i {input} -fo assembly/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}.wtdbg2",
            "wtpoa-cns -t 16 -i {output.layout} -fo {output.asm}"
        ])

rule ra:
    input:
        "scrubbing/{prefix}_{tech,[^\.]+}.{scrubber}.fasta"

    output:
        "assembly/{prefix}_{tech,[^\.]+}.{scrubber}.ra.fasta"

    benchmark:
        "benchmarks/{prefix}_{tech}.{scrubber}.ra.txt"

    params:
        tech=lambda wildcards, output: config["prefix_tech2tech_ra"][wildcards.prefix + "_" + wildcards.tech]

    shell:
        " && ".join([
            "mkdir -p ra/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}/",
            "cd ra/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}/",
            "/home/pierre.marijon/tools/ra/build/bin/ra -t 16 -x {params.tech} ../../{input} > ../../{output}"
        ])
        
rule shasta:
    input:
        "scrubbing/{prefix}_{tech,[^\.]+}.{scrubber}.fasta"

    output:
        "assembly/{prefix}_{tech,[^\.]+}.{scrubber}.shasta.fasta"

    benchmark:
        "benchmarks/{prefix}_{tech}.{scrubber}.shasta.txt"

    shell:
        " && ".join([
            "rm -rf assembly/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}.shasta/",
            "rm -rf assembly/real_reads_pb.miniscrub.cpu.shasta.fasta",
            "/home/pierre.marijon/tools/shasta/shasta --memoryMode anonymous --memoryBacking 4K --input {input} --output assembly/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}.shasta/",
            "ln -s $(readlink -f assembly/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}.shasta/Assembly.fasta) {output}"
        ])
        
