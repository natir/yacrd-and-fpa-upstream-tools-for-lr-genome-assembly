
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


genome_size={"real_reads": "5.2m", "d_melanogaster_reads": "143.7m", "c_elegans": "100.2m", "h_sapiens_chr1": "248.9m"}

prefix_tech2tech={"real_reads_ont": "ont", "real_reads_pb": "sq", "c_elegans_pb": "rs", "h_sapiens_chr1_ont": "ont", "d_melanogaster_reads_ont": "ont"}

rule wdbtg2:
    input:
        "scrubbing/{prefix}_{tech,[^\.]+}.{scrubber}.fasta"

    output:
        asm="assembly/{prefix}_{tech,[^\.]+}.{scrubber}.wtdbg2.fasta",
        layout="assembly/{prefix}_{tech,[^\.]+}.{scrubber}.wtdbg2.ctg.lay.gz"

    benchmark:
        "benchmarks/{prefix}_{tech}.{scrubber}.wdbtg2.txt",

    params:
        genome_size=lambda wildcards, output: genome_size[wildcards.prefix],
        tech=lambda wildcards, output: prefix_tech2tech[wildcards.prefix + "_" + wildcards.tech]
        
    shell:
        " && ".join([
            "wtdbg2 -t 16 -g {params.genome_size} -x {params.tech} -i {input} -fo assembly/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}.wtdbg2",
            "wtpoa-cns -t 16 -i {output.layout} -fo {output.asm}"
        ])
