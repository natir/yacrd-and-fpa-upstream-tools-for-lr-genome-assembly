rule miniasm:
    input:
        "scrubbing/{prefix}_{tech}.{scrubber}.fasta"

    output:
        ovl="assembly/{prefix}_{tech}.{scrubber}.miniasm.paf",
        asm="assembly/{prefix}_{tech}.{scrubber}.miniasm.fasta",
        graph="assembly/{prefix}_{tech}.{scrubber}.miniasm.gfa",

    benchmark:
        "benchmarks/{prefix}_{tech}.{scrubber}.miniasm.txt",
        
    shell:
        " && ".join([
            "minimap2 -t 16 -x ava-{wildcards.tech} {input} {input} > {output.ovl}",
            "miniasm -f {input} {output.ovl} > {output.graph}",
            "/home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly/script/gfaminiasm2fasta.py {output.graph} {output.asm}"
        ])


genome_size={"real_reads": "5.2m", "d_melanogaster_reads": "143.7m"}
rule wdbtg2:
    input:
        "scrubbing/{prefix}_{tech}.{scrubber}.fasta"

    output:
        asm="assembly/{prefix}_{tech}.{scrubber}.wtdbg2.fasta",
        layout="assembly/{prefix}_{tech}.{scrubber}.wtdbg2.ctg.lay.gz"

    params:
        genome_size=lambda wildcards, output: genome_size[wildcards.prefix] 
        
    shell:
        " && ".join([
            "/home/pierre.marijon/tools/wtdbg2/wtdbg2 -t 16 -g {params.genome_size} -x ont -i {input} -fo assembly/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubber}.wtdbg2",
            "/home/pierre.marijon/tools/wtdbg2/wtpoa-cns -16 -x ont -i {output.layout} -fo {output.asm}"
        ])
