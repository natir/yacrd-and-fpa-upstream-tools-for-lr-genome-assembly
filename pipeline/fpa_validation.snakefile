rule all:
    input:
        "fpa/quast/basic/report.txt",
        "fpa/quast/fpa/report.txt",

        
rule quast:
    input:
        ref="data/d_melanogaster_ref.fasta",
        asm="fpa/assembly/{prefix}.fasta",

    output:
        "fpa/quast/{prefix}/report.txt",

    shell:
        "quast -o fpa/quast/{wildcards.prefix}/ -r {input.ref} -t 8 {input.asm}"


rule gfa2fasta:
    input:
        "{prefix}.gfa",

    output:
        "{prefix}.fasta",

    benchmark:
        "fpa/benchmark/gfa2fasta_{prefix}.txt",
        
    shell:
        "./script/gfaminiasm2fasta.py {input} {output}"


rule miniasm:
    input:
        paf="fpa/mapping/{prefix}.paf",
        reads="data/d_melanogaster_reads_ont.fasta",

    output:
        "fpa/assembly/{prefix}.gfa",
        
    benchmark:
        "fpa/benchmark/miniasm_{prefix}.txt",
        
    shell:
        "miniasm -f {input.reads} {input.paf} > {output}"

        
rule minimap:
    input:
        "data/d_melanogaster_reads_ont.fasta",

    output:
        "fpa/mapping/basic.paf",

    benchmark:
        "fpa/benchmark/mapping_basic.txt",
        
    shell:
        "minimap2 -t 8 -x ava-ont {input} {input} > {output}"

        
rule minimap_fpa:
    input:
        "data/d_melanogaster_reads_ont.fasta",

    output:
        "fpa/mapping/fpa.paf"

    benchmark:
        "fpa/benchmark/mapping_fpa.txt"
        
    shell:
        "minimap2 -t 8 -x ava-ont {input} {input} | fpa -i -l 2000 > {output}"
