ref = {"real_reads_ont": "ref_e_coli_cft073.fasta", "all_real_reads_ont": "ref_e_coli_cft073.fasta","d_melanogaster_reads_ont": "d_melanogaster_ref.fasta"}

rule all:
    input:
        "fpa/quast/d_melanogaster_reads_ont/report.txt",
        "fpa/quast/d_melanogaster_reads_ont.fpa/report.txt",
        "fpa/quast/real_reads_ont/report.txt",
        "fpa/quast/real_reads_ont.fpa/report.txt",
        "fpa/quast/all_real_reads_ont/report.txt",
        "fpa/quast/all_real_reads_ont.fpa/report.txt",
        
rule quast:
    input:
        asm="fpa/assembly/{prefix}{type}.fasta",

    output:
        "fpa/quast/{prefix,[^\.]+}{type,\.?.*}/report.txt",

    params:
        ref=lambda wildcards, output: ref[wildcards.prefix]
        
    shell:
        "quast -o fpa/quast/{wildcards.prefix}/ -r {params.ref} -t 8 {input.asm}"


rule miniasm:
    input:
        reads="data/{prefix}.fasta",
        ovl="fpa/overlap/{prefix}{type}.paf"

    output:
        graph="fpa/assembly/{prefix,[^\.]+}{type,\.?.*}.gfa",
        contig="fpa/assembly/{prefix,[^\.]+}{type,\.?.*}.fasta"

    benchmark:
        "fpa/benchmark/miniasm_{prefix}{type}.txt",
        
    shell:
        " && ".join([
            "miniasm -f {input.reads} {input.ovl} > {output.graph}",
            "./script/gfaminiasm2fasta.py {output.graph} {output.contig}"
        ])

        
rule minimap:
    input:
        "data/{prefix}.fasta",

    output:
        "fpa/overlap/{prefix}.paf",

    benchmark:
        "fpa/benchmark/overlap_{prefix}.txt",
        
    shell:
        "minimap2 -t 16 -x ava-ont {input} {input} > {output}"

        
rule minimap_fpa:
    input:
        "data/{prefix}.fasta",

    output:
        "fpa/overlap/{prefix,[^\.]+}.fpa.paf"

    benchmark:
        "fpa/benchmark/overlap_{prefix}.fpa.txt"
        
    shell:
        "minimap2 -t 16 -x ava-ont {input} {input} |  /home/pierre.marijon/tools/fpa/target/release/fpa -i -l 2000 > {output}"
