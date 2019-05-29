ref = {"real_reads_ont": "ref_e_coli_cft073.fasta", "real_reads_pb": "ref_e_coli_cft073.fasta","d_melanogaster_reads_ont": "d_melanogaster_ref.fasta", "c_elegans_pb": "c_elegans_ref.fasta", "h_sapiens_chr1_ont": "h_sapiens_chr1_ref.fasta"}

def tech2tech(wildcards, output):
    if "ont" in wildcards.tech:
        return "ont"
    else:
        return "pb"

rule all:
    input:
        "fpa/quast/c_elegans_pb/report.txt",
        "fpa/quast/fpa_c_elegans_pb/report.txt",
        "fpa/quast/h_sapiens_chr1_ont/report.txt",
        "fpa/quast/fpa_h_sapiens_chr1_ont/report.txt",
        "fpa/quast/d_melanogaster_reads_ont/report.txt",
        "fpa/quast/fpa_d_melanogaster_reads_ont/report.txt",
        "fpa/quast/real_reads_pb/report.txt",
        "fpa/quast/fpa_real_reads_pb/report.txt",
        "fpa/quast/real_reads_ont/report.txt",
        "fpa/quast/fpa_real_reads_ont/report.txt",
        
        
rule quast:
    input:
        asm="fpa/assembly/{fpa}{prefix}.fasta",

    output:
        "fpa/quast/{fpa,(fpa_)?}{prefix}/report.txt",

    params:
        ref=lambda wildcards, output: ref[wildcards.prefix]
        
    shell:
        "quast -o fpa/quast/{wildcards.fpa}{wildcards.prefix}/ -r data/{params.ref} -t 16 {input.asm}  --min-identity 80.0"


rule miniasm:
    input:
        reads="data/{prefix}_{tech}.fasta",
        ovl="fpa/overlap/{fpa}{prefix}_{tech}.paf",

    output:
        asm="fpa/assembly/{fpa,(fpa_)?}{prefix}_{tech}.fasta",
        graph="fpa/assembly/{fpa,(fpa_)?}{prefix}_{tech}.gfa",

    benchmark:
        "fpa/benchmarks/{fpa}{prefix}_{tech}_miniasm.txt",

    shell:
        " && ".join([
            "miniasm -f {input.reads} {input.ovl} > {output.graph}",
            "/home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly/script/gfaminiasm2fasta.py {output.graph} {output.asm}"
        ])

        
rule minimap:
    input:
        "data/{prefix}.fasta",

    output:
        "fpa/overlap/{prefix,[^f]+}.paf",

    benchmark:
        "fpa/benchmarks/overlap_{prefix}.txt",
        
    shell:
        "minimap2 -t 16 -x ava-ont {input} {input} > {output}"

        
rule minimap_fpa:
    input:
        "data/{prefix}.fasta",

    output:
        "fpa/overlap/fpa_{prefix}.paf"

    benchmark:
        "fpa/benchmarks/overlap_fpa_{prefix}.txt"
        
    shell:
        "minimap2 -t 16 -x ava-ont {input} {input} | fpa drop -i -l 2000 > {output}"
