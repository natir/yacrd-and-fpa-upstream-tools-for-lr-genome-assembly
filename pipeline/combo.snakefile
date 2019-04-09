rule miniasm_pipeline:
    input:
        "data/{file}_{tech}.fasta"

    output:
        paf="combo/{file}_{tech}_mm.paf",
        gfa="combo/{file}_{tech}_mm.gfa",
        asm="combo/{file}_{tech}_mm.fasta"

    benchmark:
        "benchmarks/combo/{file}_{tech}_mm.txt"
        
    shell:
        " && ".join([
            "minimap2 -t16 -x ava-{wildcards.tech} {input} {input} > {output.paf}",
            "miniasm -f {input} {output.paf} > {output.gfa}",
            "/home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly/script/gfaminiasm2fasta.py {output.gfa} {output.asm}"
            ])

rule yacrd_minimap_fpa_miniasm_pipeline:
    input:
        "data/{file}_{tech}.fasta"

    output:
        paf_yacrd="combo/{file}_{tech}_yacrd.paf",
        yacrd="combo/{file}_{tech}_yacrd.report",
        scrubbed_read="combo/{file}_{tech}_yacrd.fasta",
        paf="combo/{file}_{tech}_ymfm.paf",
        gfa="combo/{file}_{tech}_ymfm.gfa",
        asm="combo/{file}_{tech}_ymfm.fasta",
        
    benchmark:
        "benchmarks/combo/{file}_{tech}_ymfm.txt"
        
    shell:
        " && ".join([
            "minimap2 -t16 -x ava-{wildcards.tech} -g 1000 -n 3 {input} {input} > {output.paf_yacrd}",
            "yacrd scrubbing -m {output.paf_yacrd} -s {input} -r {output.yacrd} -S {output.scrubbed_read} -c 4 -n 0.4",
            "minimap2 -t16 -x ava-{wildcards.tech} {input} {input} | fpa -l 2000 -i  > {output.paf}",
            "miniasm -1 -2 -f {input} {output.paf} > {output.gfa}",
            "/home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly/script/gfaminiasm2fasta.py {output.gfa} {output.asm}"
        ])

ref = {"real_reads": "ref_e_coli_cft073.fasta", "d_melanogaster_reads": "d_melanogaster_ref.fasta", "c_elegans": "c_elegans_ref.fasta", "h_sapiens_chr1": "h_sapiens_chr1_ref.fasta"}
        
rule quast:
    input:
        asm="combo/{prefix}_{tech}_{suffix}.fasta",

    output:
        "combo/quast/{prefix}_{tech}_{suffix}/report.txt"

    params:
        ref=lambda wildcards, output: ref[wildcards.prefix]    
        
    shell:
        "quast -o combo/quast/{wildcards.prefix}_{wildcards.tech}_{wildcards.suffix}/ --min-identity 80.0 -r data/{params.ref} -t 16 {input.asm}"
        
def generate_template(template, files):
    return [template.format(f) for f in files]

def generate_asm_mm(files):
    return generate_template("combo/{}_mm.fasta", files)

def generate_asm_ymfm(files):
    return generate_template("combo/{}_ymfm.fasta", files)

def generate_quast_mm(files):
    return generate_template("combo/quast/{}_mm/report.txt", files)

def generate_quast_ymfm(files):
    return generate_template("combo/quast/{}_ymfm/report.txt", files)

rule all:
    input:
        generate_asm_mm(("c_elegans_pb", "d_melanogaster_reads_ont", "h_sapiens_chr1_ont")),
        generate_asm_ymfm(("c_elegans_pb", "d_melanogaster_reads_ont", "h_sapiens_chr1_ont")),
        generate_quast_mm(("c_elegans_pb", "d_melanogaster_reads_ont", "h_sapiens_chr1_ont")),
        generate_quast_ymfm(("c_elegans_pb", "d_melanogaster_reads_ont", "h_sapiens_chr1_ont")),

        
