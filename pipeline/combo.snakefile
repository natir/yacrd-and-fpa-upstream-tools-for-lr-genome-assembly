rule miniasm_pipeline:
    input:
        "data/{file}_{tech}.fasta"

    output:
        paf="combo/basic/{file}_{tech}_mm.paf",
        gfa="combo/basic/{file}_{tech}_mm.gfa",
        asm="combo/basic/{file}_{tech}_mm.fasta"

    benchmark:
        "benchmarks/combo/{file}_{tech}_mm.txt"
        
    shell:
        " && ".join([
            "minimap2 -t16 -x ava-{wildcards.tech} {input} > {output.paf}",
            "miniasm -f {input} {output.paf} > {output.gfa}",
            "/home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly/script/gfaminiasm2fasta.py {output.gfa} {output.asm}"
            ])

rule yacrd_minimap_fpa_miniasm_pipeline:
    input:
        "data/{file}_{tech}.fasta"

    output:
        paf_yacrd="combo/basic/{file}_{tech}_yacrd.paf",
        yacrd="combo/basic/{file}_{tech}_yacrd.report",
        scrubbed_read="combo/basic/{file}_{tech}_yacrd.fasta",
        paf="combo/basic/{file}_{tech}_ymfm.paf",
        gfa="combo/basic/{file}_{tech}_ymfm.gfa",
        asm="combo/basic/{file}_{tech}_ymfm.fasta",
        
    benchmark:
        "benchmarks/combo/{file}_{tech}_ymfm.txt"
        
    shell:
        " && ".join([
            "minimap2 -t16 -x ava-{wildcards.tech} -g 1000 -n 3 {input} {input} > {output.paf_yacrd}",
            "yacrd scrubbing -m {output.paf_yacrd} -s {input} -r {output.yacrd} -S {output.scrubbed_read} -c 4 -n 0.4",
            "minimap2 -t16 -x ava-{wildcards.tech} {input} {input} | fpa -l 2000 -i  > {output.paf_yacrd}",
            "miniasm -1 -2 -f {input} {output.paf} > {output.gfa}",
            "/home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly/script/gfaminiasm2fasta.py {output.gfa} {output.asm}"
        ])

def generate_combo(template, files):
    return [template.format(f) for f in files]

def generate_combo_mm(files):
    return generate_combo("combo/basic/{}_mm.fasta", files)

def generate_combo_ymfm(files):
    return generate_combo("combo/basic/{}_ymfm.fasta", files)

rule all:
    input:
        generate_combo_mm(("c_elegans_pb", "d_melanogaster_reads_ont", "h_sapiens_chr1_ont")),
        generate_combo_ymfm(("c_elegans_pb", "d_melanogaster_reads_ont", "h_sapiens_chr1_ont"))
