ref = {"real_reads": "ref_e_coli_cft073.fasta", "d_melanogaster_reads": "d_melanogaster_ref.fasta", "c_elegans": "c_elegans_ref.fasta", "h_sapiens_chr1": "h_sapiens_chr1_ref.fasta"}

def tech2tech_bwa(wildcards, output):
    if "ont" in wildcards.tech:
        return "ont2d"
    else:
        return "pacbio"

def tech2tech_paf(wildcards, output):
    if "ont" in wildcards.tech:
        return "map-ont"
    else:
        return "map-pb"
    
rule quast:
    input:
        asm="assembly/{prefix}_{tech}.{scrubbing}.{asm}.fasta",

    output:
        "quast/{prefix}_{tech}.{scrubbing}.{asm}/report.txt"

    params:
        ref=lambda wildcards, output: ref[wildcards.prefix]    
        
    shell:
        "quast -o quast/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubbing}.{wildcards.asm}/ --min-identity 80.0 -r data/{params.ref} -t 16 {input.asm}"

rule quast_lr:
    input:
        asm="assembly/{prefix}_{tech}.{scrubbing}.{asm}.fasta",

    output:
        "quast_lr/{prefix}_{tech}.{scrubbing}.{asm}/report.txt"

    params:
        ref=lambda wildcards, output: ref[wildcards.prefix]    
        
    shell:
        "quast -o quast_lr/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubbing}.{wildcards.asm}/ --min-identity 80.0 -r data/{params.ref} -t 16 {input.asm} --extensive-mis-size 10000"


rule indexing:
    input:
        "data/{ref}.fasta",

    output:
        "data/{ref}.fasta.bwt"
        
    shell:
        "bwa index {input}"
        
rule mapping:
    input:
        reads="scrubbing/{prefix}_{tech}.{suffix}.fasta",
        ref=lambda wildcards: f"data/{ref[wildcards.prefix]}.bwt"

    output:
        "mapping/{prefix}_{tech}.{suffix}.bam"
        
    params:
        ref=lambda wildcards, output: ref[wildcards.prefix],
        tech=lambda wildcards, output: tech2tech_bwa(wildcards, output)
        
    shell:
        " && ".join([
            "bwa mem -t 16 -x {params.tech} data/{input.ref} {input.reads} | samtools sort > {output}",
            "samtools index {output}"
        ])

rule minimap2:
    input:
        "scrubbing/{prefix}_{tech}.{suffix}.fasta"
        
    output:
        "mapping/{prefix}_{tech}.{suffix}.paf"
        
    params:
        ref=lambda wildcards, output: ref[wildcards.prefix],
        tech=lambda wildcards, output: tech2tech_paf(wildcards, output)
        
    shell:
        "minimap2 -t 16 -x {params.tech} data/{params.ref} {input} > {output}"


rule porechop:
    input:
        "scrubbing/{prefix}_{tech}.{suffix}.fasta"
        
    output:
        "porechop/{prefix}_{tech}.{suffix}.out"
        
    shell:
        config['porechop_path'] + " -i {input} -o /dev/null > {output}"


rule nucmer:
    input:
        asm="assembly/{prefix}_{tech}.{scrubbing}.{asm}.fasta",

    output:
        "nucmer/{prefix}_{tech}.{scrubbing}.{asm}.delta"

    params:
        ref=lambda wildcards, output: ref[wildcards.prefix]    
        
    shell:
        "nucmer -t 16 --prefix nucmer/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubbing}.{wildcards.asm} --maxmatch -l 20 -c 500 data/{params.ref} {input.asm}"

    
