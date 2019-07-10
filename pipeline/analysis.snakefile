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

rule mapping:
    input:
        "scrubbing/{prefix}_{tech}.{suffix}.fasta"
        
    output:
        "mapping/{prefix}_{tech}.{suffix}.bam"
        
    params:
        ref=lambda wildcards, output: ref[wildcards.prefix],
        tech=lambda wildcards, output: tech2tech_bwa(wildcards, output)
        
    shell:
        " && ".join([
            "bwa mem -t 16 -x {params.tech} data/{params.ref} {input} | samtools sort > {output}",
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
        "/home/pierre.marijon/tools/Porechop/porechop-runner.py -i {input} -o /dev/null --discard_middle > {output}"
    
        
