ref = {"real_reads": "ref_e_coli_cft073.fasta", "d_melanogaster_reads": "d_melanogaster_ref.fasta"}
rule quast:
    input:
        asm="assembly/{prefix}_{tech}.{scrubbing}.{asm}.fasta",

    output:
        "quast/{prefix}_{tech}.{scrubbing}.{asm}/report.txt"

    params:
        ref=lambda wildcards, output: ref[wildcards.prefix]    
        
    shell:
        "quast -o quast/{wildcards.prefix}/ -r {params.ref} -t 16 {input.asm}"

rule mapping:
    input:
        "scrubbing/{prefix}_{tech}.{suffix}.fasta"
        
    output:
        "mapping/{prefix}_{tech}.{suffix}.bam"
        
    params:
        ref=lambda wildcards, output: ref[wildcards.prefix],
        tech=lambda wildcards, output: tech2tech(wildcards, output)
        
    shell:
        " && ".join([
            "bwa mem -t 16 -x {params.tech} data/{params.ref} {input} | samtools sort > {output}",
            "samtools index {output}"
        ])

