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

rule referenceseeker:
    input:
        asm = "assembly/{prefix}.raw.wtdbg2.fasta",
        database = "referenceseeker/bacteria/"
    output:
        result = "referenceseeker/{prefix}_possible_ref.csv",
    wildcard_constraints:
        prefix="[^\.]*"
    threads:
        16
    shell:
        "referenceseeker --unfiltered -t 16 -c 10 {input.database} {input.asm} > {output.result} --crg 10"

rule dl_reference:
    input:
        result = "referenceseeker/{prefix}_possible_ref.csv",
    output:
        ref = "references/{prefix}_ref.fasta",
    wildcard_constraints:
        prefix="[^\.]*"
    shell:
        " && ".join([
            "REF=$(cut -d$'\t' -f 1 {input.result} | head -n 2 | tail -n 1)",
            "URL=$(esearch -db nuccore -query $REF | elink -db assembly -target assembly | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq)",
            "curl $(echo $URL | sed \"s|/\([^/]*\)$|/\\1/\\1_genomic.fna.gz|g\") | seqtk seq -A - > {output.ref}"
        ])
        
rule quast:
    input:
        asm="assembly/{prefix}_{tech}.{scrubbing}.{asm}.fasta",
        ref="references/{prefix}_{tech}_ref.fasta"
    output:
        "quast/{prefix}_{tech}.{scrubbing}.{asm}/report.txt"
    wildcard_constraints:
        tech="[^\.]*"
    shell:
        "quast -o quast/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubbing}.{wildcards.asm}/ --min-identity 80.0 -r {input.ref} -t 16 {input.asm}"

rule indexing:
    input:
        "references/{ref}.fasta",

    output:
        "references/{ref}.fasta.bwt"
        
    shell:
        "bwa index {input}"
        
rule mapping:
    input:
        reads="scrubbing/{prefix}_{tech}.{suffix}.fasta",
        ref="references/{prefix}_{tech}_ref.fasta.bwt"
    wildcard_constraints:
        tech="[^\.]+"
    output:
        "mapping/{prefix}_{tech}.{suffix}.bam"
        
    params:
        tech=lambda wildcards, output: tech2tech_bwa(wildcards, output)
        
    shell:
        " && ".join([
            "bwa mem -t 16 -x {params.tech} data/{input.ref} {input.reads} | samtools sort > {output}",
            "samtools index {output}"
        ])

rule minimap2:
    input:
        reads="scrubbing/{prefix}_{tech}.{suffix}.fasta",
        ref="references/{prefix}_{tech}_ref.fasta"
    output:
        "mapping/{prefix}_{tech}.{suffix}.paf"
    wildcard_constraints:
        tech="[^\.]*"
    params:
        tech=lambda wildcards, output: tech2tech_paf(wildcards, output)
        
    shell:
        "minimap2 -t 16 -x {params.tech} {input.ref} {input.reads} > {output}"


rule porechop:
    input:
        "scrubbing/{prefix}_{tech}.{suffix}.fasta"
        
    output:
        "porechop/{prefix}_{tech}.{suffix}.out"
        
    shell:
        config['porechop_path'] + " -i {input} -o /dev/null --discard_middle > {output}"


rule nucmer:
    input:
        asm="assembly/{prefix}_{tech}.{scrubbing}.{asm}.fasta",
        ref="references/{prefix}_{tech}_ref.fasta"
    output:
        "nucmer/{prefix}_{tech}.{scrubbing}.{asm}.delta"
    wildcard_constraints:
        tech="[^\.]*"
    shell:
        "nucmer -t 16 --prefix nucmer/{wildcards.prefix}_{wildcards.tech}.{wildcards.scrubbing}.{wildcards.asm} --maxmatch -l 20 -c 500 {input.ref} {input.asm}"

    
