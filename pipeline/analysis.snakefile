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

rule referenceseeker:
    input:
        asm = "assembly/{prefix}.raw.wtdbg2.fasta",
        database = "referenceseeker/bacteria/"
    output:
        result = "referenceseeker/{prefix}_possible_ref.csv",
    threads:
        16
    shell:
        "referenceseeker --unfiltered -t 16 -c 10 {input.database} {input.asm} > {output.result} --crg 10"

rule dl_reference:
    input:
        result = "referenceseeker/{prefix}_possible_ref.csv",
    output:
        ref = "references/{prefix}_ref.fasta",
    shell:
        " && ".join([
            "REF=$(cut -d$'\t' -f 1 {input.result} | head -n 2 | tail -n 1)",
            "URL=$(esearch -db nuccore -query $REF | elink -db assembly -target assembly | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq)",
            "curl $(echo $URL | sed \"s|/\([^/]*\)$|/\\1/\\1_genomic.fna.gz|g\") | seqtk seq -A - > {output.ref}"
        ])
        
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

    
