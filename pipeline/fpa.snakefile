
def tech2tech(wildcards, output):
    if "ont" in wildcards.tech:
        return "ont"
    else:
        return "pb"

dataset = [
    "c_elegans_pb",
    "h_sapiens_chr1_ont",
    "d_melanogaster_reads_ont",
    "real_reads_ont",
    "real_reads_pb",
    "ERR2531093_pb",
    "ERR2531094_pb",
    "ERR2531095_pb",
    "ERR2531096_pb",
    "ERR2564864_pb",
    "ERR2564865_pb",
    "ERR2603033_pb",
    "ERR2615933_pb",
    "ERR2615934_pb",
    "ERR2651535_pb",
    "ERR2651536_pb",
    "ERR2632359_pb",
    "ERR2681660_pb",
    "ERR2672424_pb",
    "ERR2672425_pb", 
    "ERR3253076_pb",
    "ERR3253077_pb",
    "ERR3253078_pb",
    "ERR3253079_pb",
    "ERR3253080_pb",
    "ERR3253081_pb",
    "ERR3253082_pb",
    "ERR3253083_pb",
    "ERR3253084_pb",
    "ERR3253085_pb",
    "ERR3253086_pb",
    "ERR3253087_pb",
    "ERR3253088_pb",
    "ERR3253089_pb",
    "ERR3253090_pb",
    "ERR3253091_pb",
    "ERR3253093_pb",
    "ERR3253094_pb",
    "ERR3253095_pb",
    "ERR3253096_pb",
    "ERR3253097_pb",
    "ERR3253098_pb",
    "ERR3253099_pb",
    "ERR3253100_pb",
    "ERR3253102_pb",
    "ERR3253103_pb",
    "ERR3253104_pb",            
    "ERR3500074_pb",
    "ERR2695058_pb",
    "ERR2695057_pb",
    "SRR8494915_ont",
    "SRR8494916_ont",
    "SRR8494917_ont",
    "SRR8494918_ont",
    "SRR8494919_ont",
    "SRR8494920_ont",
    "SRR8494921_ont",
    "SRR8494922_ont",
    "SRR8494923_ont",
    "SRR8494924_ont",
    "SRR8494935_ont",
    "SRR8494936_ont",
    "SRR8494937_ont",
    "SRR8494938_ont",
    "SRR8494939_ont",
    "SRR8494941_ont",
    "SRR8494942_ont",
    "SRR8494943_ont",
    "SRR8494944_ont",
    "SRR8494905_pb",
    "SRR8494906_pb",
    "SRR8494907_pb",
    "SRR8494908_pb",
    "SRR8494909_pb",
    "SRR8494910_pb",
    "SRR8494912_pb",
    "SRR8494913_pb",
    "SRR8494914_pb",
    "SRR8494925_pb",
    "SRR8494926_pb",
    "SRR8494927_pb",
    "SRR8494928_pb",
    "SRR8494929_pb",
    "SRR8494930_pb",
    "SRR8494931_pb",
    "SRR8494932_pb",
    "SRR8494933_pb",
    "SRR8494934_pb",
]
    
rule all:
    input:
        [f"fpa/quast/{dataset}/report.txt" for dataset in dataset],
        [f"fpa/quast/fpa_{dataset}/report.txt" for dataset in dataset]

        
rule quast:
    input:
        asm="fpa/assembly/{fpa}{prefix}_{tech}.fasta",
        ref="references/{prefix}_{tech}_ref.fasta"
    output:
        "fpa/quast/{fpa,(fpa_)?}{prefix}_{tech}/report.txt"
    output:
        "quast/{prefix}_{tech}.{scrubbing}.{asm}/report.txt"
    wildcard_constraints:
        tech="[^\.]*"
    shell:
        "quast -o fpa/quast/{wildcards.fpa}{wildcards.prefix}_{wildcards.tech}/ --min-identity 80.0 -r {input.ref} -t 16 {input.asm}"
    
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
            "./script/gfaminiasm2fasta.py {output.graph} {output.asm}"
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
