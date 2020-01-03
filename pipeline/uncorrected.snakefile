def tech2tech(wildcards, output):
    if "ont" in wildcards.tech:
        return "ont"
    else:
        return "pb"

include: "scrubbing.snakefile"
include: "assembly.snakefile"
include: "analysis.snakefile"

bwa_str = "mapping/{dataset}.{scrubber}.bam"
minimap_str = "mapping/{dataset}.{scrubber}.paf"        
porechop_str = "porechop/{dataset}.{scrubber}.out"

wtdbg2_str = "assembly/{dataset}.{scrubber}.wtdbg2.fasta"
miniasm_str = "assembly/{dataset}.{scrubber}.miniasm.fasta"

quast_str = "quast/{dataset}.{scrubber}.{assembly}/report.txt"
nucmer_str = "nucmer/{dataset}.{scrubber}.{assembly}.delta"
quast_lr_str = "quast_lr/{dataset}.{scrubber}.{assembly}/report.txt"

assembly_list = [
    "miniasm",
    "wtdbg2",
#    "ra",
#    "shasta"
]

scrubber_list = [
    "raw",
    "dascrubber",
]

yacrd_p6c4 = ["g800.c4.yacrd"]
yacrd_sequel = ["g5000.c3.yacrd"]
yacrd_nanopore = ["g500.c4.yacrd"]

def c_elegans_out():
    d = "c_elegans_pb"
    for s in scrubber_list + yacrd_p6c4:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        

        for a in assembly_list:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
            yield quast_lr_str.format(dataset=d, scrubber=s, assembly=a)
            
rule c_elegans:
    input:
        c_elegans_out()
        
def h_sapiens_out():
    d = "h_sapiens_chr1_ont"
    for s in scrubber_list + yacrd_nanopore:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in assembly_list:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
            yield quast_lr_str.format(dataset=d, scrubber=s, assembly=a)
                    
rule h_sapiens:
    input:
        h_sapiens_out()

def d_melanogaster_out():
    d = "d_melanogaster_reads_ont"
    for s in scrubber_list + yacrd_nanopore:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in assembly_list:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
            yield quast_lr_str.format(dataset=d, scrubber=s, assembly=a)
            
rule d_melanogaster:
    input:
        d_melanogaster_out()

def e_coli_ont_out():
    d = "real_reads_ont"
    for s in scrubber_list + yacrd_nanopore:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in assembly_list:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
            yield quast_lr_str.format(dataset=d, scrubber=s, assembly=a)
            
def e_coli_pb_out():
    d = "real_reads_pb"
    for s in scrubber_list + yacrd_sequel:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in assembly_list:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
            yield quast_lr_str.format(dataset=d, scrubber=s, assembly=a)
            
def NCTC():
    dataset = [
        "ERR2531093_pb",
        "ERR2531094_pb",
        "ERR2531095_pb",
        "ERR2531096_pb",
        "ERR2564864_pb",
        "ERR2564865_pb",
        "ERR2603033_pb",
        "ERR2615933_pb",
        "ERR2615934_pb",
        "ERR2632359_pb",
        "ERR2651535_pb",
        "ERR2651536_pb",
        "ERR2672424_pb",
        "ERR2672425_pb",
        "ERR2681660_pb",
        "ERR2695057_pb",
        "ERR2695058_pb",
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
    ]

    dascrubber_skip = [
        "ERR2695058_pb",
        "ERR2531095_pb",
        "ERR2695057_pb",
        "ERR2672425_pb",
    ]
    for d in dataset:
        for s in scrubber_list + yacrd_sequel:
            if s == "dascrubber" and d in dascrubber_skip:
                continue
            #yield bwa_str.format(dataset=d, scrubber=s)
            #yield minimap_str.format(dataset=d, scrubber=s)
            yield porechop_str.format(dataset=d, scrubber=s)

            yield wtdbg2_str.format(dataset=d, scrubber=s)
            yield miniasm_str.format(dataset=d, scrubber=s)

rule NCTC:
    input:
        NCTC(),
            
def nanopore2pacbio():
    dataset = [
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
        #"SRR8494940_ont", real_reads_ont.fasta
        "SRR8494941_ont",
        "SRR8494942_ont",
        "SRR8494943_ont",
        "SRR8494944_ont",
    ]
    for d in dataset:
        for s in scrubber_list + yacrd_nanopore:
            #yield bwa_str.format(dataset=d, scrubber=s)
            #yield minimap_str.format(dataset=d, scrubber=s)
            yield porechop_str.format(dataset=d, scrubber=s)

            yield wtdbg2_str.format(dataset=d, scrubber=s)
            yield miniasm_str.format(dataset=d, scrubber=s)

    dataset = [
        "SRR8494905_pb",
        "SRR8494906_pb",
        "SRR8494907_pb",
        "SRR8494908_pb",
        "SRR8494909_pb",
        "SRR8494910_pb",
        #"SRR8494911_pb", real_reads_pb.fasta
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
    for d in dataset:
        for s in scrubber_list + yacrd_p6c4:
            #yield bwa_str.format(dataset=d, scrubber=s)
            #yield minimap_str.format(dataset=d, scrubber=s)
            yield porechop_str.format(dataset=d, scrubber=s)

            yield wtdbg2_str.format(dataset=d, scrubber=s)
            yield miniasm_str.format(dataset=d, scrubber=s)

rule nanopore2pacbio:
    input:
        nanopore2pacbio(),
            
rule e_coli:
    input:
        e_coli_ont_out(),
        e_coli_pb_out()
        
rule all:
    input:
        rules.c_elegans.input,
        rules.h_sapiens.input,
        rules.d_melanogaster.input,
        rules.e_coli.input,
        rules.NCTC.input,
        rules.nanopore2pacbio.input,


