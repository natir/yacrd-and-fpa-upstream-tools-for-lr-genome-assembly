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

quast_str = "quast/{dataset}.{scrubber}.{assembly}/report.txt"
nucmer_str = "nucmer/{dataset}.{scrubber}.{assembly}.delta"

def c_elegans_out():
    d = "c_elegans_pb"
    for s in ["raw", "g800.c4.yacrd", "dascrubber"]:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        

        for a in ["miniasm", "wtdbg2", "ra", "shasta"]:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)

rule c_elegans:
    input:
        c_elegans_out()

def h_sapiens_out():
    d = "h_sapiens_chr1_ont"
    for s in ["raw", "g500.c4.yacrd", "dascrubber"]:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in ["miniasm", "wtdbg2", "ra", "shasta"]:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
        
rule h_sapiens:
    input:
        h_sapiens_out()

def d_melanogaster_out():
    d = "d_melanogaster_reads_ont"
    for s in ["raw", "g500.c4.yacrd", "dascrubber"]:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in ["miniasm", "wtdbg2", "ra", "shasta"]:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
        
rule d_melanogaster:
    input:
        d_melanogaster_out()


def e_coli_ont_out():
    d = "real_reads_ont"
    for s in ["raw", "g500.c4.yacrd", "dascrubber", "miniscrub.cpu"]:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in ["miniasm", "wtdbg2", "ra", "shasta"]:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
        
def e_coli_pb_out():
    d = "real_reads_ont"
    for s in ["raw", "g5000.c3.yacrd", "dascrubber", "miniscrub.cpu"]:
        yield bwa_str.format(dataset=d, scrubber=s)
        yield minimap_str.format(dataset=d, scrubber=s)
        yield porechop_str.format(dataset=d, scrubber=s)        
        
        for a in ["miniasm", "wtdbg2", "ra", "shasta"]:
            yield quast_str.format(dataset=d, scrubber=s, assembly=a)
            yield nucmer_str.format(dataset=d, scrubber=s, assembly=a)
            
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
