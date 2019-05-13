def tech2tech(wildcards, output):
    if "ont" in wildcards.tech:
        return "ont"
    else:
        return "pb"

include: "scrubbing.snakefile"
include: "assembly.snakefile"
include: "analysis.snakefile"

scrubbing_suffix = ["raw", "4.4.yacrd", "4.4.precision.yacrd", "dascrubber"]#, "miniscrub.cpu"]

def generate_scrubb(prefix):
    return ["scrubbing/{}.{}.fasta".format(prefix, suffix) for suffix in scrubbing_suffix]


rule scrubb_worm:
    input:
        generate_scrubb("c_elegans_pb")


rule scrubb_human:
    input:
        generate_scrubb("h_sapiens_chr1_ont")


rule scrubb_droso:
    input:
        generate_scrubb("d_melanogaster_reads_ont")


rule scrubb_ont:
    input:
        generate_scrubb("real_reads_ont")


rule scrubb_pb:
    input:
        generate_scrubb("real_reads_pb")



assembly_suffix = ["miniasm", "wtdbg2"]

def generate_assembly(prefix):
    return ["scrubbing/{}.{}.{}.fasta".format(prefix, scrub, asm) for scrub in scrubbing_suffix for asm in assembly_suffix]


rule asm_worm:
    input:
        generate_assembly("c_elegans_pb")


rule asm_human:
    input:
        generate_assembly("h_sapiens_chr1_ont")


rule asm_droso:
    input:
        generate_assembly("d_melanogaster_reads_ont")


rule asm_ont:
    input:
        generate_assembly("real_reads_ont")


rule asm_pb:
    input:
        generate_assembly("real_reads_pb")


def generate_quast(prefix):
    return ["quast/{}.{}.{}/report.txt".format(prefix, scrub, asm) for scrub in scrubbing_suffix for asm in assembly_suffix]

def generate_mapping(prefix):
    return ["mapping/{}.{}.bam".format(prefix, scrub) for scrub in scrubbing_suffix]

def generate_minimap2(prefix):
    return ["mapping/{}.{}.paf".format(prefix, scrub) for scrub in scrubbing_suffix]

rule all:
    input:
        generate_quast("c_elegans_pb"),
        generate_quast("h_sapiens_chr1_ont"),
        generate_quast("d_melanogaster_reads_ont"),
        generate_quast("real_reads_ont"),
        generate_quast("real_reads_pb"),

        generate_mapping("c_elegans_pb"),
        generate_mapping("h_sapiens_chr1_ont"),
        generate_mapping("d_melanogaster_reads_ont"),
        generate_mapping("real_reads_ont"),
        generate_mapping("real_reads_pb"),

        generate_minimap2("c_elegans_pb"),
        generate_minimap2("h_sapiens_chr1_ont"),
        generate_minimap2("d_melanogaster_reads_ont"),
        generate_minimap2("real_reads_ont"),
        generate_minimap2("real_reads_pb"),
