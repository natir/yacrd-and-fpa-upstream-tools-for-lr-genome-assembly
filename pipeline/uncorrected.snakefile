def tech2tech(wildcards, output):
    if "ont" in wildcards.tech:
        return "ont"
    else:
        return "pb"

include: "scrubbing.snakefile"
include: "assembly.snakefile"
include: "analysis.snakefile"

rule scrubb_droso:
    input:
        "scrubbing/d_melanogaster_reads_ont.raw.fasta",
        "scrubbing/d_melanogaster_reads_ont.4.4.yacrd.fasta",
        "scrubbing/d_melanogaster_reads_ont.4.4.precision.yacrd.fasta",
        "scrubbing/d_melanogaster_reads_ont.dascrubber.fasta",
        "scrubbing/d_melanogaster_reads_ont.miniscrub.cpu.fasta",
    
rule scrubb_ont:
    input:
        "scrubbing/real_reads_ont.raw.fasta",
        "scrubbing/real_reads_ont.4.4.yacrd.fasta",
        "scrubbing/real_reads_ont.4.4.precision.yacrd.fasta",
        "scrubbing/real_reads_ont.dascrubber.fasta",
        "scrubbing/real_reads_ont.miniscrub.cpu.fasta",
    
rule scrubb_pb:
    input:
        "scrubbing/real_reads_pb.raw.fasta",
        "scrubbing/real_reads_pb.4.4.yacrd.fasta",
        "scrubbing/real_reads_pb.4.4.precision.yacrd.fasta",
        "scrubbing/real_reads_pb.dascrubber.fasta",
        "scrubbing/real_reads_pb.miniscrub.cpu.fasta",

rule asm_droso:
    input:
        # miniasm
        "assembly/d_melanogaster_reads_ont.raw.miniasm.fasta",
        "assembly/d_melanogaster_reads_ont.4.4.yacrd.miniasm.fasta",
        "assembly/d_melanogaster_reads_ont.4.4.precision.yacrd.miniasm.fasta",
        "assembly/d_melanogaster_reads_ont.dascrubber.miniasm.fasta",
        "assembly/d_melanogaster_reads_ont.miniscrub.cpu.miniasm.fasta",
        
        # wtdbg2
        "assembly/d_melanogaster_reads_ont.raw.wtdbg2.fasta",
        "assembly/d_melanogaster_reads_ont.4.4.yacrd.wtdbg2.fasta",
        "assembly/d_melanogaster_reads_ont.4.4.precision.yacrd.wtdbg2.fasta",
        "assembly/d_melanogaster_reads_ont.dascrubber.wtdbg2.fasta",
        "assembly/d_melanogaster_reads_ont.miniscrub.cpu.wtdbg2.fasta",

    
rule asm_ont:
    input:
        # miniasm
        "assembly/real_reads_ont.raw.miniasm.fasta",
        "assembly/real_reads_ont.4.4.yacrd.miniasm.fasta",
        "assembly/real_reads_ont.4.4.precision.yacrd.miniasm.fasta",
        "assembly/real_reads_ont.dascrubber.miniasm.fasta",
        "assembly/real_reads_ont.miniscrub.cpu.miniasm.fasta",

        # wtdbg2
        "assembly/real_reads_ont.raw.wtdbg2.fasta",
        "assembly/real_reads_ont.4.4.yacrd.wtdbg2.fasta",
        "assembly/real_reads_ont.4.4.precision.yacrd.wtdbg2.fasta",
        "assembly/real_reads_ont.dascrubber.wtdbg2.fasta",
        "assembly/real_reads_ont.miniscrub.cpu.wtdbg2.fasta",

    
rule asm_pb:
    input:
        # miniasm
        "assembly/real_reads_pb.raw.miniasm.fasta",
        "assembly/real_reads_pb.4.4.yacrd.miniasm.fasta",
        "assembly/real_reads_pb.4.4.precision.yacrd.miniasm.fasta",
        "assembly/real_reads_pb.dascrubber.miniasm.fasta",
        "assembly/real_reads_pb.miniscrub.cpu.miniasm.fasta",

        # wtdbg2
        "assembly/real_reads_pb.raw.wtdbg2.fasta",
        "assembly/real_reads_pb.4.4.yacrd.wtdbg2.fasta",
        "assembly/real_reads_pb.4.4.precision.yacrd.wtdbg2.fasta",
        "assembly/real_reads_pb.dascrubber.wtdbg2.fasta",
        "assembly/real_reads_pb.miniscrub.cpu.wtdbg2.fasta",

rule all:
    input:
        # quast
        ## droso
        ### miniasm
        "quast/d_melanogaster_reads_ont.raw.miniasm/report.txt",
        "quast/d_melanogaster_reads_ont.4.4.yacrd.miniasm/report.txt",
        "quast/d_melanogaster_reads_ont.4.4.precision.yacrd.miniasm/report.txt",
        "quast/d_melanogaster_reads_ont.dascrubber.miniasm/report.txt",
        "quast/d_melanogaster_reads_ont.miniscrub.cpu.miniasm/report.txt",
        
        ### wtdbg2
        "quast/d_melanogaster_reads_ont.raw.wtdbg2/report.txt",
        "quast/d_melanogaster_reads_ont.4.4.yacrd.wtdbg2/report.txt",
        "quast/d_melanogaster_reads_ont.4.4.precision.yacrd.wtdbg2/report.txt",
        "quast/d_melanogaster_reads_ont.dascrubber.wtdbg2/report.txt",
        "quast/d_melanogaster_reads_ont.miniscrub.cpu.wtdbg2/report.txt",

        ## ont
        ### miniasm
        "quast/real_reads_ont.raw.miniasm/report.txt",
        "quast/real_reads_ont.4.4.yacrd.miniasm/report.txt",
        "quast/real_reads_ont.4.4.precision.yacrd.miniasm/report.txt",
        "quast/real_reads_ont.dascrubber.miniasm/report.txt",
        "quast/real_reads_ont.miniscrub.cpu.miniasm/report.txt",

        ### wtdbg2
        "quast/real_reads_ont.raw.wtdbg2/report.txt",
        "quast/real_reads_ont.4.4.yacrd.wtdbg2/report.txt",
        "quast/real_reads_ont.4.4.precision.yacrd.wtdbg2/report.txt",
        "quast/real_reads_ont.dascrubber.wtdbg2/report.txt",
        "quast/real_reads_ont.miniscrub.cpu.wtdbg2/report.txt",

        ## pb
        ### miniasm
        "quast/real_reads_pb.raw.miniasm/report.txt",
        "quast/real_reads_pb.4.4.yacrd.miniasm/report.txt",
        "quast/real_reads_pb.4.4.precision.yacrd.miniasm/report.txt",
        "quast/real_reads_pb.dascrubber.miniasm/report.txt",
        "quast/real_reads_pb.miniscrub.cpu.miniasm/report.txt",

        ### wtdbg2
        "quast/real_reads_pb.raw.wtdbg2/report.txt",
        "quast/real_reads_pb.4.4.yacrd.wtdbg2/report.txt",
        "quast/real_reads_pb.4.4.precision.yacrd.wtdbg2/report.txt",
        "quast/real_reads_pb.dascrubber.wtdbg2/report.txt",
        "quast/real_reads_pb.miniscrub.cpu.wtdbg2/report.txt",

        # mapping
        ## droso
        "mapping/d_melanogaster_reads_ont.raw.bam",
        "mapping/d_melanogaster_reads_ont.4.4.yacrd.bam",
        "mapping/d_melanogaster_reads_ont.4.4.precision.yacrd.bam",
        "mapping/d_melanogaster_reads_ont.dascrubber.bam",
        "mapping/d_melanogaster_reads_ont.miniscrub.cpu.bam",
    
        ## ont
        "mapping/real_reads_ont.raw.bam",
        "mapping/real_reads_ont.4.4.yacrd.bam",
        "mapping/real_reads_ont.4.4.precision.yacrd.bam",
        "mapping/real_reads_ont.dascrubber.bam",
        "mapping/real_reads_ont.miniscrub.cpu.bam",

        ## pb
        "mapping/real_reads_pb.raw.bam",
        "mapping/real_reads_pb.4.4.yacrd.bam",
        "mapping/real_reads_pb.4.4.precision.yacrd.bam",
        "mapping/real_reads_pb.dascrubber.bam",
        "mapping/real_reads_pb.miniscrub.cpu.bam",
        
