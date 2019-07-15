def tech2tech(wildcards, output):
    if "ont" in wildcards.tech:
        return "ont"
    else:
        return "pb"

include: "scrubbing.snakefile"
include: "assembly.snakefile"
include: "analysis.snakefile"


rule c_elegans:
    input:
        "quast/c_elegans_pb.raw.miniasm/report.txt",
        "quast/c_elegans_pb.raw.wtdbg2/report.txt",
        "quast/c_elegans_pb.raw.ra/report.txt",
        "quast/c_elegans_pb.g800.c4.yacrd.miniasm/report.txt",
        "quast/c_elegans_pb.g800.c4.yacrd.wtdbg2/report.txt",
        "quast/c_elegans_pb.g800.c4.yacrd.ra/report.txt",
        "quast/c_elegans_pb.dascrubber.miniasm/report.txt",
        "quast/c_elegans_pb.dascrubber.wtdbg2/report.txt",
        "quast/c_elegans_pb.dascrubber.ra/report.txt",
        #"quast/c_elegans_pb.miniscrub.cpu.miniasm/report.txt",
        #"quast/c_elegans_pb.miniscrub.cpu.wtdbg2/report.txt",
        #"quast/c_elegans_pb.miniscrub.cpu.ra/report.txt",
        
        "mapping/c_elegans_pb.raw.bam",
        "mapping/c_elegans_pb.g800.c4.yacrd.bam",
        "mapping/c_elegans_pb.dascrubber.bam",
        #"mapping/c_elegans_pb.miniscrub.cpu.bam",

        "mapping/c_elegans_pb.raw.paf",
        "mapping/c_elegans_pb.g800.c4.yacrd.paf",
        "mapping/c_elegans_pb.dascrubber.paf",
        #"mapping/c_elegans_pb.miniscrub.cpu.paf",

    
rule h_sapiens:
    input:
        "quast/h_sapiens_chr1_ont.raw.miniasm/report.txt",
        "quast/h_sapiens_chr1_ont.raw.wtdbg2/report.txt",
        "quast/h_sapiens_chr1_ont.raw.ra/report.txt",
        "quast/h_sapiens_chr1_ont.g500.c4.yacrd.miniasm/report.txt",
        "quast/h_sapiens_chr1_ont.g500.c4.yacrd.wtdbg2/report.txt",
        "quast/h_sapiens_chr1_ont.g500.c4.yacrd.ra/report.txt",
        "quast/h_sapiens_chr1_ont.dascrubber.miniasm/report.txt",
        "quast/h_sapiens_chr1_ont.dascrubber.wtdbg2/report.txt",
        "quast/h_sapiens_chr1_ont.dascrubber.ra/report.txt",
        #"quast/h_sapiens_chr1_ont.miniscrub.cpu.miniasm/report.txt",
        #"quast/h_sapiens_chr1_ont.miniscrub.cpu.wtdbg2/report.txt",
        #"quast/h_sapiens_chr1_ont.miniscrub.cpu.ra/report.txt",
        
        "mapping/h_sapiens_chr1_ont.raw.bam",
        "mapping/h_sapiens_chr1_ont.g500.c4.yacrd.bam",
        "mapping/h_sapiens_chr1_ont.dascrubber.bam",
        #"mapping/h_sapiens_chr1_ont.miniscrub.cpu.bam",

        "mapping/h_sapiens_chr1_ont.raw.paf",
        "mapping/h_sapiens_chr1_ont.g500.c4.yacrd.paf",
        "mapping/h_sapiens_chr1_ont.dascrubber.paf",
        #"mapping/h_sapiens_chr1_ont.miniscrub.cpu.paf",

        "porechop/h_sapiens_chr1_ont.raw.out",
        "porechop/h_sapiens_chr1_ont.g500.c4.yacrd.out",
        "porechop/h_sapiens_chr1_ont.dascrubber.out",

        
rule d_melanogaster:
    input:
        "quast/d_melanogaster_reads_ont.raw.miniasm/report.txt",
        "quast/d_melanogaster_reads_ont.raw.wtdbg2/report.txt",
        "quast/d_melanogaster_reads_ont.raw.ra/report.txt",
        "quast/d_melanogaster_reads_ont.g500.c4.yacrd.miniasm/report.txt",
        "quast/d_melanogaster_reads_ont.g500.c4.yacrd.wtdbg2/report.txt",
        "quast/d_melanogaster_reads_ont.g500.c4.yacrd.ra/report.txt",
        "quast/d_melanogaster_reads_ont.dascrubber.miniasm/report.txt",
        "quast/d_melanogaster_reads_ont.dascrubber.wtdbg2/report.txt",
        "quast/d_melanogaster_reads_ont.dascrubber.ra/report.txt",
        #"quast/d_melanogaster_reads_ont.miniscrub.cpu.miniasm/report.txt",
        #"quast/d_melanogaster_reads_ont.miniscrub.cpu.wtdbg2/report.txt",
        #"quast/d_melanogaster_reads_ont.miniscrub.cpu.ra/report.txt",

        "mapping/d_melanogaster_reads_ont.raw.bam",
        "mapping/d_melanogaster_reads_ont.g500.c4.yacrd.bam",
        "mapping/d_melanogaster_reads_ont.dascrubber.bam",
        #"mapping/d_melanogaster_reads_ont.miniscrub.cpu.bam",

        "mapping/d_melanogaster_reads_ont.raw.paf",
        "mapping/d_melanogaster_reads_ont.g500.c4.yacrd.paf",
        "mapping/d_melanogaster_reads_ont.dascrubber.paf",
        #"mapping/d_melanogaster_reads_ont.miniscrub.cpu.paf",

        "porechop/d_melanogaster_reads_ont.raw.out",
        "porechop/d_melanogaster_reads_ont.g500.c4.yacrd.out",
        "porechop/d_melanogaster_reads_ont.dascrubber.out",

        
rule e_coli:
    input:
        "quast/real_reads_ont.raw.miniasm/report.txt",
        "quast/real_reads_ont.raw.wtdbg2/report.txt",
        "quast/real_reads_ont.raw.ra/report.txt",
        "quast/real_reads_ont.g500.c4.yacrd.miniasm/report.txt",
        "quast/real_reads_ont.g500.c4.yacrd.wtdbg2/report.txt",
        "quast/real_reads_ont.g500.c4.yacrd.ra/report.txt",
        "quast/real_reads_ont.dascrubber.miniasm/report.txt",
        "quast/real_reads_ont.dascrubber.wtdbg2/report.txt",
        "quast/real_reads_ont.dascrubber.ra/report.txt",
        "quast/real_reads_ont.miniscrub.cpu.miniasm/report.txt",
        "quast/real_reads_ont.miniscrub.cpu.wtdbg2/report.txt",
        "quast/real_reads_ont.miniscrub.cpu.ra/report.txt",

        "mapping/real_reads_ont.raw.bam",
        "mapping/real_reads_ont.g500.c4.yacrd.bam",
        "mapping/real_reads_ont.dascrubber.bam",
        "mapping/real_reads_ont.miniscrub.cpu.bam",

        "mapping/real_reads_ont.raw.paf",
        "mapping/real_reads_ont.g500.c4.yacrd.paf",
        "mapping/real_reads_ont.dascrubber.paf",
        "mapping/real_reads_ont.miniscrub.cpu.paf",

        "porechop/real_reads_ont.raw.out",
        "porechop/real_reads_ont.g500.c4.yacrd.out",
        "porechop/real_reads_ont.dascrubber.out",
        "porechop/real_reads_ont.miniscrub.cpu.out",
        

        "quast/real_reads_pb.raw.miniasm/report.txt",
        "quast/real_reads_pb.raw.wtdbg2/report.txt",
        "quast/real_reads_pb.raw.ra/report.txt",
        "quast/real_reads_pb.g5000.c3.yacrd.miniasm/report.txt",
        "quast/real_reads_pb.g5000.c3.yacrd.wtdbg2/report.txt",
        "quast/real_reads_pb.g5000.c3.yacrd.ra/report.txt",
        "quast/real_reads_pb.dascrubber.miniasm/report.txt",
        "quast/real_reads_pb.dascrubber.wtdbg2/report.txt",
        "quast/real_reads_pb.dascrubber.ra/report.txt",
        "quast/real_reads_pb.miniscrub.cpu.miniasm/report.txt",
        "quast/real_reads_pb.miniscrub.cpu.wtdbg2/report.txt",
        "quast/real_reads_pb.miniscrub.cpu.ra/report.txt",

        "mapping/real_reads_pb.raw.bam",
        "mapping/real_reads_pb.g5000.c3.yacrd.bam",
        "mapping/real_reads_pb.dascrubber.bam",
        "mapping/real_reads_pb.miniscrub.cpu.bam",
        
        "mapping/real_reads_pb.raw.paf",
        "mapping/real_reads_pb.g5000.c3.yacrd.paf",
        "mapping/real_reads_pb.dascrubber.paf",
        "mapping/real_reads_pb.miniscrub.cpu.paf",
    
rule all:
    input:
        rules.c_elegans.input,
        rules.h_sapiens.input,
        rules.d_melanogaster.input,
        rules.e_coli.input,
