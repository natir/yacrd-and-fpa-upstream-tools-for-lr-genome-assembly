rule mapping_scrubbing:
    input:
        "mapping/scrubbing/real_reads_pb.raw.bam",
        "mapping/scrubbing/real_reads_pb.yacrd.bam",
        "mapping/scrubbing/real_reads_pb.yacrd2.bam",
        "mapping/scrubbing/real_reads_pb.dascrubber.bam",
        "mapping/scrubbing/real_reads_pb.miniscrub.bam",
        "mapping/scrubbing/real_reads_ont.raw.bam",
        "mapping/scrubbing/real_reads_ont.yacrd.bam",
        "mapping/scrubbing/real_reads_ont.yacrd2.bam",
        "mapping/scrubbing/real_reads_ont.dascrubber.bam",
        "mapping/scrubbing/real_reads_ont.miniscrub.bam",

rule mapping_correction:
    input:
        # raw
        "mapping/correction/real_reads_pb.raw.raw.bam",
        "mapping/correction/real_reads_ont.raw.raw.bam",

        # canu
        "mapping/correction/real_reads_pb.raw.canu.bam",
        "mapping/correction/real_reads_ont.raw.canu.bam",
        "mapping/correction/real_reads_pb.yacrd.canu.bam",
        "mapping/correction/real_reads_ont.yacrd.canu.bam",
        "mapping/correction/real_reads_pb.yacrd2.canu.bam",
        "mapping/correction/real_reads_ont.yacrd2.canu.bam",
        "mapping/correction/real_reads_pb.miniscrub.canu.bam",
        "mapping/correction/real_reads_ont.miniscrub.canu.bam",
        "mapping/correction/real_reads_pb.dascrubber.canu.bam",
        "mapping/correction/real_reads_ont.dascrubber.canu.bam",

        # mecat
        "mapping/correction/real_reads_pb.raw.mecat.bam",
        "mapping/correction/real_reads_ont.raw.mecat.bam",
        "mapping/correction/real_reads_pb.yacrd.mecat.bam",
        "mapping/correction/real_reads_ont.yacrd.mecat.bam",
        "mapping/correction/real_reads_pb.yacrd2.mecat.bam",
        "mapping/correction/real_reads_ont.yacrd2.mecat.bam",
        "mapping/correction/real_reads_pb.miniscrub.mecat.bam",
        "mapping/correction/real_reads_ont.miniscrub.mecat.bam",
        "mapping/correction/real_reads_pb.dascrubber.mecat.bam",
        # "mapping/correction/real_reads_ont.dascrubber.mecat.bam",
        
        # consent
        "mapping/correction/real_reads_pb.raw.consent.bam",
        "mapping/correction/real_reads_ont.raw.consent.bam",
        "mapping/correction/real_reads_pb.yacrd.consent.bam",
        "mapping/correction/real_reads_ont.yacrd.consent.bam",
        "mapping/correction/real_reads_pb.yacrd2.consent.bam",
        "mapping/correction/real_reads_ont.yacrd2.consent.bam",
        "mapping/correction/real_reads_pb.miniscrub.consent.bam",
        "mapping/correction/real_reads_ont.miniscrub.consent.bam",
        "mapping/correction/real_reads_pb.dascrubber.consent.bam",
        "mapping/correction/real_reads_ont.dascrubber.consent.bam",
        
        
rule mapping_ref_pb:
    input:
        "data/ref_e_coli_cft073.fasta.sa",
        "data/ref_e_coli_cft073.fasta.amb",
        "data/ref_e_coli_cft073.fasta.ann",
        "data/ref_e_coli_cft073.fasta.bwt",
        "data/ref_e_coli_cft073.fasta.pac",
        reads = "{prefix}pb{suffix}.fasta",
        
    output:
        "mapping/{prefix}pb{suffix}.bam"        
    shell:
        " && ".join([
            "bwa mem -t 8 -x pacbio data/ref_e_coli_cft073.fasta {input.reads} | samtools sort > {output}",
            "samtools index {output}"
        ])

rule mapping_ref_ont:
    input:
        "data/ref_e_coli_cft073.fasta.sa",
        "data/ref_e_coli_cft073.fasta.amb",
        "data/ref_e_coli_cft073.fasta.ann",
        "data/ref_e_coli_cft073.fasta.bwt",
        "data/ref_e_coli_cft073.fasta.pac",
        reads = "{prefix}ont{suffix}.fasta",
        
    output:
        "mapping/{prefix}ont{suffix}.bam"
        
    shell:
        " && ".join([
            "bwa mem -t 8 -x ont2d data/ref_e_coli_cft073.fasta {input.reads} | samtools sort > {output}",
            "samtools index {output}"
        ])
