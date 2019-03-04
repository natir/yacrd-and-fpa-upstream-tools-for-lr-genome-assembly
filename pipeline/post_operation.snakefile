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

rule quast:
    input:
        ref="data/ref_e_coli_cft073.fasta",
        asm="assembly/{prefix}.fasta",

    output:
        "quast/{prefix}/report.tsv"

    shell:
        "quast -o quast/{wildcards.prefix}/ -r {input.ref} -t 8 {input.asm}"

rule sub_sample:
    input:
        "data/real_reads_pb.25.fasta",
        "data/real_reads_pb.50.fasta",
        
        
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

rule mapping_corrected:
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

rule evaluate_assembly:
    input:
        # miniasm
        ## raw
        "quast/real_reads_pb.raw.raw.miniasm/report.tsv",
        "quast/real_reads_ont.raw.raw.miniasm/report.tsv",
        "quast/real_reads_pb.raw.canu.miniasm/report.tsv",
        "quast/real_reads_ont.raw.canu.miniasm/report.tsv",
        "quast/real_reads_pb.raw.consent.miniasm/report.tsv",
        "quast/real_reads_ont.raw.consent.miniasm/report.tsv",
        "quast/real_reads_pb.raw.mecat.miniasm/report.tsv",
        "quast/real_reads_ont.raw.mecat.miniasm/report.tsv",
        ## yacrd
	# "quast/real_reads_pb.yacrd.raw.miniasm/report.tsv",
        # "quast/real_reads_ont.yacrd.raw.miniasm/report.tsv",
        "quast/real_reads_pb.yacrd.canu.miniasm/report.tsv",
        "quast/real_reads_ont.yacrd.canu.miniasm/report.tsv",
        "quast/real_reads_pb.yacrd.consent.miniasm/report.tsv",
        "quast/real_reads_ont.yacrd.consent.miniasm/report.tsv",
        "quast/real_reads_pb.yacrd.mecat.miniasm/report.tsv",
        "quast/real_reads_ont.yacrd.mecat.miniasm/report.tsv",
        "quast/real_reads_pb.yacrd2.canu.miniasm/report.tsv",
        "quast/real_reads_ont.yacrd2.canu.miniasm/report.tsv",
        "quast/real_reads_pb.yacrd2.consent.miniasm/report.tsv",
        "quast/real_reads_ont.yacrd2.consent.miniasm/report.tsv",
        "quast/real_reads_pb.yacrd2.mecat.miniasm/report.tsv",
        "quast/real_reads_ont.yacrd2.mecat.miniasm/report.tsv",
        ## dascrubber
	# "quast/real_reads_pb.dascrubber.raw.miniasm/report.tsv",
        # "quast/real_reads_ont.dascrubber.raw.miniasm/report.tsv",
        "quast/real_reads_pb.dascrubber.canu.miniasm/report.tsv",
        "quast/real_reads_ont.dascrubber.canu.miniasm/report.tsv",
        "quast/real_reads_pb.dascrubber.consent.miniasm/report.tsv",
        "quast/real_reads_ont.dascrubber.consent.miniasm/report.tsv",
        "quast/real_reads_pb.dascrubber.mecat.miniasm/report.tsv",
        # "quast/real_reads_ont.dascrubber.mecat.miniasm/report.tsv", # correction not work
        ## miniscrub
	# "quast/real_reads_pb.miniscrub.raw.miniasm/report.tsv",
        # "quast/real_reads_ont.miniscrub.raw.miniasm/report.tsv",
        "quast/real_reads_pb.miniscrub.canu.miniasm/report.tsv",
        "quast/real_reads_ont.miniscrub.canu.miniasm/report.tsv",
        "quast/real_reads_pb.miniscrub.consent.miniasm/report.tsv",
        "quast/real_reads_ont.miniscrub.consent.miniasm/report.tsv",
        "quast/real_reads_pb.miniscrub.mecat.miniasm/report.tsv",
        "quast/real_reads_ont.miniscrub.mecat.miniasm/report.tsv",
        
        # canu
	## raw
	# "quast/real_reads_pb.raw.raw.canu/report.tsv",
        # "quast/real_reads_ont.raw.raw.canu/report.tsv",
        "quast/real_reads_pb.raw.canu.canu/report.tsv",
        "quast/real_reads_ont.raw.canu.canu/report.tsv",
        "quast/real_reads_pb.raw.consent.canu/report.tsv",
        "quast/real_reads_ont.raw.consent.canu/report.tsv",
        "quast/real_reads_pb.raw.mecat.canu/report.tsv",
        "quast/real_reads_ont.raw.mecat.canu/report.tsv",
        ## yacrd
	# "quast/real_reads_pb.yacrd.raw.canu/report.tsv",
        # "quast/real_reads_ont.yacrd.raw.canu/report.tsv",
        "quast/real_reads_pb.yacrd.canu.canu/report.tsv",
        "quast/real_reads_ont.yacrd.canu.canu/report.tsv",
        "quast/real_reads_pb.yacrd.consent.canu/report.tsv",
        "quast/real_reads_ont.yacrd.consent.canu/report.tsv",
        "quast/real_reads_pb.yacrd.mecat.canu/report.tsv",
        "quast/real_reads_ont.yacrd.mecat.canu/report.tsv",
        "quast/real_reads_pb.yacrd2.canu.canu/report.tsv",
        "quast/real_reads_ont.yacrd2.canu.canu/report.tsv",
        "quast/real_reads_pb.yacrd2.consent.canu/report.tsv",
        "quast/real_reads_ont.yacrd2.consent.canu/report.tsv",
        "quast/real_reads_pb.yacrd2.mecat.canu/report.tsv",
        "quast/real_reads_ont.yacrd2.mecat.canu/report.tsv",        ## dascrubber
	# "quast/real_reads_pb.dascrubber.raw.canu/report.tsv",
        # "quast/real_reads_ont.dascrubber.raw.canu/report.tsv",
        "quast/real_reads_pb.dascrubber.canu.canu/report.tsv",
        "quast/real_reads_ont.dascrubber.canu.canu/report.tsv",
        "quast/real_reads_pb.dascrubber.consent.canu/report.tsv",
        "quast/real_reads_ont.dascrubber.consent.canu/report.tsv",
        "quast/real_reads_pb.dascrubber.mecat.canu/report.tsv",
        # "quast/real_reads_ont.dascrubber.mecat.canu/report.tsv", # correction not work
        ## miniscrub
	# "quast/real_reads_pb.miniscrub.raw.canu/report.tsv",
        # "quast/real_reads_ont.miniscrub.raw.canu/report.tsv",
        "quast/real_reads_pb.miniscrub.canu.canu/report.tsv",
        "quast/real_reads_ont.miniscrub.canu.canu/report.tsv",
        # "quast/real_reads_pb.miniscrub.consent.canu/report.tsv",
        "quast/real_reads_ont.miniscrub.consent.canu/report.tsv",
        # "quast/real_reads_pb.miniscrub.mecat.canu/report.tsv",
        "quast/real_reads_ont.miniscrub.mecat.canu/report.tsv",
        
