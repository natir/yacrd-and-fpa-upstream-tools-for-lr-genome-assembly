rule all:
    input:
        # raw
        "correction/real_reads_pb.raw.raw.fasta",
        "correction/real_reads_ont.raw.raw.fasta",

        # canu
        "correction/real_reads_pb.raw.canu.fasta",
        "correction/real_reads_ont.raw.canu.fasta",
        "correction/real_reads_pb.yacrd.canu.fasta",
        "correction/real_reads_ont.yacrd.canu.fasta",
        "correction/real_reads_pb.miniscrub.canu.fasta",
        "correction/real_reads_ont.miniscrub.canu.fasta",
        "correction/real_reads_pb.dascrubber.canu.fasta",
        "correction/real_reads_ont.dascrubber.canu.fasta",

        # consent
        "correction/real_reads_pb.raw.consent.fasta",
        "correction/real_reads_ont.raw.consent.fasta",
        "correction/real_reads_pb.yacrd.consent.fasta",
        "correction/real_reads_ont.yacrd.consent.fasta",
        "correction/real_reads_pb.miniscrub.consent.fasta",
        "correction/real_reads_ont.miniscrub.consent.fasta",
        "correction/real_reads_pb.dascrubber.consent.fasta",
        "correction/real_reads_ont.dascrubber.consent.fasta",
        
rule raw:
    input:
        "scrubbing/{file}.{ext}"

    output:
        "correction/{file}.raw.{ext}"

    shell:
        "ln -s {input} {output}"
        
rule canu_correct_pb:
    input:
        "scrubbing/{prefix}_pb.{scrubber}.fasta"

    output:
        "correction/{prefix}_pb.{scrubber}.canu.fasta"

    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.canu.txt",
        
    shell:
        " && ".join([
            "canu -trim -p canu -d canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber} genomeSize=5.2g -pacbio-raw {input} gnuplotTested=true executiveThread=8",
            "gzip -d -c canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber}/canu.trimmedReads.fasta.gz > {output}"
            ])

rule canu_correct_ont:
    input:
        "scrubbing/{prefix}_ont.{scrubber}.fasta"

    output:
        "correction/{prefix}_ont.{scrubber}.canu.fasta"

    benchmark:
        "benchmarks/{prefix}_ont.{scrubber}.canu.txt",
        
    shell:
        " && ".join([
            "canu -trim -p canu -d canu_asm/canu_{wildcards.prefix}_ont_{wildcards.scrubber} genomeSize=5.2g -nanopore-raw {input} gnuplotTested=true executiveThread=8",
            "gzip -d -c canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber}/canu.trimmedReads.fasta.gz > {output}"
            ])
        
rule consent_pb:
    input:
        "scrubbing/{prefix}_pb.{scrubber}.fasta"

    output:
        "correction/{prefix}_pb.{scrubber}.consent.fasta"

    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.consent.txt"

    shell:
        "CONSENT-correct --in {input} --out {output} --type PB -j 8"

rule consent_ont:
    input:
        "scrubbing/{prefix}_ont.{scrubber}.fasta"

    output:
        "correction/{prefix}_ont.{scrubber}.consent.fasta"

    benchmark:
        "benchmarks/{prefix}_ont.{scrubber}.consent.txt",
        
    shell:
        "CONSENT-correct --in {input} --out {output} --type ONT -j 8"

