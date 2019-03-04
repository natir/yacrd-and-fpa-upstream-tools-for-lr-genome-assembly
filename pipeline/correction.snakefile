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
        "correction/real_reads_pb.yacrd2.canu.fasta",
        "correction/real_reads_ont.yacrd2.canu.fasta",
        "correction/real_reads_pb.miniscrub.canu.fasta",
        "correction/real_reads_ont.miniscrub.canu.fasta",
        "correction/real_reads_pb.dascrubber.canu.fasta",
        "correction/real_reads_ont.dascrubber.canu.fasta",

        # mecat
        "correction/real_reads_pb.raw.mecat.fasta",
        "correction/real_reads_ont.raw.mecat.fasta",
        "correction/real_reads_pb.yacrd.mecat.fasta",
        "correction/real_reads_ont.yacrd.mecat.fasta",
        "correction/real_reads_pb.yacrd2.mecat.fasta",
        "correction/real_reads_ont.yacrd2.mecat.fasta",
        "correction/real_reads_pb.miniscrub.mecat.fasta",
        "correction/real_reads_ont.miniscrub.mecat.fasta",
        "correction/real_reads_pb.dascrubber.mecat.fasta",
        "correction/real_reads_ont.dascrubber.mecat.fasta",
        
        # consent
        "correction/real_reads_pb.raw.consent.fasta",
        "correction/real_reads_ont.raw.consent.fasta",
        "correction/real_reads_pb.yacrd.consent.fasta",
        "correction/real_reads_ont.yacrd.consent.fasta",
        "correction/real_reads_pb.yacrd2.consent.fasta",
        "correction/real_reads_ont.yacrd2.consent.fasta",
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
        "ln -s $(readlink -f {input}) {output}"
        
rule canu_correct_pb:
    input:
        "scrubbing/{prefix}_pb.{scrubber}.fasta"

    output:
        "correction/{prefix}_pb.{scrubber}.canu.fasta"

    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.canu.txt",
        
    shell:
        " && ".join([
            "canu -correct -p canu -d canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber} genomeSize=5.2m -pacbio-raw {input} executiveThreads=8 batMemory=160 batThreads=8 maxThreads=8",
        "canu -trim -p canu -d canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber} genomeSize=5.2m -pacbio-raw {input} executiveThreads=8 batMemory=160 batThreads=8 maxThreads=8",
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
            "canu -correct -p canu -d canu_asm/{wildcards.prefix}_ont_{wildcards.scrubber} genomeSize=5.2m -nanopore-raw {input} executiveThreads=8 batMemory=160 batThreads=8 maxThreads=8",
            "canu -trim -p canu -d canu_asm/{wildcards.prefix}_ont_{wildcards.scrubber} genomeSize=5.2m -nanopore-raw {input} executiveThreads=8 batMemory=160 batThreads=8 maxThreads=8",
            "gzip -d -c canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber}/canu.trimmedReads.fasta.gz > {output}"
            ])

rule mecat_pb:
    input:
        "scrubbing/{prefix}_pb.{scrubber}.fasta"

    output:
        corrected = "correction/{prefix}_pb.{scrubber}.mecat.fasta",
        work_dir = "mecat/{prefix}_pb_{scrubber}/fileindex.txt"
        
    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.mecat.txt"

    shell:
        " && ".join([
            "mecat2pw -j 0 -t 8 -d {input} -o correction/{wildcards.prefix}_pb.{wildcards.scrubber}.pm.can -w mecat/{wildcards.prefix}_pb_{wildcards.scrubber}",
            "mecat2cns -i 0 -t 8 correction/{wildcards.prefix}_pb.{wildcards.scrubber}.pm.can {input} {output.corrected}",
            ])

rule mecat_ont:
    input:
        "scrubbing/{prefix}_ont.{scrubber}.fasta"

    output:
        corrected = "correction/{prefix}_ont.{scrubber}.mecat.fasta",
        work_dir = "mecat/{prefix}_ont_{scrubber}/fileindex.txt"
        
    benchmark:
        "benchmarks/{prefix}_ont.{scrubber}.mecat.txt"

    shell:
        " && ".join([
            "mecat2pw -j 0 -t 8 -x 1 -d {input} -o correction/{wildcards.prefix}_ont.{wildcards.scrubber}.candidatex.txt -w mecat/{wildcards.prefix}_ont_{wildcards.scrubber}",
            "mecat2cns -i 0 -t 8 -x 1 correction/{wildcards.prefix}_ont.{wildcards.scrubber}.candidatex.txt {input} {output.corrected}"
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

