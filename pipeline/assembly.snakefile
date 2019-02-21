rule all:
    input:
        # raw
        "correction/real_reads_pb.raw.raw.canu.fasta"
        "correction/real_reads_ont.raw.raw.canu.fasta"
        "correction/real_reads_pb.raw.raw.miniasm.fasta"
        "correction/real_reads_ont.raw.raw.miniasm.fasta"
        
        # canu
        "correction/real_reads_pb.raw.canu.canu.fasta"
        "correction/real_reads_ont.raw.canu.canu.fasta"
        "correction/real_reads_pb.yacrd.canu.canu.fasta"
        "correction/real_reads_ont.yacrd.canu.canu.fasta"
        "correction/real_reads_pb.miniscrub.canu.canu.fasta"
        "correction/real_reads_ont.miniscrub.canu.canu.fasta"
        "correction/real_reads_pb.dascrubber.canu.canu.fasta"
        "correction/real_reads_ont.dascrubber.canu.canu.fasta"
        "correction/real_reads_pb.raw.consent.canu.fasta"
        "correction/real_reads_ont.raw.consent.canu.fasta"
        "correction/real_reads_pb.yacrd.consent.canu.fasta"
        "correction/real_reads_ont.yacrd.consent.canu.fasta"
        "correction/real_reads_pb.miniscrub.consent.canu.fasta"
        "correction/real_reads_ont.miniscrub.consent.canu.fasta"
        "correction/real_reads_pb.dascrubber.consent.canu.fasta"
        "correction/real_reads_ont.dascrubber.consent.canu.fasta"

        # miniasm
        "correction/real_reads_pb.raw.canu.miniasm.fasta"
        "correction/real_reads_ont.raw.canu.miniasm.fasta"
        "correction/real_reads_pb.yacrd.canu.miniasm.fasta"
        "correction/real_reads_ont.yacrd.canu.miniasm.fasta"
        "correction/real_reads_pb.miniscrub.canu.miniasm.fasta"
        "correction/real_reads_ont.miniscrub.canu.miniasm.fasta"
        "correction/real_reads_pb.dascrubber.canu.miniasm.fasta"
        "correction/real_reads_ont.dascrubber.canu.miniasm.fasta"
        "correction/real_reads_pb.raw.consent.miniasm.fasta"
        "correction/real_reads_ont.raw.consent.miniasm.fasta"
        "correction/real_reads_pb.yacrd.consent.miniasm.fasta"
        "correction/real_reads_ont.yacrd.consent.miniasm.fasta"
        "correction/real_reads_pb.miniscrub.consent.miniasm.fasta"
        "correction/real_reads_ont.miniscrub.consent.miniasm.fasta"
        "correction/real_reads_pb.dascrubber.consent.miniasm.fasta"
        "correction/real_reads_ont.dascrubber.consent.miniasm.fasta"

rule canu_pb:
    input:
        "correction/{prefix}_pb.{scrubber}.{corrector}.{ext}"

    output:
        "assembly/{prefix}_pb.{scrubber}.{corrector}.gfa"
        "assembly/{prefix}_pb.{scrubber}.{corrector}.fasta"
        
    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.{corrector}.canu.txt",
        
    shell:
        "canu -p canu -d canu_asm/{prefix}_pb_{scrubber}_corrector genomeSize=4.6g -pacbio-corrected {input} executiveThread=8"

rule canu_ont:
    input:
        "correction/{prefix}_ont.{scrubber}.{corrector}.{ext}"

    output:
        "assembly/{prefix}_ont.{scrubber}.{corrector}.gfa"
        "assembly/{prefix}_ont.{scrubber}.{corrector}.fasta"
        
    benchmark:
        "benchmarks/{prefix}_ont.{scrubber}.{corrector}.canu.txt",
        
    shell:
        "canu -p canu -d canu_asm/{prefix}_ont_{scrubber}_corrector genomeSize=4.6g -nanopore-corrected {input} executiveThread=8"

rule miniasm_pb:
    input:
        "correction/{prefix}_pb.{scrubber}.{corrector}.{ext}"

    output:
        "assembly/{prefix}_pb.{scrubber}.{corrector}.gfa"
        "assembly/{prefix}_pb.{scrubber}.{corrector}.fasta"
        
    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.{corrector}.miniasm.txt",
        
    shell:
        " && ".join([
            "minimap2 -t 8 -x ava-pb {input} {input} > assembly/{prefix}_pb.{scrubber}.{corrector}.paf",
            "miniasm -t 8 -f {input} assembly/{prefix}_pb.{scrubber}.{corrector}.paf > assembly/{prefix}_pb.{scrubber}.{corrector}.gfa",
            "./script/gfaminiasm2fasta.py assembly/{prefix}_pb.{scrubber}.{corrector}.gfa assembly/{prefix}_pb.{scrubber}.{corrector}.fasta"
        ])

rule miniasm_ont:
    input:
        "correction/{prefix}_ont.{scrubber}.{corrector}.{ext}"

    output:
        "assembly/{prefix}_ont.{scrubber}.{corrector}.gfa"
        "assembly/{prefix}_ont.{scrubber}.{corrector}.fasta"
        
    benchmark:
        "benchmarks/{prefix}_ont.{scrubber}.{corrector}.miniasm.txt",
        
    shell:
        " && ".join([
            "minimap2 -t 8 -x ava-ont {input} {input} > assembly/{prefix}_ont.{scrubber}.{corrector}.paf",
            "miniasm -t 8 -f {input} assembly/{prefix}_ont.{scrubber}.{corrector}.paf > assembly/{prefix}_ont.{scrubber}.{corrector}.gfa",
            "./script/gfaminiasm2fasta.py assembly/{prefix}_ont.{scrubber}.{corrector}.gfa assembly/{prefix}_ont.{scrubber}.{corrector}.fasta"
        ])
