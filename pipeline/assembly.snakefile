rule all:
    input:
	# miniasm
        ## raw
        "assembly/real_reads_pb.raw.raw.miniasm.fasta",
        "assembly/real_reads_ont.raw.raw.miniasm.fasta",
        "assembly/real_reads_pb.raw.canu.miniasm.fasta",
        "assembly/real_reads_ont.raw.canu.miniasm.fasta",
        "assembly/real_reads_pb.raw.consent.miniasm.fasta",
        "assembly/real_reads_ont.raw.consent.miniasm.fasta",
        "assembly/real_reads_pb.raw.mecat.miniasm.fasta",
        "assembly/real_reads_ont.raw.mecat.miniasm.fasta",
        ## yacrd
	# "assembly/real_reads_pb.yacrd.raw.miniasm.fasta",
        # "assembly/real_reads_ont.yacrd.raw.miniasm.fasta",
        "assembly/real_reads_pb.yacrd.canu.miniasm.fasta",
        "assembly/real_reads_ont.yacrd.canu.miniasm.fasta",
        "assembly/real_reads_pb.yacrd.consent.miniasm.fasta",
        "assembly/real_reads_ont.yacrd.consent.miniasm.fasta",
        "assembly/real_reads_pb.yacrd.mecat.miniasm.fasta",
        "assembly/real_reads_ont.yacrd.mecat.miniasm.fasta",
        "assembly/real_reads_pb.yacrd2.canu.miniasm.fasta",
        "assembly/real_reads_ont.yacrd2.canu.miniasm.fasta",
        "assembly/real_reads_pb.yacrd2.consent.miniasm.fasta",
        "assembly/real_reads_ont.yacrd2.consent.miniasm.fasta",
        "assembly/real_reads_pb.yacrd2.mecat.miniasm.fasta",
        "assembly/real_reads_ont.yacrd2.mecat.miniasm.fasta",
        ## dascrubber
	# "assembly/real_reads_pb.dascrubber.raw.miniasm.fasta",
        # "assembly/real_reads_ont.dascrubber.raw.miniasm.fasta",
        "assembly/real_reads_pb.dascrubber.canu.miniasm.fasta",
        "assembly/real_reads_ont.dascrubber.canu.miniasm.fasta",
        "assembly/real_reads_pb.dascrubber.consent.miniasm.fasta",
        "assembly/real_reads_ont.dascrubber.consent.miniasm.fasta",
        "assembly/real_reads_pb.dascrubber.mecat.miniasm.fasta",
        # "assembly/real_reads_ont.dascrubber.mecat.miniasm.fasta", # correction not work
        ## miniscrub
	# "assembly/real_reads_pb.miniscrub.raw.miniasm.fasta",
        # "assembly/real_reads_ont.miniscrub.raw.miniasm.fasta",
        "assembly/real_reads_pb.miniscrub.canu.miniasm.fasta",
        "assembly/real_reads_ont.miniscrub.canu.miniasm.fasta",
        "assembly/real_reads_pb.miniscrub.consent.miniasm.fasta",
        "assembly/real_reads_ont.miniscrub.consent.miniasm.fasta",
        "assembly/real_reads_pb.miniscrub.mecat.miniasm.fasta",
        "assembly/real_reads_ont.miniscrub.mecat.miniasm.fasta",
        
        # canu
	## raw
	"assembly/real_reads_pb.raw.raw.canu.fasta",
        "assembly/real_reads_ont.raw.raw.canu.fasta",
        "assembly/real_reads_pb.raw.canu.canu.fasta",
        "assembly/real_reads_ont.raw.canu.canu.fasta",
        "assembly/real_reads_pb.raw.consent.canu.fasta",
        "assembly/real_reads_ont.raw.consent.canu.fasta",
        "assembly/real_reads_pb.raw.mecat.canu.fasta",
        "assembly/real_reads_ont.raw.mecat.canu.fasta",
        ## yacrd
	# "assembly/real_reads_pb.yacrd.raw.canu.fasta",
        # "assembly/real_reads_ont.yacrd.raw.canu.fasta",
        "assembly/real_reads_pb.yacrd.canu.canu.fasta",
        "assembly/real_reads_ont.yacrd.canu.canu.fasta",
        "assembly/real_reads_pb.yacrd.consent.canu.fasta",
        "assembly/real_reads_ont.yacrd.consent.canu.fasta",
        "assembly/real_reads_pb.yacrd.mecat.canu.fasta",
        "assembly/real_reads_ont.yacrd.mecat.canu.fasta",
        "assembly/real_reads_pb.yacrd2.canu.canu.fasta",
        "assembly/real_reads_ont.yacrd2.canu.canu.fasta",
        "assembly/real_reads_pb.yacrd2.consent.canu.fasta",
        "assembly/real_reads_ont.yacrd2.consent.canu.fasta",
        "assembly/real_reads_pb.yacrd2.mecat.canu.fasta",
        "assembly/real_reads_ont.yacrd2.mecat.canu.fasta",        ## dascrubber
	# "assembly/real_reads_pb.dascrubber.raw.canu.fasta",
        # "assembly/real_reads_ont.dascrubber.raw.canu.fasta",
        "assembly/real_reads_pb.dascrubber.canu.canu.fasta",
        "assembly/real_reads_ont.dascrubber.canu.canu.fasta",
        "assembly/real_reads_pb.dascrubber.consent.canu.fasta",
        "assembly/real_reads_ont.dascrubber.consent.canu.fasta",
        "assembly/real_reads_pb.dascrubber.mecat.canu.fasta",
        # "assembly/real_reads_ont.dascrubber.mecat.canu.fasta", # correction not work
        ## miniscrub
	# "assembly/real_reads_pb.miniscrub.raw.canu.fasta",
        # "assembly/real_reads_ont.miniscrub.raw.canu.fasta",
        "assembly/real_reads_pb.miniscrub.canu.canu.fasta",
        "assembly/real_reads_ont.miniscrub.canu.canu.fasta",
        "assembly/real_reads_pb.miniscrub.consent.canu.fasta",
        "assembly/real_reads_ont.miniscrub.consent.canu.fasta",
        "assembly/real_reads_pb.miniscrub.mecat.canu.fasta",
        "assembly/real_reads_ont.miniscrub.mecat.canu.fasta",

rule canu_pb:
    input:
        "correction/{prefix}_pb.{scrubber}.{corrector}.fasta"

    output:
        graph = "assembly/{prefix}_pb.{scrubber}.{corrector}.canu.gfa",
        contigs = "assembly/{prefix}_pb.{scrubber}.{corrector}.canu.fasta"
        
    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.{corrector}.canu.txt",
        
    shell:
        " && ".join([
            "rm -rf canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber}_{wildcards.corrector}"
            "canu -p canu -d canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber}_{wildcards.corrector} genomeSize=5.2m -pacbio-corrected {input} executiveThreads=8",
            "cp canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber}_{wildcards.corrector}/canu.contigs.gfa {output.graph}",
            "cp canu_asm/{wildcards.prefix}_pb_{wildcards.scrubber}_{wildcards.corrector}/canu.contigs.fasta {output.contigs}"
        ])

rule canu_ont:
    input:
        "correction/{prefix}_ont.{scrubber}.{corrector}.fasta"

    output:
        graph = "assembly/{prefix}_ont.{scrubber}.{corrector}.canu.gfa",
        contigs = "assembly/{prefix}_ont.{scrubber}.{corrector}.canu.fasta"
        
    benchmark:
        "benchmarks/{prefix}_ont.{scrubber}.{corrector}.canu.txt",
        
    shell:
        " && ".join([
            "rm canu_asm/{wildcards.prefix}_ont_{wildcards.scrubber}_{wildcards.corrector}"
            "canu -p canu -d canu_asm/{wildcards.prefix}_ont_{wildcards.scrubber}_{wildcards.corrector} genomeSize=5.2m -nanopore-corrected {input} executiveThreads=8",
            "cp canu_asm/{wildcards.prefix}_ont_{wildcards.scrubber}_{wildcards.corrector}/canu.contigs.gfa {output.graph}",
            "cp canu_asm/{wildcards.prefix}_ont_{wildcards.scrubber}_{wildcards.corrector}/canu.contigs.fasta {output.contigs}"
        ])
        
rule miniasm_pb:
    input:
        "correction/{prefix}_pb.{scrubber}.{corrector}.fasta"

    output:
        "assembly/{prefix}_pb.{scrubber}.{corrector}.miniasm.gfa",
        "assembly/{prefix}_pb.{scrubber}.{corrector}.miniasm.fasta"
        
    benchmark:
        "benchmarks/{prefix}_pb.{scrubber}.{corrector}.miniasm.txt",
        
    shell:
        " && ".join([
            "minimap2 -t 8 -x ava-pb {input} {input} > assembly/{wildcards.prefix}_pb.{wildcards.scrubber}.{wildcards.corrector}.miniasm.paf",
            "miniasm -f {input} assembly/{wildcards.prefix}_pb.{wildcards.scrubber}.{wildcards.corrector}.miniasm.paf > assembly/{wildcards.prefix}_pb.{wildcards.scrubber}.{wildcards.corrector}.miniasm.gfa",
            "./script/gfaminiasm2fasta.py assembly/{wildcards.prefix}_pb.{wildcards.scrubber}.{wildcards.corrector}.miniasm.gfa assembly/{wildcards.prefix}_pb.{wildcards.scrubber}.{wildcards.corrector}.miniasm.fasta"
        ])

rule miniasm_ont:
    input:
        "correction/{prefix}_ont.{scrubber}.{corrector}.fasta"

    output:
        "assembly/{prefix}_ont.{scrubber}.{corrector}.miniasm.gfa",
        "assembly/{prefix}_ont.{scrubber}.{corrector}.miniasm.fasta"
        
    benchmark:
        "benchmarks/{prefix}_ont.{scrubber}.{corrector}.miniasm.txt",
        
    shell:
        " && ".join([
            "minimap2 -t 8 -x ava-ont {input} {input} > assembly/{wildcards.prefix}_ont.{wildcards.scrubber}.{wildcards.corrector}.miniasm.paf",
            "miniasm -f {input} assembly/{wildcards.prefix}_ont.{wildcards.scrubber}.{wildcards.corrector}.miniasm.paf > assembly/{wildcards.prefix}_ont.{wildcards.scrubber}.{wildcards.corrector}.miniasm.gfa",
            "./script/gfaminiasm2fasta.py assembly/{wildcards.prefix}_ont.{wildcards.scrubber}.{wildcards.corrector}.miniasm.gfa assembly/{wildcards.prefix}_ont.{wildcards.scrubber}.{wildcards.corrector}.miniasm.fasta"
        ])
