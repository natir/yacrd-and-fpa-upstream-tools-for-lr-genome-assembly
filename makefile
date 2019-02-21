
# Make configuration 
.PHONY: false_negative synthetic_dataset real_dataset clean
.PRECIOUS: %.yacrd %.paf %.fasta %.fastq %.gfa # avoid intermediate remove

# Major rules
false_negative: data/t_roseus.pb.fasta data/t_roseus.pb.yacrd data/t_roseus.pb.dascrubber.fasta data/t_roseus.pb.miniscrub.fastq 
	@echo "yacrd"
	./script/stats_yacrd.py -c data/t_roseus.pb.fasta -y data/t_roseus.pb.yacrd
	@echo "dascrubber"
	./script/stats_dascrubber.py -c data/t_roseus.pb.fasta -d data/t_roseus.pb.dascrubber.fasta

synthetic_dataset: generate/synthetic_chimera.pb.fasta generate/synthetic_chimera.pb.yacrd generate/synthetic_chimera.pb.dascrubber.fasta generate/synthetic_chimera.pb.miniscrub.fastq
	@echo "yacrd"
	./script/stats_yacrd.py -c generate/synthetic_chimera.pb.fasta -y generate/synthetic_chimera.pb.yacrd
	@echo "dascrubber"
	./script/stats_dascrubber.py -c generate/synthetic_chimera.pb.fasta -d generate/synthetic_chimera.pb.dascrubber.fasta

real_dataset: data/LC.ont.assembly.canu.fasta data/LC.yacrd.ont.assembly.canu.fasta data/LC.ont.dascrubber.assembly.canu.fasta data/LC.ont.miniscrub.assembly.canu.fasta data/LC.ont.assembly.miniasm.fasta data/LC.yacrd.ont.assembly.miniasm.fasta data/LC.ont.dascrubber.assembly.miniasm.fasta data/LC.ont.miniscrub.assembly.miniasm.fasta

# Assembly rules

## Canu
%.assembly.canu.gfa %.assembly.canu.fasta: %.fasta
	@echo "assembly canu "$*
	canu -p canu -nanopore-raw $< -d canu_$* genomeSize=4.6M #correctedErrorRate=0.12
	cp canu_$*/canu.contigs.fasta ./$*.assembly.canu.fasta
	cp canu_$*/canu.contigs.gfa ./$*.assembly.canu.gfa

## Miniasm 
%.assembly.miniasm.gfa: %.fasta
	@echo "assembly miniasm "$*
	minimap2 -t 8 -x ava-ont $< $< > $*.assembly.miniasm.paf
	miniasm -f $< $*.assembly.miniasm.paf > $*.assembly.miniasm.gfa

%.assembly.miniasm.gfa: %.fastq
	@echo "assembly miniasm "$*
	minimap2 -t 8 -x ava-ont $< $< > $*.assembly.miniasm.paf
	miniasm -f $< $*.assembly.miniasm.paf > $*.assembly.miniasm.gfa

%.assembly.miniasm.fasta: %.assembly.miniasm.gfa
	./script/gfaminiasm2fasta.py $< $@

# Download data
data/LC.ont.fastq:
	curl https://portal.nersc.gov/archive/home/r/regan/www/X0124/nanopore03_jgi_psf_org_20170608_FNFAH07719_MN17641_sequencing_run_170608_1002000055_001-combined.pass-1D.fastq.gz | seqtk sample -s 42 /dev/stdin 0.2 > data/LC.ont.fastq

# Generic rules
%.fasta: %.fastq
	sed -n '1~4s/^@/>/p;2~4p' $*.fastq > $@

#%.fastq: %.fasta
#	seqtk seq -F I $< > $@

%.ont.paf: %.ont.fasta
	minimap2 -t 8 -x ava-ont $< $< > $@

%.pb.paf: %.pb.fasta
	minimap2 -t 8 -x ava-pb $< $< > $@

%.ont.yacrd %.yacrd.ont.fasta: %.ont.fasta
	minimap2 -t 8 -x ava-ont $< $< | fpa -l 500 -i > $*.ont.paf
	yacrd -c 1 -i $*.ont.paf -o $*.ont.yacrd -s $*.ont.fasta --splited-suffix .yacrd

%.pb.yacrd %.yacrd.pb.fasta: %.pb.fasta
	minimap2 -t 8 -x ava-pb $< $< | fpa -l 500 -i > $*.pb.paf
	yacrd -c 1 -i $*.pb.paf -o $*.pb.yacrd -s $*.pb.fasta --splited-suffix .yacrd

%.dascrubber.fasta: %.fasta
	dascrubber_wrapper.py -i $< -g 5.5M > $@

%.miniscrub.fastq: %.fastq
	run_miniscrub.sh --verbose --output $@ $<

clean:
	rm -rf generate/*
	rm -rf data/*.paf
	rm -rf data/*.yacrd
	rm -rf data/LC.ont.fasta
	rm -rf data/LC.ont.*assembly*
	rm -rf data/t_roseus.pb.fastq
	rm -rf data/t_roseus.yacrd.pb.fasta
	rm -rf dascrubber_*

# Other rules 
generate/synthetic_chimera.pb.fasta:
	./script/generation_merge.py -i data/t_roseus.pb.fasta -o generate/synthetic_chimera.pb.fasta -n 2535

