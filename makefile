
# Make configuration 
.PHONY: false_negative synthetic_dataset real_dataset clean
.PRECIOUS: %.yacrd %.paf %.fasta %.fastq # avoid intermediate remove

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

real_dataset: data/LC.ont.assembly.fasta data/LC.ont.yacrd.assembly.fasta data/LC.ont.dascrubber.assembly.fasta data/LC.ont.miniscrub.assembly.fastq


# Assembly rules
%.assembly.fasta: %.fasta
	@echo "assembly "$^
	minimap2 -t 8 -x ava-ont $< $< > $^.assembly.paf
	miniasm -f $< $^.assembly.paf > $^.assembly.gfa

%.assembly.fasta: %.assembly.gfa
	./script/gfaminiasm2fasta.py $< $@

# Download data
data/LC.ont.fastq.gz:
	curl https://portal.nersc.gov/archive/home/r/regan/www/X0124/nanopore03_jgi_psf_org_20170608_FNFAH07719_MN17641_sequencing_run_170608_1002000055_001-combined.pass-1D.fastq.gz > data/LC.ont.fastq.gz

# Generic rules
%.fasta: %.fastq.gz
	gunzip -k $<
	sed -n '1~4s/^@/>/p;2~4p' $*.fastq > $@

%.fastq: %.fasta
	seqtk seq -F I $< > $@

%.ont.paf: %.ont.fasta
	minimap2 -t 8 -x ava-ont $< $< > $@

%.pb.paf: %.pb.fasta
	minimap2 -t 8 -x ava-pb $< $< > $@

%.yacrd %.yacrd.fasta: %.paf
	yacrd -i $< -o $*.yacrd -s $*.fasta --extracted-suffix .yacrd

%.dascrubber.fasta: %.fasta
	dascrubber_wrapper.py -i $< -g 5.5M > $@

%.miniscrub.fastq: %.fastq
	run_miniscrub.sh --verbose --output $@ $<

clean:
	rm -rf generate/*
	rm -rf data/*.paf
	rm -rf data/*.yacrd
	rm -rf data/LC.fasta
	rm -rf data/LC.fastq
	rm -rf data/t_roseus.yacrd.fasta
	rm -rf data/t_roseus.dascrubber.fasta
	rm -rf dascrubber_*

# Other rules 
generate/synthetic_chimera.pb.fasta:
	./script/generation_merge.py -i data/t_roseus.pb.fasta -o generate/synthetic_chimera.pb.fasta -n 2535

