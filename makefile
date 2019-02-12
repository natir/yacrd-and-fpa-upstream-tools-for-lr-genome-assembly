
# Make configuration 
.PHONY: false_negative synthetic_dataset real_dataset clean
.PRECIOUS: %.yacrd %.paf %.fasta # avoid intermediate remove

# Major rules
false_negative: data/t_roseus.fasta data/t_roseus.yacrd data/t_roseus.dascrubber.fasta 
	@echo "yacrd"
	./script/stats_yacrd.py -c data/t_roseus.fasta -y data/t_roseus.yacrd
	@echo "dascrubber"
	./script/stats_dascrubber.py -c data/t_roseus.fasta -d data/t_roseus.dascrubber.fasta

synthetic_dataset: generate/synthetic_chimera.fasta generate/synthetic_chimera.yacrd generate/synthetic_chimera.dascrubber.fasta
	@echo "yacrd"
	./script/stats_yacrd.py -c generate/synthetic_chimera.fasta -y generate/synthetic_chimera.yacrd
	@echo "dascrubber"
	./script/stats_dascrubber.py -c generate/synthetic_chimera.fasta -d generate/synthetic_chimera.dascrubber.fasta

real_dataset: data/LC.fastq.gz data/LC.yacrd.fasta data/LC.dascrubber.fasta
	@echo "assembly yacrd"
	minimap2 -x ava-pb data/LC.yacrd.fasta data/LC.yacrd.fasta > data/LC.yacrd.paf
	miniasm -f data/LC.yacrd.fasta > data/LC.yacrd.gfa

# Download data
LC.fastq.gz:
	curl https://portal.nersc.gov/archive/home/r/regan/www/X0124/nanopore03_jgi_psf_org_20170608_FNFAH07719_MN17641_sequencing_run_170608_1002000055_001-combined.pass-1D.fastq.gz > LC.fastq.gz

# Generic rules
%: %.gz
	gunzip $<

%.fasta: %.fastq
	sed -n '1~4s/^@/>/p;2~4p' $< > $@

%.paf: %.fasta
	minimap2 -t 4 -x ava-pb $< $< > $@

%.yacrd.fasta: %.paf
	yacrd -i $< -o /dev/null -s $*.fasta --extracted-suffix .yacrd

%.dascrubber.fasta: %.fastq.gz
	dascrubber_wrapper.py -i $< -g 5.5M > $@

%.dascrubber.fasta: %.fasta
	dascrubber_wrapper.py -i $< -g 5.5M > $@

%.miniscrub.fasta: %.fasta
	/home/pierre/dev/jgi-miniscrub/venv/bin/python /home/pierre/dev/jgi-miniscrub/miniscrub.py

clean:
	rm -rf generate/*
	rm -rf data/*.paf
	rm -rf data/*.yacrd


# Other rules 
generate/synthetic_chimera.fasta:
	./script/generation_merge.py -i data/t_roseus.fasta -o generate/synthetic_chimera.fasta -n 2535

