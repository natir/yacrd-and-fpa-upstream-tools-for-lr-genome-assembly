# How to test yacrd and fpa

If you want run yacrd and fpa speedly you can run :

```
./script/small_test.sh
```

This script download the *E. coli* Nanopore dataset run a subsampling on it and run yacrd, fpa and a combination of this tools on this dataset.

# Requirements

This tools need to be avaible in your path :

- [seqtk](https://github.com/lh3/seqtk) 1.3-r106
- [fpa](https://gitlab.inria.fr/pmarijon/fpa) 0.5
- [yacrd](https://gitlab.inria.fr/pmarijon/yacrd) 0.5.1
- [dascrubber](https://github.com/thegenemyers/DASCRUBBER/) commit 0e90524 you can follow [dascrubber-wrapper](https://github.com/rrwick/DASCRUBBER-wrapper) instruction to install all dascrubber requirements
- [snakemake](https://snakemake.readthedocs.io/en/stable/) 5.4.3
- [wtdbg2](https://github.com/ruanjue/wtdbg2) 2.3
- [miniasm](https://github.com/lh3/miniasm) 0.3-r179
- [quast](http://bioinf.spbau.ru/quast) v5.0.2
- [nucmer](http://mummer.sourceforge.net/) 4.0.0beta2
- [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) 0.7.17
- [samtools](https://samtools.github.io/) 1.9

You need change path of this tools in snakemake pipeline file:

- [miniscrub](https://bitbucket.org/berkeleylab/jgi-miniscrub) commit 3d11d3e
- [ra](https://github.com/lbcb-sci/ra) commit 07364a1 
- [porechop](https://github.com/rrwick/Porechop/) v0.2.3-C++11
- [shasta](https://github.com/chanzuckerberg/shasta/) 0.1.0

Update miniscrub path in file `pipeline/scrubbing.snakefile` line 136.
Update ra path in file `pipeline/assembly.snakefile` line 69.
Update shasta path in `pipeline/assembly.snakefile` line 85.
Update porechop path in `pipeline/analysis.snakefile` line 68.

If you execute `conda env create -f conda_env.yml` conda create environment `yacrd_fpa` with all dependency except dascrubber, miniscrub, ra, porechop and shasta.

# Dataset

- Reference:
  * [E. coli CFT073](https://www.uniprot.org/taxonomy/199310)  5.231428 Mb
  * [D. melanogaster](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4) 143.726002 Mb
  * [C. elegans](ftp://ftp.ensembl.org/pub/release-95/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz) 100.2 Mb
  * [H. sapiens chr1](ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz) 248.9 Mb
  
- Reads:
  * E. coli:
	+ [Pacbio](https://www.ebi.ac.uk/ena/data/view/SRX5299472)
	+ [Oxford nanopore](https://www.ebi.ac.uk/ena/data/view/SRR8494940)
  * [Oxford nanopore D melanogaster](https://www.ebi.ac.uk/ena/data/view/SRX3676783)
  * [Oxford nanopore H sapiens chr1](http://s3.amazonaws.com/nanopore-human-wgs/chr1.sorted.bam)
  * [Pacbio RS P6-C4 C elegans](http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/list.html)

## Build data directory

Download all data and build a data script like this:
```
all_real_reads_ont.fasta         -> E. coli Nanopore reads 
all_real_reads_ont.fastq         -> fastq of previous reads
all_real_reads_pb.fasta          -> E. coli Pacbio reads
all_real_reads_pb.fastq          -> fastq of previous reads
c_elegans_pb.fasta               -> C. elegans pacbio reads
c_elegans_ont.fasta               -> C. elegans nanopore reads
d_melanogaster_reads_ont.fasta   -> D. melanogaster reads
h_sapiens_chr1_ont.fasta         -> H. Sapiens maps against chromosomes 1 reads 

c_elegans_ref.fasta              -> C. elegans reference
d_melanogaster_ref.fasta         -> D. melanogaster reference
h_sapiens_chr1_ref.fasta         -> H. sapiens reference
ref_e_coli_cft073.fasta          -> E. coli reference
```

Run :
```
seqtk -s 42 data/all_real_reads_ont.fasta 0.1618 > data/real_reads_ont.fasta
seqtk -s 42 data/all_real_reads_ont.fastq 0.1618 > data/real_reads_ont.fastq
seqtk -s 42 data/all_real_reads_pb.fasta 0.1838 > data/real_reads_pb.fasta
seqtk -s 42 data/all_real_reads_pb.fastq 0.1838 > data/real_reads_pb.fastq
```

# Rerun analysis

- Run scrubbing+assembly+analysis `snakemake --snakefile pipeline/uncorrected.snakefile all`
- Run fpa+assembly+analysis `snakemake --snakefile pipeline/fpa.snakefile all`
- Run comparaison against minimap+miniasm and yacrd+minimap+fpa+miniasm pipeline `snakemake --snakefile pipeline/combo.snakefile all`

## Analysis script

Get number of chimera per file in csv format:
```
for i in \$(ls mapping/*.paf)
do
./script/found_chimera.py \$i
done
```

Get information one read mapped on reference (this script required pysam): 
```
./script/get_mapping_info.py mapping/*.bam
```

Get information about assembly in csv format:
```
./script/simplify_quast_result.py quast/*/report.tsv       # assembly stats of scrubbed read
./script/simplify_quast_result.py fpa/quast/*/report.tsv   # assembly stats of after fpa
./script/simplify_quast_result.py combo/quast/*/report.tsv # assembly stats of modified miniasm pipeline assembly with yacrd and fpa
```
