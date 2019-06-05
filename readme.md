# Requirements

- [seqtk](https://github.com/lh3/seqtk) 1.3-r106
- [fpa](https://gitlab.inria.fr/pmarijon/fpa) 0.5
- [yacrd](https://gitlab.inria.fr/pmarijon/yacrd) 0.5.1
- [dascrubber](https://github.com/thegenemyers/DASCRUBBER/) commit 0e90524 you can follow [dascrubber-wrapper](https://github.com/rrwick/DASCRUBBER-wrapper) instruction to install all dascrubber requirements
- [miniscrub](https://bitbucket.org/berkeleylab/jgi-miniscrub) commit 3d11d3e
- [snakemake](https://snakemake.readthedocs.io/en/stable/) 5.4.3
- [wtdbg2](https://github.com/ruanjue/wtdbg2) commit 8908a31
- [miniasm](https://github.com/lh3/miniasm) 3d11d3e
- [quast](http://bioinf.spbau.ru/quast) v5.0.2

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
all_real_reads_pb.fasta          -> E. coli Pacbio reads
c_elegans_pb.fasta               -> C. elegans reads
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
seqtk -s 42 data/all_real_reads_pb.fasta 0.1838 > data/real_reads_pb.fasta
```

# Rerun analysis

Update miniscrub path in file `pipeline/scrubbing.snakefile` line 136 for miniscrub in gpu mode and line 154 for miniscrub in cpu mode

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

Get information one read mapped on reference:
```
./script/get_mapping_info.py mapping/*.bam
```

Get information about assembly in csv format:
```
./script/simplify_quast_result.py quast/*/report.tsv       # assembly stats of scrubbed read
./script/simplify_quast_result.py fpa/quast/*/report.tsv   # assembly stats of after fpa
./script/simplify_quast_result.py combo/quast/*/report.tsv # assembly stats of modified miniasm pipeline assembly with yacrd and fpa
```
