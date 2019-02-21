# Dataset

- Reference [Escherichia coli str. K-12 substr. MDS42](https://www.ncbi.nlm.nih.gov/nuccore/AP012306.1?report=fasta)
- Reads [Nanopore 1D](https://portal.nersc.gov/archive/home/r/regan/www/X0124/nanopore03_jgi_psf_org_20170608_FNFAH07719_MN17641_sequencing_run_170608_1002000055_001-combined.pass-1D.fastq.gz)

# Tools

scrubber:
- [yacrd](https://gitlab.inria.fr/pmarijon/yacrd)
- [dascrubber](https://github.com/rrwick/DASCRUBBER-wrapper)
- [miniscrub](https://bitbucket.org/berkeleylab/jgi-miniscrub)

corrector:
- canu 1.8
- [MECAT](https://github.com/xiaochuanle/MECAT)
- [CONSENT](https://github.com/morispi/CONSENT) 1.1.1

assembly tools:
- canu 1.8
- miniasm 0.3-r179

# Analysis

## direct scrubber output

- venn diagrame of modified read
- \# of read
- \# of base
- \# reads length histogram
- % base removed
- \# mapped read ?
- \# mismatch ?
- timing

## compare correction output

- \# of read
- \# of base
- \# reads length histogram
- \# mapped read ?
- \# mismatch ?
- timing
- other stat provides in corrector publication

Anaylsis if time computation reduction is correlate to read reduction

## compare contig output

- quast info
- other stat provides in assembly publication

# Result

|                   | raw  | yacrd | dascrubber | miniscrub |
| ----------------  | ---- | ----- | ---------- | --------- |
| # of read         |      |       |            |           |
| # of base         |      |       |            |           | 
| % of base removed |      |       |            |           | 
| % mapped read     |      |       |            |           |
| % mapped base     |      |       |            |           |
| # mismatch        |      |       |            |           |
| time              |      |       |            |           |
| memory            |      |       |            |           |

# How to run

1. Install all tools in your PATH

2. Run:
   - download data `snakemake --snakefiles pipeline/download.snakefile all`
   - run scrubbing `snakemake --snakefiles pipeline/scrubbing.snakefile all`
   - run correction `snakemake --snakefiles pipeline/correction.snakefile all`
   - run assembly `snakemake --snakefiles pipeline/assembly.snakefile all`

3. Analysis:
   TODO
