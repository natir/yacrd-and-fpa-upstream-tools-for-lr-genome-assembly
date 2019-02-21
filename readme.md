# Dataset

- Reference [E. cole CFT073](https://www.uniprot.org/taxonomy/199310)
- Reads:
  * [Pacbio](https://www.ebi.ac.uk/ena/data/view/SRX5299472)
  * [Oxford nanopore](https://www.ebi.ac.uk/ena/data/view/SRR8494940)

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

## Scrubbing

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

## Correction

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

## Correction

|                   | raw  | yacrd | dascrubber | miniscrub |
| ----------------  | ---- | ----- | ---------- | --------- |
| # of contig       |      |       |            |           |
| # of base         |      |       |            |           | 
| % assembly mapped |      |       |            |           |
| # mismatch        |      |       |            |           |
| time              |      |       |            |           |
| memory            |      |       |            |           |

# How to run

1. Install all tools in your PATH, and python3

2. Install python dependency `pip install -r requirements.txt`

3. Run:
   - download data `snakemake --snakefile pipeline/download.snakefile all`
   - run scrubbing `snakemake --snakefile pipeline/scrubbing.snakefile all`
   - run correction `snakemake --snakefile pipeline/correction.snakefile all`
   - run assembly `snakemake --snakefile pipeline/assembly.snakefile all`

4. Analysis:
   TODO
