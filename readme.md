# Dataset

- Reference [E. coli CFT073](https://www.uniprot.org/taxonomy/199310)  5.231428 Mb
- Reads:
  * [Pacbio](https://www.ebi.ac.uk/ena/data/view/SRX5299472)
  * [Oxford nanopore](https://www.ebi.ac.uk/ena/data/view/SRR8494940)


|                   | Pacbio         | Nanopore       | Pacbio subsample     | Nanopore subsample  |
| ----------------  | --------------:| --------------:| --------------------:| -------------------:|
| # of sequences    | 207069         | 158590         | 37404         (18 %) | 25469        (16 %) |
| Total length      | 1425.446392 Mb | 1621.000527 Mb | 257.882884 Mb (18 %) | 257508441 Mb (16 %) |
| Longest sequence  | 41.631 Kb      | 164.088 Kb     | 38.331 Kb            | 137.142 Kb          |
| Shortest sequence | 35 b           | 88 b           | 35 b                 | 152 b               |
| Mean sequences    | 6.883 Kb       | 10.221 Kb      | 6.894 Kb             | 10.11 Kb            |
| Median Length     | 6.679 Kb       | 5.591 Kb       | 6.672 Kb             | 5.515 Kb            |
| N10               | 7467           | 2531           | 1354                 | 400                 |
| N50               | 58081          | 23903          | 10502                | 3807                |
| N90               | 142591         | 86781          | 25787                | 13969               |
| L10               | 15.631 Kb      | 50.77 Kb       | 15.591 Kb            | 51.316 Kb           |
| L50               | 9.052 Kb       | 20.189 Kb      | 9.064 Kb             | 20.073 Kb           |
| L90               | 4.191 Kb       | 4.763 Kb       | 4.218 Kb             | 4.701 Kb            |
| Coverage          | 272x           | 309x           | 49x                  | 49x                 |



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

### Pacbio

|                        | raw        | yacrd      | dascrubber | miniscrub  |
|:---------------------- | ----------:| ----------:| ----------:| ----------:|
| # of read              | 37404      | 38890      | 33926      | 59174      |
| # of base (Mb)         | 257.882884 | 224.056471 | 177.025153 | 79.459162  |
| % of base removed (Mb) |            | 33.826413  | 80.857731  | 178.423722 |
| coverage               |            | 43x        | 34x        | 15x        |
| % mapped read          |            |            |            |            |
| % mapped base          |            |            |            |            |
| # mismatch             |            |            |            |            |
| time (s)               |            | 21.2511    | 381.8210   | 19995.2767 |
| memory (Mo)            |            | 2258.18    | 13832.30   | 17234.49   |

### Nanopore

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
