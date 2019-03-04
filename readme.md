# Dataset

- Reference [E. coli CFT073](https://www.uniprot.org/taxonomy/199310)  5.231428 Mb
- Reads:
  * [Pacbio](https://www.ebi.ac.uk/ena/data/view/SRX5299472)
  * [Oxford nanopore](https://www.ebi.ac.uk/ena/data/view/SRR8494940)


|                   | Pacbio         | Nanopore       | Pacbio subsample     | Nanopore subsample  |
| ----------------  | --------------:| --------------:| --------------------:| -------------------:|
| # of sequences    | 207069         | 158590         | 37404         (18 %) | 25469        (16 %) |
| Total length      | 1425.446392 Mb | 1621.000527 Mb | 257.882884 Mb (18 %) | 257.508441 Mb (16 %) |
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

|                | pb.raw | pb.yacrd | pb.yacrd2 | pb.dascrubber | pb.miniscrub |
| - | -:|  -:|  -:|  -:|  -:| 
| \# of read     | 37404 | 38890 | 38068 | 33926 | 59174 |
| \# of base     | 257882884 | 224056471 | 209162570 | 177025153 | 79459162 |
| N10            | 1354 | 1465 | 1371 | 1287 | 1349 |
| L10            | 15590 | 12712 | 12630 | 11714 | 4625 |
| N50            | 10502 | 10727 | 10077 | 9092 | 13609 |
| L50            | 10502 | 10727 | 10077 | 9092 | 13609 |
| N90            | 25787 | 26767 | 25133 | 23020 | 44810 |
| L90            | 4218 | 3333 | 3311 | 2900 | 652 |
| % base removed | 0.00 | 13.12 | 18.89 | 31.35 | 69.19 |
| \# mapped read | 36206 | 37632 | 35465 | 33750 | 59171 |
| \# mismatch    | 29863768 | 29254044 | 27275079 | 21534830 | 7642568 |
| time           | 0 | 21.2511 | 24.5972 | 381.8210 | 19995.2767 |
| memory         | 0 | 2258.18 | 2276.36 | 13832.30 | 17234.49 |

### Nanopore

|                | ont.raw | ont.yacrd | ont.yacrd2 | ont.dascrubber | ont.miniscrub |
| - | -:|  -:|  -:|  -:|  -:| 
| \# of read     | 25469 | 25660 | 25590 | 22538 | 62915 |
| \# of base     | 257508441 | 252085244 | 250396634 | 236474486 | 180605742 |
| N10            | 400 | 396 | 393 | 366 | 772 |
| L10            | 51296 | 50779 | 50835 | 51571 | 17727 |
| N50            | 3807 | 3765 | 3743 | 3469 | 8711 |
| L50            | 3807 | 3765 | 3743 | 3469 | 8711 |
| N90            | 13969 | 13868 | 13774 | 12588 | 38354 |
| L90            | 4700 | 4585 | 4592 | 4806 | 1093 |
| % base removed | 0.00 | 2.11 | 2.76 | 8.17 | 29.86 |
| \# mapped read | 25393 | 25575 | 25418 | 22536 | 62905 |
| \# mismatch    | 37383344 | 37109735 | 36828994 | 29638330 | 21041230 |
| time           | 0 | 28.4964 | 30.9635 | 526.5879 | 21352.3889 |
| memory         | 0 | 2103.60 | 2104.96 | 27901.59 | 18240.23 |

## Correction

### Pacbio

|                | pb.raw.raw | pb.raw.canu | pb.raw.consent | pb.raw.mecat | pb.yacrd.canu | pb.yacrd.consent | pb.yacrd.mecat | pb.yacrd2.canu | pb.yacrd2.consent | pb.yacrd2.mecat | pb.dascrubber.canu | pb.dascrubber.consent | pb.dascrubber.mecat | pb.miniscrub.canu | pb.miniscrub.consent | pb.miniscrub.mecat |
| - | -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:| 
| \# of read     | 37404 | 24439 | 33642 | 10905 | 25569 | 35021 | 12150 | 24657 | 32818 | 11817 | 28518 | 32775 | 14102 | 26254 | 49976 | 0 |
| \# of base     | 257882884 | 147344360 | 188870837 | 81199441 | 149811236 | 183105081 | 90110565 | 145211033 | 174306872 | 87363457 | 161809806 | 168948669 | 105499016 | 54369709 | 69959182 | 0 |
| N10            | 1354 | 1088 | 1331 | 662 | 1113 | 1357 | 728 | 1078 | 1291 | 708 | 1204 | 1258 | 836 | 868 | 1190 | 0 |
| L10            | 15590 | 11636 | 12012 | 10501 | 11537 | 11528 | 10582 | 11543 | 11535 | 10569 | 11486 | 11450 | 10778 | 5097 | 4655 | 0 |
| N50            | 10502 | 7562 | 9521 | 4464 | 7768 | 9590 | 4951 | 7522 | 9102 | 4816 | 8461 | 8903 | 5713 | 7677 | 11718 | 0 |
| L50            | 10502 | 7562 | 9521 | 4464 | 7768 | 9590 | 4951 | 7522 | 9102 | 4816 | 8461 | 8903 | 5713 | 7677 | 11718 | 0 |
| N90            | 25787 | 18148 | 23550 | 9375 | 18858 | 24219 | 10450 | 18209 | 22786 | 10167 | 20829 | 22522 | 12118 | 21232 | 37850 | 0 |
| L90            | 4218 | 3538 | 3240 | 5625 | 3384 | 2915 | 5610 | 3410 | 3000 | 5586 | 3209 | 2855 | 5618 | 1179 | 674 | 0 |
| % base removed | 0.00 | 42.86 | 26.76 | 68.51 | 33.14 | 18.28 | 59.78 | 30.58 | 16.66 | 58.23 | 8.60 | 4.56 | 40.40 | 31.58 | 11.96 | 100.00 |
| \# mapped read | 36206 | 24437 | 33640 | 10905 | 25569 | 35017 | 12150 | 24657 | 32813 | 11817 | 28518 | 32775 | 14102 | 26252 | 49976 | 0 |
| \# mismatch    | 29863768 | 1103030 | 4084253 | 642480 | 1090079 | 3794849 | 696601 | 1047965 | 3592823 | 675251 | 924368 | 1886849 | 706539 | 322402 | 604671 | 0 |
| time           | 0 | 3103.2063 | 2118.4074 | 162.7338 | 2743.7773 | 2067.4703 | 189.3940 | 2707.6767 | 1984.7657 | 187.9458 | 2429.9100 | 1782.7734 | 151.4521 | 793.1057 | 661.6553 | 17.1732 |
| memory         | 0 | 3357.47 | 4279.85 | 246.45 | 4208.07 | 4054.53 | 1882.32 | 3320.95 | 3187.07 | 1808.07 | 4421.41 | 2511.33 | 1489.20 | 2891.51 | 1301.25 | 1165.45 |


### Nanopore

|                | ont.raw.raw | ont.raw.canu | ont.raw.consent | ont.raw.mecat | ont.yacrd.canu | ont.yacrd.consent | ont.yacrd.mecat | ont.yacrd2.canu | ont.yacrd2.consent | ont.yacrd2.mecat | ont.dascrubber.canu | ont.dascrubber.consent | ont.miniscrub.canu | ont.miniscrub.consent | ont.miniscrub.mecat |
| - | -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:|  -:| 
| \# of read     | 25469 | 24439 | 24208 | 19972 | 25569 | 24322 | 19835 | 24657 | 24208 | 19734 | 28518 | 22302 | 26254 | 57412 | 23761 |
| \# of base     | 257508441 | 147344360 | 252343673 | 239684772 | 149811236 | 251645063 | 238062600 | 145211033 | 250468181 | 236971782 | 161809806 | 237219041 | 54369709 | 175260949 | 138989242 |
| N10            | 400 | 1088 | 390 | 472 | 1113 | 390 | 469 | 1078 | 388 | 467 | 1204 | 365 | 868 | 740 | 556 |
| L10            | 51296 | 11636 | 51674 | 45872 | 11537 | 51349 | 45843 | 11543 | 51380 | 45832 | 11486 | 51885 | 5097 | 18078 | 19610 |
| N50            | 3807 | 7562 | 3698 | 3832 | 7768 | 3707 | 3805 | 7522 | 3694 | 3786 | 8461 | 3452 | 7677 | 8243 | 5527 |
| L50            | 3807 | 7562 | 3698 | 3832 | 7768 | 3707 | 3805 | 7522 | 3694 | 3786 | 8461 | 3452 | 7677 | 8243 | 5527 |
| N90            | 13969 | 18148 | 13464 | 13032 | 18858 | 13511 | 12937 | 18209 | 13461 | 12870 | 20829 | 12506 | 21232 | 35375 | 17878 |
| L90            | 4700 | 3538 | 4821 | 5301 | 3384 | 4783 | 5290 | 3410 | 4776 | 5292 | 3209 | 4894 | 1179 | 1183 | 2800 |
| % base removed | 0.00 | 42.78 | 2.01 | 6.92 | 40.57 | 0.17 | 5.56 | 42.01 | -0.03 | 5.36 | 31.57 | -0.31 | 69.90 | 2.96 | 23.04 |
| \# mapped read | 25393 | 24437 | 24205 | 19972 | 25569 | 24319 | 19835 | 24657 | 24206 | 19734 | 28518 | 22302 | 26252 | 57412 | 23761 |
| \# mismatch    | 37383344 | 1102903 | 4334858 | 7278735 | 1089999 | 4324461 | 7133121 | 1047670 | 4295612 | 7097879 | 924362 | 3063901 | 322402 | 1905281 | 2604923 |
| time           | 0 | 13676.9856 | 3970.4994 | 658.4983 | 4.0438 | 3935.5846 | 639.6991 | 5.3595 | 4169.6171 | 667.4647 | 4.3111 | 3609.1146 | 1.9398 | 2265.2232 | 234.6207 |
| memory         | 0 | 4428.85 | 3693.43 | 2045.31 | 44.34 | 3663.08 | 2022.72 | 29.71 | 3440.56 | 2012.88 | 44.27 | 3608.59 | 21.70 | 2655.66 | 1665.74 |

# How to run

1. Install all tools in your PATH, and python3

2. Install python dependency `pip install -r requirements.txt`

3. Run:
   - download data `snakemake --snakefile pipeline/download.snakefile all`
   - run scrubbing `snakemake --snakefile pipeline/scrubbing.snakefile all`
   - run correction `snakemake --snakefile pipeline/correction.snakefile all`
   - run assembly `snakemake --snakefile pipeline/assembly.snakefile all`

4. Analysis:
   - map scrubbed 
   - map corrected
   - run quast on assembly
