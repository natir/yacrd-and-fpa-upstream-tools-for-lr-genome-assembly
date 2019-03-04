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

## Read type 

### Pacbio

Count:

|           | yacrd | dascrubber | miniscrub |
| --------- | -----:| ----------:| ---------:|
| discard   |  1842 |      13400 |      5928 |
| splited   |  1628 |      15708 |      2208 |
| trimmed   | 28583 |       6743 |     29268 |
| nmodified |  5351 |       1553 |        73 |

Jacard distance discard:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.07 |            |           |
| miniscrub  |  0.09 |       0.40 |           |

Jacard distance splited:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.33 |            |           |
| miniscrub  |  0.05 |       0.05 |           |

Jacard distance trimmed:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.87 |            |           |
| miniscrub  |  0.20 |       0.20 |           |

Jacard distance not modified:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.00 |            |           |
| miniscrub  |  0.04 |       0.00 |           |

### Nanopore

Count:

|           | yacrd | dascrubber | miniscrub |
| --------- | -----:| ----------:| ---------:|
| discard   |   138 |       3293 |      3096 |
| splited   |   184 |      10611 |       158 |
| trimmed   | 22590 |       2892 |     22215 |
| nmodified |  2557 |       8673 |         0 |

Jacard distance discard:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.03 |            |           |
| miniscrub  |  0.03 |       0.46 |           |

Jacard distance splited:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.34 |            |           |
| miniscrub  |  0.01 |       0.01 |           |

Jacard distance trimmed:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.94 |            |           |
| miniscrub  |  0.12 |       0.12 |           |

Jacard distance not modified:

|            | yacrd | dascrubber | miniscrub |
| ---------- | -----:| ----------:| ---------:|
| yacrd      |       |            |           |
| dascrubber |  0.00 |            |           |
| miniscrub  |  0.12 |       0.00 |           |

## Scrubbing

### Pacbio

|                | pb.raw | pb.yacrd | pb.yacrd2 | pb.dascrubber | pb.miniscrub |
| - | -:| -:| -:| -:| -:|
| \# of read     | 37404 | 38890 | 38068 | 33926 | 59174 |
| \# of base     | 257882884 | 224056471 | 209162570 | 177025153 | 79459162 |
| coverage       | 49.29x | 42.83x | 39.98x | 33.84x | 15.19x |
| N10            | 1354 | 1465 | 1371 | 1287 | 1349 |
| L10            | 15590 | 12712 | 12630 | 11714 | 4625 |
| N50            | 10502 | 10727 | 10077 | 9092 | 13609 |
| L50            | 9064 | 7951 | 7905 | 7500 | 1672 |
| N90            | 25787 | 26767 | 25133 | 23020 | 44810 |
| L90            | 4218 | 3333 | 3311 | 2900 | 652 |
| % base removed | 0.00 | 13.12 | 18.89 | 31.35 | 69.19 |
| \# mapped read | 36206 | 37632 | 35465 | 33750 | 59171 |
| \# mismatch    | 29863768 | 29254044 | 27275079 | 21534830 | 7642568 |
| time           | 0 | 21.2511 | 24.5972 | 381.8210 | 19995.2767 |
| memory         | 0 | 2258.18 | 2276.36 | 13832.30 | 17234.49 |

### Nanopore

|                | ont.raw | ont.yacrd | ont.yacrd2 | ont.dascrubber | ont.miniscrub |
| - | -:| -:| -:| -:| -:|
| \# of read     | 25469 | 25660 | 25590 | 22538 | 62915 |
| \# of base     | 257508441 | 252085244 | 250396634 | 236474486 | 180605742 |
| coverage       | 49.22x | 48.19x | 47.86x | 45.20x | 34.52x |
| N10            | 400 | 396 | 393 | 366 | 772 |
| L10            | 51296 | 50779 | 50835 | 51571 | 17727 |
| N50            | 3807 | 3765 | 3743 | 3469 | 8711 |
| L50            | 20071 | 19856 | 19856 | 20333 | 5447 |
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
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of read     | 37404 | 24439 | 33642 | 10905 | 25569 | 35021 | 12150 | 24657 | 32818 | 11817 | 28518 | 32775 | 14102 | 26254 | 49976 | 0 |
| \# of base     | 257882884 | 147344360 | 188870837 | 81199441 | 149811236 | 183105081 | 90110565 | 145211033 | 174306872 | 87363457 | 161809806 | 168948669 | 105499016 | 54369709 | 69959182 | 0 |
| coverage       | 49.29x | 28.17x | 36.10x | 15.52x | 28.64x | 35.00x | 17.22x | 27.76x | 33.32x | 16.70x | 30.93x | 32.29x | 20.17x | 10.39x | 13.37x | 0.00x |
| N10            | 1354 | 1088 | 1331 | 662 | 1113 | 1357 | 728 | 1078 | 1291 | 708 | 1204 | 1258 | 836 | 868 | 1190 | 0 |
| L10            | 15590 | 11636 | 12012 | 10501 | 11537 | 11528 | 10582 | 11543 | 11535 | 10569 | 11486 | 11450 | 10778 | 5097 | 4655 | 0 |
| N50            | 10502 | 7562 | 9521 | 4464 | 7768 | 9590 | 4951 | 7522 | 9102 | 4816 | 8461 | 8903 | 5713 | 7677 | 11718 | 0 |
| L50            | 9064 | 7631 | 7618 | 7514 | 7514 | 7342 | 7472 | 7526 | 7386 | 7450 | 7404 | 7305 | 7538 | 2260 | 1743 | 0 |
| N90            | 25787 | 18148 | 23550 | 9375 | 18858 | 24219 | 10450 | 18209 | 22786 | 10167 | 20829 | 22522 | 12118 | 21232 | 37850 | 0 |
| L90            | 4218 | 3538 | 3240 | 5625 | 3384 | 2915 | 5610 | 3410 | 3000 | 5586 | 3209 | 2855 | 5618 | 1179 | 674 | 0 |
| % base removed | 0.00 | 42.86 | 26.76 | 68.51 | 33.14 | 18.28 | 59.78 | 30.58 | 16.66 | 58.23 | 8.60 | 4.56 | 40.40 | 31.58 | 11.96 | 100.00 |
| \# mapped read | 36206 | 24437 | 33640 | 10905 | 25569 | 35017 | 12150 | 24657 | 32813 | 11817 | 28518 | 32775 | 14102 | 26252 | 49976 | 0 |
| \# mismatch    | 29863768 | 1103030 | 4084253 | 642480 | 1090079 | 3794849 | 696601 | 1047965 | 3592823 | 675251 | 924368 | 1886849 | 706539 | 322402 | 604671 | 0 |
| time           | 0 | 3103.2063 | 2118.4074 | 162.7338 | 2743.7773 | 2067.4703 | 189.3940 | 2707.6767 | 1984.7657 | 187.9458 | 2429.9100 | 1782.7734 | 151.4521 | 793.1057 | 661.6553 | 17.1732 |
| memory         | 0 | 3357.47 | 4279.85 | 246.45 | 4208.07 | 4054.53 | 1882.32 | 3320.95 | 3187.07 | 1808.07 | 4421.41 | 2511.33 | 1489.20 | 2891.51 | 1301.25 | 1165.45 |

### Nanopore

|                | ont.raw.raw | ont.raw.canu | ont.raw.consent | ont.raw.mecat | ont.yacrd.consent | ont.yacrd.mecat | ont.yacrd2.canu | ont.yacrd2.consent | ont.yacrd2.mecat | ont.dascrubber.canu | ont.dascrubber.consent | ont.miniscrub.canu | ont.miniscrub.consent | ont.miniscrub.mecat |
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of read     | 25469 | 24439 | 24208 | 19972 | 24322 | 19835 | 24657 | 24208 | 19734 | 28518 | 22302 | 26254 | 57412 | 23761 |
| \# of base     | 257508441 | 147344360 | 252343673 | 239684772 | 251645063 | 238062600 | 145211033 | 250468181 | 236971782 | 161809806 | 237219041 | 54369709 | 175260949 | 138989242 |
| coverage       | 49.22x | 28.17x | 48.24x | 45.82x | 48.10x | 45.51x | 27.76x | 47.88x | 45.30x | 30.93x | 45.34x | 10.39x | 33.50x | 26.57x |
| N10            | 400 | 1088 | 390 | 472 | 390 | 469 | 1078 | 388 | 467 | 1204 | 365 | 868 | 740 | 556 |
| L10            | 51296 | 11636 | 51674 | 45872 | 51349 | 45843 | 11543 | 51380 | 45832 | 11486 | 51885 | 5097 | 18078 | 19610 |
| N50            | 3807 | 7562 | 3698 | 3832 | 3707 | 3805 | 7522 | 3694 | 3786 | 8461 | 3452 | 7677 | 8243 | 5527 |
| L50            | 20071 | 7631 | 20272 | 19356 | 20203 | 19385 | 7526 | 20164 | 19428 | 7404 | 20469 | 2260 | 5640 | 7368 |
| N90            | 13969 | 18148 | 13464 | 13032 | 13511 | 12937 | 18209 | 13461 | 12870 | 20829 | 12506 | 21232 | 35375 | 17878 |
| L90            | 4700 | 3538 | 4821 | 5301 | 4783 | 5290 | 3410 | 4776 | 5292 | 3209 | 4894 | 1179 | 1183 | 2800 |
| % base removed | 0.00 | 42.78 | 2.01 | 6.92 | 0.17 | 5.56 | 42.01 | -0.03 | 5.36 | 31.57 | -0.31 | 69.90 | 2.96 | 23.04 |
| \# mapped read | 25393 | 24437 | 24205 | 19972 | 24319 | 19835 | 24657 | 24206 | 19734 | 28518 | 22302 | 26252 | 57412 | 23761 |
| \# mismatch    | 37383344 | 1102903 | 4334858 | 7278735 | 4324461 | 7133121 | 1047670 | 4295612 | 7097879 | 924362 | 3063901 | 322402 | 1905281 | 2604923 |
| time           | 0 | 13676.9856 | 3970.4994 | 658.4983 | 3935.5846 | 639.6991 | 5.3595 | 4169.6171 | 667.4647 | 4.3111 | 3609.1146 | 1.9398 | 2265.2232 | 234.6207 |
| memory         | 0 | 4428.85 | 3693.43 | 2045.31 | 3663.08 | 2022.72 | 29.71 | 3440.56 | 2012.88 | 44.27 | 3608.59 | 21.70 | 2655.66 | 1665.74 |

## Assembly

### Pacbio

|                             | pb.raw.raw.miniasm | pb.raw.canu.miniasm | pb.raw.canu.canu | pb.raw.consent.miniasm | pb.raw.consent.canu | pb.raw.mecat.canu | pb.yacrd.canu.miniasm | pb.yacrd.canu.canu | pb.yacrd.consent.miniasm | pb.yacrd.consent.canu | pb.yacrd.mecat.canu | pb.yacrd2.canu.miniasm | pb.yacrd2.canu.canu | pb.yacrd2.consent.miniasm | pb.yacrd2.consent.canu | pb.yacrd2.mecat.canu | pb.dascrubber.canu.miniasm | pb.dascrubber.canu.canu | pb.dascrubber.consent.miniasm | pb.dascrubber.consent.canu | pb.dascrubber.mecat.canu |
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of contig                | 4 | 14 | 4 | 12 | 3 | 76 | 13 | 3 | 10 | 2 | 66 | 13 | 3 | 10 | 10 | 71 | 10 | 3 | 10 | 14 | 59 |
| \# of base                  | 5417095 | 5291308 | 5266503 | 5311045 | 5243442 | 5115332 | 5291304 | 5256701 | 5303726 | 5253082 | 5144689 | 5275568 | 5256427 | 5298556 | 5232660 | 5129748 | 5275853 | 5258099 | 5288414 | 5280884 | 5181184 |
| N10                         | 1 | 1 | 5267 | 1 | 5244 | 5116 | 1 | 5257 | 1 | 5254 | 5145 | 1 | 5257 | 1 | 5233 | 5130 | 1 | 5259 | 1 | 5281 | 5182 |
| L10                         | 1908669 | 767710 | 100 | 926662 | 100 | 100 | 842515 | 100 | 925237 | 100 | 100 | 842137 | 100 | 924945 | 100 | 100 | 1596591 | 100 | 1596472 | 100 | 100 |
| N50                         | 2 | 1 | 26333 | 2 | 26218 | 25577 | 2 | 26284 | 2 | 26266 | 25724 | 2 | 26283 | 2 | 26164 | 25649 | 1 | 26291 | 1 | 26405 | 25906 |
| L50                         | 1131354 | 767710 | 100 | 919875 | 100 | 100 | 767633 | 100 | 922177 | 100 | 100 | 652264 | 100 | 920774 | 100 | 100 | 1596591 | 100 | 1596472 | 100 | 100 |
| N90                         | 3 | 5 | 47399 | 5 | 47191 | 46038 | 5 | 47311 | 5 | 47278 | 46303 | 7 | 47308 | 5 | 47094 | 46168 | 3 | 47323 | 3 | 47528 | 46631 |
| L90                         | 402183 | 122965 | 100 | 171101 | 100 | 100 | 124240 | 100 | 171812 | 100 | 100 | 122211 | 100 | 171012 | 100 | 100 | 123792 | 100 | 171241 | 100 | 100 |
| genome fraction             | 1.791 | 99.884 | 99.976 | 98.987 | 99.809 | 96.170 | 99.901 | 99.983 | 99.108 | 99.984 | 96.565 | 99.892 | 99.983 | 99.073 | 99.736 | 96.549 | 99.889 | 99.984 | 99.865 | 99.723 | 97.346 |
| unaligned length            | 5323396 | 1978 | 1969 | 59155 | 1969 | 1825 | 1974 | 1969 | 53291 | 2906 | 1983 | 4992 | 1969 | 53300 | 1969 | 1983 | 1630 | 1897 | 5456 | 1969 | 1825 |
| misassemblies               | 0 | 3 | 4 | 5 | 5 | 5 | 3 | 3 | 5 | 3 | 3 | 3 | 3 | 7 | 3 | 4 | 3 | 3 | 2 | 4 | 5 |
| misassembled contigs length | 0 | 3349610 | 4222263 | 4607332 | 4544185 | 492065 | 1553459 | 5245605 | 1920816 | 5239708 | 567359 | 1235716 | 5245336 | 3999292 | 1644489 | 376646 | 4552591 | 5247333 | 2956300 | 2641331 | 676631 |
| time                        | 25.7627 | 47.2837 | 490.4215 | 68.3803 | 622.8198 | 374.2829 | 48.7166 | 569.5921 | 65.2462 | 708.1638 | 393.4832 | 46.0358 | 515.3420 | 60.9253 | 559.6089 | 383.3671 | 56.7287 | 558.8035 | 59.8127 | 569.9428 | 408.1536 |
| memory                      | 2489.82 | 1203.89 | 4282.86 | 1762.65 | 5799.73 | 696.57 | 1342.86 | 4326.30 | 1735.23 | 4842.86 | 2746.73 | 1284.82 | 4252.66 | 1681.21 | 4524.48 | 731.81 | 1547.96 | 4210.20 | 1471.61 | 4236.09 | 3433.34 |

### Nanopore

|                             | ont.raw.canu.miniasm | ont.raw.canu.canu | ont.raw.consent.miniasm | ont.raw.consent.canu | ont.raw.mecat.canu | ont.yacrd.canu.miniasm | ont.yacrd.canu.canu | ont.yacrd.consent.miniasm | ont.yacrd.consent.canu | ont.yacrd.mecat.canu | ont.yacrd2.canu.miniasm | ont.yacrd2.canu.canu | ont.yacrd2.consent.miniasm | ont.yacrd2.consent.canu | ont.yacrd2.mecat.canu | ont.dascrubber.canu.miniasm | ont.dascrubber.canu.canu | ont.dascrubber.consent.miniasm | ont.dascrubber.consent.canu | ont.miniscrub.consent.canu | ont.miniscrub.mecat.canu |
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of contig                | 14 | 10 | 1 | 3 | 1 | 11 | 12 | 1 | 1 | 1 | 14 | 12 | 1 | 1 | 1 | 11 | 6 | 1 | 3 | 5 | 12 |
| \# of base                  | 5257542 | 5311197 | 5194071 | 5275063 | 5225082 | 5272269 | 5331390 | 5198237 | 5257416 | 5220515 | 5276320 | 5331611 | 5194621 | 5262107 | 5220532 | 5274623 | 5283771 | 5208605 | 5253766 | 5244047 | 5249717 |
| N10                         | 1 | 5312 | 0 | 5276 | 5226 | 1 | 5332 | 0 | 5258 | 5221 | 1 | 5332 | 0 | 5263 | 5221 | 1 | 5284 | 0 | 5254 | 5245 | 5250 |
| L10                         | 1492921 | 100 | 0 | 100 | 100 | 886172 | 100 | 0 | 100 | 100 | 842131 | 100 | 0 | 100 | 100 | 1591843 | 100 | 0 | 100 | 100 | 100 |
| N50                         | 2 | 26556 | 0 | 26376 | 26126 | 2 | 26657 | 0 | 26288 | 26103 | 2 | 26659 | 0 | 26311 | 26103 | 2 | 26419 | 0 | 26269 | 26221 | 26249 |
| L50                         | 766050 | 100 | 0 | 100 | 100 | 842531 | 100 | 0 | 100 | 100 | 652278 | 100 | 0 | 100 | 100 | 924381 | 100 | 0 | 100 | 100 | 100 |
| N90                         | 5 | 47801 | 0 | 47476 | 47026 | 5 | 47983 | 0 | 47317 | 46985 | 7 | 47985 | 0 | 47359 | 46985 | 4 | 47554 | 0 | 47284 | 47197 | 47248 |
| L90                         | 122985 | 100 | 0 | 100 | 100 | 136628 | 100 | 0 | 100 | 100 | 122234 | 100 | 0 | 100 | 100 | 127117 | 100 | 0 | 100 | 100 | 100 |
| genome fraction             | 99.253 | 99.984 | 99.433 | 99.983 | 99.976 | 99.745 | 99.984 | 99.536 | 99.975 | 99.965 | 99.815 | 99.983 | 99.225 | 99.979 | 99.966 | 99.738 | 99.984 | 99.966 | 99.968 | 99.927 | 99.863 |
| unaligned length            | 1969 | 1969 | 25993 | 3149 | 3549 | 1969 | 1969 | 20617 | 3262 | 3445 | 4471 | 1969 | 33968 | 3043 | 3219 | 1621 | 1897 | 2670 | 2811 | 2043 | 3830 |
| misassemblies               | 5 | 4 | 3 | 5 | 3 | 4 | 4 | 3 | 3 | 3 | 3 | 4 | 5 | 3 | 4 | 5 | 4 | 3 | 3 | 3 | 3 |
| misassembled contigs length | 2239884 | 2779335 | 5194071 | 5275063 | 5225082 | 1562836 | 2779259 | 5198237 | 5257416 | 5220515 | 1232167 | 2779135 | 5194621 | 5262107 | 5220532 | 2535565 | 2783164 | 5208605 | 5216952 | 4823962 | 915334 |
| time                        | 60.5426 | 700.0259 | 157.2728 | 11009.5875 | 28701.7988 | 62.0122 | 702.6146 | 157.0813 | 9769.6387 | 28579.6023 | 59.1903 | 650.0218 | 163.3751 | 8089.7369 | 28569.4188 | 72.1994 | 716.8257 | 147.5832 | 5088.9658 | 767.2521 | 972.6400 |
| memory                      | 1777.24 | 4407.60 | 3703.86 | 8522.98 | 7609.25 | 1794.04 | 4437.46 | 3819.04 | 7945.97 | 7254.09 | 1759.77 | 4359.96 | 3746.79 | 7998.99 | 7181.24 | 1861.88 | 4285.22 | 3628.72 | 7287.12 | 4287.68 | 4306.69 |

# How to run

1. Install all tools in your PATH, and python3

2. Install python dependency `pip install -r requirements.txt`

3. Run:
   - download data `snakemake --snakefile pipeline/download.snakefile all`
   - run scrubbing `snakemake --snakefile pipeline/scrubbing.snakefile all`
   - run correction `snakemake --snakefile pipeline/correction.snakefile all`
   - run assembly `snakemake --snakefile pipeline/assembly.snakefile all`

4. Analysis:
   - map scrubbed `snakemake --snakefile pipeline/post_operation.snakefile mapping_scrubbing `
   - map corrected `snakemake --snakefile pipeline/post_operation.snakefile mapping_corrected`
   - run quast on assembly `snakemake --snakefile pipeline/post_operation.snakefile quast`
   - summarize info `./script/get_info.py -t [ont|pb] -s (raw) (yacrd) (yacrd2) (dascrubber) (miniscrub) -c (canu) (consent) (mecat) -a (miniasm) (canu)`
