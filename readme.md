# Dataset

- Reference:
  * [E. coli CFT073](https://www.uniprot.org/taxonomy/199310)  5.231428 Mb
  * [D. melanogaster] 143.726002 Mb
  * [C. elegans](ftp://ftp.ensembl.org/pub/release-95/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz) 100.2 Mb
  * [H. sapiens chr1](ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz) 248.9 Mb
- Reads:
  * E. coli:
	+ [Pacbio](https://www.ebi.ac.uk/ena/data/view/SRX5299472)
	+ [Oxford nanopore](https://www.ebi.ac.uk/ena/data/view/SRR8494940)
  * [Oxford nanopore D melanogaster](https://www.ebi.ac.uk/ena/data/view/SRX3676783)
  * [Oxford nanopore H sapiens chr1](http://s3.amazonaws.com/nanopore-human-wgs/chr1.sorted.bam)
  * [Pacbio RS P6-C4 C elegans](http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/list.html)


|                   | Pacbio         | Nanopore       | Pacbio subsample     | Nanopore subsample   | D melanogaster |
| ----------------  | --------------:| --------------:| --------------------:| --------------------:| --------------:|
| # of sequences    | 207069         | 158590         | 37404         (18 %) | 25469        (16 %)  | 1327569		 | 
| Total length      | 1425.446392 Mb | 1621.000527 Mb | 257.882884 Mb (18 %) | 257.508441 Mb (16 %) | 9064.470438 Mb | 
| Longest sequence  | 41.631 Kb      | 164.088 Kb     | 38.331 Kb            | 137.142 Kb           | 446.05 kb	     | 
| Shortest sequence | 35 b           | 88 b           | 35 b                 | 152 b                | 5 b			 |
| Mean Length       | 6.883 Kb       | 10.221 Kb      | 6.894 Kb             | 10.11 Kb             | 6.827 kb	     |
| Median Length     | 6.679 Kb       | 5.591 Kb       | 6.672 Kb             | 5.515 Kb             | 4.568 kb	     |
| N10               | 7467           | 2531           | 1354                 | 400                  | 29049		     |
| N50               | 58081          | 23903          | 10502                | 3807                 | 243356		 |
| N90               | 142591         | 86781          | 25787                | 13969                | 779045 		 |
| L10               | 15.631 Kb      | 50.77 Kb       | 15.591 Kb            | 51.316 Kb            | 25.964 kb	     |
| L50               | 9.052 Kb       | 20.189 Kb      | 9.064 Kb             | 20.073 Kb            | 11.853 kb	     |
| L90               | 4.191 Kb       | 4.763 Kb       | 4.218 Kb             | 4.701 Kb             | 3.533 kb       |
| Coverage          | 272x           | 309x           | 49x                  | 49x                  | 63x            |

# Tools

- [fpa](https://gitlab.inria.fr/pmarijon/fpa)

scrubber:
- [yacrd](https://gitlab.inria.fr/pmarijon/yacrd)
- [dascrubber](https://github.com/rrwick/DASCRUBBER-wrapper)
- [miniscrub](https://bitbucket.org/berkeleylab/jgi-miniscrub)

assembly tools:
- wtdbg2 2.4
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

Discard {intersection} {union} {jacard}: 

|            | yacrd           | dascrubber      | miniscrub |
| ---------- | ---------------:| ---------------:| ---------:|
| yacrd      |                 |                 |           |
| dascrubber | 487 7283 0.07   |                 |           |
| miniscrub  | 1247 13995 0.09 | 5554 13774 0.40 |           |

Splited {intersection} {union} {jacard}: 

|            | yacrd          | dascrubber     | miniscrub |
| ---------- | --------------:| --------------:| ---------:|
| yacrd      |                |                |           |
| dascrubber | 943 2893 0.33  |                |           |
| miniscrub  | 816 16520 0.05 | 920 16996 0.05 |           |

Trimmed {intersection} {union} {jacard}:

|            | yacrd            | dascrubber      | miniscrub |
| ---------- | ----------------:| ---------------:| ---------:|
| yacrd      |                  |                 |           |
| dascrubber | 26874 30977 0.87 |                 |           |
| miniscrub  | 5790 29536 0.20  | 6010 30001 0.20 |           |

Not modified {intersection} {union} {jacard}:

|            | yacrd         | dascrubber     | miniscrub |
| ---------- | -------------:| --------------:| ---------:|
| yacrd      |               |                |           |
| dascrubber | 0 5424 0.00   |                |           |
| miniscrub  | 276 6628 0.04 | 6 1620 0.00    |           |

### Nanopore

Count:

|           | yacrd | dascrubber | miniscrub |
| --------- | -----:| ----------:| ---------:|
| discard   |   138 |       3293 |      3096 |
| splited   |   184 |      10611 |       158 |
| trimmed   | 22590 |       2892 |     22215 |
| nmodified |  2557 |       8673 |         0 |

Discard {intersection} {union} {jacard}:

|            | yacrd        | dascrubber     | miniscrub |
| ---------- | ------------:| --------------:| ---------:|
| yacrd      |              |                |           |
| dascrubber | 81 3153 0.03 |                |           |
| miniscrub  | 90 3341 0.03 | 2005 4384 0.46 |           |


Splited {intersection} {union} {jacard}:

|            | yacrd          | dascrubber      | miniscrub |
| ---------- | --------------:| ---------------:| ---------:|
| yacrd      |                |                 |           |
| dascrubber | 87 255 0.34    |                 |           |
| miniscrub  | 130 10665 0.01 | 130 10639 0.01  |           |


Trimmed {intersection} {union} {jacard}:

|            | yacrd            | dascrubbe r     | miniscrub |
| ---------- | ----------------:| ---------------:| ---------:|
| yacrd      |                  |                 |           |
| dascrubber | 21679 23126 0.94 |                 |           |
| miniscrub  | 2718 22764 0.12  | 2782 22325 0.12 |           |

Not modified {intersection} {union} {jacard}:

|            | yacrd           | dascrubber     | miniscrub |
| ---------- | ---------------:| --------------:| ---------:|
| yacrd      |                 |                |           |
| dascrubber | 0 2557 0.00     |                |           |
| miniscrub  | 1186 10044 0.12 | 0 8673 0.00    |           |

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
| time           | 0 | 29.0794 | 29.4860 | 381.8210 | 19995.2767 |
| memory         | 0 | 2257.67 | 2229.11 | 13832.30 | 17234.49 |

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
| time           | 0 | 36.3609 | 36.6981 | 526.5879 | 21352.3889 |
| memory         | 0 | 2057.56 | 2090.84 | 27901.59 | 18240.23 |

## Correction

### Pacbio

|                | pb.raw.canu | pb.raw.consent | pb.raw.mecat | pb.yacrd.canu | pb.yacrd.consent | pb.yacrd.mecat | pb.yacrd2.canu | pb.yacrd2.consent | pb.yacrd2.mecat | pb.dascrubber.canu | pb.dascrubber.consent | pb.dascrubber.mecat | pb.miniscrub.canu | pb.miniscrub.consent | pb.miniscrub.mecat |
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of read     | 24439 | 33642 | 10905 | 25569 | 35021 | 12150 | 24657 | 32818 | 11817 | 28518 | 32775 | 14102 | 26254 | 49976 | 0 |
| \# of base     | 147344360 | 188870837 | 81199441 | 149811236 | 183105081 | 90110565 | 145211033 | 174306872 | 87363457 | 161809806 | 168948669 | 105499016 | 54369709 | 69959182 | 0 |
| coverage       | 28.17x | 36.10x | 15.52x | 28.64x | 35.00x | 17.22x | 27.76x | 33.32x | 16.70x | 30.93x | 32.29x | 20.17x | 10.39x | 13.37x | 0.00x |
| N10            | 1088 | 1331 | 662 | 1113 | 1357 | 728 | 1078 | 1291 | 708 | 1204 | 1258 | 836 | 868 | 1190 | 0 |
| L10            | 11636 | 12012 | 10501 | 11537 | 11528 | 10582 | 11543 | 11535 | 10569 | 11486 | 11450 | 10778 | 5097 | 4655 | 0 |
| N50            | 7562 | 9521 | 4464 | 7768 | 9590 | 4951 | 7522 | 9102 | 4816 | 8461 | 8903 | 5713 | 7677 | 11718 | 0 |
| L50            | 7631 | 7618 | 7514 | 7514 | 7342 | 7472 | 7526 | 7386 | 7450 | 7404 | 7305 | 7538 | 2260 | 1743 | 0 |
| N90            | 18148 | 23550 | 9375 | 18858 | 24219 | 10450 | 18209 | 22786 | 10167 | 20829 | 22522 | 12118 | 21232 | 37850 | 0 |
| L90            | 3538 | 3240 | 5625 | 3384 | 2915 | 5610 | 3410 | 3000 | 5586 | 3209 | 2855 | 5618 | 1179 | 674 | 0 |
| % base removed | 42.86 | 26.76 | 68.51 | 33.14 | 18.28 | 59.78 | 30.58 | 16.66 | 58.23 | 8.60 | 4.56 | 40.40 | 31.58 | 11.96 | 100.00 |
| \# mapped read | 24437 | 33640 | 10905 | 25569 | 35017 | 12150 | 24657 | 32813 | 11817 | 28518 | 32775 | 14102 | 26252 | 49976 | 0 |
| \# mismatch    | 1103030 | 4084253 | 642480 | 1090079 | 3794849 | 696601 | 1047965 | 3592823 | 675251 | 924368 | 1886849 | 706539 | 322402 | 604671 | 0 |
| time           | 3103.2063 | 2118.4074 | 162.7338 | 2883.5208 | 2207.0397 | 150.1093 | 2736.1575 | 2101.9235 | 140.8604 | 2429.9100 | 1782.7734 | 151.4521 | 793.1057 | 661.6553 | 17.1732 |
| memory         | 3357.47 | 4279.85 | 246.45 | 3279.80 | 4050.11 | 385.61 | 4201.26 | 3191.35 | 485.85 | 4421.41 | 2511.33 | 1489.20 | 2891.51 | 1301.25 | 1165.45 |

### Nanopore

|                | ont.raw.canu | ont.raw.consent | ont.raw.mecat | ont.yacrd.canu | ont.yacrd.consent | ont.yacrd.mecat | ont.yacrd2.canu | ont.yacrd2.consent | ont.yacrd2.mecat | ont.dascrubber.canu | ont.dascrubber.consent | ont.miniscrub.canu | ont.miniscrub.consent | ont.miniscrub.mecat |
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of read     | 24439 | 24208 | 19972 | 25569 | 24322 | 19835 | 24657 | 24208 | 19734 | 28518 | 22302 | 26254 | 57412 | 23761 |
| \# of base     | 147344360 | 252343673 | 239684772 | 149811236 | 251645063 | 238062600 | 145211033 | 250468181 | 236971782 | 161809806 | 237219041 | 54369709 | 175260949 | 138989242 |
| coverage       | 28.17x | 48.24x | 45.82x | 28.64x | 48.10x | 45.51x | 27.76x | 47.88x | 45.30x | 30.93x | 45.34x | 10.39x | 33.50x | 26.57x |
| N10            | 1088 | 390 | 472 | 1113 | 390 | 469 | 1078 | 388 | 467 | 1204 | 365 | 868 | 740 | 556 |
| L10            | 11636 | 51674 | 45872 | 11537 | 51349 | 45843 | 11543 | 51380 | 45832 | 11486 | 51885 | 5097 | 18078 | 19610 |
| N50            | 7562 | 3698 | 3832 | 7768 | 3707 | 3805 | 7522 | 3694 | 3786 | 8461 | 3452 | 7677 | 8243 | 5527 |
| L50            | 7631 | 20272 | 19356 | 7514 | 20203 | 19385 | 7526 | 20164 | 19428 | 7404 | 20469 | 2260 | 5640 | 7368 |
| N90            | 18148 | 13464 | 13032 | 18858 | 13511 | 12937 | 18209 | 13461 | 12870 | 20829 | 12506 | 21232 | 35375 | 17878 |
| L90            | 3538 | 4821 | 5301 | 3384 | 4783 | 5290 | 3410 | 4776 | 5292 | 3209 | 4894 | 1179 | 1183 | 2800 |
| % base removed | 42.78 | 2.01 | 6.92 | 40.57 | 0.17 | 5.56 | 42.01 | -0.03 | 5.36 | 31.57 | -0.31 | 69.90 | 2.96 | 23.04 |
| \# mapped read | 24437 | 24205 | 19972 | 25569 | 24319 | 19835 | 24657 | 24206 | 19734 | 28518 | 22302 | 26252 | 57412 | 23761 |
| \# mismatch    | 1102903 | 4334858 | 7278735 | 1089999 | 4324461 | 7133121 | 1047670 | 4295612 | 7097879 | 924362 | 3063901 | 322402 | 1905281 | 2604923 |
| time           | 13676.9856 | 3970.4994 | 658.4983 | 4.2001 | 4224.7420 | 603.6600 | 4.0139 | 4304.9277 | 627.2861 | 4.3111 | 3609.1146 | 1.9398 | 2265.2232 | 234.6207 |
| memory         | 4428.85 | 3693.43 | 2045.31 | 20.23 | 3641.38 | 507.81 | 43.94 | 3646.96 | 506.05 | 44.27 | 3608.59 | 21.70 | 2655.66 | 1665.74 |

## Assembly

### Pacbio

|                             | pb.raw.canu.miniasm | pb.raw.canu.canu | pb.raw.consent.miniasm | pb.raw.consent.canu | pb.raw.mecat.miniasm | pb.raw.mecat.canu | pb.yacrd.canu.miniasm | pb.yacrd.canu.canu | pb.yacrd.consent.miniasm | pb.yacrd.consent.canu | pb.yacrd.mecat.miniasm | pb.yacrd.mecat.canu | pb.yacrd2.canu.miniasm | pb.yacrd2.canu.canu | pb.yacrd2.consent.miniasm | pb.yacrd2.consent.canu | pb.yacrd2.mecat.miniasm | pb.yacrd2.mecat.canu | pb.dascrubber.canu.miniasm | pb.dascrubber.canu.canu | pb.dascrubber.consent.miniasm | pb.dascrubber.consent.canu | pb.dascrubber.mecat.miniasm | pb.dascrubber.mecat.canu | pb.miniscrub.canu.miniasm | pb.miniscrub.canu.canu | pb.miniscrub.consent.miniasm |
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of contig                | 14 | 4 | 12 | 3 | 134 | 76 | 13 | 3 | 10 | 2 | 129 | 66 | 13 | 3 | 10 | 10 | 129 | 71 | 10 | 3 | 10 | 14 | 123 | 59 | 9 | 182 | 6 |
| \# of base                  | 5291308 | 5266503 | 5311045 | 5243442 | 4143043 | 5115332 | 5291304 | 5256701 | 5303726 | 5253082 | 4410612 | 5144689 | 5275568 | 5256427 | 5298556 | 5232660 | 4216954 | 5129748 | 5275853 | 5258099 | 5288414 | 5280884 | 4728589 | 5181184 | 76245 | 4875069 | 51551 |
| N10                         | 1 | 5267 | 1 | 5244 | 5 | 5116 | 1 | 5257 | 1 | 5254 | 5 | 5145 | 1 | 5257 | 1 | 5233 | 5 | 5130 | 1 | 5259 | 1 | 5281 | 3 | 5182 | 1 | 4876 | 1 |
| L10                         | 767710 | 100 | 926662 | 100 | 76598 | 100 | 842515 | 100 | 925237 | 100 | 86026 | 100 | 842137 | 100 | 924945 | 100 | 82282 | 100 | 1596591 | 100 | 1596472 | 100 | 122136 | 100 | 10814 | 100 | 10763 |
| N50                         | 1 | 26333 | 2 | 26218 | 37 | 25577 | 2 | 26284 | 2 | 26266 | 35 | 25724 | 2 | 26283 | 2 | 26164 | 35 | 25649 | 1 | 26291 | 1 | 26405 | 31 | 25906 | 4 | 24376 | 3 |
| L50                         | 767710 | 100 | 919875 | 100 | 35865 | 100 | 767633 | 100 | 922177 | 100 | 42987 | 100 | 652264 | 100 | 920774 | 100 | 38279 | 100 | 1596591 | 100 | 1596472 | 100 | 44550 | 100 | 7736 | 100 | 7304 |
| N90                         | 5 | 47399 | 5 | 47191 | 104 | 46038 | 5 | 47311 | 5 | 47278 | 98 | 46303 | 7 | 47308 | 5 | 47094 | 100 | 46168 | 3 | 47323 | 3 | 47528 | 92 | 46631 | 8 | 43876 | 0 |
| L90                         | 122965 | 100 | 171101 | 100 | 16556 | 100 | 124240 | 100 | 171812 | 100 | 17509 | 100 | 122211 | 100 | 171012 | 100 | 17292 | 100 | 123792 | 100 | 171241 | 100 | 19066 | 100 | 5818 | 100 | 0 |
| genome fraction             | 99.884 | 99.976 | 98.987 | 99.809 | 78.704 | 96.170 | 99.901 | 99.983 | 99.108 | 99.984 | 83.700 | 96.565 | 99.892 | 99.983 | 99.073 | 99.736 | 79.946 | 96.549 | 99.889 | 99.984 | 99.865 | 99.723 | 89.352 | 97.346 | 1.460 | 92.977 | 0.985 |
| unaligned length            | 1978 | 1969 | 59155 | 1969 | 0 | 1825 | 1974 | 1969 | 53291 | 2906 | 1275 | 1983 | 4992 | 1969 | 53300 | 1969 | 1865 | 1983 | 1630 | 1897 | 5456 | 1969 | 1966 | 1825 | 0 | 0 | 0 |
| misassemblies               | 3 | 4 | 5 | 5 | 7 | 5 | 3 | 3 | 5 | 3 | 7 | 3 | 3 | 3 | 7 | 3 | 5 | 4 | 3 | 3 | 2 | 4 | 9 | 5 | 1 | 5 | 0 |
| misassembled contigs length | 3349610 | 4222263 | 4607332 | 4544185 | 344743 | 492065 | 1553459 | 5245605 | 1920816 | 5239708 | 370620 | 567359 | 1235716 | 5245336 | 3999292 | 1644489 | 228108 | 376646 | 4552591 | 5247333 | 2956300 | 2641331 | 458225 | 676631 | 6201 | 246671 | 0 |
| time                        | 47.2837 | 490.4215 | 68.3803 | 622.8198 | 18.2072 | 374.2829 | 48.7166 | 569.5921 | 65.2462 | 708.1638 | 21.9514 | 393.4832 | 46.0358 | 515.3420 | 60.9253 | 559.6089 | 20.9881 | 383.3671 | 56.7287 | 558.8035 | 59.8127 | 569.9428 | 27.4294 | 408.1536 | 9.6378 | 200.8158 | 15.6122 |
| memory                      | 1203.89 | 4282.86 | 1762.65 | 5799.73 | 784.16 | 696.57 | 1342.86 | 4326.30 | 1735.23 | 4842.86 | 851.81 | 2746.73 | 1284.82 | 4252.66 | 1681.21 | 4524.48 | 827.32 | 731.81 | 1547.96 | 4210.20 | 1471.61 | 4236.09 | 924.59 | 3433.34 | 489.35 | 1458.04 | 653.78 |

### Nanopore

|                             | ont.raw.canu.miniasm | ont.raw.canu.canu | ont.raw.consent.miniasm | ont.raw.consent.canu | ont.raw.mecat.miniasm | ont.raw.mecat.canu | ont.yacrd.canu.miniasm | ont.yacrd.canu.canu | ont.yacrd.consent.miniasm | ont.yacrd.consent.canu | ont.yacrd.mecat.miniasm | ont.yacrd.mecat.canu | ont.yacrd2.canu.miniasm | ont.yacrd2.canu.canu | ont.yacrd2.consent.miniasm | ont.yacrd2.consent.canu | ont.yacrd2.mecat.miniasm | ont.yacrd2.mecat.canu | ont.dascrubber.canu.miniasm | ont.dascrubber.canu.canu | ont.dascrubber.consent.miniasm | ont.dascrubber.consent.canu | ont.miniscrub.canu.miniasm | ont.miniscrub.canu.canu | ont.miniscrub.consent.miniasm | ont.miniscrub.consent.canu | ont.miniscrub.mecat.miniasm | ont.miniscrub.mecat.canu |
| - | -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:| -:|
| \# of contig                | 14 | 10 | 1 | 3 | 1 | 1 | 11 | 12 | 1 | 1 | 1 | 1 | 14 | 12 | 1 | 1 | 1 | 1 | 11 | 6 | 1 | 3 | 9 | 178 | 6 | 5 | 9 | 12 |
| \# of base                  | 5257542 | 5311197 | 5194071 | 5275063 | 5098951 | 5225082 | 5272269 | 5331390 | 5198237 | 5257416 | 5101041 | 5220515 | 5276320 | 5331611 | 5194621 | 5262107 | 5099858 | 5220532 | 5274623 | 5283771 | 5208605 | 5253766 | 81541 | 4869118 | 5068024 | 5244047 | 5126504 | 5249717 |
| N10                         | 1 | 5312 | 0 | 5276 | 0 | 5226 | 1 | 5332 | 0 | 5258 | 0 | 5221 | 1 | 5332 | 0 | 5263 | 0 | 5221 | 1 | 5284 | 0 | 5254 | 1 | 4870 | 1 | 5245 | 1 | 5250 |
| L10                         | 1492921 | 100 | 0 | 100 | 0 | 100 | 886172 | 100 | 0 | 100 | 0 | 100 | 842131 | 100 | 0 | 100 | 0 | 100 | 1591843 | 100 | 0 | 100 | 11418 | 100 | 604843 | 100 | 646536 | 100 |
| N50                         | 2 | 26556 | 0 | 26376 | 0 | 26126 | 2 | 26657 | 0 | 26288 | 0 | 26103 | 2 | 26659 | 0 | 26311 | 0 | 26103 | 2 | 26419 | 0 | 26269 | 4 | 24346 | 1 | 26221 | 2 | 26249 |
| L50                         | 766050 | 100 | 0 | 100 | 0 | 100 | 842531 | 100 | 0 | 100 | 0 | 100 | 652278 | 100 | 0 | 100 | 0 | 100 | 924381 | 100 | 0 | 100 | 9584 | 100 | 604843 | 100 | 592502 | 100 |
| N90                         | 5 | 47801 | 0 | 47476 | 0 | 47026 | 5 | 47983 | 0 | 47317 | 0 | 46985 | 7 | 47985 | 0 | 47359 | 0 | 46985 | 4 | 47554 | 0 | 47284 | 8 | 43823 | 4 | 47197 | 6 | 47248 |
| L90                         | 122985 | 100 | 0 | 100 | 0 | 100 | 136628 | 100 | 0 | 100 | 0 | 100 | 122234 | 100 | 0 | 100 | 0 | 100 | 127117 | 100 | 0 | 100 | 5805 | 100 | 339331 | 100 | 113483 | 100 |
| genome fraction             | 99.253 | 99.984 | 99.433 | 99.983 | 99.903 | 99.976 | 99.745 | 99.984 | 99.536 | 99.975 | 99.918 | 99.965 | 99.815 | 99.983 | 99.225 | 99.979 | 99.941 | 99.966 | 99.738 | 99.984 | 99.966 | 99.968 | 1.537 | 92.916 | 97.209 | 99.927 | 99.190 | 99.863 |
| unaligned length            | 1969 | 1969 | 25993 | 3149 | 3710 | 3549 | 1969 | 1969 | 20617 | 3262 | 3219 | 3445 | 4471 | 1969 | 33968 | 3043 | 3677 | 3219 | 1621 | 1897 | 2670 | 2811 | 1224 | 0 | 1751 | 2043 | 1920 | 3830 |
| misassemblies               | 5 | 4 | 3 | 5 | 2 | 3 | 4 | 4 | 3 | 3 | 5 | 3 | 3 | 4 | 5 | 3 | 3 | 4 | 5 | 4 | 3 | 3 | 0 | 7 | 3 | 3 | 3 | 3 |
| misassembled contigs length | 2239884 | 2779335 | 5194071 | 5275063 | 5098951 | 5225082 | 1562836 | 2779259 | 5198237 | 5257416 | 5101041 | 5220515 | 1232167 | 2779135 | 5194621 | 5262107 | 5099858 | 5220532 | 2535565 | 2783164 | 5208605 | 5216952 | 0 | 269916 | 890726 | 4823962 | 909581 | 915334 |
| time                        | 60.5426 | 700.0259 | 157.2728 | 11009.5875 | 133.3762 | 28701.7988 | 62.0122 | 702.6146 | 157.0813 | 9769.6387 | 134.2463 | 28579.6023 | 59.1903 | 650.0218 | 163.3751 | 8089.7369 | 134.0414 | 28569.4188 | 72.1994 | 716.8257 | 147.5832 | 5088.9658 | 12.7905 | 209.7028 | 89.0136 | 767.2521 | 53.9216 | 972.6400 |
| memory                      | 1777.24 | 4407.60 | 3703.86 | 8522.98 | 3406.86 | 7609.25 | 1794.04 | 4437.46 | 3819.04 | 7945.97 | 3414.40 | 7254.09 | 1759.77 | 4359.96 | 3746.79 | 7998.99 | 3390.58 | 7181.24 | 1861.88 | 4285.22 | 3628.72 | 7287.12 | 687.51 | 1880.11 | 2138.55 | 4287.68 | 1765.57 | 4306.69 |

# How to run

1. Install all tools in your PATH, and python3

2. Install python dependency `pip install -r requirements.txt`

3. Run:
   - download data `snakemake --snakefile pipeline/download.snakefile all`
   - run yacrd analysis `snakemake --snakefile pipeline/uncorrected.snakefile --directory /home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly -j 999 --cluster-config config/cluster.json -c "sbatch --gres={cluster.gres} --mem={cluster.mem} --nodes={cluster.nodes} --mail-type={cluster.mail-type} --time={cluster.time} --ntasks-per-node={cluster.ntasks-per-node} --mail-user={cluster.mail-user}" -r -p all`
   - run fpa analysis `snakemake --snakefile pipeline/fpa.snakefile --directory /home/pierre.marijon/data/optimizing-early-steps-of-lr-assembly -j 999 --cluster-config config/cluster.json -c "sbatch --gres={cluster.gres} --mem={cluster.mem} --nodes={cluster.nodes} --mail-type={cluster.mail-type} --time={cluster.time} --ntasks-per-node={cluster.ntasks-per-node} --mail-user={cluster.mail-user}" -r -p all`

4. Get result:
   - summarize info `./script/get_info.py -t [ont|pb] -s (raw) (yacrd) (yacrd2) (dascrubber) (miniscrub) -c (canu) (consent) (mecat) -a (miniasm) (canu)`


