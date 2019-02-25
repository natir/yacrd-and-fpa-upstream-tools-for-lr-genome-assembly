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

|                        | raw        | yacrd      | dascrubber | miniscrub  |
|:---------------------- | ----------:| ----------:| ----------:| ----------:|
| # of read              | 25469      | 25660      | 22538      | 62915      |
| # of base (Mb)         | 257.508441 | 252.085244 | 236.474486 | 180.605742 |
| % of base removed (Mb) |            | 5.423197   | 21.033955  | 76.902699  |
| coverage               |            | 48x        | 45x        | 34x        |
| % mapped read          |            |            |            |            |
| % mapped base          |            |            |            |            |
| # mismatch             |            |            |            |            |
| time (s)               |            | 28.4964    | 526.5879   | 21352.3889 |
| memory (Mo)            |            | 2103.60    | 27901.59   | 18240.23   |

## Correction

### Pacbio

#### Raw

|        | raw.canu	 | raw.consent | raw.mecat |
| ------ | ---------:| -----------:| ---------:|
| time   | 3103.2063 | 2118.4074   | 162.7338  |
| memory | 3357.47   | 4279.85	   | 246.45    |

#### yacrd

|        | yacrd.canu |	yacrd.consent | yacrd.mecat |
| ------ | ----------:| -------------:| -----------:|
| time   | 2743.7773  | 2067.4703     | 189.394     |
| memory | 4208.07    | 4054.53	      | 1882.32     |

#### dascrubber

|        | dascrubber.canu	| dascrubber.consent	| dascrubber.mecat |
| ------ | ----------------:| ---------------------:| ----------------:|
| time   | 2429.91			| 1782.7734			    | 151.4521         |
| memory | 4421.41			| 2511.33				| 1489.2           |

#### mecat

|        | miniscrub.canu	| miniscrub.consent	| miniscrub.mecat |
| ------ | ----------------:| -----------------:| ---------------:|
| time   | 793.1057		    | 661.6553			| 17.1732         |
| memory | 2891.51			| 1301.25			| 1165.45         |


### Nanopore

#### raw

| raw.canu	    | raw.mecat | raw.consent |
| 13676.9856	| 658.4983	| 3970.4994   |
| 4428.85		| 2045.31	| 3693.43     |

#### yacrd

|	yacrd.consent	| yacrd.mecat	| yacrd.canu	|
|	3935.5846		| 639.6991		| 4.0438        |
|	3663.08			| 2022.72		| 44.34         |


#### dascrubber

|	dascrubber.consent	| dascrubber.canu |
|	3609.1146			| 4.3111          |
|	3608.59				| 44.27           |

#### miniscrub

|	miniscrub.consent	| miniscrub.mecat	| miniscrub.canu	|
|	2265.2232			| 234.6207			| 1.9398			|
|	2655.66				| 1665.74			| 21.7				|

## Correction

### Pacbio

#### raw 

|        | raw.raw.miniasm	| raw.consent.miniasm	| raw.consent.canu	| raw.canu.canu	| raw.canu.miniasm	| raw.raw.canu	| raw.mecat.canu | raw.mecat.miniasm |
| ------ | ----------------:| ---------------------:| -----------------:| -------------:| -----------------:| -------------:| --------------:| -----------------:|
| time   | 25.7627			| 68.3803				| 0.3703			| 0.3781		| 47.2837			| 22.8496		| 0.3645		 | 18.2072           |
| memory | 2489.82			| 1762.65				| 2.88				| 3	            | 1203.89	 	    | 952.37		| 2.87			 | 784.16            |

#### yacrd

|        | yacrd.consent.miniasm | yacrd.consent.canu	| yacrd.canu.canu	| yacrd.canu.miniasm	| yacrd.mecat.miniasm	| yacrd.mecat.canu |
| ------ | ---------------------:| --------------------:| -----------------:| ---------------------:| ---------------------:| ----------------:|
| time   | 65.2462	             | 0.3688	            | 0.3674    	    | 48.7166	            | 21.9514            	| 0.3704           |
| memory | 1735.23	             | 2.92	                | 2.95	            | 1342.86	            | 851.81	            | 2.92             | 

#### dascrubber

| dascrubber.consent.canu | dascrubber.consent.miniasm | dascrubber.canu.miniasm | dascrubber.canu.canu | dascrubber.mecat.canu | dascrubber.mecat.miniasm |
| -----------------------:| --------------------------:| -----------------------:| --------------------:| ---------------------:| ------------------------:|
| 0.371                   | 59.8127                    | 56.7287                 | 0.3825               | 0.3665                | 27.4294                  |
| 2.94                    | 1471.61                    | 1547.96                 | 2.97                 | 2.95                  | 924.59                   |

#### miniscrub

| miniscrub.mecat.miniasm | miniscrub.consent.miniasm | miniscrub.canu.canu	| miniscrub.canu.miniasm	| miniscrub.consent.canu	| miniscrub.mecat.canu |
| -----------------------:| -------------------------:| -------------------:| -------------------------:| -------------------------:| --------------------:|
| 0.0709				  | 15.6122					  | 0.3663				| 9.6378					| 134.5584					| 0.3913               |
| 2.88					  | 653.78					  | 2.94				| 489.35					| 1413.37					| 2.92                 |

### Nanopore

#### raw

|	raw.raw.miniasm | raw.consent.miniasm	| raw.consent.canu	| raw.mecat.canu	| raw.mecat.miniasm | raw.canu.miniasm	| raw.canu.canu | raw.raw.canu |
|	40.2725			| 157.2728				| 0.4003			| 0.4384			| 133.3762			| 60.5426			| 0.3701		| 51.6651      |
|	3489.91			| 3703.86				| 2.94				| 3					| 3406.86			| 1777.24			| 2.95			| 410.22       |

#### yacrd

|	yacrd.consent.miniasm	| yacrd.consent.canu	| yacrd.mecat.canu	| yacrd.mecat.miniasm	| yacrd.canu.miniasm	| yacrd.canu.canu |
|	157.0813				| 0.3906				| 0.3787			| 134.2463				| 62.0122				| 0.3665          |
|	3819.04					| 2.95					| 2.95				| 3414.4				| 1794.04				| 2.95            |

#### dascrubber

|	dascrubber.consent.miniasm	| dascrubber.consent.canu	| dascrubber.canu.miniasm	| dascrubber.canu.canu |
|	147.5832					| 0.3816					| 72.1994					| 0.3677               |
|	3628.72						| 3							| 1861.88					| 2.94                 |


#### miniscrub

|	miniscrub.mecat.miniasm | miniscrub.consent.miniasm | miniscrub.consent.canu	| miniscrub.mecat.canu	| miniscrub.canu.miniasm	| miniscrub.canu.canu |
|	53.9216					| 89.0136					| 0.3934					| 0.3715				| 12.7905					| 0.3716              |
|	1765.57					| 2138.55					| 2.88						| 2.95					| 687.51					| 3                   |


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
