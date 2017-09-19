# TGIRT smRNA correction #

This package used a linear Ridge model to correct count for TGIRT-seq data. A linear model follows:

![](https://raw.githubusercontent.com/wckdouglas/tgirt_smRNA_correction/master/img/linear_model.png?token=AGGQnE6rCxzOmhjEbdxH0VEkaRLlYCcSks5ZysD1wA%3D%3D)

where Y is the $$\Delta$$CPM between true CPM and observed CPM, X are the positional nucleotides while $$\hat{a}$$ are the influencys of each position nucleotides.


Workflow:

1. build index
2. Run correction 

## Build index ##

This module build an index for correction from (1) TGIRT count on a set of known transcript; (2) small RNA transcriptome fasta file; (3) ground truth gene count of each known transcript


* file needed:
	* fasta file (see **test/mir9_2.fa**)
	* TGIRT experimental count (see **test/tgirt_count.csv**)

* file created:
  * Index file ($OUTPUT_PREFIX$\_index.pkl)

```
usage: build_index.py [-h] -f FASTA [-n NUCLEOTIDE] [-o OUTPUT_PREFIX]
                      [-e EXPECTED_COUNT] -c EXPERIMENTAL_COUNT

Building nucleotide table for each transcripts

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Small RNA fasta file
  -n NUCLEOTIDE, --nucleotide NUCLEOTIDE
                        How many nucleotide to look at from both end?
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        How many nucleotide to look at from both end?
  -e EXPECTED_COUNT, --expected_count EXPECTED_COUNT
                        Expected count (comma delimintaed: seq,count)
                        (default: all equal)
  -c EXPERIMENTAL_COUNT, --experimental_count EXPERIMENTAL_COUNT
                        Experimental count (comma delimintaed: seq,count)
```


## Add correction factor to fragment BED file ##

This module use the index from last step to add a column indicating the correction factor to a fragment bed file.

* file needed:
  * genome fasta file
  * TGIRT-seq bed file containing fragment coordinate (see **test/test.bed**; bam_to_bed.py from [tgirt-seq-tools](https://github.com/wckdouglas/tgirt_seq_tools))
  * index file 

* file created:
  * bed file with a new column storing correction factor  (see: **test/out.bed**)

```
usage: tgirt_correction.py [-h] -f FASTA [-b BED] [-i INDEX]
                           [-o OUTPUT_PREFIX]

Building nucleotide table for each transcripts

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Genome fasta file
  -b BED, --bed BED     Input fragment bed file (default: -)
  -i INDEX, --index INDEX
                        Index
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
```

