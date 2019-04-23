# TGIRT smRNA correction #

This package used a read reweighing scheme, as detailed in [Hansen *et al*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896536/). This package is designed for use in TGIRT-se paired-end data.

![](https://github.com/wckdouglas/tgirt_smRNA_correction/blob/master/img/reweighing.png?raw=true)

where *w(s1, s2)* is the computed weight for a read pair having trinucleotides *s1* as the start of read 1 and *s2* as the start of read 2. This weight is the geometric mean of weight of seeing *s1* in the start of read1 and *s2* in the start of read2. The weight for *s* in each read is calculated by the dividing the probability of seeing *s* in the middle of the read (any position *i* > 3) by the probability of seeing *s* in the start of the read.


Workflow:

1. build index
2. Run correction 

---

## Train ##
This function trains a reweighing model using first three nucleotides from each reads.



* file needed:
	* paired-end bam-file (can be sorted or unsorted)

* file created:
  * Index file (e.g. weights.pkl)

```
usage: tgirt_correction.py train [-h] -i INBAM [-x WEIGHT_INDEX] [-c ITER]

optional arguments:
  -h, --help            show this help message and exit
  -i INBAM, --inbam INBAM
                        Input bam file
  -x WEIGHT_INDEX, --weight_index WEIGHT_INDEX
                        Output weight index
  -c ITER, --iter ITER  How many reads to analyze for each end
```


## Add weight to each alignment in bam file ##

This module use the index from last step to add a tag (**ZW**) indicating the weight to a bam file.

* file needed:
  * bam file (name-sorted!!)
  * index file 

* file created:
  * bam file with AS tag added

```
usage: tgirt_correction.py correct [-h] -i INF [-x INDEX] [-o OUTF] [--bed]
                                   [-f FASTA]

optional arguments:
  -h, --help            show this help message and exit
  -i INF, --inf INF     Input fragment (name sorted bam file!!)
  -x INDEX, --index INDEX
                        Weight index trained by "train" command (default: /sto
                        r/work/Lambowitz/cdw2854/src/miniconda3/lib/python3.6/
                        site-
                        packages/tgirt_smRNA_correction/model/weights.pkl)
  -o OUTF, --outf OUTF  output file (default: -)
  --bed                 input and output files are bed files, otherwise bam
                        files
  -f FASTA, --fasta FASTA
                        Genome fasta file, for correction with bed only

```

