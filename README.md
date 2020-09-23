# CDCscreen
A computational pipeline for **C**as13**d**-mediated **c**ircRNA **screen** (CDCscreen) to identify negatively selected functional circular RNAs.


## Schema
![image](https://github.com/xueweireally/CDCscreen/blob/master/doc/CDCscreen_pipeline.jpg)

## Installation requirements

* Software
perl (version 5.26.2)
bowtie (version 1.1.2)
cutadapt (version 1.18)
samtools (version: 1.9)
MAGeCK (version 0.5.9.2)
R (version 3.5.1)

## Installation
```bash
git clone https://github.com/xueweireally/CDCscreen
```

## Usage
Usage: sh run_CDCscreen.sh ref_gRNA_seq.fa FPBcirc.txt D1_rep1_R1.fq D1_rep1_R2.fq D1_rep2_R1.fq D1_rep2_R2.fq D30_rep1_R1.fq D30_rep1_R2.fq D30_rep2_R1.fq D30_rep2_R2.fq

### Example
sh run_CDCscreen.sh test_data/ref_gRNA_seq.fa test_data/FPBcirc.txt D1_rep1_R1.fq D1_rep1_R2.fq D1_rep2_R1.fq D1_rep2_R2.fq D30_rep1_R1.fq D30_rep1_R2.fq D30_rep2_R1.fq D30_rep2_R2.fq

* D1_rep1_R1.fq, D1_rep1_R2.fq, D1_rep2_R1.fq, D1_rep2_R2.fq, D30_rep1_R1.fq, D30_rep1_R2.fq, D30_rep2_R1.fq and D30_rep2_R2.fq were downloaded from NCBI GEO dataset
