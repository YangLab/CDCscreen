# CDCscreen
A computational pipeline for **C**as13**d**-mediated **c**ircRNA **screen** (CDCscreen) to identify negatively selected functional circular RNAs.


## Schema
![image](https://github.com/xueweireally/CDCscreen/blob/master/doc/CDCscreen_pipeline.jpg)

## Installation requirements

* Software
- perl (version 5.26.2)
- bowtie (version 1.1.2)
- cutadapt (version 1.18)
- samtools (version: 1.9)
- MAGeCK (version 0.5.9.2)
- R (version 3.5.1)

## Date requirements
* Cas13d BSJ-gRNA reference sequences
* Expression (FPBcirc) of circRNAs in examined cells
* FASTQ files (Paired-End, 2 biology replicates of control [Day 1] and treatment [Day 30])
- Day 1 of biology replicate 1
- Day 1 of biology replicate 2
- Day 30 of biology replicate 1
- Day 30 of biology replicate 2

## Installation
```bash
git clone https://github.com/xueweireally/CDCscreen
```

## Usage
```bash
Usage: sh run_CDCscreen.sh ref_gRNA_seq.fa FPBcirc.txt D1_rep1_R1.fq D1_rep1_R2.fq D1_rep2_R1.fq D1_rep2_R2.fq D30_rep1_R1.fq D30_rep1_R2.fq D30_rep2_R1.fq D30_rep2_R2.fq
```

### Example
```bash
sh run_CDCscreen.sh test_data/ref_gRNA_seq.fa test_data/FPBcirc.txt D1_rep1_R1.fq D1_rep1_R2.fq D1_rep2_R1.fq D1_rep2_R2.fq D30_rep1_R1.fq D30_rep1_R2.fq D30_rep2_R1.fq D30_rep2_R2.fq
```

* Test Cas13d BSJ-gRNA reference sequences file was in directory [test_data]
- test_data/ref_gRNA_seq.fa

* Test expression (FPBcirc) of circRNAs of 293FT cells was in directory [test_data]
- test_data/FPBcirc.txt

* Test data (D1 and D30 of 293FT cells) were downloaded from NCBI GEO dataset (GES:xxxxxxx)
- D1_rep1_R1.fq, Fastq file of day 1 biology replicate 1 R1
- D1_rep1_R2.fq, Fastq file of day 1 biology replicate 1 R2
- D1_rep2_R1.fq, Fastq file of day 1 biology replicate 2 R1
- D1_rep2_R2.fq, Fastq file of day 1 biology replicate 2 R2
- D30_rep1_R1.fq, Fastq file of day 30 biology replicate 1 R1
- D30_rep1_R2.fq, Fastq file of day 30 biology replicate 1 R2
- D30_rep2_R1.fq, Fastq file of day 30 biology replicate 2 R1
- D30_rep2_R2.fq, Fastq file of day 30 biology replicate 2 R2



