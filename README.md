# CDCscreen
A computational pipeline for **C**as13**d**-mediated **c**ircRNA **screen** (CDCscreen) to identify negatively selected functional circular RNAs.


## Schema
![image](https://github.com/xueweireally/CDCscreen/blob/master/doc/CDCscreen_pipeline.jpg)

Authors: Wei Xue (xuewei@picb.ac.cn), Li Yang (liyang@picb.ac.cn)

Maintainer: Wei Xue (xuewei@picb.ac.cn)

## Installation requirements

* Software
    - perl (version 5.26.2)
    - bowtie (version 1.1.2)
    - cutadapt (version 1.18)
    - samtools (version: 1.9)
    - MAGeCK (version 0.5.9.2)
    - R (version 3.5.1)

## Date requirements
* 1. Cas13d BSJ-gRNA reference sequences
    - ref_gRNA_seq.fa
* 2. Expression (FPBcirc) of circRNAs in examined cells
    - FPBcirc.txt
* 3. Raw FASTQ files (Paired-End, 2 biology replicates of control [Day 1] and treatment [Day 30])
    - Day 1 of biology replicate 1, D1_rep1_R1.fq and D1_rep1_R2.fq
    - Day 1 of biology replicate 2, D1_rep2_R1.fq and D1_rep2_R2.fq
    - Day 30 of biology replicate 1, D30_rep1_R1.fq and D30_rep1_R2.fq
    - Day 30 of biology replicate 2, D30_rep2_R1.fq and D30_rep2_R2.fq

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
sh run_CDCscreen.sh test_data/ref_gRNA_seq.fa test_data/FPBcirc_293FT.txt D1_rep1_R1.fq D1_rep1_R2.fq D1_rep2_R1.fq D1_rep2_R2.fq D30_rep1_R1.fq D30_rep1_R2.fq D30_rep2_R1.fq D30_rep2_R2.fq
```

### Input
* Test Cas13d BSJ-gRNA reference sequences file was in directory [test_data]
    - test_data/ref_gRNA_seq.fa

* Test expression (FPBcirc) of circRNAs of 293FT cells was in directory [test_data]
    - test_data/FPBcirc_293FT.txt

* Test raw data (D1 and D30 of 293FT cells) were downloaded from NCBI GEO dataset (GES:xxxxxxx)
    - D1_rep1_R1.fq, Fastq file of day 1 biology replicate 1 R1
    - D1_rep1_R2.fq, Fastq file of day 1 biology replicate 1 R2
    - D1_rep2_R1.fq, Fastq file of day 1 biology replicate 2 R1
    - D1_rep2_R2.fq, Fastq file of day 1 biology replicate 2 R2
    - D30_rep1_R1.fq, Fastq file of day 30 biology replicate 1 R1
    - D30_rep1_R2.fq, Fastq file of day 30 biology replicate 1 R2
    - D30_rep2_R1.fq, Fastq file of day 30 biology replicate 2 R1
    - D30_rep2_R2.fq, Fastq file of day 30 biology replicate 2 R2

### Output
* 05_CDCscreen/06_CDCscreen_circRNA.txt

| Field       | Description                      |
| :---------- | :--------------------------------|
| chrom       | Chromosome                       |
| start       | Start of circRNA                 |
| end         | End of circRNA                   |
| geneName    | Gene symbol of circRNA           |
| FPBcirc     | Expression of circRNA            |
| CDCscreen   | CDCscreen score of circRNA       |


## Citation
Siqi Li#, Xiang Li#, Wei Xue#, Lin Zhang#, Liang-Zhong Yang, Shi-Meng Cao, Yun-Ni Lei, Chu-Xiao Liu, Si-Kun Guo, Lin Shan, Man Wu, Xiao Tao, Jia-Lin Zhang, Xiang Gao, Jun Zhang, Jia Wei, Jinsong Li\*, Li Yang\*, Ling-Ling Chen\*. Screening for functional circular RNAs using the CRISPR-Cas13 system. 2020, xxxxxx


## License
Copyright (C) 2019 YangLab. Licensed GPLv3 for open source use or contact YangLab (yanglab@@picb.ac.cn) for commercial use.
