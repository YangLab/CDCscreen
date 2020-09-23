#!/usr/bin/bash
if [ ! $# = 8 ]; then
    echo "Usage: `basename $0` D1_rep1_R1.fq D1_rep1_R2.fq D1_rep2_R1.fq D1_rep2_R2.fq D30_rep1_R1.fq D30_rep1_R2.fq D30_rep2_R1.fq D30_rep2_R2.fq";
    exit 0;
fi

# 0. bowtie-build index gRNA reference sequences (bowtie version 1.1.2)
mkdir 00_ref_gRNA_seq
cd 00_ref_gRNA_seq
cp ../ref_gRNA_seq.fa .
bowtie-build ref_gRNA_seq.fa ref_gRNA_seq

cd ..


# 1.1 ln -s D1 rep1 and rep2, and D30 rep1 and rep2 fastq(or fastq.gz)
mkdir 01_gRNA_DNA_seq
cd 01_gRNA_DNA_seq
ln -s $1 01_D1_rep1_R1.fastq.gz
ln -s $2 01_D1_rep1_R2.fastq.gz
ln -s $3 01_D1_rep2_R1.fastq.gz
ln -s $4 01_D1_rep2_R2.fastq.gz

ln -s $5 01_D30_rep1_R1.fastq.gz
ln -s $6 01_D30_rep1_R2.fastq.gz
ln -s $7 01_D30_rep2_R1.fastq.gz
ln -s $8 01_D30_rep2_R2.fastq.gz

# 1.2 remove adapter sequences (cutadapt version 1.18)
# For R1 3' (TTTTTTAAGCTTGGCGTAACTAGATCT), 5' (CCCTACCAACTGGTCGGGGTTTGAAAC)
# For R2 3' (GTTTCAAACCCC or GTTTCAAACCCCGACCAGTTGGTAGGG, which one was selected by your R2 sequences), 5' (AGATCTAGTTACGCCAAGCTTAAAAAA)
# D1 rep1 R1 rm3p
cutadapt -a TTTTTTAAGCTTGGCGTAACTAGATCT -m 15 -o 02_D1_rep1_R1_rm3p.fq 01_D1_rep1_R1.fastq.gz
# D1 rep1 R1 rm5p
cutadapt -g CCCTACCAACTGGTCGGGGTTTGAAAC -m 15 -o 02_D1_rep1_R1_trimmed.fq 02_D1_rep1_R1_rm3p.fq
# D1 rep1 R2 rm3p
cutadapt -a GTTTCAAACCCC -m 15 -o 02_D1_rep1_R2_rm3p.fq 01_D1_rep1_R2.fastq.gz
## cutadapt -a GTTTCAAACCCCGACCAGTTGGTAGGG -m 15 -o 02_D1_rep1_R2_rm3p.fq 01_D1_rep1_R2.fastq.gz
# D1 rep1 R2 rm5p
cutadapt -g AGATCTAGTTACGCCAAGCTTAAAAAA -m 15 -o 02_D1_rep1_R2_trimmed.fq 02_D1_rep1_R2_rm3p.fq

# D1 rep2 R1 rm3p
cutadapt -a TTTTTTAAGCTTGGCGTAACTAGATCT -m 15 -o 02_D1_rep2_R1_rm3p.fq 01_D1_rep2_R1.fastq.gz
# D1 rep2 R1 rm5p
cutadapt -g CCCTACCAACTGGTCGGGGTTTGAAAC -m 15 -o 02_D1_rep2_R1_trimmed.fq 02_D1_rep2_R1_rm3p.fq
# D1 rep2 R2 rm3p
cutadapt -a GTTTCAAACCCC -m 15 -o 02_D1_rep2_R2_rm3p.fq 01_D1_rep2_R2.fastq.gz
## cutadapt -a GTTTCAAACCCCGACCAGTTGGTAGGG -m 15 -o 02_D1_rep2_R2_rm3p.fq 01_D1_rep2_R2.fastq.gz
# D1 rep2 R2 rm5p
cutadapt -g AGATCTAGTTACGCCAAGCTTAAAAAA -m 15 -o 02_D1_rep2_R2_trimmed.fq 02_D1_rep2_R2_rm3p.fq

# D30 rep1 R1 rm3p
cutadapt -a TTTTTTAAGCTTGGCGTAACTAGATCT -m 15 -o 02_D30_rep1_R1_rm3p.fq 01_D30_rep1_R1.fastq.gz
# D30 rep1 R1 rm5p
cutadapt -g CCCTACCAACTGGTCGGGGTTTGAAAC -m 15 -o 02_D30_rep1_R1_trimmed.fq 02_D30_rep1_R1_rm3p.fq
# D30 rep1 R2 rm3p
cutadapt -a GTTTCAAACCCC -m 15 -o 02_D30_rep1_R2_rm3p.fq 01_D30_rep1_R2.fastq.gz
## cutadapt -a GTTTCAAACCCCGACCAGTTGGTAGGG -m 15 -o 02_D30_rep1_R2_rm3p.fq 01_D30_rep1_R2.fastq.gz
# D30 rep1 R2 rm5p
cutadapt -g AGATCTAGTTACGCCAAGCTTAAAAAA -m 15 -o 02_D30_rep1_R2_trimmed.fq 02_D30_rep1_R2_rm3p.fq

# D30 rep2 R1 rm3p
cutadapt -a TTTTTTAAGCTTGGCGTAACTAGATCT -m 15 -o 02_D30_rep2_R1_rm3p.fq 01_D30_rep2_R1.fastq.gz
# D30 rep2 R1 rm5p
cutadapt -g CCCTACCAACTGGTCGGGGTTTGAAAC -m 15 -o 02_D30_rep2_R1_trimmed.fq 02_D30_rep2_R1_rm3p.fq
# D30 rep2 R2 rm3p
cutadapt -a GTTTCAAACCCC -m 15 -o 02_D30_rep2_R2_rm3p.fq 01_D30_rep2_R2.fastq.gz
## cutadapt -a GTTTCAAACCCCGACCAGTTGGTAGGG -m 15 -o 02_D30_rep2_R2_rm3p.fq 01_D30_rep2_R2.fastq.gz
# D30 rep2 R2 rm5p
cutadapt -g AGATCTAGTTACGCCAAGCTTAAAAAA -m 15 -o 02_D30_rep2_R2_trimmed.fq 02_D30_rep2_R2_rm3p.fq

cd ..


# 2.1 align to gRNA reference seuqnces with bowtie
mkdir 02_bowtie_align
cd 02_bowtie_align
ls ../01_gRNA_DNA_seq/ |grep "trimmed" |awk -F"." '{print "bowtie ../00_ref_gRNA_seq/ref_gRNA_seq ../01_gRNA_DNA_seq/"$0" -S 01_"$1".sam -p 12 -v 3 -m 1 -k 1 2>log_01_"$1".txt"}' |sed 's/01_02/01/g' |sh

# 2.2 sam to bam (samtools Version: 1.9)
ls |grep "sam" |awk -F"." '{print "samtools view -bh -F 4 "$1".sam |samtools sort -o 02_"$1".bam"}' |sed 's/02_01/02/g' |sh

# 2.3 select read ID
ls |grep "sam" |awk -F"." '{print "samtools view 02_"$1".bam |awk '\''{print $3""\"\\t""\"$1}'\'' > 03_"$1".txt"}' |sed 's/02_01/02/g' |sed 's/03_01/03/g' |sh

# 2.4 remove duplication reads
ls |grep "sam" |awk -F"." '{print "sort -u 03_"$1".txt > 04_"$1".txt"}' |sed 's/03_01/03/g' |sed 's/04_01/04/g' |sh

# 2.5 count sgRNA mapped reads
ls |grep "sam" |awk -F"." '{print "cut -f 1 04_"$1".txt |sort |uniq -c |awk '\''{print $2""\"\\t""\"$1}'\'' > 05_"$1".txt"}' |sed 's/04_01/04/g' |sed 's/05_01/05/g' |sh

cd ..


# 3.1 calculate total gRNA reads
mkdir 03_gRNA_reads
cd 03_gRNA_reads
cat ../02_bowtie_align/05_D1_rep1_R1_trimmed.txt ../02_bowtie_align/05_D1_rep1_R2_trimmed.txt |sort -k1,1n |awk '{a[$1]+=$2;b[$1]++}END{for(n in a)print n"\t"a[n]}' |sort -k1,1n > 01_D1_rep1.txt
cat ../02_bowtie_align/05_D1_rep2_R1_trimmed.txt ../02_bowtie_align/05_D1_rep2_R2_trimmed.txt |sort -k1,1n |awk '{a[$1]+=$2;b[$1]++}END{for(n in a)print n"\t"a[n]}' |sort -k1,1n > 01_D1_rep2.txt
cat ../02_bowtie_align/05_D30_rep1_R1_trimmed.txt ../02_bowtie_align/05_D30_rep1_R2_trimmed.txt |sort -k1,1n |awk '{a[$1]+=$2;b[$1]++}END{for(n in a)print n"\t"a[n]}' |sort -k1,1n > 01_D30_rep1.txt
cat ../02_bowtie_align/05_D30_rep2_R1_trimmed.txt ../02_bowtie_align/05_D30_rep2_R2_trimmed.txt |sort -k1,1n |awk '{a[$1]+=$2;b[$1]++}END{for(n in a)print n"\t"a[n]}' |sort -k1,1n > 01_D30_rep2.txt

# 3.2 combine gRNA reads
perl ../combine_diff.pl -in "01_D1_rep1.txt,01_D1_rep2.txt,01_D30_rep1.txt,01_D30_rep2.txt" -out 02_mapped_reads.tmp -suff bed -fill none
tail -n +2 02_mapped_reads.tmp |sed 's/none/0\t0/g' |sort -k1,1n > 02_mapped_reads.txt
grep "gRNA" 02_mapped_reads.txt > 02_gRNA_reads.txt

# 3.3 normalized by total gRNA reads
perl -alne 'next if $F[0]=~/#/; next if $F[0] eq "Geneid"; $rpk=0; $rpk=$F[1] if $F[1]>0; $totalrpk += $rpk; push @gene,$F[0];push @generpk,$rpk; END{ $scale=$totalrpk/100000; for($i=0;$i<=$#gene;$i++){$,="\t";print $gene[$i],$generpk[$i]/$scale+1;} }' 02_gRNA_reads.txt > 03_gRNA_D0_rep1.txt
perl -alne 'next if $F[0]=~/#/; next if $F[0] eq "Geneid"; $rpk=0; $rpk=$F[2] if $F[2]>0; $totalrpk += $rpk; push @gene,$F[0];push @generpk,$rpk; END{ $scale=$totalrpk/100000; for($i=0;$i<=$#gene;$i++){$,="\t";print $gene[$i],$generpk[$i]/$scale+1;} }' 02_gRNA_reads.txt > 03_gRNA_D0_rep2.txt
perl -alne 'next if $F[0]=~/#/; next if $F[0] eq "Geneid"; $rpk=0; $rpk=$F[3] if $F[3]>0; $totalrpk += $rpk; push @gene,$F[0];push @generpk,$rpk; END{ $scale=$totalrpk/100000; for($i=0;$i<=$#gene;$i++){$,="\t";print $gene[$i],$generpk[$i]/$scale+1;} }' 02_gRNA_reads.txt > 03_gRNA_D30_rep1.txt
perl -alne 'next if $F[0]=~/#/; next if $F[0] eq "Geneid"; $rpk=0; $rpk=$F[4] if $F[4]>0; $totalrpk += $rpk; push @gene,$F[0];push @generpk,$rpk; END{ $scale=$totalrpk/100000; for($i=0;$i<=$#gene;$i++){$,="\t";print $gene[$i],$generpk[$i]/$scale+1;} }' 02_gRNA_reads.txt > 03_gRNA_D30_rep2.txt
paste 03_gRNA_D0_rep1.txt 03_gRNA_D0_rep2.txt 03_gRNA_D30_rep1.txt 03_gRNA_D30_rep2.txt |cut -f 1,2,4,6,8 > 03_normalized_gRNA_reads.txt

# 3.4 mean of normalized reads in rep1 and rep2
awk '{print $0"\t"($2+$3)/2"\t"($4+$5)/2}' 03_normalized_gRNA_reads.txt > 04_normalized_gRNA_reads_mean.txt

# 3.5 Fold change of rep1, rep2, and mean of rep1 + rep2
awk '{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,($4+0.01)/($2+0.01),($5+0.01)/($3+0.01),($7+0.01)/($6+0.01)}' 04_normalized_gRNA_reads_mean.txt |sort -k1,1 > 05_normalized_gRNA_FC.txt

cd ..


# 4. MAGeCK
mkdir 04_MAGeCK
cd 04_MAGeCK

# 4.1 normalized reads
cut -f 1,6,7 ../03_gRNA_reads/04_normalized_gRNA_reads_mean.txt > 01_gRNA_normalized_mean.txt

# 4.2 generate input data
awk -F"|" '{print $1"|"$2"\t"$0}' 01_gRNA_normalized_mean.txt |awk '{print $2"\t"$1"\t"$3"\t"$4}' > 02_all_count.txt
awk '{print "gRNA\tgene\tD1\tD30"}' 02_all_count.txt |sed -n '1p' > 02_all_count.header

# 4.3 MAGeCK input gRNA normalized reads
cat 02_all_count.header 02_all_count.txt > 03_gRNA_all_MAGeCK.txt

# 4.4 run MAGeCK
mageck test -k 03_gRNA_all_MAGeCK.txt -t D30 -c D1 -n 04_mageck --keep-tmp

cd ..


# 5. CDCscreen
mkdir 05_CDCscreen
cd 05_CDCscreen
# 5.1 circRNA expression
cp ../FPBcirc.txt 01_FPBcirc.txt

# 5.2 circRNA permutation test P value
cp ../04_MAGeCK/04_mageck.gene_summary.txt 02_gene_summary.txt
sed '1d' 02_gene_summary.txt |cut -f 1,4 > 02_negP.txt

# 5.3 Fold change
sed 's/|gRNA/\tgRNA/g' ../03_gRNA_reads/05_normalized_gRNA_FC.txt |awk '{if($11<1)print $0}' |awk '{a[$1]+=$11;b[$1]++}END{for(n in a)print n"\t"(a[n]/b[n])}' > 03_FC.txt

# 5.4 joing P value, Fold change and FPBcirc
perl ../join_ID.pl 03_FC.txt 02_negP.txt 1 1 |cut -f 1-2,4 > 04_join_FC_negP.txt
perl ../join_ID.pl 04_join_FC_negP.txt 01_FPBcirc.txt 1 1 |awk '{if($5>0)print $1"\t"$2"\t"$3}' > 04_join_FC_negP_FPBcirc.txt

# 5.5 Rscript CDCscreen (R version 3.5.1)
Rscript ../CDCscreen_scale.R

cd ..

