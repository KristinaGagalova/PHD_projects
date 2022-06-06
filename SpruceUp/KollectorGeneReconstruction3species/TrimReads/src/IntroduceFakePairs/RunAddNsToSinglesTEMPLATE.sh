#!/bin/bash

#add here the path to software
#SBATCH --job-name=AddNsSort
#SBATCH --partition=all
#SBATCH --ntasks=12
#SBATCH --mem=100G

export LC_ALL=C

singles_file=_RePlAcE_

tmp_dirSort=/projects/btl/scratch/kgagalova/Reads/ReadsPG29singles/tmp/$RANDOM
mkdir -p $tmp_dirSort

echo "start time: "
date;

zgrep "^@" $singles_file | grep "/1$\|/2$" > tmp.singles

sed '1~3 s/$/\nNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' <(sed 's/1$/2tmp/;s/2$/1/;s/2tmp/2/' tmp.singles | awk -v n=1 '1; NR % n == 0 {print "+"}' | awk -v n=2 '1; NR % n == 0 {print "FAKEFAKEFAKEFAKEFAKEFAKEFAKE"}') > tmp.reads

nam=$(echo $singles_file | rev | cut -d/ -f1 | rev | sed 's/_R[12]trim.fq.gz_singles.fastq.gz//' )

zcat $singles_file | cat - tmp.reads > tmp.reads2

./deinterleave_fastq.sh < tmp.reads2 ${nam}_R1trimSingl.fastq ${nam}_R2trimSingl.fastq

##sort the reads

cat ${nam}_R1trimSingl.fastq | sed 's/^+H.*/+/' | sed 's/^+M.*/+/' | paste - - - - | sort -T $tmp_dirSort -k1,1 | tr "\t" "\n" | gzip -9c > ${nam}_R1trimSingl_sort.fastq.gz
cat ${nam}_R2trimSingl.fastq | sed 's/^+H.*/+/' | sed 's/^+M.*/+/' | paste - - - - | sort -T $tmp_dirSort -k1,1 | tr "\t" "\n" | gzip -9c > ${nam}_R2trimSingl_sort.fastq.gz

echo ""
echo "---------------------------"
echo "Job successfully finished!"
echo "End time: "
date;

rm tmp.reads tmp.reads2 tmp.singles *.fastq

