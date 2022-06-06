#!/bin/bash

#add here the path to software
#SBATCH --job-name=AddNsSort
#SBATCH --partition=all
#SBATCH --ntasks=12
#SBATCH --mem=100G

fastq=_RePlAcE_

nam=$(echo $fastq | rev | cut -d/ -f1 | rev | sed 's/trimSingl_sort.fastq.gz//' )

entries=$(zcat $fastq | wc -l | cut -f1 -d" ")
reads=`echo $entries/4 | bc -l`
div=$(zcat $fastq | grep "^+$" | wc -l)
printf "%b\t%0.1f\t%0.1f\t%0.1f\n" $nam $entries $reads $div >> LinesReadsDiv.txt
sort -k,k1 LinesReadsDiv.txt > LinesReadsDivSort.txt
rm LinesReadsDiv.txt
