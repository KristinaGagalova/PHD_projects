#!/bin/bash

#$ -S /bin/bash
#$ -N TrimReads
#$ -q all.q,mpi.q
#$ -pe ncpus 6
#$ -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G
#$ -R y
#$ -j y
#$ -cwd
#$ -V


######################################
#######Kristina Gagalova##############
######################################
##########1 Aug 2017##################
######################################

prinseq=/home/kgagalova/bin/prinseq-lite-0.20.4/prinseq-lite.pl

fq=$1

namfq=$(echo $fq | rev | cut -d/ -f1 | rev | sed 's/.fq.gz//')

zcat $fq | paste - - - - | sort -k1,1 -t " " |tr "\t" "\n" | perl $prinseq -noniupac -ns_max_p 95 -lc_method dust -lc_threshold 80 -trim_qual_right 20 -min_len 25 -fastq stdin -out_good stdout -out_bad null | gzip  > ${namfq}trim.fq.gz
