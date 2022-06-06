#!/bin/bash

#$ -S /bin/bash
#$ -N SortFiles
#$ -q all.q,mpi.q
#$ -pe ncpus 4
#$ -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G
#$ -R y
#$ -j y

export LC_ALL=C

f1=_RePlAcE_
tmp_dirSort=/extscratch/spruceup/interior_spruce/PG29/data/reads/TrimmesFastq/SolitaryReads/tmp/$RANDOM

mkdir -p $tmp_dirSort

zcat $f1 | paste - - - - | sort -T $tmp_dirSort -k1,1 | tr "\t" "\n" | gzip -9c > ${f1}_sort.gz
