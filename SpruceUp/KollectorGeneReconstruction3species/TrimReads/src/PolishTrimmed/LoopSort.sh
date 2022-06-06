#!/bin/bash

fastq_list=/extscratch/spruceup/interior_spruce/PG29/data/reads/TrimmesFastq/SolitaryReads/all_fastqTrimmed.in

while read fastq
do
nam=$(echo $fastq | rev | cut -d/ -f1 | rev | sed 's/.gz_pairs_R1.fastq.gz//' |sed 's/.gz_pairs_R2.fastq.gz//')
mkdir -p $nam && cd $nam
cp ../SortFastqTEMPLATE.sh SortFastq$nam.sh
sed -i -e "s#_RePlAcE_#${fastq}#g" SortFastq$nam.sh
qsub ./SortFastq$nam.sh
cd ..
done < $fastq_list

