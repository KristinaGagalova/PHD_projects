#!/bin/bash

singles_list=/extscratch/spruceup/interior_spruce/PG29/data/reads/TrimmesFastq/SolitaryReads/SinglesList.in

while read fastq
do
nam=$(echo $fastq | rev | cut -d/ -f1 | rev | sed 's/_R1trim.fq.gz_singles.fastq.gz//')
mkdir -p $nam && cd $nam
cp ../deinterleave_fastq.sh ./ && cp ../RunAddNsToSinglesTEMPLATE.sh RunAddNsToSingles$nam.sh
sed -i -e "s#_RePlAcE_#${fastq}#g" RunAddNsToSingles$nam.sh
qsub ./RunAddNsToSingles$nam.sh
cd ..
done < $singles_list

