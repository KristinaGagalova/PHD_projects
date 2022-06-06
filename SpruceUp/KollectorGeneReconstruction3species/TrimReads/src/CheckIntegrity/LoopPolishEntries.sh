#!/bin/bash

fastq_list=/projects/spruceup/scratch/dev/kollector_comparativeGen/Reads/TrimmedReads/PG29/SolitaryReads/list_solitary.in

cd /projects/spruceup/scratch/dev/kollector_comparativeGen/Reads/TrimmedReads/PG29/SolitaryReads/IntegrityCheck

while read fastq
do
nam=$(echo $fastq | rev | cut -d/ -f1 | rev | sed 's/trimSingl_sort.fastq.gz//')
cp ./CheckIntegrityTEMPLATE.sh CheckIntegrity$nam.sh
sed -i -e "s#_RePlAcE_#${fastq}#g" CheckIntegrity$nam.sh
sleep 20
sbatch ./CheckIntegrity$nam.sh
done < $fastq_list
