#!/usr/bin/env bash

fasta_batch=/projects/spruceup_scratch/dev/MakerKollectorContigs/WS77111/RepeatMasker/BatchesFasta/list_batches.in

while read fasta
do
	nam=$(echo $fasta | rev | cut -d/ -f1 | rev | cut -f1 -d. )
	mkdir -p $nam && cd $nam
	cp -p ../RunRepeatMaskerTEMPLATE RunRepeatMasker_$nam
	sed -i -e "s#_RePlaCe_#${fasta}#g" RunRepeatMasker_$nam #file path
	sed -i -e "s#_RePlaCeNam_#${nam}#g" RunRepeatMasker_$nam #output name
	sbatch ./RunRepeatMasker_$nam
	cd ..
done < $fasta_batch
