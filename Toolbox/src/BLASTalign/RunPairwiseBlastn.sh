#!/usr/bin/env bash

###############################
#####Kristina Gagalova#########
###############################
########June 14 2017###########
###############################

#Description: given a list of pairs of fasta sequence names, do pairwise alignment and output to unique file
## How to run script: ./RunPairwiseBlastn.sh <list_pair_names> <A1_fasta> <A2_fasta> <---- A1,A2 may be the same is contigs are listed in the same file

LISTA_TARGETS=$1
A1_fasta=$2
A2_fasta=$3

while read target
do
a1=$(cut -d" " -f1 <<< $target)
a2=$(cut -d" " -f2 <<< $target)
samtools faidx $A1_fasta $a1 > a1tmp.fa
samtools faidx $A2_fasta $a2 > a2tmp.fa
blastn -query a1tmp.fa -subject a2tmp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' >> Abyss_alignments_blastn.txt
done < $LISTA_TARGETS

rm *tmp*
