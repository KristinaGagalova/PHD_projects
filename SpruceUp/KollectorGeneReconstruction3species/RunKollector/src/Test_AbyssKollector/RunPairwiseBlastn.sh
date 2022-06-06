#!/usr/bin/env bash

###############################
#####Kristina Gagalova#########
###############################
########June 14 2017###########
###############################

#Description: given a list of pairs of fasta headers, do pairwise alignment and output to unique file



LISTA_TARGETS=/projects/btl/kgagalova/compare_kollector/MergedCommonPolished.txt
A1_fasta=/projects/btl/kgagalova/compare_kollector/contigs_succeed_1.fa
A2_fasta=/projects/btl/kgagalova/compare_kollector/contigs_succeed_2.fa

while read target
do
a1=$(cut -d" " -f3 <<< $target)
a2=$(cut -d" " -f5 <<< $target)
samtools faidx $A1_fasta $a1 > a1tmp.fa
samtools faidx $A2_fasta $a2 > a2tmp.fa
blastn -query a1tmp.fa -subject a2tmp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' >> Abyss_alignments_blastn.txt
done < $LISTA_TARGETS

rm *tmp*
