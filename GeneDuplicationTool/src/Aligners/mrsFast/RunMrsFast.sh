#!/bin/bash

#######################################
###########Kristina Gagalova###########
#######################################
###############Jun 2018################
#######################################


#Run mrsFast on whole genome, make index and align reads

#https://github.com/sfu-compbio/mrsfast

########IMPORTANT: the genome needs to be hard masked for repeats - otherwise the alignment outputs all the repeats alignments. MrsFasta aligns multimapping reads (plus SNPs)
fasta_file=MyFasta

mrsfast=/home/kgagalova/src/mrsfast/mrsfast

###make index
$mrsfast --index $fasta_file

##Align paired end reads

out_nam=NA12878_S1_L001_R1_001

##############################
## search: indexed fasta file
## min/max = fragment size for pe
## --pe: --seq1 + --seq2
## --seqcomp: compressed input file
## --disable-nohits: disable output of non-mapped reads
###############################

$mrsfast --search $fasta_file --pe --seq1 $pe1 --seq2 $pe2 --seqcomp --disable-nohits \
        --min 300 --max 600 --threads 32 -o ${out_nam}.sam

#Convert to bam and sort, index
samtools view -bS ${out_nam}.sam > ${out_nam}.bam
samtools sort ${out_nam}.bam -o ${out_nam}Sorted.bam
samtools index -o ${out_nam}Sorted.bam

rm ${out_nam}.bam ${out_nam}.sam
