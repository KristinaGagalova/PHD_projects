#!/bin/bash

#$ -S /bin/bash
#$ -N gmapAlign
#$ -q all.q,mpi.q
#$ -pe ncpus 1
#$ -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G
#$ -R y
#$ -j y

################################
######Kristina Gagalova#########
################################
##########8 May#################
################################

#Description: for a given list of seed and reconstructed genes, output GMAP alignment in gff format

export PATH+=:~/.linuxbrew/bin

namegmap=$3
contigs=$2
seed=$1
WORK_DIR=gmap_align

mkdir -p $WORK_DIR
gmap_build -D $WORK_DIR -d $namegmap $contigs
gmap -D $WORK_DIR -d $namegmap -f gff3_gene -O -B 4 $seed >> $namegmap.gff
gmap -D $WORK_DIR -d $namegmap -f psl -O -B 4  $seed >> $namegmap.psl
gmap -D $WORK_DIR -d $namegmap -f gff3_gene -n 1 -A -O -B 4  $seed >> $namegmap.align
                      
 
