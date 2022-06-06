#!/usr/bin/env bash

############################
#####Kristina Gagalova######
############################

#Description: Run EvidentialGene for the selection of complete transcripts

export PATH=$PATH:/home/kgagalova/src/EvidentialGene/evigene/scripts/prot/
export PATH=$PATH:/home/shammond/src/ncbi-blast-2.2.28+/bin:/gsc/btl/linuxbrew/bin

dataset=$1

tr2aacds.pl -mrnaseq $dataset -NCPU=12 -MAXMEM=100000

