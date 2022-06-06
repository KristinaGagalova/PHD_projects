#!/bin/bash

#############################
######Kristina Gagalova######
#############################
########23 Dec 2017##########
#############################

#Description: extract 2nt upstream first exon and 2nt last exon - this can be used to check if there is a splicing motif which may indicate incomplete sequence

export PATH=$PATH:/gsc/btl/linuxbrew/bin

#fasta sequences
fasta=/projects/spruceup/interior_spruce/PG29/annotation/PG29v3-renamedID_500nt.fa
bed_all=/projects/spruceup/Collaborators_spruce_resources/CombinedResources/Kollector_targetSets/cdhit/cdhit-output-4genes/HCGmaker.bed

#get length seqs
bioawk -c fastx -F $'\t' '{ print $name,length($seq)}' $fasta > PG29v3.chromsizes

#get genes from gff
grep gene $bed_all > HCGmakerGenes.bed

#####################################Splicing motifs
##upstream 2 nt
bedtools flank -i HCGmakerGenes.bed -g PG29v3.chromsizes -l 2 -r 0 -s > HCGgenes.PG29.2ntUpstream.bed
#downstream 2 nt
bedtools flank -i HCGmakerGenes.bed -g PG29v3.chromsizes -l 0 -r 2 -s > HCGgenes.PG29.2ntDownstream.bed

#---------------------------------------------------
#LEGEND:
#bedtools flank - extracts the regions before or after a coordinate
#-l - left
#-r - right
#-s - consider the strand orientation when extracting
#---------------------------------------------------

#extract sequences - exclude the ones with 0 length
bedtools getfasta -s -fi $fasta -bed <(awk '$1 > 0 && $2 > 0 {print $0}' HCGgenes.PG29.2ntUpstream.bed) -fo HCGgenes.PG29.2ntUpstream.fa
bedtools getfasta -s -fi $fasta -bed <(awk '$1 > 0 && $2 > 0 {print $0}' HCGgenes.PG29.2ntDownstream.bed) -fo HCGgenes.PG29.2ntDownstream.fa 

#output stats of the frequencies
grep -v ">" HCGgenes.PG29.2ntUpstream.fa | sort | uniq -c | sort -nrk1,1 >  FreqUpstreamHCGPG29.txt
grep -v ">" HCGgenes.PG29.2ntDownstream.fa | sort | uniq -c | sort -nrk1,1 >  FreqDownstreamHCGPG29.txt

