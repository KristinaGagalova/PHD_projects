#!/usr/bin/bash

#The names in the concatenated file nee to be different
#list_RNAseq.in is a list of fasta files

#in this case the name is parsed from a directory name, it an ba also the name of the fasta files. 
#parse from here: nam=$(echo $line | rev | cut -d"/" -f3 | rev )

while read line; do nam=$(echo $line | rev | cut -d"/" -f3 | rev ) && sed "s/>/>${nam}_/g" $line  > ${nam}.fasta; done < list_RNAseq.in
cat *fasta > RNAseqAll_pooled.fa
rm *fasta
