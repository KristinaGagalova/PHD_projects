#!/bin/bash

annotation=/projects/btl/datasets/hsapiens/hg38.p10/annotation/Homo_sapiens.GRCh38.90.gff3
nam=
genome_genes=


#################Get genes
#get coding regions
nam_ann=$(echo $annotation | rev | cut -d"/" -f1 | rev | sed 's/.gff3//;s/.gff//')
awk '$3 == "mRNA" {print $0}' ${annotation} > ${nam_ann}.gff
#add flanking regions
awk -v OFS='\t' '{print $1,$2,$3,$4-100,$5+100,$6,$7,$8,$9}' ${nam_ann}.gff | sort -k1,1 -k4,4n -k5,5n > ${nam_ann}flanking.gff
bedtools merge -i ${nam_ann}flanking.gff > ${nam_ann}flankingMerged.bed
#fasta sequences of the intervals
bedtools getfasta -fi $genome -bed ${nam_ann}flankingMerged.bed -fo ${nam_ann}flankingMerged.fa


###############Create GMAPindex to determine exons of the genes
mkdir -p GMAPindex
gmap_build -D GMAPindex -d ${nam_ann}flankingMerged ${nam_ann}flankingMerged.fa


################align with gmap
gmap -d ${nam} -D ./${nam} -f 2 -O $genome_genes >> ${nam}.gff

#Discard duplicated
paste <(cat ${nam}.gff) <(awk '{print $9}' ${nam}.gff | awk -F  "coverage=" '{print $2}' | cut -d";" -f1) \
	<(awk '{print $9}' ${nam}.gff | awk -F  "identity=" '{print $2}' | cut -d";" -f1)  \
	<(awk '{print $9}' ${nam}.gff | awk -F  "Name=" '{print $2}' | cut -d";" -f1) | awk '$3 == "mRNA"' | awk '$10 >= 90' | sort -k12,10rn -k11rn > ${nam}uniq.gff
#get the couples of alignments that are the highest
awk '{print $1,$12}' ${nam}uniq.gff > list_genes.in
#get exons only
awk '$3 == "exon"' ${nam}.gff > ${nam}exons.gff
#extract the couple of alignments
while read line; do contig=$(echo $line | cut -d" " -f1 ) && gene=$(echo $line | cut -d" " -f2 ) && grep $contig ${nam}exons.gff | grep $gene ; done < list_genes.in >> sel${nam}.gff

#paste ()
