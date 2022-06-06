#!/usr/bin/bash

export PATH=$PATH:/gsc/btl/linuxbrew/bin


gff_annotation=$1
og_all=$2
genome=$3
prefix=$4

while read line
do
nam=$(echo $line | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
echo $nam
grep -F -f <(cat $line | grep ">" | sed 's/>//g' | cut -d"-" -f1) <(cat $gff_annotation | awk '$3 == "gene" {print $0}') > tmp${nam}.gff
if [ -s tmp${nam}.gff ]
then
	paste <(cat tmp${nam}.gff | awk -F$'\t' -v OFS='\t' '{print $1,$4 - 1,$5}') \
	<(cat tmp${nam}.gff | awk -F$'\t' -v OFS='\t' '{print $9}' | cut -d";" -f1 | sed 's/ID=//') \
	<(cat tmp${nam}.gff | awk -F$'\t' -v OFS='\t' '{print $6,$7}') | tr ' ' '\t' > ${prefix}_${nam}.bed  
else
	continue
fi
bedtools getfasta -s -name -fi $genome -bed ${prefix}_${nam}.bed > ${prefix}_${nam}Extract.fa
rm tmp${nam}.gff
done < $og_all
