#!/usr/bin/bash

bbmap=/home/kgagalova/src/bbmap/fuse.sh
base=/projects/spruceup_scratch/dev/SprucePaper2018/GeneFamilies/miBF/ExtractAndProcessGeneFamilies/*/GeneSeq

list_og=$1

while read line
do
cat $base/*${line}Extract.fa > tmp_file${line}.fasta
sed -i "s#^>.*#>$line#" tmp_file${line}.fasta
$bbmap -in=tmp_file${line}.fasta -out=tmp_file_out${line}.fasta -pad=100
cat tmp_file_out${line}.fasta >> GeneFamilies3sp.fasta
rm tmp_file_out*
done < $list_og
