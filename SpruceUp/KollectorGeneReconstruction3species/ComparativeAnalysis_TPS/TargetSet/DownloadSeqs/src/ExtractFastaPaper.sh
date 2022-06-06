#!/bin/bash

list_ids=/projects/btl/kgagalova/AssembleIndividualTargets/TargetSeqs/DownloadSeqs/list_genbankids.in
out_fastaProt=TPS_pubProt.fa
out_fastaNucl=TPS_pubNucl.fa

#extract prot seqs
while read id
do
efetch -db protein -id $id -format fasta >> $out_fastaProt
done < $list_ids

#convert to accession - genebank nucleotides
while read id
do
efetch -db protein -id $id -format gp | grep "DBSOURCE" | sed 's/accession/?/' | cut -d'?' -f2 >> list_nucleotideids.in
done < $list_ids

#extract metadata
while read id
do
efetch -db protein -id $id -format gp >> metadata_nucleotides.in
done < $list_ids


#extract nucleotides
while read id
do
efetch -db nucleotide -id $id -format fasta >> $out_fastaNucl
done < list_nucleotideids.in

#reshape entries
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <(sed -r 's/\s+//g' $out_fastaProt)  > TPS_pubProtOneline.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <(sed -r 's/\s+//g' $out_fastaNucl)  > TPS_pubNuclOneline.fa

#clean
rm $out_fastaProt $out_fastaNucl
