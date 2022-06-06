#!/bin/bash

list_ids=/projects/btl/kgagalova/AssembleIndividualTargets/TargetSeqs/SerachNCBI/NCBIesearchProt.in
out_fastaProt=TPS_pubNCBIProt.fa
out_fastaNucl=TPS_pubNCBINucl.fa

#extract prot seqs
while read id
do
efetch -db protein -id $id -format fasta >> $out_fastaProt
done < $list_ids

#extract etadata
while read id
do
efetch -db protein -id $id -format gp  >> metadata_proteins.in
done < $list_ids

#convert to accession - genebank nucleotides
grep "DBSOURCE" metadata_proteins.in | awk '!/UniProtKB/ && !/pdb/' | sed 's/accession/?/' | cut -d'?' -f2 >> list_nucleotidesGenBank.in
sed -e '/^\s*$/d' list_nucleotidesGenBank.in 

#convert to accession - genebank nucleotides
grep "xrefs:" metadata_proteins.in | sed 's/xrefs://' |  cut -d',' -f1 | sed 's/^[ \t]*//' >> list_nucleotidesUniProt.in

#extract nucleotides
while read id
do
efetch -db nucleotide -id $id -format fasta >> $out_fastaNucl
done < list_nucleotidesGenBank.in

#extract nucleotides uniprot
while read id
do
efetch -db nucleotide -id $id -format fasta >> $out_fastaNucl
done < list_nucleotidesUniProt.in

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <(sed -r 's/\s+/_/g' $out_fastaProt)  > TPS_pubNCBIProtOneline.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <(sed -r 's/\s+/_/g' $out_fastaNucl)  > TPS_pubNCBINuclOneline.fa

rm $out_fastaProt $out_fastaNucl
