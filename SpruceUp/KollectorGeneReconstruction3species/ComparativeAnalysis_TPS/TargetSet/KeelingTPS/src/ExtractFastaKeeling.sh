#!/bin/bash

list_pep=/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/ComparativeAnalysis_TPS/TargetSet/KeelingTPS/NCBIesearchProt.in
list_nucl=/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/ComparativeAnalysis_TPS/TargetSet/KeelingTPS/NCBIesearchNucl.in
out_fastaProt=TPS_pubProt.fa
out_fastaNucl=TPS_pubNucl.fa

#extract prot seqs
while read id
do
efetch -db protein -id $id -format fasta >> $out_fastaProt
done < $list_pep

#extract nucl seqs
while read id
do
efetch -db nucleotide -id $id -format fasta >> $out_fastaNucl
done < $list_nucl

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <(sed -r 's/\s+//g' $out_fastaProt)  > TPS_KeelingProtOneline.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <(sed -r 's/\s+//g' $out_fastaNucl)  > TPS_KeelingNuclOneline.fa

#clean
rm $out_fastaProt $out_fastaNucl
