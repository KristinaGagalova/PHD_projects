#!/usr/bin/env bash

###############################
#####Kristina Gagalova#########
###############################
########June 14 2017###########
###############################

#Description: given a list of pairs of fasta headers, do pairwise alignment and output to unique file. Align longer sequence to the shorter. 
## How to run script: ./RunPairwiseBlastn.sh <list_pair_names> <A1_fasta> <A2_fasta> <---- A1,A2 may be the same is contigs are listed in the same file

LISTA_TARGETS=$1
A1_fasta=$2
A2_fasta=$3

while read target
do
a1=$(cut -d" " -f1 <<< $target)
a2=$(cut -d" " -f2 <<< $target)
samtools faidx $A1_fasta $a1 > a1tmp.fa
samtools faidx $A2_fasta $a2 > a2tmp.fa
lena1=$(wc -c a1tmp.fa | awk '{print $1}')
lena2=$(wc -c a2tmp.fa | awk '{print $1}')
if (("$lena1" > "$lena2"))
then
echo "1" >> Order_align.txt
blastn -query a1tmp.fa -subject a2tmp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' >> Shuf_alignments_blastn.txt
else
echo "2" >> Order_align.txt
blastn -query a2tmp.fa -subject a1tmp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' >> Shuf_alignments_blastn.txt
fi
done < $LISTA_TARGETS

paste Shuf_alignments_blastn.txt Order_align.txt > Shuf_alignments_blastn2.txt
rm *tmp* Order_align.txt Shuf_alignments_blastn.txt

