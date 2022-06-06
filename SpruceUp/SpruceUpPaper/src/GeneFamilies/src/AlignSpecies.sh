#!/usr/bin/bash

########################################
#########Kristina Gagalova##############
############Mar 21 2017#################
########################################

# Description: given a fasta file with Kollector seed sequences, Kollector output against it
# BLAST documentation: https://www.ncbi.nlm.nih.gov/books/NBK279675/

#############################################################
#------------------------------------------------------------

#blast program version
BLAST=/gsc/btl/linuxbrew/Cellar/blast/2.2.31/bin

list_og=$1
fasta_firstSp=$2
fasta_secondSp=$3
firstSpPrefix=$4
secondSpPrefix=$5

#out names
OUT_DB=${firstSpPrefix}_vs_${secondSpPrefix}
BASE_OUT=/projects/spruceup_scratch/dev/SprucePaper2018/GeneFamilies/OthofinderFinalRuns/AllConifers_noDF1/Results_Out_6Species_pine_spruce/Orthogroups/ExtractGenes_singleCopy/3speciesSpruce/${OUT_DB}
BLAST_OUT=$OUT_DB.allvsseed.tsv

##############################################################
#     Program here
##############################################################

#-------------------------------------------------------------
#create DB: SEED
#-------------------------------------------------------------
mkdir -p $BASE_OUT
cd $BASE_OUT
$BLAST/makeblastdb -in $fasta_firstSp -dbtype prot -out $OUT_DB

#-------------------------------------------------------------
#blast SEQ sequences with BLASTN
#-------------------------------------------------------------

$BLAST/blastp -db $OUT_DB -query $fasta_secondSp -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' -out seq.${BLAST_OUT} -num_threads 24

#filter data set
sort -k1,1 -k2,2 -k13,13nr -k3,3nr seq.${BLAST_OUT} > seq.sort_${BLAST_OUT}
#awk '$11 < 1e-5 {print $0}' seq.sort_${BLAST_OUT} > seq.sort2_${BLAST_OUT}
awk -F'\t' '!x[$1,$2]++' seq.sort_${BLAST_OUT} > seq.sortSel_${BLAST_OUT}

awk -F'\t' 'NR==FNR{e[$1,$2]=1;next};e[$1,$2]' <(awk -F'\t' -v OFS='\t' '{print $2,$3}' $list_og) seq.sortSel_${BLAST_OUT} > seq.sortOut_${BLAST_OUT}

rm *.phr *.pin *.psq

