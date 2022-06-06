#!/usr/bin/env bash

########################################
#########Kristina Gagalova##############
############Feb 03 2017#################
########################################

# Description: given a fasta file with nucleotide sequences, blast all against all
# BLAST documentation: https://www.ncbi.nlm.nih.gov/books/NBK279675/

#############################################################
#------------------------------------------------------------
#blast program version
BLAST=/home/shammond/src/ncbi-blast-2.2.31+/bin

#seq file that need to be blasted
SEQS=/projects/spruceup/pglauca/WS77111/assemblies/kollector/target-sequences/evaluation/cdhit-output-1.fa
SEQS=/projects/spruceup/pglauca/WS77111/assemblies/kollector/target-sequences/evaluation/cdhit-output-2.fa
SEQS=/projects/spruceup/pglauca/WS77111/assemblies/kollector/target-sequences/evaluation/cdhit-output-3.fa
SEQS=/projects/spruceup/pglauca/WS77111/assemblies/kollector/target-sequences/evaluation/cdhit-output-3-500plus.fa

#out names
BASE_OUT=cd_hit_3_500plus
OUT_DB=$BASE_OUT.db
BLAST_OUT=$BASE_OUT.allvsall.tsv

##############################################################
#     Program here
##############################################################

#-------------------------------------------------------------
#create DB from the file
#-------------------------------------------------------------

$BLAST/makeblastdb -in $SEQS -dbtype nucl -out $OUT_DB > logDB

#-------------------------------------------------------------
#blast sequences with BLASTN
#-------------------------------------------------------------

$BLAST/blastn -db $OUT_DB -query $SEQS -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' -out $BLAST_OUT -num_threads 4
 
