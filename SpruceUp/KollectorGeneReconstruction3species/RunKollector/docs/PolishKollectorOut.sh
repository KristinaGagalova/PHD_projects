#!/bin/bash

##############################
######Kristina Gagalova#######
##############################
#######June 16 2017###########
##############################

#Description: remove intermediate files from Kollector pipeline

# Usage: ./PolishKollectorOut.sh <dir kollector out> 

koll_dir=$1

cd $koll_dir

echo "STARTING FILES DELETION"

echo "intermediate files..."

find -name "allcigarseq.txt" -exec rm {} \;
find -name "allcigars.txt" -exec rm {} \;
find -name "allcoverages.txt" -exec rm {} \;
find -name "allmapped.sam" -exec rm {} \;
find -name "allnames.txt" -exec rm {} \;
find -name "allseq.txt" -exec rm {} \;
find -name "alltranscripts.txt" -exec rm {} \;
find -name "failedtranscripts.fa" -exec rm {} \;
find -name "hitlist.txt" -exec rm {} \;
find -name "*.transcript2assembly.sam" -exec rm {} \;

echo "fastq files..."

find -name "*recruited_pe.fastq" -exec rm {} \;

echo "Abyss files..."

find -name "*.dot" -exec rm {} \;
find -name "*.dot1" -exec rm {} \;
find -name "*.fai" -exec rm {} \;
find -name "*amb" -exec rm {} \;
find -name "*ann" -exec rm {} \;
find -name "*bwt" -exec rm {} \;
find -name "*pac" -exec rm {} \;
find -name "*sa" -exec rm {} \;

echo "GMAP indexes..."

find -type d -name "*.gmap" -print0 | xargs -0 rm -r

echo "Remove empty files..."

find . -size 0 -exec rm {} \;

echo "Kollector cleaning finished!!"
