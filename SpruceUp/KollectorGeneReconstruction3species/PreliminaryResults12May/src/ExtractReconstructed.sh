#!/bin/bash

################################
######Kristina Gagalova#########
################################
##########1 May#################
################################

#Description: Extract the results from assembled transcripts

find /extscratch/btl/kgagalova/EasterRun -name "allcoveredmapped.sam" -exec awk -v OFS='\t' '{ if ($1>=0.90) print $1,$2,$4;}' {} > Out.results.WS77111 \;
sort -rnk1 Out.results.WS77111 | awk '!x[$2]++' > Out.results.WS77111.top

find /extscratch/btl/bullfrog/genome/tga/spruce/myseq*.fa.r112 -name "allcoveredmapped.sam" -exec awk -v OFS='\t' '{ if ($1>=0.90) print $1,$2,$4;}' {} > Out.results.PG29  \;
sort -rnk1 Out.results.PG29 | awk '!x[$2]++' > Out.results.PG29.top
