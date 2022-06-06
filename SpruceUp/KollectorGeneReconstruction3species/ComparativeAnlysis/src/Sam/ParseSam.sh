#!/usr/bin/env bash

###############################
#####Kristina Gagalova#########
###############################
########June 27 2017###########
###############################

#Description: given a sam file (with or without header doesn't matter), parse each line and output info in tabular format
#Output info: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, NM, MD

#Reference here: https://samtools.github.io/hts-specs/SAMv1.pdf
#Other links: http://broadinstitute.github.io/picard/explain-flags.html

SAM=$1

grep -v "^@" $SAM > sam.tmp

grep -o 'NM:[^[:space:]]*' sam.tmp > nm.tmp
grep -o 'MD:[^[:space:]]*' sam.tmp > md.tmp
awk '{print $1,$2,$3,$4,$5,$6}' sam.tmp > Out.results.tmp

paste Out.results.tmp md.tmp nm.tmp > ParsedSam.txt

rm *tmp

