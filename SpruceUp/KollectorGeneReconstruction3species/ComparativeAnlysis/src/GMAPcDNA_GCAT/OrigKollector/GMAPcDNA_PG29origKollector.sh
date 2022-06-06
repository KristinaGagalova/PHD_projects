#!/bin/bash 

#Run on hpce706
#PG29v4

export PATH=$PATH:/gsc/btl/linuxbrew/bin

GMAP_index=/projects/spruceup/interior_spruce/PG29/assemblies/PG29-v4/GMAPindex/
query=/projects/spruceup/Collaborators_spruce_resources/CombinedResources/Kollector_targetSets/cdhit/cdhit-output-4noSpace.fa

gmapl -D $GMAP_index -d PG29v4 -S -O -B 4 -t 4 $query >> PG29.transcript2assembly.align
gmapl -D $GMAP_index -d PG29v4 -f samse -O -B 4 -t 4 $query  >> PG29.transcript2assembly.sam 

#parse align GMAP
paste <(grep "Percent identity:\|^>" PG29.transcript2assembly.align | awk '{$1=$1};1' |awk '/^>/{s=$0;next}{print s$0}' | sed 's/Percent identity://;s/>//' | cut -d" " -f1,2) <(grep "Coverage:\|^>" PG29.transcript2assembly.align |awk '{$1=$1};1' |awk '/^>/{s=$0;next}{print $0}'| sed 's/Coverage://' | cut -d" " -f2) > CoverageIdentityPG29targets.txt

