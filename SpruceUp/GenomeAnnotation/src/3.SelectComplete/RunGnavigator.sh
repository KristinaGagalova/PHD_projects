#!/bin/bash

#####################################
########Kristina Gagalova############
#####################################
############22 Feb 2018##############
#####################################

#Description: Run Gnavigator to extract contgs with complete 

export PATH=/home/shammond/src/gmap-20171115-build/bin:/projects/spruceup_scratch/dev/AnalysisDownstream/Gnavigator/gnavigator:$PATH

genome=$1
#in case of existing gmap index
path_gmapIndex=$2
prefix_gmapIndex=$3

gnavigator=/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Step0annoatation/SelectContigs/SrcGnavigatorLowerIdent/gnavigator.py

cDNAass=/projects/spruceup/pglauca/WS77111/annotation/genome-annotation/WS77111v2/Maker/MakerDependencies/GCAT_WS-3.3.cluseq_completeFinal.fa

$gnavigator -p ${prefix_gmapIndex}_out -d $path_gmapIndex -n $prefix_gmapIndex -t 6 $cDNAass $genome

