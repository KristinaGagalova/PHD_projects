#!/bin/bash

##########LIBS
export PATH=$PATH:/gsc/btl/linuxbrew/bin
export AUGUSTUS_CONFIG_PATH=/home/kgagalova/src/augustus.2.5.5/config
export PERL5LIB=$PERL5LIB:/home/kgagalova/localperl/lib/site_perl/5.20.1/Parallel
###requires: faSomeRecords

scripts_augustus=/home/kgagalova/src/augustus.2.5.5
genome=/Genome/with/gene/models/genome.fasta
gff_genes=/file/with/gene/coord/genes.gff
nam=GenesMaker


echo "preparing files for Augustus"
perl $scripts_augustus/scripts/gff2gbSmallDNA.pl $gff_genes $genome 100 ${nam}.gb
${scripts_augustus}/bin/etraining --species=spruce17 ${nam}.gb

