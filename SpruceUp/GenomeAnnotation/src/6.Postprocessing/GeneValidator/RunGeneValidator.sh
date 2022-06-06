#!/usr/bin/env bash

#Description: run GeneValidator on batches of proteins

fasta=/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneValidator/Batch1/WSgenome_annotation_0.fa
blast_db=/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneValidator/BlastDBSprot/SelFastaSprot.fa

genevalidator -d $blast_db -n 24 -b /gsc/btl/linuxbrew/Cellar/blast/2.7.1/bin/ -b /home/kgagalova/bin $fasta
