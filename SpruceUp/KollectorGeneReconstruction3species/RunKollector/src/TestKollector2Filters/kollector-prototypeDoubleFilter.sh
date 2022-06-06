#!/bin/bash

#$ -S /bin/bash
#$ -N KollectorPrototype
#$ -q all.q,mpi.q
#$ -pe ncpus 6
#$ -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G,excl=true
#$ -R y
#$ -j y

################################################
##########Kristina Gagalova#####################
################################################
##############24 Jul 2017#######################
################################################

#Description: 2 step filtering -B, filtered reads

#add here the path to software
PATH+=:/home/kgagalova/bin/biobloom_develop2working/bin
PATH+=:/home/kgagalova/src/kollector_develop
PATH+=:/home/shammond/src/abyss-2.0.1-build/k256/bin
PATH+=:/usr/mpi/gcc/openmpi-1.4.1/bin
PATH+=:~/.linuxbrew/bin

UNITIGS=/extscratch/btl/kgagalova/trial_kollectorAbyssSparamNorm/trial_kollector.abyss/trial_kollector-1.fa
TARGETS=/extscratch/btl/kgagalova/trial_kollectorAbyssSparam/seed_noSpaces/cdhit-output-4_1.fasta
READS_B0=/extscratch/btl/kgagalova/trial_kollectorAbyssSparam/trial_kollector.recruited_pe.fastq
READS_B1_R1=/extscratch/btl/kgagalova/trial_kollectorAbyssSparamNormTEST/RecrutedReads/read_1.fastq
READS_B1_R2=/extscratch/btl/kgagalova/trial_kollectorAbyssSparamNormTEST/RecrutedReads/read_2.fastq

base_dir=/extscratch/btl/kgagalova/trial_kollectorImproved/FilteredReadsDoubleFilter
abyss_dir=kollector.abyss
rep_fil=/extscratch/spruceup/interior_spruce/PG29/assemblies/kollector/PG29alt/RepFilters/dsk25_100.bf
nam_bbt="PrototypeB0"
K=25
k=109
#1300 X 1000 = 1300000 -> ~1500000
n0=1500000
n1=8000000
s=0.7
j=6
r=0.9

cd $base_dir

biobloommaker -i -k $K -p $nam_bbt -f 0.001 -t $j -n $n0 -s $rep_fil $TARGETS 

biobloomcategorizer -t $j -d $nam_bbt -f ${nam_bbt}.bf -s $s $UNITIGS >> ${nam_bbt}.fa

cat ${nam_bbt}.fa $TARGETS > ${nam_bbt}tar.fa

biobloommaker -i -k $K -p ${nam_bbt}tar -f 0.001 -t $j -n $n1 -s $rep_fil ${nam_bbt}tar.fa

biobloomcategorizer -i -t $j -d ${nam_bbt}tar -f ${nam_bbt}tar.bf -s $s $READS_B1_R1 $READS_B1_R2 >> ${nam_bbt}tar.fastq

mkdir -p $abyss_dir
cd $abyss_dir
abyss-pe -C $abyss_dir v=-v k=$k name=$nam_bbt np=$j lib='pet' pet=../${nam_bbt}tar.fastq long='longlib' longlib=$TARGETS B=18G H=4 kc=3

abyss_fa=$abyss_dir/$nam_bbt-10.fa
kollector-extract.sh $nam_bbt $abyss_fa $TARGETS $abyss_dir $j

cut -f1 -d " " hitlist.txt|sort|uniq > succeedtranscripts.txt
