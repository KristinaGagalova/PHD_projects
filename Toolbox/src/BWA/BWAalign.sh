#!/bin/bash

#Philip Richmond - Genome analysis introduction workshop

## Define variables to be used in the commands below
WORKING_DIR='/path/to/mydir/GAGALOVA_KRISTINA'
GENOME='/path/to/genome/dir/my_genome.fa'
FASTQR1='/path/to/fastq/dir/myfile1.fastq'
FASTQR2='/path/to/fastq/dir/myfile2.fastq'
SORTED_BAM='my_bam_sorted.bam'
BAM='my_bam.bam'
SAM='my_sam.sam'

## Let's go to directory we want to work in
cd $WORKING_DIR

## Index fasta file
ln -s $GENOME ./
bwa index $GENOME

## Map with BWA
### bwa mem -t <numProcs> -T <minScore> -k <kmerSize> <genome_index> <R1.fastq> <R2.fastq>  >  <out.sam>
### it can also use fasta files as input (note that the won't be the alignment quality)
bwa mem -t 4  $GENOME \
$FASTQR1 $FASTQR2 \
> $WORKING_DIR/$SAM

#convert sam to bam
samtools view  -b $WORKING_DIR/$SAM  -o $WORKING_DIR/$BAM 

# Sort the bam file using samtools sort
samtools sort $WORKING_DIR/$BAM -o $WORKING_DIR/$SORTED_BAM

## Index the sorted bam 
samtools index $WORKING_DIR/$SORTED_BAM

#remove intermediate
rm $WORKING_DIR/$BAM $WORKING_DIR/$SAM
