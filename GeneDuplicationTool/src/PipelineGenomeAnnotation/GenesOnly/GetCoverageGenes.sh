#!/bin/bash

#SBATCH --job-name=mrfast_human
#SBATCH --partition=all
#SBATCH --ntasks=32
#SBATCH --mem=200G

#https://github.com/sfu-compbio/mrsfast

#################Software paths - to substitute with source
export PATH=$PATH:$HOME/arch/filer/bin/biobloom_ntHashBF3/bin:/gsc/btl/linuxbrew/bin
mrsfast=/home/kgagalova/src/mrsfast/mrsfast
##################

##################input files
list_pe=/projects/btl_scratch/kgagalova/PipelineCNV/ProofOfConcept/GenomeAnnotationsGenesOnly/H.sapiens/mrsFast/list_pe.in
fasta_file=/projects/btl_scratch/kgagalova/PipelineCNV/ProofOfConcept/GenomeAnnotationsGenesOnly/H.sapiens/mrsFast/Index/mRNAHomo_sapiens.GRCh38.90flankingMerged.fa
annotation=/projects/btl/datasets/hsapiens/hg38.p10/annotation/Homo_sapiens.GRCh38.90.gff3
genome=/projects/btl/datasets/hsapiens/hg38.p10/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
###################

################Filter reads
#get the reads
#nam=$(echo ${nam}flankingMerged.fa | sed 's/.fa//')
biobloommaker -p $nam -f 0.0001 ${nam_ann}flankingMerged.fa
biobloomcategorizer -t 24 -d ${nam_ann} -f ${nam_ann}.bf -s 175 -e -l $list_pe | gzip > ${nam_ann}_175.fastq

################Align reads
$mrsfast --search $fasta_file --pe --seq ${nam_ann}_175.fastq.gz --seqcomp \
        --min 300 --max 600 --threads 12 -o {nam_ann}.sam
samtools view -bS ${nam_ann}.sam > ${nam_ann}.bam
samtools sort ${nam_ann}.bam -o ${nam_ann}Sorted.bam
samtools index ${nam_ann}Sorted.bam


