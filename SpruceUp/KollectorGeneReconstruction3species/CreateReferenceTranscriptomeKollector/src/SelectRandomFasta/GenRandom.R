#!/gsc/software/linux-x86_64-centos6/R-3.3.2/bin/Rscript


####################################
#########Kristina Gagalova##########
####################################
############23/02/2017##############
####################################

#Description: Sample from fasta file randomly

library(Biostrings)

#options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE) 
# in this case input file, output file and number of seqs

#------------------------------------------------

fasta = readDNAStringSet(args[1])

#------------------------------------------------

new_fasta = sample(1:length(fasta),args[3],replace=F)

writeXStringSet(fasta[new_fasta], file=args[2])