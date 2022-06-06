################################################
##########Kristina Gagalova#####################
################################################
##############1 Mar 2017########################
################################################

#Description: Run kollector on sitka spruce reads

export SHELL=bash #-o pipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %

#--------------------------------------------------------------------------------------------
# Kollector options
j=12
r=0.9
s=0.7
k=32
K=25
n=50000000
a=30000M
o=trial_kollectorHalf_RepFil
p=/extscratch/btl/bullfrog/genome/tga/repeat_filters/spruce_250.bf

#---------------------------------------------------------------------------------------------
#input files
seed=/extscratch/btl/kgagalova/trial_kollector/seed_500/cdhit-output-4_1.fasta
pe1=/extscratch/btl/kgagalova/trial_kollectorHalf500/WS77111-reads-1pet.in
pe2=/extscratch/btl/kgagalova/trial_kollectorHalf500/WS77111-reads-2pet.in

list_r1=$(shell sed "s/ /?/g" <(cat $(pe1) | tr '\n' ' '))
list_r2=$(shell sed "s/ /?/g" <(cat $(pe2) | tr '\n' ' '))

# Report Run time and memory usage
#Path to Kollector executables

.DELETE_ON_ERROR:
.SECONDARY:

SOFTWARE_DIR = $PATH

all: assembledtargets.fa
#all: check_versions

#auxiliary function to pass values with space
#check "makefile, passing arguments with white spaces": http://stackoverflow.com/questions/9838384/can-gnu-make-handle-filenames-with-spaces
+s = $(subst ?,\ ,$1)

assembledtargets.fa: check_versions
	kollector.sh -j$j -r$r -s$s -k$k -K$K -n$n -o$o -a$a -p$p $(seed) $(call +s,$(list_r1)) $(call +s,$(list_r2)) 2>&1 | tee trial_kollectorHalfRepFil.log
	echo 'Kollector completed.'

check_versions:
	echo $(value SOFTWARE_DIR)
	biobloommaker --version
	abyss-fac --version
	gsnap --version
	samtools --version
