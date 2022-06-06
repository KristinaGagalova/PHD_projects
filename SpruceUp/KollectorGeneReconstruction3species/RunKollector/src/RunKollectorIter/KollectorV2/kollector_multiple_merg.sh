#!/bin/bash

#------------------------------------------------------------
# Usage
#------------------------------------------------------------

PROGRAM=$(basename $0)
read -r -d '' USAGE <<HEREDOC
Usage: $PROGRAM [options] <seed> <merg> <pe> 

Description:

Do a targeted assembly  using ABySS. The input files are  
PET sequencing reads which must be a FASTA/FASTQ pair and a 
seed sequence FASTA file to recruit reads. The input files may be gzipped.


AbySS(1.5+),BioBloom Tools and GMAP should be in your path.

Options:

    
    -h        show this help message
    -j N      threads [1]
    -r N      min match length for tagging  reads [0.7]
    -s N      min match length for recruiting reads [0.50]
    -k N      k-mer size for ABySS contig assembly [32]
    -K N      k-mer size for read overlap detection [25]
    -n N      max k-mers to recruit in total [10000]
    -o FILE   output file prefix ['kollector']
    -p FILE   Bloom filter containing repeat k-mers for
              exclusion from scoring calculations; must match
              k-mer size selected with -K opt [disabled]
    -max_iterations N  number of iterations to be performed [5]
    -decrement N       decrement of the r parameter in each iteration [0.1]


HEREDOC

set -eu -o pipefail

#------------------------------------------------------------
# Parse command line opts
#------------------------------------------------------------

# default values for options
a=0
align=0
evaluate=0
j=1
k=32
K=25
r=0.7
s=0.5
o=kollector
n=10000
help=0
d=0.10
max_iterations=5
e=1

#parse command line options
while getopts a:A:d:e:g:hH:Cj:k:K:r:s:m:n:o:p: opt; do
	case $opt in
		a) a=$OPTARG;;
		A) abyss_opt="$OPTARG";;
		C) clean=0;;
		d) decrement=$OPTARG;;
		e) e=$OPTARG;;
		g) ref=$OPTARG;;
		h) help=1;;
		H) num_hash=$OPTARG;;
		j) j=$OPTARG;;
		k) k=$OPTARG;;
		K) K=$OPTARG;;
		r) r=$OPTARG;;
		s) s=$OPTARG;;
		m) max_iterations=$OPTARG;;
		n) n=$OPTARG;;
		o) o=$OPTARG;;
		p) p=$OPTARG;;

		\?) echo "$PROGRAM: invalid option: $OPTARG"; exit 1;;
	esac
done
shift $((OPTIND-1))

# -h for help message
if [ $help -ne 0 ]; then
	echo "$USAGE"
	exit 0;
fi

# we expect 3 file arguments
if [ $# -lt 3  ]; then
    echo "Error: number of file args must be  3" >&2
	echo "$USAGE" >&2
	exit 1
fi

seed=$1; shift;
merg=$1; shift;
pe=$1; shift;

for i in $(seq 2 $e)
do
if [ "$i" = 2 ]
then
mkdir -p iterationMerg.$i
cd iterationMerg.$i
shuffled=$(echo $pe | rev | cut -d/ -f1 | rev | cut -f1 -d. )
shuf $pe --output=./${shuffled}.$i
kollector.sh -j$j -d$d -k$k -K$K -r$r -s$s -p$p -n$n -o$o -a$a -e$i $seed $merg ./${shuffled}.$i
cut -f1 -d " " hitlist.txt|sort|uniq > succeedtranscripts.txt
grep ">" $seed | sed 's/>//g' >alltranscripts.txt
awk '($1<0.90) && ($1>=0.80)  { print $2,$4;}' allcoveredmapped.sam > hitlist_partial.txt
cut -f2 -d " " hitlist_partial.txt|sort|uniq|xargs samtools faidx $o.abyss/$o-10.fa > assembledtargets_partial.fa
grep -w -v -f succeedtranscripts.txt alltranscripts.txt|xargs samtools faidx $seed > failedtranscripts.fa
cat assembledtargets_partial.fa failedtranscripts.fa > failedtranscripts_all.fa
cd ..
else
preve=$(($i-1))
mkdir -p iterationMerg.$i
cd iterationMerg.$i
#check if there is an existing file (to avoid makefile re-run targets when error)
if [ ! -f prevfailed.fa ]
then
cp ../iterationMerg.$preve/failedtranscripts_all.fa prevfailed.fa
fi
if [ ! -f prevfailed.fa.fai ]
then
samtools faidx prevfailed.fa
fi
shuffled=$(echo $pe | rev | cut -d/ -f1 | rev | cut -f1 -d. )
shuf $pe --output=./${shuffled}.$i
seed_new="$(pwd)"/prevfailed.fa
kollector.sh -j$j -d$d -k$k -K$K -r$r -s$s -p$p -n$n -o$o -a$a -e$i $seed_new $merg ./${shuffled}.$i
cut -f1 -d " " hitlist.txt|sort|uniq > succeedtranscripts.txt
grep ">" prevfailed.fa | sed 's/>//g' >alltranscripts.txt
awk '($1<0.90) && ($1>=0.80)  { print $2,$4;}' allcoveredmapped.sam > hitlist_partial.txt
cut -f2 -d " " hitlist_partial.txt|sort|uniq|xargs samtools faidx $o.abyss/$o-10.fa > assembledtargets_partial.fa
#re-run same targets in next iteration if no assembled targets
if [ ! -s succeedtranscripts.txt ]
then
cp prevfailed.fa failedtranscripts.fa
else
grep -w -v -f succeedtranscripts.txt alltranscripts.txt|xargs samtools faidx prevfailed.fa > failedtranscripts.fa
fi
cat assembledtargets_partial.fa failedtranscripts.fa > failedtranscripts_all.fa
cd ..
fi
done
