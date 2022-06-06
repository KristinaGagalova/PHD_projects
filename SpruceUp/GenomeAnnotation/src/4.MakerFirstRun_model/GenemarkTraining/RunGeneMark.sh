#!/bin/bash

##########LIBS
export PATH=$PATH:/gsc/btl/linuxbrew/bin
export PERL5LIB=$PERL5LIB:/home/kgagalova/localperl/lib/site_perl/5.20.1

genome=/path/to/the/genome/assembly/genome.fasta
genemark=/home/kgagalova/src/gm_et_linux_64/gmes_petap/gmes_petap.pl

/gsc/software/linux-x86_64-centos7/perl-5-5.20.3/bin/perl $genemark --cores 32 --ES --sequence $genome

