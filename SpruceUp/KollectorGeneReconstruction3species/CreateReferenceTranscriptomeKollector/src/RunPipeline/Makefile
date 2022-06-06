###########################################################
#####################Kristina Gagalova##################### 
###########################################################
####################21.02.2017#############################
###########################################################

#Description: evaluate multeplicity from cdhit output
#------------------------------------------------------------------------------

#-----------------------------------------------------------
#Required sofware
#-----------------------------------------------------------
#Python 2.7
#awk
#
#ParseBLASTv2.1.py (dep: csv, getopt) 
#RunParseBlastv2.1.sh

BLAST := /home/shammond/src/ncbi-blast-2.2.31+/bin

BASE := input
IN := $(shell ls $(BASE)/*.fa | sed 's/.fa//' | rev | cut -f1 -d/ | rev)
OUT := $(addsuffix .blast.filt.out, $(IN))  

OUT_DIR := out

.SECONDARY: *.blast.filt.out

all: counts clean

clean:
	rm *.nin *.nhr *.nsq


#----------------------------------------
# in: .blast.filt.out
# out: counts of aligned seq per query
#----------------------------------------

counts: $(OUT)
	bash RunParseBlastv2.1.sh .blast.filt.out $(OUT_DIR)
	#mkdir -p $(OUT_DIR)	
	#cp count_* $(OUT_DIR)
	#rm $(OUT)
	
#----------------------------------------
# in: .blast.out
# out: filtered blast (remove self hits and reverse compl)
#----------------------------------------

%.blast.filt.out: %.blast.out
	awk '$$1 != $$2 && $$7 < $$8' $< > $@


#----------------------------------------
# in: .fa file from cdhit
# out: blasted seq
#----------------------------------------

%.blast.out: $(BASE)/%.fa
	echo $<
	$(BLAST)/makeblastdb -in $< -dbtype nucl -out $*.db > $*.logDB
	$(BLAST)/blastn -db $*.db -query $< -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' -out $@ -num_threads 4
