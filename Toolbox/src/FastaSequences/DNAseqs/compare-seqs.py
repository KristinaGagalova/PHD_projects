#!/usr/bin/env python

import sys

from fastatools import fasta_iter


oldsqn = sys.argv[1]
newsqn = sys.argv[2]

# read in old scaffs with sqn as key and id as value
# then read in new scaffs and compare sqn to keys and write out old and new ids upon match

seqdict = {}

collisions = open("collisions.txt","w")

with open(oldsqn,"r") as infile:
    for rec in fasta_iter(infile):
        seqid = rec[0]
        sqn = rec[1]
        if sqn in seqdict:
            print >> collisions, "Collision between " + seqid + " AND previously added " + seqdict[sqn]
#            sys.exit()
        else:
            seqdict[sqn] = seqid

print "OldID\tNewID"

with open(newsqn,"r") as infile:
    for rec in fasta_iter(infile):
        seqid = rec[0]
        sqn = rec[1]
        if sqn in seqdict:
            print seqdict[sqn] + "\t" + seqid

collisions.close()

### EOF ###
