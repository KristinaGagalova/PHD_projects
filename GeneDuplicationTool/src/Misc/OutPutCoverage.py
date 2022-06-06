#!/usr/bin/env python

#Description: https://github.com/alimanfoo/pysamstats
#Bined options: https://github.com/alimanfoo/pysamstats/blob/master/pysamstats/binned.py
#Python 2.7

import sys
import pysam
import pysamstats
import numpy
import matplotlib.pyplot as plt


START=1
BIN_SIZE=1000

fafile = "/media/kgagalova/SD2562/Documents/Spruce/GeneDuplicationPipeline/1.Regions/TestDataSets/ENSG00000005889.fa.masked"
mybam = pysam.AlignmentFile('/media/kgagalova/SD2562/Documents/Spruce/GeneDuplicationPipeline/1.Regions/TestDataSets/NA12878_S1_L001_R1_001.bbt175Sorted.bam')

#numpy version
a = pysamstats.load_coverage(pysamstats.stat_coverage_binned(mybam, fafile= fafile,start=100,end=10000, window_size=100, chrom='ENSG00000005889'))
print a
#plt.plot(a.pos, a.reads_all)
#plt.show()

#looping version
#for rec in pysamstats.stat_coverage_binned(mybam, fafile= fafile,start=100,end=10000, window_size=100, chrom='ENSG00000005889'):
	#print(rec)    
	#print rec['chrom'], rec['pos'], rec['reads_all'], rec['reads_pp']
    
#a = pysamstats.load_coverage(mybam, chrom='ENSG00000005889')
#print(a.pos)
#plt.plot(a.pos, a.reads_all)
#plt.show()
