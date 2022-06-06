#!/usr/bin/env python

import sys
#geneCoord_file = "/projects/btl_scratch/kgagalova/PipelineCNV/ProofOfConcept/GenomeAnnotationsMasked/A.thaliana/BWA/Genes_rpkm/Test.gff"

geneCoord_file = sys.argv[1]
nameGene_coord = geneCoord_file.strip(".gff")

segm_range = 1000
#left_over_max = 800

f= open(nameGene_coord + "Segm" + ".gff","w+")
with open(geneCoord_file) as gff3_file:
    for line in gff3_file:
        if not line.startswith("#"):
            chrom,source,type,start,stop,score, strand, frame, attribute = line.rstrip(" ").rstrip('\n').split('\t')
            start = int(start)
            stop = int(stop)
            if int(stop) - int(start) >= segm_range:
                ranges_elem = list(range(start,stop+segm_range,segm_range))
                for elem in range(0,len(ranges_elem)-1):
                    f.write(str(chrom) + '\t' + source + '\t' + type + '\t' + str(ranges_elem[elem]) + '\t' + str(ranges_elem[elem+1]) + '\t' + strand + '\t' + frame + '\t' + attribute + '\n')
            else:
                print("The gene range is shorter than the segmentation interval: gene on " + str(chrom) + '_____' + str(start) + '_____' + str(stop) + '_____' + "won't be output")
f.close()