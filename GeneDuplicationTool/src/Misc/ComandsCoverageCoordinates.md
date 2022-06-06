## Use simple command line tools to calculate coverage in intervals


**pysamstats**

Create a bedgraph for binned coverage based on bam file     

Example: calculated binned for window size 100 and symmetric offset 50.     
You will need the sequence fast aand the corresponding bam file. with ```-d``` all the positions are included, also those with zero coverage.
```
pysamstats --type coverage_binned --window-size 100 --window-offset 50 -d -f ENSG00000005889.fa NA12878_S1_L001_R1_001.bbt175unmaskedSorted.bam | \ 
awk -v OFS='\t' '{print $1,$2,$2 + 100,$4}' > ENSG00000005889binned100.bedgraph
```

Example: calculated binned for window size 50 and symmetric offset 25.     
```
pysamstats --type coverage_binned --window-size 50 --window-offset 25 -d -f ENSG00000005889.fa NA12878_S1_L001_R1_001.bbt175unmaskedSorted.bam | \
awk -v OFS='\t' '{print $1,$2,$2 + 50,$4}' > ENSG00000005889binned50.bedgraph
```

**bedtools map**      
This function allows to calculate the average of a given set of intervals in a bed file.     
Input: bedfile (or gff file), begraph file from pysamstats    
Output: bedfile with average coverage in the last column

-c needs to be set to 4, the operation performed by bedtools id ```-o mean```
```
bedtools map -c 4 -a ../ENST00000379188.gff  -b ENSG00000005889binned50.bedgraph -o mean 
```
