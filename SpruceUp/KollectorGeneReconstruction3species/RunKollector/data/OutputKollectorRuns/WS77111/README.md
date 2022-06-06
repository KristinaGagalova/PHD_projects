## Assembled contigs and coverages; 5 iterations Kollector run + 2 extra iterations

allcoveredmapped_WS77111: coverege percentage after alignment contig-target      

allassembledtargetsWS77111_succeedOneline: contig sequences, one per line fasta format       

**NOTE THAT SOME CONTIGS ALIGN TO MULTIPLE TARGETS, usually those genes are close into space**    

```
[kgagalova@kgagalova01]$ wc -l allcoveredmapped_WS77111.out
12377 allcoveredmapped_WS77111.out
[kgagalova@kgagalova01]$ grep ">" allassembledtargetsWS77111_succeedOneline.fa | wc -l
12276
```

blacklistWS77111 - targets assembled by multiple contigs      

The files are compressed in one unique tar.gz
