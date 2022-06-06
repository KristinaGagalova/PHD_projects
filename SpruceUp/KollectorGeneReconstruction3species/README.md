# Create reference transcriptome for Kollector runs (CreateReferenceTranscriptomeKollector)

[Kollector](https://github.com/bcgsc/kollector) - targeted *de novo* genome assembly.    

Given several resources (cDNA from EST, high confidence gene data set and genomic map) from 2 spruce strains, define the most representative data set to input in Kollector.     

In order to create the optimal transcriptome data set it is necessary to reduce the redundancy and to increase the comprehensiveness. A suitable threshold needs to be set.  

The transcriptome data set will be used to assemble the gene regions for comparative genomics.      


## Useful links

- [SpruceUp genome assembly](https://www.bcgsc.ca/wiki/pages/viewpage.action?pageId=22643711)

- [Info about PG mapping data](https://www.bcgsc.ca/jira/browse/BTL-236)

- [Exploratory analysis of Genetic Map and Sequence capture mapping](https://www.bcgsc.ca/jira/browse/BTL-253)


## Previous work

S. Austin Hammond: run Cdhit over several data sets. Please find the ticket [here](https://www.bcgsc.ca/jira/browse/BTL-820)


# Optimize Kollector for spruce (RunKollector)

Test Kollector on Genesis cluster and optimize parameters for targeted gene assembly.

# PreliminaryResults12May

Preliminary comparison between assembled WS77111 genes and previously assembled PG29 genes. Analysis is performed for Spruce meeting on the 12th of May.   

# PreliminaryResults1June
Preliminary comparison between assembled WS77111, Q903 genes and previously assembled PG29 genes. Analysis is performed for Spruce meeting on the 1st of June.  


## References

Limin Fu, Beifang Niu, Zhengwei Zhuy, Sitao Wu and Weizhong Li (2012) CD-HIT: accelerated for clustering the next-generation sequencing data
