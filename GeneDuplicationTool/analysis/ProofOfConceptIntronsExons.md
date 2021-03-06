Plot coverage exons and introns
================

Load data
---------

``` r
library( dplyr )
library( ggplot2 )
library( tidyr )
library( data.table )

nobins <- read.delim("/projects/btl/scratch/kgagalova/PipelineCNV/ProofOfConcept/CheckCopies/HumanReference/CheckDifferentCopies/mrFastAligner/MultipleCopies/CalculateAverage/BinnedIntervals/NoBinnedENSG00000005889.masked.bedgraph1", header=FALSE)
nobins$len = nobins$V5 - nobins$V4
nobins2 = nobins[,c("V1","V3","V10","len")]
nobins2 = nobins2[3:nrow(nobins2),]
nobins2$V3 = as.factor(nobins2$V3)

binned50 <- read.delim("/projects/btl/scratch/kgagalova/PipelineCNV/ProofOfConcept/CheckCopies/HumanReference/CheckDifferentCopies/mrFastAligner/MultipleCopies/CalculateAverage/BinnedIntervals/BinnedENSG00000005889_50.masked.bedgraph1", header=FALSE)
binned50$len = binned50$V5 - binned50$V4
binned502 = binned50[,c("V1","V3","V10","len")]
binned502 = binned502[3:nrow(binned502),]
binned502$V3 = as.factor(binned502$V3)

binned100 <- read.delim("/projects/btl/scratch/kgagalova/PipelineCNV/ProofOfConcept/CheckCopies/HumanReference/CheckDifferentCopies/mrFastAligner/MultipleCopies/CalculateAverage/BinnedIntervals/BinnedENSG00000005889_100.masked.bedgraph1", header=FALSE)
binned100$len = binned100$V5 - binned100$V4
binned1002 = binned100[,c("V1","V3","V10","len")]
binned1002 = binned1002[3:nrow(binned1002),]
binned1002$V3 = as.factor(binned1002$V3)

testType <- read.delim("/projects/btl/scratch/kgagalova/PipelineCNV/ProofOfConcept/CheckCopies/HumanReference/CheckDifferentCopies/mrFastAligner/MultipleCopies/CalculateAverage/Intervals.all.bedgraph", header=FALSE)
testType$len = testType$V3 - testType$V2
testType$V3 = as.factor(testType$V5)
```

Plot coverage
-------------

``` r
ggplot(nobins2, aes(x=V3, y=V10,color=V3)) + 
  geom_point(aes(size=len)) + 
  xlab("Type") +
  ylab("Mean coverage") +
  ggtitle("No binning") + 
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/exons_intronsCov-1.png)

``` r
ggplot(binned502, aes(x=V3, y=V10,color=V3)) + 
  geom_point(aes(size=len)) + 
  xlab("Type") +
  ylab("Mean coverage") +
  ggtitle("Binned 50") + 
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/exons_intronsCov-2.png)

``` r
ggplot(binned1002, aes(x=V3, y=V10,color=V3)) + 
  geom_point(aes(size=len)) + 
  xlab("Type") +
  ylab("Mean coverage") +
  ggtitle("Binned 100") + 
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/exons_intronsCov-3.png)

``` r
ggplot(testType, aes(x=V5, y=V4,color=V5)) + 
  geom_point(aes(size=len)) + 
  xlab("Type") +
  ylab("Mean averge") +
  ggtitle("Binned type") + 
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/exons_intronsCov-4.png)
