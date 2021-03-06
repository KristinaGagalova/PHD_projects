Plot coverage introns-exons for different mrsFast parameters
================

Load data
---------

``` r
library(reshape)
library( dplyr )
library( ggplot2 )
library( tidyr )
library( data.table )
library(plyr)

dataPath="/projects/btl/kgagalova/PHD_projects2/GeneDuplicationTool/data/CheckBinningAlignment"
allFiles <- list.files( path = dataPath, pattern = ".bedgraph2", full.names = TRUE )

l <- lapply( allFiles, function( fn ){
  d <- read.table( fn, header = F );
  d$fileName <- fn;
  d
  } );

d <- bind_rows( l );
dim(d)
```

    ## [1] 114  11

``` r
d$class = sapply(strsplit(sapply(strsplit(d$fileName,"/"),tail,1),"\\."),"[[",1)
d$params = as.factor(ifelse(grepl("lowere", d$class),"lower","default"))
d$len = d$V5 - d$V4
d$bins = as.factor(c(rep(100,38),rep(50,38),rep(1,38)))
```

Including Plots
---------------

``` r
dvals = subset(d,d$V3 == "intron" | d$V3 == "exon")

ggplot(dvals, aes(x=V3, y=V10,color=V3)) + 
  geom_point(aes(size=len)) + 
  xlab("Type") +
  ylab("Mean coverage") +
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(face='bold',size=14,color='black')) + facet_grid( params ~ bins  )
```

![](images/exons_intronsDiffParams-1.png)
