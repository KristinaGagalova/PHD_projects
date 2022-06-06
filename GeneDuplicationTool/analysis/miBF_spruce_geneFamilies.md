miBF gene family
================

Load the matches
----------------

``` r
library(ggplot2)

dataHits <- read.delim("/projects/btl/kgagalova/PHD_projects2/GeneDuplicationTool/data/miBF/WS77111geneFamilies/all_summaryHits.tsv", header=FALSE)

rownames(dataHits) = dataHits[,1]
dataHits = dataHits[,2:ncol(dataHits)]

#dataHits<- sapply(sapply(dataHits, as.character),as.numeric)
#dataHits

dataHits$median = apply(dataHits,1,median)
dataHits$mad = apply(dataHits,1,mad)

summary(dataHits$mad)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##     0.00     8.90    16.31   182.74    37.06 75296.06

``` r
summary(dataHits$median)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0      44      91    1284     273  429232

``` r
ggplot(dataHits, aes(x=row.names(dataHits), y=median)) + geom_point() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
```

![](images/loadETL-1.png)

``` r
ggplot(dataHits, aes(x=median, y=mad)) + geom_point() #+ 
```

![](images/loadETL-2.png)

``` r
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
```

    ## List of 3
    ##  $ axis.title.x: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.text.x : list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  $ axis.ticks.x: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

``` r
######remove zero values
dataHitsNoZero = subset(dataHits, dataHits$median > 0 )
dim(dataHits)
```

    ## [1] 13035    19

``` r
dim(dataHitsNoZero)
```

    ## [1] 13030    19

``` r
summary(dataHits$median)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0      44      91    1284     273  429232

``` r
summary(dataHitsNoZero$median)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       6      44      91    1285     273  429232
