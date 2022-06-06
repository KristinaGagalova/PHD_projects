GMAP alignments for 3 species - cdhit4
================

Upload the results from the alignment - coverage/identity and length of the targets
-----------------------------------------------------------------------------------

``` r
library( dplyr )
library( ggplot2 )
library( tidyr )
library( data.table )
library(plyr)
library(splitstackshape)

dataPath=c("/projects/btl/kgagalova/PHD_projects2/SpruceUp/KollectorGeneReconstruction3species/ComparativeAnlysis/data/GMAPcDNA_GCAT")

allFiles <- list.files( path = dataPath, pattern = "Coverage", full.names = TRUE )

l <- lapply( allFiles, function( fn ){
  d <- read.table( fn, header = F)
  d$fileName <- fn;
  d
  } );

allGMAP <- bind_rows( l );
dim(allGMAP)
```

    ## [1] 98772     4

``` r
colnames(allGMAP)[1:3] = c("query","identity","coverage")

#load the len of the sequences
lenSeqs <- read.delim("/projects/btl/kgagalova/PHD_projects2/SpruceUp/KollectorGeneReconstruction3species/ComparativeAnlysis/data/GMAPcDNA_GCAT/lenSeqs.txt", header=FALSE)
colnames(lenSeqs) = c("query","len")
allGMAPm = merge(allGMAP,lenSeqs,by="query",all.x=T)

allGMAPm$species = gsub("CoverageIdentity||targetsTargets.txt||targets.txt","",sapply(strsplit(allGMAPm$fileName,"/"),"[[",11))

countsPaths = as.data.frame(table(allGMAPm[,c("query","species")]))

allGMAPm = merge(allGMAPm,countsPaths,by.x=c("query","species"),by.y=c("query","species"))

allGMAPm$type = ifelse(grepl("augustus|maker|snap|genemark",allGMAPm$query), "Maker", "GCAT")#number of total alignments

#number of unique alignments - targets
detach("package:plyr", unload=TRUE) 
test = allGMAPm %>% group_by(query,species,type) %>% summarize(count=n())

table(subset(test,test$species=="PG29")$type)
```

    ## 
    ##  GCAT Maker 
    ##  7269 10785

``` r
table(subset(test,test$species=="Q903")$type)
```

    ## 
    ##  GCAT Maker 
    ##  7398 10760

``` r
table(subset(test,test$species=="WS77111")$type)
```

    ## 
    ##  GCAT Maker 
    ##  7273 10759

Plot the coverage and identity for the GMAP alignments
------------------------------------------------------

``` r
ggplot(allGMAPm, aes(x=coverage, y=identity)) + 
    geom_point(aes(shape=Freq > 1,colour = Freq > 1)) + facet_wrap( ~ species ) 
```

![](images/plots_cov_identity-1.png)

``` r
ggplot(allGMAPm, aes(x=coverage, y=identity)) + 
    geom_point(aes(shape=Freq > 1,colour = Freq > 1)) + facet_wrap( ~ species ) + 
  coord_cartesian(xlim = c(90,100),ylim= c(95,100))
```

![](images/plots_cov_identity-2.png)

``` r
#select only the ones that are >90% cov and 95% identity
allGMAPmSel = subset(allGMAPm,allGMAPm$coverage>=90 & allGMAPm$identity>=95)

tmp = allGMAPmSel %>% group_by(query,species,type) %>% summarize(count=n())
table(tmp[,c("species","type")])
```

    ##          type
    ## species    GCAT Maker
    ##   PG29     3768 10765
    ##   Q903     2809  7407
    ##   WS77111  3843  8589

``` r
ggplot(allGMAPmSel, aes(x=coverage, y=identity)) + 
    geom_point(aes(shape=Freq > 1,colour = Freq > 1),size=3) + facet_wrap( ~ species )
```

![](images/plots_cov_identity-3.png)

``` r
#Quantify exactly for the 3 species
table(subset(allGMAPmSel, allGMAPmSel$Freq == 1)$species)
```

    ## 
    ##    PG29    Q903 WS77111 
    ##   10389    7081    7527

``` r
#divide per type
ggplot(allGMAPm, aes(x=coverage, y=identity)) + 
    geom_point(aes(shape=type,colour = type),size=3) + facet_wrap( ~ species )
```

![](images/plots_cov_identity-4.png)

``` r
ggplot(allGMAPmSel, aes(x=coverage, y=identity)) + 
    geom_point(aes(shape=type,colour = Freq > 1),size=3) + facet_wrap( ~ species ) + 
    scale_shape_manual(values=c(3, 16)) + 
    scale_color_manual(values=c('#999999','#E69F00'))
```

![](images/plots_cov_identity-5.png)

``` r
#GCAT only
allGMAPmSelGCAT = subset(allGMAPmSel,allGMAPmSel$type == "GCAT")
allGMAPmSelMaker = subset(allGMAPmSel,allGMAPmSel$type == "Maker")

ggplot(allGMAPmSelGCAT, aes(x=coverage, y=identity)) + 
    geom_point(aes(shape=type,colour = Freq > 1),size=3) + facet_wrap( ~ species ) + 
    scale_shape_manual(values=c(3)) + 
    scale_color_manual(values=c('#999999','#E69F00'))
```

![](images/plots_cov_identity-6.png)

``` r
ggplot(allGMAPmSelMaker, aes(x=coverage, y=identity)) + 
    geom_point(aes(shape=type,colour = Freq > 1),size=3) + facet_wrap( ~ species ) + 
    scale_shape_manual(values=c(16)) + 
    scale_color_manual(values=c('#999999','#E69F00'))
```

![](images/plots_cov_identity-7.png)

``` r
allGMAPmSelGCATfreq = subset(melt(table(allGMAPmSelGCAT[,c("query","species")]),id.var=query),melt(table(allGMAPmSelGCAT[,c("query","species")]),id.var=query)$value>0)
allGMAPmSelMakerfreq = subset(melt(table(allGMAPmSelMaker[,c("query","species")]),id.var=query),melt(table(allGMAPmSelMaker[,c("query","species")]),id.var=query)$value>0)


allGMAPmSelGCATfreq[allGMAPmSelGCATfreq$value == 1,] %>% group_by(species) %>% summarize(count=n())
```

    ## # A tibble: 3 x 2
    ##   species count
    ##    <fctr> <int>
    ## 1    PG29  3315
    ## 2    Q903  2396
    ## 3 WS77111  2896

``` r
allGMAPmSelMakerfreq[allGMAPmSelMakerfreq$value == 1,] %>% group_by(species) %>% summarize(count=n())
```

    ## # A tibble: 3 x 2
    ##   species count
    ##    <fctr> <int>
    ## 1    PG29  8700
    ## 2    Q903  6278
    ## 3 WS77111  6217

``` r
allGMAPmSelGCATfreq[allGMAPmSelGCATfreq$value == 2,] %>% group_by(species) %>% summarize(count=n())
```

    ## # A tibble: 3 x 2
    ##   species count
    ##    <fctr> <int>
    ## 1    PG29   373
    ## 2    Q903   331
    ## 3 WS77111   705

``` r
allGMAPmSelMakerfreq[allGMAPmSelMakerfreq$value == 2,] %>% group_by(species) %>% summarize(count=n())
```

    ## # A tibble: 3 x 2
    ##   species count
    ##    <fctr> <int>
    ## 1    PG29  1671
    ## 2    Q903   778
    ## 3 WS77111  1584

``` r
allGMAPmSelGCATfreq[allGMAPmSelGCATfreq$value > 2,] %>% group_by(species) %>% summarize(count=n())
```

    ## # A tibble: 3 x 2
    ##   species count
    ##    <fctr> <int>
    ## 1    PG29    80
    ## 2    Q903    82
    ## 3 WS77111   242

``` r
allGMAPmSelMakerfreq[allGMAPmSelMakerfreq$value > 2,] %>% group_by(species) %>% summarize(count=n())
```

    ## # A tibble: 3 x 2
    ##   species count
    ##    <fctr> <int>
    ## 1    PG29   394
    ## 2    Q903   351
    ## 3 WS77111   788

``` r
#plot the graphs
ggplot(allGMAPmSelGCATfreq, aes(x=value, fill=species)) + geom_density(alpha=.3) +
  ggtitle("GCAT targets") +
    xlab("Number of complete aligmnetsper target")
```

![](images/plots_multiple_alignments-1.png)

``` r
#ggplot(df, aes(x=dose, y=, fill=species)) 
ggplot(allGMAPmSelGCATfreq, aes(x=value, fill=species)) + geom_histogram(alpha=.3) +
    ggtitle("GCAT targets") +
      xlab("Number of complete aligmnets per target") + 
          facet_wrap( ~ species )
```

![](images/plots_multiple_alignments-2.png)

``` r
ggplot(allGMAPmSelMakerfreq, aes(x=value, fill=species)) + geom_density(alpha=.3) +
  ggtitle("Maker targets") +
    xlab("Number of complete aligmnets per target")
```

![](images/plots_multiple_alignments-3.png)

``` r
ggplot(allGMAPmSelMakerfreq, aes(x=value, fill=species)) + geom_histogram(alpha=.3) +
    ggtitle("Maker targets") +
      xlab("Number of complete aligmnets per target") + 
                facet_wrap( ~ species )
```

![](images/plots_multiple_alignments-4.png)

Intersect with the common Kollector
-----------------------------------

``` r
Common <- read.table("/projects/btl/kgagalova/PHD_projects2/SpruceUp/KollectorGeneReconstruction3species/ComparativeAnlysis/data/GMAPcDNA_GCAT/common_intersection3.txt", quote="\"", comment.char="")
table(as.character(Common$V1) %in% unique(as.character(allGMAPmSelMakerfreq$query))) 
```

    ## 
    ## FALSE  TRUE 
    ##  2288  4794

``` r
namsCorr = gsub("_status","status",gsub("_cluster","cluster", gsub("_clone","clone", unique(as.character(Common$V1)))))

#check how many are reconstructed as complete
allGMAPmSelMakerfreq[as.character(allGMAPmSelMakerfreq$query) %in% as.character(Common$V1),] %>% group_by(species) %>% count(value)
```

    ## # A tibble: 15 x 3
    ## # Groups:   species [3]
    ##    species value     n
    ##     <fctr> <int> <int>
    ##  1    PG29     1  3936
    ##  2    PG29     2   684
    ##  3    PG29     3   135
    ##  4    PG29     4    30
    ##  5    PG29     5     8
    ##  6    Q903     1  3572
    ##  7    Q903     2   324
    ##  8    Q903     3    88
    ##  9    Q903     4    23
    ## 10    Q903     5    17
    ## 11 WS77111     1  3265
    ## 12 WS77111     2   767
    ## 13 WS77111     3   309
    ## 14 WS77111     4    71
    ## 15 WS77111     5    11

``` r
allGMAPmSelGCATfreq[as.character(allGMAPmSelGCATfreq$query) %in% as.character(namsCorr),] %>% group_by(species) %>% count(value)
```

    ## # A tibble: 12 x 3
    ## # Groups:   species [3]
    ##    species value     n
    ##     <fctr> <int> <int>
    ##  1    PG29     1   199
    ##  2    PG29     2    28
    ##  3    PG29     3     5
    ##  4    PG29     4     4
    ##  5    Q903     1   193
    ##  6    Q903     2    22
    ##  7    Q903     3    14
    ##  8    Q903     4     2
    ##  9 WS77111     1   151
    ## 10 WS77111     2    62
    ## 11 WS77111     3    22
    ## 12 WS77111     4     6

``` r
ggplot(allGMAPmSelMakerfreq[as.character(allGMAPmSelMakerfreq$query) %in% as.character(Common$V1),] %>% group_by(species) %>% count(value),aes(x=value,y=n,colour=species)) + geom_bar(stat="identity") +  facet_wrap( ~ species ) + 
ggtitle("Maker targets") +
xlab("Number of complete aligmnets per target") + ylab("Count")
```

![](images/intersect_KollectorCommon-1.png)

``` r
ggplot(allGMAPmSelGCATfreq[as.character(allGMAPmSelGCATfreq$query) %in% as.character(namsCorr),] %>% group_by(species) %>% count(value),aes(x=value,y=n,colour=species)) + geom_bar(stat="identity") +  facet_wrap( ~ species ) + 
ggtitle("GCAT targets") +
xlab("Number of complete aligmnets per target") + ylab("Count")
```

![](images/intersect_KollectorCommon-2.png)
