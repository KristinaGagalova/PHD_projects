Introns size from WS77111 and PG29
================

Plot the intron length
----------------------

``` r
library(ggplot2)
library(cowplot)
library(reshape)

WS_int <- read.table("/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/PreliminaryResults12May/data/NumLengthIntrons/WS77111mutipleIntrons.txt", header=TRUE, quote="\"")
PG_int <- read.table("/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/PreliminaryResults12May/data/NumLengthIntrons/PG29mutipleIntrons.txt", header=TRUE, quote="\"")

size_WS = as.numeric(unlist(strsplit(as.character(WS_int$introns_len), ",")))
size_PG = as.numeric(unlist(strsplit(as.character(PG_int$introns_len), ",")))

#size_all = as.numeric(unlist(strsplit(as.character(PG_WS$introns_len), ",")))

##plot the 2 distributions
## see if there is a distribution of 3, alignment artifact
ggplot(as.data.frame(size_WS), aes(x=size_WS)) +
    geom_histogram(binwidth=.5, position="dodge") + 
    coord_cartesian(ylim = c(0, 60),xlim= c(100,120))
```

![](images/non_overlapping_transcripts-1.png)

``` r
#plot together
WS_plot = ggplot(as.data.frame(size_WS), aes(x=size_WS)) +
    geom_histogram(binwidth=.5, position="dodge",colour="red") + 
    coord_cartesian(ylim = c(0, 60),xlim= c(0,200)) + 
    labs(x = "intron length")
PG_plot = ggplot(as.data.frame(size_PG), aes(x=size_PG)) +
    geom_histogram(binwidth=.5, position="dodge",colour="darkblue") + 
    coord_cartesian(ylim = c(0, 60),xlim= c(0,200)) + 
    labs(x = "intron length")
plot_grid(WS_plot, PG_plot, labels=c("WS77111", "PG29"), ncol = 1, nrow = 2)
```

![](images/non_overlapping_transcripts-2.png)

``` r
#plot the full length
WS_plotFull = ggplot(as.data.frame(size_WS), aes(x=size_WS)) +
    geom_histogram(binwidth=.5, position="dodge",colour="red") +
    labs(x = "intron length")  +
    geom_vline(xintercept = mean(size_WS),colour="black",linetype="dashed") + 
    coord_cartesian(xlim= c(0,max(size_WS)))
PG_plotFull = ggplot(as.data.frame(size_PG), aes(x=size_PG)) +
    geom_histogram(binwidth=.5, position="dodge",colour="darkblue") +
    labs(x = "intron length") + 
    geom_vline(xintercept = mean(size_PG),colour="black",linetype="dashed") + 
    coord_cartesian(xlim= c(0,max(size_WS)))
plot_grid(WS_plotFull, PG_plotFull, labels=c("WS77111", "PG29"), ncol = 1, nrow = 2)
```

![](images/non_overlapping_transcripts-3.png)

Only for overlapping introns (where the comparison - difference was made)
-------------------------------------------------------------------------

``` r
PG_WS <- read.table("/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/PreliminaryResults12May/data/NumLengthIntrons/OverlappingMultipleIntrons.txt", header=TRUE, quote="\"")

size_WScomm = as.numeric(unlist(strsplit(as.character(PG_WS$introns_len.WS), ",")))
size_PGcomm = as.numeric(unlist(strsplit(as.character(PG_WS$introns_len.PG29), ",")))


WS_plot1 = ggplot(as.data.frame(size_WScomm), aes(x=size_WScomm)) +
    geom_histogram(binwidth=.5, position="dodge",colour="red") + 
    coord_cartesian(ylim = c(0, 60),xlim= c(0,200)) + 
    labs(x = "intron length")
PG_plot1 = ggplot(as.data.frame(size_PGcomm), aes(x=size_PGcomm)) +
    geom_histogram(binwidth=.5, position="dodge",colour="darkblue") + 
    coord_cartesian(ylim = c(0, 60),xlim= c(0,200)) + 
    labs(x = "intron length")
plot_grid(WS_plot1, PG_plot1, labels=c("WS77111", "PG29"), ncol = 1, nrow = 2)
```

![](images/overlapping_transcripts-1.png)

``` r
WS_plotComm = ggplot(as.data.frame(size_WScomm), aes(x=size_WScomm)) +
    geom_histogram(binwidth=.5, position="dodge",colour="red") +
    labs(x = "intron length")  +
    geom_vline(xintercept = mean(size_WScomm),colour="black",linetype="dashed") + 
    coord_cartesian(xlim= c(0,max(size_PGcomm)))
PG_plotComm = ggplot(as.data.frame(size_PGcomm), aes(x=size_PGcomm)) +
    geom_histogram(binwidth=.5, position="dodge",colour="darkblue") +
    labs(x = "intron length") + 
    geom_vline(xintercept = mean(size_PGcomm),colour="black",linetype="dashed") + 
    coord_cartesian(xlim= c(0,max(size_PGcomm)))

plot_grid(WS_plotComm, PG_plotComm, labels=c("WS77111", "PG29"), ncol = 1, nrow = 2)
```

![](images/overlapping_transcripts-2.png)

Summary
-------

-   The histogram plot does not seem very different for the 2 species
-   When considering all the intron lengths, WS has a higher mean value. WHne considering only the overlapping transcripts, PG29 has longer intron size mean
-   No irregular mapping detected from the plots
