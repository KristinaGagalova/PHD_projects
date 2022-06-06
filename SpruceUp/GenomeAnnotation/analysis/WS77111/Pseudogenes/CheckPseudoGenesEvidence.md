Analysis of pseudogenes WS77111 - Maker-P - in Evidence based annotation
================

Load and process data
---------------------

``` r
library(ggplot2)

d = read.table(gzfile("/projects/btl/kgagalova/PHD_projects2/SpruceUp/GenomeAnnotation/WS77111/data/Pseudogenes/OutForR_EvidenceOnly.txt.gz"),header=F)
nams=c("trans","contig","coord","eval","stop_high","stop_low","frame_high","farme_low")
colnames(d) = nams
#convert to char
d$trans = as.character(d$trans)
d$contig = as.character(d$contig)
d$coord = as.character(d$coord)
d$start_match_prot= sapply(strsplit(sapply(strsplit(d$coord,":"),"[[",1),"-"),"[[",1)
d$stop_match_prot= sapply(strsplit(sapply(strsplit(d$coord,":"),"[[",1),"-"),"[[",2)
d$start_match_prot = as.numeric(d$start_match_prot)
d$stop_match_prot = as.numeric(d$stop_match_prot)
d$len_align = d$stop_match_prot - d$start_match_prot
d$prot_nam = sapply(strsplit(d$trans,"_protein"),"[[",1)

#load prot lens
lens = read.table("/projects/btl/kgagalova/PHD_projects2/SpruceUp/GenomeAnnotation/WS77111/data/Pseudogenes/proteinsPolishedEvidence.len")
colnames(lens) = c("prot_nam","len")
lens$prot_nam = as.character(lens$prot_nam)
d1 = merge(d,lens,by="prot_nam")
d1$perc_align = d1$len_align / d1$len
d1$contig1 = sapply(strsplit(d1$contig,":"),"[[",1) 
d1$start_coord = sapply(strsplit(sapply(strsplit(d1$contig,"\\|"),"[[",2),"-"),"[[",1)
d1$stop_coord = sapply(strsplit(sapply(strsplit(d1$contig,"\\|"),"[[",2),"-"),"[[",2)
d1$contig2 = paste(d1$contig1,d1$start_coord,d1$stop_coord,sep="_")

stop_c = ifelse(d1$stop_high > 0, "stop_codon", "none")
framesh = ifelse(d1$frame_high > 0, "frameshift", "none")

chech_mut <-function(stopc,fr){
  type_mut = c()
  for(i in 1:length(stopc)){
    si = stopc[i]
    fi = fr[i]
    if (si == "stop_codon" && fi == "none") {
      type_mut = c(type_mut,"stop_codon")
    } else if(si == "none" && fi == "frameshift") {
      type_mut = c(type_mut,"frameshift")
    } else if(si == "stop_codon" && fi == "frameshift") {
      type_mut = c(type_mut,"stop_codon/frameshift")
    } else 
      type_mut = c(type_mut,"none")
  }
  return(type_mut)
  }

type_mut = chech_mut(stop_c, framesh)
d1$type_mut = factor(type_mut,levels = c("none", "stop_codon","frameshift","stop_codon/frameshift"))
```

``` r
ggplot(lens, aes(x=len)) + geom_histogram(color="darkgreen", fill="white", binwidth=2) + 
        theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
          axis.title.y = element_text(face='bold',size=16,vjust=1),
          axis.text.x = element_text(face='bold',size=14,color='black'),
          axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/genes_len-1.png)

``` r
summary(lens$len)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    32.0   216.0   338.0   388.6   488.0  3677.0

Plot distributions
------------------

``` r
d1_counts = as.data.frame(table(d1$prot_nam))
d2_counts = merge(d1_counts,lens, by.x="Var1",by.y="prot_nam")
d2_counts = d2_counts[with(d2_counts, order(len)),]

ggplot(d2_counts, aes(len,Freq)) + geom_point(colour="darkgreen",alpha=0.2) + ylab("# of pseudogenes") + xlab("Peptide length (aa)") +
    geom_density2d(colour="black") +   
    theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
          axis.title.y = element_text(face='bold',size=16,vjust=1),
          axis.text.x = element_text(face='bold',size=14,color='black'),
          axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/psedogenes_len-1.png)

``` r
summary(d2_counts$Freq)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     1.0     4.0     9.0    17.8    20.0   400.0

``` r
ggplot(d1, aes(len,perc_align,colour=type_mut)) + geom_point(size = 3) + ylab("Pseudogene len - coverage orig pept") + xlab("Peptide length (aa)") + 
      theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
          axis.title.y = element_text(face='bold',size=16,vjust=1),
          axis.text.x = element_text(face='bold',size=14,color='black'),
          axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/plot_cov_individual_pseudo-1.png)

``` r
ggplot(d1, aes(len,perc_align)) + geom_point(colour="darkgreen", alpha=0.2) + ylab("Pseudogene len - coverage orig pept") + xlab("Peptide length (aa)") + 
     geom_density2d(colour="black") + 
      theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
          axis.title.y = element_text(face='bold',size=16,vjust=1),
          axis.text.x = element_text(face='bold',size=14,color='black'),
          axis.text.y = element_text(face='bold',size=14,color='black'),
          strip.text = element_text(size=15) ) + 
      facet_wrap( ~ type_mut )
```

![](images/plot_cov_individual_pseudo-2.png)

``` r
table(d1$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                181092                128939                  9258 
    ## stop_codon/frameshift 
    ##                 23773

Quantification
--------------

``` r
#matching bp with proteins and pseudogenes
sum(abs(as.numeric(sapply(strsplit(sapply(strsplit(d1$coord, ":"),"[[",2),"-"),"[[",1)))- abs(as.numeric(sapply(strsplit(sapply(strsplit(d1$coord, ":"),"[[",2),"\\-"),"[[",2))))
```

    ## [1] 270557

``` r
#genomic region covered by pseudogenes
sum(as.numeric(d1$stop_coord) - as.numeric(d1$start_coord))
```

    ## [1] 1710639349

``` r
no_mutd1noreps = subset(d1,d1$type_mut == "none")
sum(as.numeric(no_mutd1noreps$stop_coord) - as.numeric(no_mutd1noreps$start_coord))
```

    ## [1] 443399735

``` r
stopd1noreps = subset(d1,d1$type_mut == "stop_codon")
framed1noreps = subset(d1,d1$type_mut == "frameshift")
stopframed1noreps = subset(d1,d1$type_mut == "stop_codon/frameshift")

sum(as.numeric(stopd1noreps$stop_coord) - as.numeric(stopd1noreps$start_coord))
```

    ## [1] 754768796

``` r
sum(as.numeric(framed1noreps$stop_coord) - as.numeric(framed1noreps$start_coord))
```

    ## [1] 84645435

``` r
sum(as.numeric(stopframed1noreps$stop_coord) - as.numeric(stopframed1noreps$start_coord))
```

    ## [1] 427825383

``` r
#filter by the % aligned
table(subset(d1,d1$perc_align > 0.5)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                 21101                 25479                  2124 
    ## stop_codon/frameshift 
    ##                  8259

``` r
table(subset(d1,d1$perc_align > 0.7)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                  8870                 10744                  1125 
    ## stop_codon/frameshift 
    ##                  4353

``` r
#filter by e-value
table(subset(d1,d1$eval < 1e-50)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                 19650                 28190                  2121 
    ## stop_codon/frameshift 
    ##                  9037

``` r
table(subset(d1,d1$eval < 1e-100)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                  3558                  6375                   524 
    ## stop_codon/frameshift 
    ##                  2731
