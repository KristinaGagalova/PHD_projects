Analysis of pseudogenes WS77111 - Maker-P - Evidence and prediction based
================

Load and process data
---------------------

``` r
library(ggplot2)

d = read.table(gzfile("/projects/btl/kgagalova/PHD_projects2/SpruceUp/GenomeAnnotation/WS77111/data/Pseudogenes/OutForR_EvidencePred.txt.gz"),header=F)
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
lens = read.table("/projects/btl/kgagalova/PHD_projects2/SpruceUp/GenomeAnnotation/WS77111/data/Pseudogenes/proteinsPolishedEvidencePred.len")
colnames(lens) = c("prot_nam","len")
lens$prot_nam = as.character(lens$prot_nam)
d1 = merge(d,lens,by="prot_nam")
d1$perc_align = d1$len_align / d1$len
d1$contig1 = sapply(strsplit(d1$contig,":"),"[[",1)
d1$start_coord = as.numeric(sapply(strsplit(sapply(strsplit(d1$contig,"\\|"),"[[",2),"-"),"[[",1))
d1$stop_coord = as.numeric(sapply(strsplit(sapply(strsplit(d1$contig,"\\|"),"[[",2),"-"),"[[",2))
d1$contig2 = paste(d1$contig1,d1$start_coord,d1$stop_coord,sep="_")

stop_c = ifelse(d1$stop_high > 0, "stop_codon", "none")
framesh = ifelse(d1$frame_high > 0, "frameshift", "none")

chech_mut <-function(stopc,fr){
  type_mut = c()
  for(i in 1:length(stopc)){
    si = stopc[i]
    fi = fr[i]
    #print(si)
    #print(fi)
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

![](images/genes_lenEvidencePred-1.png)

``` r
summary(lens$len)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##      33     198     319     346     444    3493

Plot distributions
------------------

``` r
d1_counts = as.data.frame(table(d1$prot_nam))
d2_counts = merge(d1_counts,lens, by.x="Var1",by.y="prot_nam")
d2_counts = d2_counts[with(d2_counts, order(len)),]

ggplot(d2_counts, aes(len,Freq)) + geom_point(colour="darkorange",alpha=0.2) + ylab("# of pseudogenes") + xlab("Peptide length (aa)") +
  geom_density2d(colour="black") +     
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
          axis.title.y = element_text(face='bold',size=16,vjust=1),
          axis.text.x = element_text(face='bold',size=14,color='black'),
          axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/psedogenes_lenEvidencePred-1.png)

``` r
summary(d2_counts$Freq)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    1.00    2.00    5.00   11.39   12.00  305.00

``` r
ggplot(d1, aes(len,perc_align,colour=type_mut)) + geom_point(size = 3) + ylab("Pseudogene len - coverage orig pept") + xlab("Peptide length (aa)") + 
      theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
          axis.title.y = element_text(face='bold',size=16,vjust=1),
          axis.text.x = element_text(face='bold',size=14,color='black'),
          axis.text.y = element_text(face='bold',size=14,color='black'))
```

![](images/plot_cov_individual_pseudoEvidencePred-1.png)

``` r
ggplot(d1, aes(len,perc_align,colour=type_mut)) + geom_point(colour="darkgreen", alpha=0.2) + ylab("Pseudogene len - coverage orig pept") + xlab("Peptide length (aa)") + 
      geom_density2d(colour="black") + 
      theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
          axis.title.y = element_text(face='bold',size=16,vjust=1),
          axis.text.x = element_text(face='bold',size=14,color='black'),
          axis.text.y = element_text(face='bold',size=14,color='black'),
          strip.text = element_text(size=15) ) + 
      facet_wrap( ~ type_mut )
```

![](images/plot_cov_individual_pseudoEvidencePred-2.png)

``` r
table(d1$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                 34599                 24510                  1772 
    ## stop_codon/frameshift 
    ##                  4757

``` r
#matching bp with proteins and pseudogenes
sum(abs(as.numeric(sapply(strsplit(sapply(strsplit(d1$coord, ":"),"[[",2),"-"),"[[",1)))- abs(as.numeric(sapply(strsplit(sapply(strsplit(d1$coord, ":"),"[[",2),"\\-"),"[[",2))))
```

    ## [1] 139079

``` r
#genomic region covered by pseudogenes
sum(as.numeric(d1$stop_coord) - as.numeric(d1$start_coord))
```

    ## [1] 288953589

``` r
no_mutd1noreps = subset(d1,d1$type_mut == "none")
sum(as.numeric(no_mutd1noreps$stop_coord) - as.numeric(no_mutd1noreps$start_coord))
```

    ## [1] 75834990

``` r
stopd1noreps = subset(d1,d1$type_mut == "stop_codon")
framed1noreps = subset(d1,d1$type_mut == "frameshift")
stopframed1noreps = subset(d1,d1$type_mut == "stop_codon/frameshift")

sum(as.numeric(stopd1noreps$stop_coord) - as.numeric(stopd1noreps$start_coord))
```

    ## [1] 129884452

``` r
sum(as.numeric(framed1noreps$stop_coord) - as.numeric(framed1noreps$start_coord))
```

    ## [1] 13309328

``` r
sum(as.numeric(stopframed1noreps$stop_coord) - as.numeric(stopframed1noreps$start_coord))
```

    ## [1] 69924819

``` r
#filter by the % aligned
table(subset(d1,d1$perc_align > 0.5)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                  3788                  5109                   332 
    ## stop_codon/frameshift 
    ##                  1532

``` r
table(subset(d1,d1$perc_align > 0.7)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                  1449                  2513                   132 
    ## stop_codon/frameshift 
    ##                   764

``` r
#filter by e-value
table(subset(d1,d1$eval < 1e-50)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                  3481                  4277                   358 
    ## stop_codon/frameshift 
    ##                  1606

``` r
table(subset(d1,d1$eval < 1e-100)$type_mut)
```

    ## 
    ##                  none            stop_codon            frameshift 
    ##                   732                  1004                    88 
    ## stop_codon/frameshift 
    ##                   431
