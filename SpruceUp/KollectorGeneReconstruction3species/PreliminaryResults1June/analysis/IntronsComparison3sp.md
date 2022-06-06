Introns comparison PG29, WS77111 and Q903
================

Given the reconstructed gene aligned to the transcript, extract intron info from the CIGAR - GMAP alignment.
Filter any alignmnet which:
- Has soft/hard clipping - Has deletions or insertions

Load the dataset and process the CIGAR strings
----------------------------------------------

``` r
library( dplyr )
library( ggplot2 )
library( tidyr )
library( data.table )
library(qdapRegex)

dataPath="/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/PreliminaryResults1June/data"
allFiles <- list.files( path = dataPath, pattern = "_sam.txt", full.names = TRUE )

l <- lapply( allFiles, function( fn ){
  d <- read.table( fn, header = F );
  d$fileName <- fn;
  d
  } );

d <- bind_rows( l );
d$species = gsub("_sam.txt","",gsub("OverlapTranscripts","",(sapply(strsplit(d$fileName,"/"),"[[",10))))

names(d)[1:8] = c("bin","coverage","transcript","strand","gene","cigar","MD","NM")

#1) soft/hard clipped 
d_sh=d[sapply(regmatches(d$cigar, gregexpr("S|H",d$cigar)), length) == 0,]
#2) exclude deletions insertions from the alignment 
d_di=d_sh[sapply(regmatches(d_sh$cigar, gregexpr("D|I",d_sh$cigar)), length) == 0,]

#get number of exons
d_exons=sapply(regmatches(d_di$cigar, gregexpr("M",d_di$cigar)), length)
# single exons genes
d_single=d_di[d_exons==1,]
table(d_single$species)
```

    ## 
    ##            PG29            Q903 WS77111extraAll 
    ##            2183            1656            1914

``` r
# multiple exons genes
d_multi=d_di[d_exons>1,]
table(d_multi$species)
```

    ## 
    ##            PG29            Q903 WS77111extraAll 
    ##            2589            1520            2521

``` r
#extract intron lengths
get_intron_len <- function(mycigar) {
  int_introns = paste(rm_between(mycigar, 'M', 'N', extract=TRUE)[[1]],collapse=",")
  return(int_introns)
}


d_di$introns_len=sapply(d_di$cigar, get_intron_len)
d_di$introns_num=sapply(gregexpr(",", d_di$introns_len, fixed = TRUE), function(x) sum(x > -1)) + 1
d_di$introns_num[d_di$introns_len=="NA"] <- 0

#functon to swap strands when antisense == 16
#enter here the two columns of the data frame, one with the strand infor and the other one with the introns length
swap_same <-function(df){
  rev_int=c()
  for(row in 1:nrow(df)){
    cur_row = df[row,]
    if(cur_row[1]==0){
      rev_int = c(rev_int,unname(cur_row[2]))
    }else{
      rev_int=c(rev_int, paste(rev(as.numeric(strsplit(as.character(cur_row[2]),",")[[1]])),collapse=","))
  }
 }
 rev_int = sapply(rev_int,"[[",1)
 return(rev_int)
}

d_di$swap.introns_len = swap_same(d_di[,c("strand","introns_len")])
#export to file
write.table(d_di,"AllIntrons3sp.txt",row.names=F,quote=F)

#get only the overlapping in all the species
df_d <- as.data.frame.matrix(table(d_di$transcript,d_di$species))
overlap_all <- row.names(df_d[which(apply(df_d,1,sum)==3),])

df_d1 <- d_di[d_di$transcript %in% overlap_all,]
#change introns number to 0 whether NA
write.table(df_d1,"AllIntrons3spOverlap.txt",row.names=F,quote=F)
```
