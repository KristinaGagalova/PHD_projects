library( dplyr )
library( ggplot2 )
library( tidyr )
library( data.table )
library( reshape )
library( UpSetR )

dataPath="/home/kgagalova/MissingBusco"
allFiles <- list.files( path = dataPath, pattern = ".in", full.names = TRUE )

l <- lapply( allFiles, function( fn ){
  d <- read.table( fn, header = F );
  d$fileName <- fn;
  d
} );

missing <- bind_rows( l );
dim(missing)

missing$nam = gsub("missing.in","",gsub("missing_busco_list_","",gsub("1.tsv","",sapply(strsplit(sapply(strsplit(missing$fileName,"/"),"[[",5),"\\="),"[[",1),"_")))

missing_long = unstack(missing[,c("V1","nam")])

upset(fromList(missing_long),order.by = "freq",nsets=8)

#movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), 
#                    header=T, sep=";" )
