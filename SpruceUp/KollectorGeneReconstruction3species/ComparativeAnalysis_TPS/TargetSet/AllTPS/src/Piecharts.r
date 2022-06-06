#!/usr/bin/Rscript

#################################
####Kristina Gagalova############
#################################
########15 Oct 2017##############
#################################

#Description: print coverage identi

### Piechart

nams <- read.delim("/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/ComparativeAnalysis_TPS/TargetSet/AllTPS/data/list_nams.in", header=FALSE)
nams$V1 <- as.character(nams$V1)

##Species
WS = grep("1PiceaglaucacloneWS|1_Picea_glauca_clone_",nams$V1)
PG29 = grep("PiceaengelmanniixPiceaglaucaclone|maker|snap|TSA:_Picea_glauca_|_Picea_engelmannii_x_Picea_glauca_clone",nams$V1)
Q903 = grep("Picea_sitchensis|Piceasitchensisclone",nams$V1)
Pinus = grep("_Pinus",nams$V1)
Abies = grep("Abies|Picea_abies",nams$V1)
Other = grep("Taxus|_Chamaecyparis|Pseudotsuga|Ginkgo|Pseudolarix",nams$V1)
all = c(WS, PG29, Q903, Pinus, Abies, Other)
nams$V1[-all]

species <- c("PG29", "WS", "Q903", "Pinus", "Abies", "Other picea")
numb <- c(107, 12, 41, 64, 29, 26)
pct <- round(numb/sum(numb)*100)
lbls <- paste(species, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
png("SpeciesPicea.png")
pie(numb,labels = lbls, col=rainbow(length(lbls)),
    main="Picea species")
dev.off()

##Sources
Trinity = grep("Trinity",nams$V1)
Abyss = grep("_k60_",nams$V1)
ESTcomplete = grep("_complete_cds|completecds",nams$V1)
#ESTpartial = grep("_partial_cds",nams$V1)
HCG = grep("maker|snap",nams$V1)
Other = grep("_partial_cds|(CPR1_gene)|_exons_1-13",nams$V1)

allSource = c(Trinity,Abyss,ESTcomplete,HCG, Other)
nams$V1[-allSource]

species <- c("Trinity", "Abyss", "FLcDNA", "HCG", "Other")
numb <- c(72, 18, 178, 9, 3)
pct <- round(numb/sum(numb)*100)
lbls <- paste(species, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
png("SourcePicea.png")
pie(numb,labels = lbls, col=rainbow(length(lbls)),
    main="Source targets")
dev.off()

##Target types
d2 <- read.csv("/projects/btl/kgagalova/PHD_projects/SpruceUp/KollectorGeneReconstruction3species/ComparativeAnalysis_TPS/TargetSet/AllTPS/data/MetadataEnzymesCurated.txt", header=FALSE)
colnames(d2) = c("id","enzyme","type")
d2$collapsed <- apply( d2[,c("enzyme","type")] , 1 , paste , collapse = "_" )
mono = grep("monoterpenes|apin1",d2$collapsed)
di = grep("diterpene",d2$collapsed)
sesqui = grep("sesquiterpene", d2$collapsed)
cyp = grep("P450|p450|CYP7",d2$collapsed)
others = grep("prenyltransferase|geranylgeranylpyrophosphatesynthase|tax|Taxus|gibberellin|13-alpha-hydroxylase",d2$collapsed)

allType = c(mono,di,sesqui,cyp, others)
d2$enzyme[-allType]

types <- c("monoterpene", "diterpene", "sesquiterpene", "CYP450", "Other")
numb <- c(76, 26, 24, 130, 22)
pct <- round(numb/sum(numb)*100)
lbls <- paste(types, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
png("TypeEnz.png")
pie(numb,labels = lbls, col=rainbow(length(lbls)),
    main="TPS groups")
dev.off()
