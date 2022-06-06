#!/usr/bin/Rscript

#################################
####Kristina Gagalova############
#################################
########15 Oct 2017##############
#################################

#Description: print coverage identity

#to find in data - in git
d <- read_delim("/media/kgagalova/SD2561/Documents/Spruce/TPS_kollector/PreliminaryStats/PG29/CovIdentityPG29.txt", 
                "\t", escape_double = FALSE, col_names = FALSE, 
                trim_ws = TRUE)

d1 = d[order(d$target, d$coverage, decreasing=T), ]
d1 = d1[ !duplicated(d1$target), ]  

commonTheme = list(labs(color="Density",fill="Density",
                        x="Coverage",
                        y="Identity"),
                   theme_bw(),
                   theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
                         axis.title.y = element_text(face='bold',size=16,vjust=1),
                         axis.text.x = element_text(face='bold',size=14,color='black'),
                         axis.text.y = element_text(face='bold',size=14,color='black'),
                         legend.position="none",
                         legend.justification=c(0,1)))

ggplot(data=d,aes(coverage,identity, fill = species)) + 
  geom_density2d(aes(colour=..level..)) + 
  scale_colour_gradient(low="green",high="red") + 
  geom_point(size=4) + commonTheme

d$species = rep("Other",nrow(d))
d$species[grep("PiceaglaucacloneWS",d$target)] = "WS"
d$species[grep("PiceaengelmanniixPiceaglaucaclone",d$target)] = "PG29"
d$species[grep("Piceasitchensisclone",d$target)] = "Q903"
d$species = as.factor(d$species)
d$species = ordered(d$species, levels = c("PG29", "WS", "Q903","Other"))

png("BoxplotCoverage.png")
ggplot(d1, aes(x=species, y=coverage, fill=species)) +
  geom_boxplot() + 
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(face='bold',size=14,color='black'),
        legend.justification=c(0,0.5)) +
  geom_point(size=4,alpha=0.4)
dev.off()

png("ScatterCovIdentity.png")
ggplot(d1, aes(x=coverage, y=identity, color=species)) +
  geom_point(size=4) + 
  theme(axis.title.x = element_text(face='bold',size=16,hjust=0.5),
        axis.title.y = element_text(face='bold',size=16,vjust=1),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(face='bold',size=14,color='black'),
        legend.justification=c(0,0.5))
dev.off()
