Check overlap targets - 3 species
================

Calculate the overlap between the common targets reconstructed in the 3 species
-------------------------------------------------------------------------------

``` r
library(Biostrings)
library(VennDiagram)

WS <- read.table("/projects/spruceup/pglauca/WS77111/assemblies/comparative_genomics/kollector/kollector_runs/kollector_cdhit4AllTargets/assembledtargets.txt", quote="\"")
PG29 <- read.table("/projects/spruceup/interior_spruce/PG29/assemblies/comparative_genomics/kollector/kollector_runs/kollector_cdhit4AllTargets/assembledtargets.txt", quote="\"")
Q903 <- read.table("/projects/spruceup/psitchensis/Q903/assembly/comparative_genomics/kollector/kollector_runs/kollector_cdhit4AllTargets/assembledtargets.txt", quote="\"")
venn.plot <- draw.triple.venn(area1           = length(unique(as.character(PG29[,1]))),
                              area2           = length(unique(as.character(WS[,1]))),
                              area3           = length(unique(as.character(Q903[,1]))),
                              n12             = length(intersect(unique(as.character(PG29[,1])), unique(as.character(WS[,1])))),
                              n23             = length(intersect(unique(as.character(WS[,1])), unique(as.character(Q903[,1])))),
                              n13             = length(intersect(unique(as.character(PG29[,1])), unique(as.character(Q903[,1])))),
                              n123            = length(Reduce(intersect, list(as.character(PG29[,1]),as.character(WS[,1]),as.character(Q903[,1])))),
                              category        = c("PG29", "WS77111","Q903"),
                              fill            = c("red", "blue","green"),
                              lty             = "blank",
                              cex             = 2,
                              cat.cex         = 2,
                              euler.d         = TRUE,
                              scaled          = TRUE, 
                              #cat.pos         = c(-110, 100),
                              #cat.dist        = c(0.2,0.05),
                              #cat.just        = list(c(-1, -1), c(1, 1)),
                              #ext.pos         = c(20,20),
                              #ext.dist        = c(1,1),
                              #ext.length      = 0.85,
                              #ext.line.lwd    = 2,
                              ext.line.lty    = "dashed"
)
```

![](images/all_targets-1.png)

Excluding gaps
--------------

``` r
PG29g <- read.table("/projects/spruceup/interior_spruce/PG29/assemblies/comparative_genomics/kollector/kollector_runs/kollector_cdhit4AllTargets/Sealer/assembledtargets_nogapsPG29.txt", quote="\"")
WSg <- read.table("/projects/spruceup/pglauca/WS77111/assemblies/comparative_genomics/kollector/kollector_runs/kollector_cdhit4AllTargets/Sealer/assembledtargets_nogapsWS.txt", quote="\"")
Q903g <- read.table("/projects/spruceup/psitchensis/Q903/assembly/comparative_genomics/kollector/kollector_runs/kollector_cdhit4AllTargets/Sealer/assembledtargets_nogapsQ903.txt", quote="\"")

venn.plot <- draw.triple.venn(area1           = length(unique(as.character(PG29g[,1]))),
                              area2           = length(unique(as.character(WSg[,1]))),
                              area3           = length(unique(as.character(Q903g[,1]))),
                              n12             = length(intersect(unique(as.character(PG29g[,1])), unique(as.character(WSg[,1])))),
                              n23             = length(intersect(unique(as.character(WSg[,1])), unique(as.character(Q903g[,1])))),
                              n13             = length(intersect(unique(as.character(PG29g[,1])), unique(as.character(Q903g[,1])))),
                              n123            = length(Reduce(intersect, list(as.character(PG29g[,1]),as.character(WSg[,1]),as.character(Q903g[,1])))),
                              category        = c("PG29", "WS77111","Q903"),
                              fill            = c("red", "blue","green"),
                              lty             = "blank",
                              cex             = 2,
                              cat.cex         = 2,
                              euler.d         = TRUE,
                              scaled          = TRUE, 
                              #cat.pos         = c(-110, 100),
                              #cat.dist        = c(0.2,0.05),
                              #cat.just        = list(c(-1, -1), c(1, 1)),
                              #ext.pos         = c(20,20),
                              #ext.dist        = c(1,1),
                              #ext.length      = 0.85,
                              #ext.line.lwd    = 2,
                              ext.line.lty    = "dashed"
)
```

![](images/all_targets_nogaps-1.png)
