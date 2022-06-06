annotationsPG29 <- read.csv("/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/ValidatePredictedGenesSecondStep/AED.remove_intron.CDS.codons/genome.proteins1.all.maker.transcripts.summary.tsv", header=T, sep = "\t")
annotationsWS <- read.csv("/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneAndModelsPredFirstStepMaker/ValidatePredictedGenesThirdStep/AED.remove_intron.CDS.codons2/genome.proteins1.all.maker.transcripts.summary.tsv", header=T, sep = "\t")
annotationsQ903 <- read.csv("/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/ValidatePredictedGenesSecondStepONT2/AED.remove_intron.CDS.codons/genome.proteins1.all.maker.transcripts.summary.tsv", header=T, sep = "\t")

annotationsPG29$genome <- rep("PG29", nrow(annotationsPG29)) 
annotationsWS$genome <- rep("WS77111", nrow(annotationsWS)) 
annotationsQ903$genome <- rep("Q903", nrow(annotationsQ903))

annotations_all <- rbind(annotationsPG29,annotationsWS,annotationsQ903)
nonZero_eAED <- subset(annotations_all, annotations_all$eAED < 1)

ggplot(nonZero_eAED, aes(x=eAED,col=genome)) +
  geom_line(stat="bin", binwidth = 0.01, size=1.5) +
  theme_bw() + xlab("eAED score") + ylab("Number of predicted transcripts") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=20)) 


###########################
annotationsWSrnaseq <- read.csv("/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneAndModelsPredFirstStepMaker/ValidatePredictedGenesThirdStepRNAseq/AED.remove_intron.CDS.codons/WS77111-v2_Maker.all.maker.transcripts.summaryGenVal.tsv", header=T, sep = "\t")
annotationsWS <- read.csv("/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneAndModelsPredFirstStepMaker/ValidatePredictedGenesThirdStep/AED.remove_intron.CDS.codons2/genome.proteins1.all.maker.transcripts.summary.tsv", header=T, sep = "\t")
annotationsWSrnaseq$genome <- rep("WSrnaseq", nrow(annotationsWSrnaseq)) 
annotationsWS$genome <- rep("WS", nrow(annotationsWS))

annotations_all <- rbind(annotationsWSrnaseq,annotationsWS)
nonZero_eAED <- subset(annotations_all, annotations_all$eAED < 1)

ggplot(nonZero_eAED, aes(x=eAED,col=genome)) +
  geom_line(stat="bin", binwidth = 0.01, size=1.5) +
  theme_bw() + xlab("eAED score") + ylab("Number of predicted transcripts") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=20)) 

######################
annotationsQ903 <- read.csv("/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/ValidatePredictedGenesSecondStepONT2/AED.remove_intron.CDS.codons/genome.proteins1.all.maker.transcripts.summary.tsv", header=T, sep = "\t")
annotationsQ903gff <- read.csv("/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/SecondStepONT2_alignGff/OutputMaker/LoweAED/Q903compl.all.maker.transcriptsLoweAED.summary.tsv", header=T, sep = "\t")
annotationsQ903$genome <- rep("Q903blast", nrow(annotationsQ903)) 
annotationsQ903gff$genome <- rep("Q903gff", nrow(annotationsQ903gff))

annotations_all <- rbind(annotationsQ903,annotationsQ903gff)
nonZero_eAED <- subset(annotations_all, annotations_all$eAED < 1)

ggplot(nonZero_eAED, aes(x=eAED,col=genome)) +
  geom_line(stat="bin", binwidth = 0.01, size=1.5) +
  theme_bw() + xlab("GAG - eAED score") + ylab("Number of predicted transcripts") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=20)) 


ggplot(nonZero_eAED, aes(x=AED,col=genome)) +
  geom_line(stat="bin", binwidth = 0.01, size=1.5) +
  theme_bw() + xlab("GAG - AED score") + ylab("Number of predicted transcripts") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=20)) 
