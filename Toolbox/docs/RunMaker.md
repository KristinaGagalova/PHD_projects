## How to efficiently and accurately run Maker for genome annotation

### Resources

[Tutorial Make](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018)      
[Maker wiki](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Main_Page)

### Procedure

#### Format

- Fasta sequences need to have short header, this is a limitation given by RepeatMasker

#### Repeat masking with repeated files

In case you may need to re-run the same annotation several times, it is better to generate the repeat masking only once. After it's generated, parse the gff files for the repeatmakser sequences and add the new gff to this option in maker_opts.ctl
```
rm_gff= #pre-identified repeat elements from an external GFF3 file
```

#### Filtering step is better before starting Maker

Mail from Austin

```
I get the impression from maker that it assumes both a highly contiguous genome assembly and that the transcript evidence that you give it represents complete sequences. 
I think that in some cases, giving it potentially incomplete protein-coding transcripts or lncRNAs leads it astray, so it predicts a (protein coding) gene where there shouldn't be one, or a partial one. 
It would be worthwhile to filter transcript and protein evidence before giving it to maker. 
I'll remember this for our next whole-genome annotation run.
```

One way is to map the transripts to the assembled genome and select only those that have a good match.

#### Selection of output predictions

Annotatios that begin with "maker-" are more reliablesione supported by high evidence. eAED can be used to rank the predictions, but alignments to a high quality set can be more informative than filtering on eAED alone.

Types of output predictions, more details can be found [here](https://groups.google.com/forum/#!topic/maker-devel/hxETPzlZoIk)

#### How I run Maker on the Kollector reconstructed targets

- I have added a protein data base. I have used all type of sequences, predicted and annotated. Consider tha the db query can be quite slow if the db is large, filter the sequences that you don't consider relevant for the annotation

```
#example- proteins from only selected species, I have avioded bacterial and viral proteins that make aroiund 75% of all the protein sequences
protein=/projects/spruceup_scratch/dev/MakerKollectorContigs/ProteinLibraries/SelectedProt/Seluniprot_tremblsprotOneline1.fasta  #protein sequence file in fasta format (i.e. from mutiple oransisms)
```

- I have used annotation for arabidopsis

```
augustus_species=arabidopsis
```

- Used protein and EST evidence

```
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
```

- Look for high quality genes with start and stop

```
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
```

The output was parsed for only maker evidence and valifated with [genevalidator](https://wurmlab.github.io/tools/genevalidator/)
