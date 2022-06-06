Protein sequences from NCBI serach 
```
esearch -db protein -query "Transcriptome mining, functional characterization, and phylogeny of a large terpene synthase gene family in spruce (Picea spp.)" | efetch -format acc > NCBIesearchProt.in
```

Nucleotide sequences from NCBI
```
esearch -db nucleotide -query "Transcriptome mining, functional characterization, and phylogeny of a large terpene synthase gene family in spruce (Picea spp.)" | efetch -format acc > NCBIesearchNucl.in
```

* 23 hits proteins
* 22 hits nucleotides

