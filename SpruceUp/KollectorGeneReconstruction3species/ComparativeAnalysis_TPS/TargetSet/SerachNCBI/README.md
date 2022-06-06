## Protein sequences from NCBI serach 

The search is performed with the following parameters:

```
esearch -db protein -query "[prenyl transferase OR monoterpene synthase OR sesquiterpene synthase OR diterpene synthase OR P450] AND pinales NOT putative NOT hypothetical NOT partial NOT unknown NOT BAC" | efetch -format acc > NCBIesearchProt.in
```

283 hits

## Search refinement

pdb entries are removed from the nucleotide fasta sequences - final sequences 273

