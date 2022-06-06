# Directories description

- BLAST - contains the self all-vs-all BLAST. Generates a database and blasts the sequences against it    

- ParseBLAST - re-shapes the BLAST output. Outputs each query sequence in a tab delimited format against all the subject sequences (and the corresponding identity %, e-values and query coverages).

- RunPipeline - BLAST and ParseBLAST in unique pipeline format

- SelectRandomFasta - select randomly sequences from fasta file. Output from this script be compared to the cdhit file output file in order to evaluate whether cdhit really reduces the redundancy or whether there is only a number effect

- CreateBins - divide the fasta sequences in bins of 1000 elements, after that remove spaces in the header
