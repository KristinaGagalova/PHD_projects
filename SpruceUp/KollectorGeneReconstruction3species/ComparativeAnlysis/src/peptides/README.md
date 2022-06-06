## Run ORFfinder on a list of cDNAs

ORF finder searches for open reading frames (ORFs) in the DNA sequence you enter. The program returns the range of each ORF, along with its protein translation. Use ORF finder to search newly sequenced DNA for potential protein encoding segments, verify predicted protein using newly developed SMART BLAST or regular BLASTP.     

ORFfinder version [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/)

```
USAGE
  ORFfinder [-h] [-help] [-xmlhelp] [-in Input_File] [-id Accession_GI]
    [-b begin] [-e end] [-c circular] [-g Genetic_code] [-s Start_codon]
    [-ml minimal_length] [-n nested_ORFs] [-strand Strand] [-out Output_File]
    [-outfmt output_format] [-logfile File_Name] [-conffile File_Name]
    [-version] [-version-full] [-dryrun]

DESCRIPTION
   Searching open reading frames in a sequence

OPTIONAL ARGUMENTS
 -h
   Print USAGE and DESCRIPTION;  ignore all other parameters
 -help
   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
 -xmlhelp
   Print USAGE, DESCRIPTION and ARGUMENTS in XML format; ignore all other
   parameters
 -logfile <File_Out>
   File to which the program log should be redirected
 -conffile <File_In>
   Program's configuration (registry) data file
 -version
   Print version number;  ignore other arguments
 -version-full
   Print extended version data;  ignore other arguments
 -dryrun
   Dry run the application: do nothing, only test all preconditions

 *** Input query options (one of them has to be provided):
 -in <File_In>
   name of file with the nucleotide sequence in FASTA format
   (more than one sequence is allowed)
   Default = `'
 -id <String>
   Accession or gi number of the nucleotide sequence
   (ignored, if the file name is provided)
   Default = `'

 *** Query sequence details:
 -b <Integer>
   Start address of sequence fragment to be processed
   Default = `1'
 -e <Integer>
   Stop address of sequence fragment to be processed (0 - to the end of the
   sequence)
   Default = `0'
 -c <Boolean>
   Is the sequence circular? (t/f) *** Under development
   Default = `false'

 *** Search parameters:
 -g <Integer>
   Genetic code to use (1-31)
   see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details
   Default = `1'
 -s <Integer>
   ORF start codon to use:
       0 = "ATG" only
       1 = "ATG" and alternative initiation codons
       2 = any sense codon
   Default = `1'
 -ml <Integer>
   Minimal length of the ORF (nt)
   Value less than 30 is automatically changed by 30.
   Default = `75'
 -n <Boolean>
   Ignore nested ORFs (completely placed within another)
   Default = `false'
 -strand <String>
   Output ORFs on specified strand only (both|plus|minus)
   Default = `both'

 *** Output options:
 -out <File_Out>
   Output file name
 -outfmt <Integer>
   Output options:
       0 = list of ORFs in FASTA format
       1 = CDS in FASTA format
       2 = Text ASN.1
       3 = Feature table
   Default = `0'
```
