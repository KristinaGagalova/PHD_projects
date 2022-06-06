## Parse Yaml format (my brain is exploding!)

ParseYaml.py script, uses **yaml** and **csv** libraries     

The script parses one type element per line, similar to gff file. Th eelement can be an exon, intron, intron? and gap. If they belong to the same alignment, the first 12 columns of the table are the same.      
Note that in case of dictionaries or list for each element in the 'matchings', those are compacted together.

**Output type:** 
```
bin - file name 
id - alignment id
target - protein name
prot_len - target len
contig - cDNA aligned
contig_len - contig len	
strand
score - Scipio score	
status - alignment outcome
reason	- more info about alignments
seq_upstream - seq before protein 
seq_downstream	- seq after protein
############################################
#elements from 'matchings' entry - Scipio
dna_end	
dna_start	
inframe_stopcodons	
mismatchlist	
nucl_end	
nucl_start	
overlap	
prot_end	
prot_start	
seq	
seqshifts	
translation	
type	
undeterminedlist
```

More infor about Scipio output are [here](http://www.webscipio.org/help/scipio)

**Example script run**
```
#ParseYaml.py <file_name.yaml> <output_dir>
python ParseYaml.py Scipio_WS_EasterRun_cdhit-output-4_8_iteration.5.yaml ./
```



