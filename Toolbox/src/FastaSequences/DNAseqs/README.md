# Basic toolkit for DNA sequences

## Austin Hammond

To use with Python 2.7 or higher      

Example of use
```
-bash-3.2$ FASTA=/projects/bullfrog_assembly/genome/annotation/maker/genomic-fastas/RC-genome-V2-20160223.submission.fa 
-bash-3.2$ grep -A1 -m1 "Rc-01r160223s0597824" $FASTA|seqselect 536 999 | revcomp | dna2aa -f "+1"
>Rc-01r160223s0597824_536-999_revcomp_frame=+1
EEHIFRQCTCLVSTY--ESSVP-HETHQNASTGYSIRDSTLWPPVVRFSFSMLLWS-NNSHCVLKIHIFVTSEIQFLCNNQCFLCLLPGEMPMKLKK-KEVSWISSRIWGRRLLGICWIR-NVQLEHVHPPPESEVI-YGISAKSLMSY-N-TY
```
**Tips:** use single line entry in fasta.    

How to change fasta multiple lines entry to single line:

```
-bash-3.2$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < file.fa | sed '/^$/d' > out.fa
```

Comparing sequences to find exact, full-length matches
compare-seqs.py reads in old scaffolds to dictionary (hash table) with sqn as key and id as value.
Then reads in new scaffoldss and compare sqn to keys and write out old and new ids upon match.
Also writes out file collisions.txt if it experiences a collision while populating the hash table.
Example
```
$ compare-seqs.py old.fa new.fa > seq-compare.tsv
```
