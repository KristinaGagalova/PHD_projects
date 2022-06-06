## Polisg trimmed reads output


1) fastqCombinePairedEnd.py - divides the trimmed/filtered reads in 3 files, correctly paired R1 and R2 + singletons  

- Runs with Python 2.7 =< 

```
Resynchronize 2 fastq or fastq.gz files (R1 and R2) after they have been
trimmed and cleaned

WARNING! This program assumes that the fastq file uses EXACTLY four lines per
    sequence

Three output files are generated. The first two files contain the reads of the
    pairs that match and the third contains the solitary reads.

Usage:
    python fastqCombinePairedEnd.py input1 input2 separator

input1 = LEFT  fastq or fastq.gz file (R1)
input2 = RIGHT fastq or fastq.gz file (R2)
separator = character that separates the name of the read from the part that
    describes if it goes on the left or right, usually with characters '1' or
    '2'.  The separator is often a space, but could be another character. A
    space is used by default. If the sequence names do not contain two parts
    and you want to use the full name info to pair your sequences, use 'None'
    (as text) for the separator. Eg:
        python fastqCombinePairedEnd.py input1 input2 None
```

2) SortFastqTEMPLATE.sh - sorts the fastq file

- Important: in order to be read by BBT maker, the reads need to be sorted in the same way. The script contains this variable ```export LC_ALL=C``` which sorts correctly. 

```
#incorrect sort below
## R1
READNEM00001
READNAME00011
...
## R2
READNAME00011
READNAME0001
...
```

- sort uses a temporary directory which may not be sufficient for large files. set ```-T /my/new/tmp/dir``` to a directory with enough space 
