# Cd-hit-est iterative

## Motivation

Cd-hit is a greedy algorithm for redundancy removal in a given data set.     
The sequences are first sorted by length and the algorithm starts with the longest sequence. In case the next sequence has lower similarity than the threshold, it creates a separate cluster: in case similarity is higher, the sequence is clustered together. 
There are several limitatins introduced by the greedy incremental clustering. 

Let's say there are 3 custers

```
1) 
A =========================== * (representative)
X ==========================
Y ==========================

2) 
W ========================== * (representative)

3) 
B ========================== *	(representative) 
Z ==========================

```
The problem in this case is that even if Y is more similar to B than to A, it happen to be in cluster 1 because Y first hits A during clustering.    

This can be improved by multi-step clustering.  

## Script - RunCdHit

The script runs multi-step iterative *cd-hit-est*.     

```
Usage: ./RunCdHit <list_param> <inputfile>
```

**n** and **c** parameters are given by the user in <list_param>. Those are written one per line, giving the sequence similarity threshold (*c*) and word size (*n*) for each iterations separated by space. An example can be found in *list_iter.in*    

The input file is given as second argument to the script.      

Please check cd-hit documentation [here](http://www.bioinformatics.org/cd-hit/cd-hit-user-guide.pdf) 
