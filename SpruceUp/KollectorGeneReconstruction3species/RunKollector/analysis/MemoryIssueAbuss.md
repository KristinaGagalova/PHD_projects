
### Error message

```
Extended 62180307 of 407277000 reads (15.3%), assembled 3010529774 bp so far
Extended 62180372 of 407278000 reads (15.3%), assembled 3010530490 bp so far
Extended 62180427 of 407279000 reads (15.3%), assembled 3010530914 bp so far
Extended 62180443 of 407279252 reads (15.3%), assembled 3010531236 bp so far
Assembly complete
AdjList -v   -k32 -m50 --dot trial_kollectorHalf-1.fa >trial_kollectorHalf-1.dot
Reading `trial_kollectorHalf-1.fa'...
terminate called after throwing an instance of 'std::bad_alloc'
  what():  St9bad_alloc
/bin/bash: line 1: 16099 Aborted                 (core dumped) AdjList -v -k32 -m50 --dot trial_kollectorHalf-1.fa > trial_kollectorHalf-1.dot
make[1]: *** [trial_kollectorHalf-1.dot] Error 134
make[1]: *** Deleting file `trial_kollectorHalf-1.dot'
make[1]: Leaving directory `/extscratch/btl/kgagalova/trial_kollectorHalf1000/trial_kollectorHalf.abyss'
Command exited with non-zero status 2
```


## Summary table: memory issue Abyss (1)

| Name   | Reads type | Tot reads | Bin size | Reads in Abyss | Kollector parameters             | abyss-bloom-dbg | Status           | Time |
|--------|------------|-----------|----------|----------------|----------------------------------|-----------------|------------------|------|
| Run #6 | PET        | 3.1G	  | 1000     | 400M           | **r=0.7,s=0.5**,k=32,K=25,n=50000000 | b30G,H4,kc=3    | 'std::bad_alloc' |82.9h |
| Run #7 | PET        | 3.1G	  | 1000     | 76M            | **r=0.9,s=0.7**,k=32,K=25,n=50000000 | b30G,H4,kc=3    |                  |12.5h |



**abyss-fac**

| n	  | n:500 | L50  | min | N80 | N50 | N20  | E-size | max  | sum     | name                     |
|---------|-------|------|-----|-----|-----|------|--------|------|---------|--------------------------|
| 73.58e6 | 5002  | 1689 | 500 | 564 | 758 | 1263 | 991    | 4869 | 3923458 | trial_kollectorHalf-1.fa |
| 9800609 | 3233  | 1021 | 500 | 610 | 874 | 1575 | 1164   | 6172 | 2843619 | trial_kollectorHalf-1.fa |


**seq length statistics**

| Len contig | 32   | 33-35 | 36-45 | 46-65 | 66-120 | 121-180 |
|------------|------|-------|-------|-------|--------|---------|
| %          | 37.8 | 20.9  | 21.3  | 15.4  | 3      | 6.3     |

Large amount of contigs with length equal to kmer

**REMOVE CONTIGS WITH SIZE == 32 (KMER SIZE): no memory issues**         


## Summary table: memory issue Abyss (2)

| Name   | Reads type | Tot reads | Bin size | Reads in Abyss | Kollector parameters             | abyss-bloom-dbg | Status           | Time |
|--------|------------|-----------|----------|----------------|----------------------------------|-----------------|------------------|------|
| Run #8 | PET        | 3.1G	  | 1000     | 400M           | r=0.7,s=0.5,k=116,K=25,n=50000000 | b7000M,H4,kc=3    | 'std::bad_alloc' |82.9h |
| Run #9 | PET        | 3.1G	  | 1000     | 76M            | r=0.9,s=0.7,k=116,K=25,n=50000000 | b2000M,H4,kc=3    |                  |75h |

**abyss-fac**

| n       | n:500 | L50   | min | N80 | N50 | N20  | E-size | max  | sum     | name                     |
|---------|-------|-------|-----|-----|-----|------|--------|------|---------|--------------------------|
| 31.15e6 | 54382 | 20363 | 500 | 549 | 673 | 987  | 843    | 6160 | 38.65e6 | trial_kollectorHalf-1.fa |
| 3492169 | 8738  | 2671  | 500 | 591 | 831 | 1618 | 1251   | 7847 | 7536026 | trial_kollectorHalf-1.fa |



**seq length statistics**

| Len contig | 116  | 117-119 | 120-130 | 131-150 | 151-300 | 301-450 |
|------------|------|---------|---------|---------|---------|---------|
| %          | 13.8 | 12      | 21.5    | 21.8    | 29.2    | 1.1     |


**REMOVE CONTIGS WITH SIZE == 116 (KMER SIZE): memory issues**       
**REMOVE CONTIGS WITH SIZE == 116..119 (25.8% of contigs removed): memory issues**     
**REMOVE CONTIGS WITH SIZE == 116..130 (47.3% of contigs removed): memory issues**      

## Observations

- contigs removal is not always a solution. In the case of kmer size == 32 it worked fine but for a more heterogenoues landscape of contigs it is not solving the memory issues (as in the case of kmer == 116)
