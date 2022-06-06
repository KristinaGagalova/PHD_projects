# Genesis test runs - WS77111

## Run #0 13 Mar 2017

Small subset of reads (2 PET) and small bin

- Directory: ```/extscratch/btl/kgagalova/trial_kollector```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Read files location:
```
/extscratch/spruce_assembly/PG/PG77111/genome/data/gsc/PET/IX2298-1..1_lane1_ATCACG_L001_R[12]_001.fastq.gz 
/extscratch/spruce_assembly/PG/PG77111/genome/data/gsc/PET/IX2298-1..1_lane1_CGATGT_L001_R[12]_001.fastq.gz
```
- Type of reads: PET
- Number of files: 2
- Number of reads: 94935192 (94M)

- Seed bin size: 500
- Seed Seq statistics:
```
n       n:500   L50     min     N80     N50     N20     E-size  max     sum     name
500     491     182     500     1039    1437    1836    1468    2850    649937  /extscratch/btl/kgagalova/trial_kollector/seed_500/cdhit-output-4_1

```

- Kollector options: j=12, r=0.7, s=0.5, k=32, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b2G -H3 -j12 --kc=3
- Reads input in abyss 2.0.2: 19544
- Bloom filter FPR: 5.45e-13%

- Repeat filter: NONE

- time: **kollector-recruit.mk** - 1803.66 (18min), **abyss-pe** - 30.88, **kollector-extract.sh** - 74.26

- **Status**: no output from kollector, coverage too low after alignment with GMAP. Amount of sequences too low after Bloom filter


## Run #1 - 7 Mar 2017

"Big" run: full data set of reads into Kollector as crash test    

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorFull1000```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G    

- Files location:   
```
/extscratch/spruceup/pglauca/WS77111/data/reads/WS77111-reads-1.in    
/extscratch/spruceup/pglauca/WS77111/data/reads/WS77111-reads-2.in     
```
- Type of reads: PET, MiSeq, Chromium    
- Number of files: 54    
- Number of reads: 6178931058 (6.1B) 
 
- Seed bin size: 1000   
- Seed Seq statistics:
```
n	n:500	L50	min	N80	N50	N20	E-size	max	sum	name
1000	986	365	500	1056	1452	1873	1492	3009	1321600	/extscratch/btl/kgagalova/trial_kollector/seed_1000/cdhit-output-4_1.fasta
```

- Kollector options: j=12, r=0.7, s=0.5, k=32, K=25, n=50000000       

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b2G -H3 -j12 --kc=3     
- Reads input in abyss 2.0.2: 1750672376 (1.7B)    
- Bloom filter FPR: 97.6%     

- Repeat filter: NONE    

- **Status**: job killed, memory requirements too high     

- time: **kollector-recruit.mk** - 62298.92 (17.3h) , **abyss-pe** - , **kollector-extract.sh** - 

## Run #2 - 13 Mar 2017

Same as Run #1 but with smaller seed bin size (half the previous).

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorFull500```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/spruceup/pglauca/WS77111/data/reads/WS77111-reads-1.in
/extscratch/spruceup/pglauca/WS77111/data/reads/WS77111-reads-2.in
```
- Type of reads: PET, MiSeq, Chromium
- Number of files: 54
- Number of reads: 6178931058 (6.1B)

- Seed bin size: 500
- Seed Seq statistics:
```
n	n:500	L50	min	N80	N50	N20	E-size	max	sum	name
500	491	182	500	1039	1437	1836	1468	2850	649937	/extscratch/btl/kgagalova/trial_kollector/seed_500/cdhit-output-4_1
```

- Kollector options: j=12, r=0.7, s=0.5, k=32, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b30G -H4 -j12 --kc=3
- Reads input in abyss 2.0.2: 1856425512 (1.8B)
- Bloom filter FPR: 0.0255%

- Repeat filter: NONE

- time: **kollector-recruit.mk** - 63671.29 (17.7h) , **abyss-pe** - , **kollector-extract.sh** -

- **Status**: job killed, memory requirements too high 

## Run #3 - 14 Mar 2017

Run Kollector only on PET reads

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorHalf500```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/btl/kgagalova/trial_kollectorHalf500/WS77111-reads-1pet.in
/extscratch/btl/kgagalova/trial_kollectorHalf500/WS77111-reads-2pet.in
```
- Type of reads: PET
- Number of files: 38
- Number of reads: 3198709675 (3.1G)

- Seed bin size: 500
- Seed Seq statistics:
```
n	n:500   L50     min     N80     N50     N20     E-size  max     sum     name
500     491     182     500     1039    1437    1836    1468    2850    649937  /extscratch/btl/kgagalova/trial_kollector/seed_500/cdhit-output-4_1
```

- Kollector options: j=12, r=0.7, s=0.5, k=32, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b30G -H4 -j12 --kc=3
- Reads input in abyss 2.0.2: 998545322 (0.9G)
- Bloom filter FPR: 0.00587%

- Repeat filter: NONE

- time: **kollector-recruit.mk** - 37056.27 (10.3h) , **abyss-pe** - , **kollector-extract.sh** -

## Run #4 - 14 Mar 2017

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorHalf500_1```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/btl/kgagalova/Reads/WS77111-reads-1pet.in
/extscratch/btl/kgagalova/Reads/WS77111-reads-2pet.in
```
- Type of reads: PET
- Number of files: 38
- Number of reads: 3198709675 (3.1G)

- Seed bin size: 500
- Seed Seq statistics:
```
n	n:500   L50     min     N80     N50     N20     E-size  max     sum     name
500     491     182     500     1039    1437    1836    1468    2850    649937  /extscratch/btl/kgagalova/trial_kollector/seed_500/cdhit-output-4_1
```

- Kollector options: j=12, r=0.7, s=0.5, k=32, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b30G -H4 -j12 --kc=3
- Reads input in abyss 2.0.2: 468026002 (0.48B) 
- Bloom filter FPR: 0.000957%

- Repeat filter: Yes - jellyfish, from Erdi: ```/extscratch/btl/bullfrog/genome/tga/repeat_filters/spruce_250.bf```, Rep filter = 1% rep kmers, kmer=25

## Run #5 - 17 Mar

Kollector with repeat filter and higher *r* and *s* values.   

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorHalf500_2```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/btl/kgagalova/Reads/WS77111-reads-1pet.in
/extscratch/btl/kgagalova/Reads/WS77111-reads-2pet.in
```
- Type of reads: PET
- Number of files: 38
- Number of reads: 3198709675 (3.1G)

- Seed bin size: 500
- Seed Seq statistics:
```
n	n:500   L50     min     N80     N50     N20     E-size  max     sum     name
500     491     182     500     1039    1437    1836    1468    2850    649937  /extscratch/btl/kgagalova/trial_kollector/seed_500/cdhit-output-4_1
```

- Kollector options: j=12, r=0.9, s=0.7, k=32, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b30G -H4 -j12 --kc=3
- Reads input in abyss 2.0.2: 39367320 (39M)
- Bloom filter FPR: 4.84e-08%

- Repeat filter: Yes - jellyfish, from Erdi: ```/extscratch/btl/bullfrog/genome/tga/repeat_filters/spruce_250.bf```, Rep filter = 1% rep kmers, kmer=25

- time: **kollector-recruit.mk** - 52988.63 (14h) , **abyss-pe** - 19741.69 (5.5h) , **kollector-extract.sh** - 437.41 (0.12h)

## Run #6 - 17 Mar

Kollector with repeat filter, larger seed

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorHalf1000```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/btl/kgagalova/Reads/WS77111-reads-1pet.in
/extscratch/btl/kgagalova/Reads/WS77111-reads-2pet.in
```
- Type of reads: PET
- Number of files: 38
- Number of reads: 3198709675 (3.1G)

- Seed bin size: 1000
- Seed Seq statistics:
```
n	n:500   L50     min     N80     N50     N20     E-size  max     sum     name
1000    986     365     500     1056    1452    1873    1492    3009    1321600 /extscratch/btl/kgagalova/trial_kollector/seed_1000/cdhit-output-4_1.fasta
```

- Kollector options: j=12, r=0.7, s=0.5, k=32, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b30G -H4 -j12 --kc=3
- Reads input in abyss 2.0.2: 407279252 (0.4G)
- Bloom filter FPR: 0.000633%

- Repeat filter: Yes - jellyfish, from Erdi: ```/extscratch/btl/bullfrog/genome/tga/repeat_filters/spruce_250.bf```, Rep filter = 1% rep kmers, kmer=25

- time: **kollector-recruit.mk** - 42681 (11.8h) , **abyss-pe** -  , **kollector-extract.sh** - 


## Run #7 - 17 Mar

Kollector with repeat filter and higher *r* and *s* values, larger bin

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorHalf1000_1```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/btl/kgagalova/Reads/WS77111-reads-1pet.in
/extscratch/btl/kgagalova/Reads/WS77111-reads-2pet.in
```
- Type of reads: PET
- Number of files: 38
- Number of reads: 3198709675 (3.1G)

- Seed bin size: 1000
- Seed Seq statistics:
```
n       n:500   L50     min     N80     N50     N20     E-size  max     sum     name
1000    986     365     500     1056    1452    1873    1492    3009    1321600 /extscratch/btl/kgagalova/trial_kollector/seed_1000/cdhit-output-4_1.fasta
```

- Kollector options: j=12, r=0.9, s=0.7, k=32, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k32 -q3 -v -b30G -H4 -j12 --kc=3
- Reads input in abyss 2.0.2:  76468502 (76M)
- Bloom filter FPR: 7.97e-07%

- Repeat filter: Yes - jellyfish, from Erdi: ```/extscratch/btl/bullfrog/genome/tga/repeat_filters/spruce_250.bf```, Rep filter = 1% rep kmers, kmer=25

- time: **kollector-recruit.mk** - 51827 (14.4h) , **abyss-pe** - 45767.95 (12.7h) , **kollector-extract.sh** - 919.76 (0.25h)


## Run #8

Use k value for Abyss equel to the WS77111 assembly and lower down the size of the bloom filter

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorHalf1000_2```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/btl/kgagalova/Reads/WS77111-reads-1pet.in
/extscratch/btl/kgagalova/Reads/WS77111-reads-2pet.in
```
- Type of reads: PET
- Number of files: 38
- Number of reads: 3198709675 (3.1G)

- Seed bin size: 1000
- Seed Seq statistics:
```
n       n:500   L50     min     N80     N50     N20     E-size  max     sum     name
1000    986     365     500     1056    1452    1873    1492    3009    1321600 /extscratch/btl/kgagalova/trial_kollector/seed_1000/cdhit-output-4_1.fasta
```

- Kollector options: j=12, r=0.7, s=0.5, k=116, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k116 -q3 -v -b7000M -H4 -j12 --kc=3
- Reads input in abyss 2.0.2: 407279252  
- Bloom filter FPR: 0.586%

- Repeat filter: Yes - jellyfish, from Erdi: ```/extscratch/btl/bullfrog/genome/tga/repeat_filters/spruce_250.bf```, Rep filter = 1% rep kmers, kmer=25

## Run #9

Same as Run #8 but with higher s and r Kollector values

- Directory: ```/extscratch/btl/kgagalova/trial_kollectorHalf1000_3```

- Memory settings: -pe ncpus 12, -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G

- Files location:
```
/extscratch/btl/kgagalova/Reads/WS77111-reads-1pet.in
/extscratch/btl/kgagalova/Reads/WS77111-reads-2pet.in
```
- Type of reads: PET
- Number of files: 38
- Number of reads: 3198709675 (3.1G)

- Seed bin size: 1000
- Seed Seq statistics:
```
n       n:500   L50     min     N80     N50     N20     E-size  max     sum     name
1000    986     365     500     1056    1452    1873    1492    3009    1321600 /extscratch/btl/kgagalova/trial_kollector/seed_1000/cdhit-output-4_1.fasta
```

- Kollector options: j=12, r=0.9, s=0.7, k=116, K=25, n=50000000

- Abyss abyss-bloom-dbg settings: -k116 -q3 -v -b7000M -H4 -j12 --kc=3
- Reads input in abyss 2.0.2: 76468502
- Bloom filter FPR: 0.0155%

- Repeat filter: Yes - jellyfish, from Erdi: ```/extscratch/btl/bullfrog/genome/tga/repeat_filters/spruce_250.bf```, Rep filter = 1% rep kmers, kmer=25

- time: **kollector-recruit.mk** - same as #7 , **abyss-pe** - 270309.42 (75h) , **kollector-extract.sh** - 758.38 (0.21h)

## Summary table

|       Run        |      Type reads      | # reads tot.  | # transcripts (seed) | # reads in Abyss |         Kollector param.         | abyss-bloom-dbg | FPR (%) |    Repeat filter    | Hall of shame |
|:----------------:|:--------------------:|:-------------:|:--------------------:|:----------------:|:--------------------------------:|:---------------:|:-------:|:-------------------:|:-------------:|
|  #1 crash test   | PET, Miseq, Chromium |      6.1G     |         1000         |       1700M       | r=0.7,s=0.5,k=32,K=25,n=50000000 |   b2G,H3,kc=3   |   97.6  |         NONE        |	qdel	 |
| #2  smaller seed | PET, Miseq, Chromium |      6.1G     |          500         |       1800M       | r=0.7,s=0.5,k=32,K=25,n=50000000 |   b30G,H4,kc=3  |  0.0255 |         NONE        |	qdel 	 |
|  #3 subset reads |          PET         |      3.1G     |          500         |       900M       | r=0.7,s=0.5,k=32,K=25,n=50000000 |   b30G,H4,kc=3  | 0.00587 |         NONE        |	abyss: terminate called after throwing an instance of 'std::bad_alloc'	 |
| #4 repeat filter |          PET         |      3.1G     |          500         |       480M      | r=0.7,s=0.5,k=32,K=25,n=50000000 |   b30G,H4,kc=3  | 0.000957|  K=25, 250X rep kmers |abyss: terminate called after throwing an instance of 'std::bad_alloc' |
| #5 higher r and s|          PET	  |	 3.1G     |	     500	 | 	 39M        | r=0.9,s=0.7,k=32,K=25,n=50000000 |   b30G,H4,kc=3  | 4.84e-08|  K=25, 250X rep kmers |               |
| #6 larger seed   |	      PET         |      3.1G	  |         1000	 | 	 400M 	    | r=0.7,s=0.5,k=32,K=25,n=50000000 |   b30G,H4,kc=3  | 0.000633|  K=25, 250X rep kmers |abyss: terminate called after throwing an instance of 'std::bad_alloc' |
| #7 higher r and s|	      PET      	  |    	 3.1G     |    	    1000         | 	 76M 	    | r=0.9,s=0.7,k=32,K=25,n=50000000 |   b30G,H4,kc=3	 | 7.97e-07|  K=25, 250X rep kmers |		 |
| #8 Abyss kmer	   | 	      PET 	  | 	 3.1G	  | 	    1000	 |       400M	    | r=0.7,s=0.5,k=116,K=25,n=50000000|   b7000M,H4,kc=3| 0.586   |  K=25, 250X rep kmers |	abyss: terminate called after throwing an instance of 'std::bad_alloc'	 |	
| #9 AByss kmer,r s| 	      PET      	  |	 3.1G     |	    1000         |  	 76M	    | r=0.9,s=0.7,k=116,K=25,n=50000000|   b2000M,H4,kc=3| 0.0155  |  K=25, 250X rep kmers |		 | 	 
