# Test different repeat filters applied to Kollector

## Different threshold repeat filters

**Total reads**: 3198709675 (3.1G)

|   Run   | Rep filter coverage | Threshold reached at reads # | Reads in Abyss | Size (G) |
|:-------:|:-------------------:|-----------------------:|--------------:|:--------:|
| Rep 250 |         250         |        1137685572       |    473224372   |    173   |
| Rep 100 |         100         |        1547829737       |    484741554   |    177   |
|  Rep 50 |          50         |        2087547853       |    275194930   |    103   |


## Test runs - repeat filters

**Total reads**: 3198709675 (3.1G)

parameters kollector: r=0.7, s=0.5, k=116, n=7000000

| kmer size | Reads read  (bbmaker) | Hit cap |Reads in Abyss | Success rate (%) | time BBT | time Abyss|
|:---------:|:---------------------:|:----:|--------------:|:------------:|:--------:|:---------:|
|     25    |       0      | yes |   1433124     |   20.9           | 6.7h    |  13.5h   |
|     30    |       0      | yes |   1241282    |    20.7          | 7h       |  15.1h         |
|     35    |       0      | yes |   588218     |   20           | 7.3h    |  0.049h    |
|     40    |       0      | yes |   546706     |   18.2           | 7h         | 0.04h          |

parameters kollector: r=0.9, s=0.7, k=116, n=7000000

| kmer size | Reads read  (bbmaker) | Hit cap | Reads in Abyss | Success rate (%) | time BBT | time Abyss|
|:---------:|:---------------------:|:----:|--------------:|:------------:|:--------:|:---------:|
|     25    |       0      |  yes |  577550  |    17          | 6.7h    |  3.2h   |
|     30    |       0      |  yes |  532874  |   16.7         | 6.8h       | 0.94h  |
|     35    |       0      |  yes |  332248  |   14.3          | 6.8h    | 0.03h    |
|     40    |       0      |  yes | 320860   |   13.8         | 6.75h         |   0.02h     |


_________________

parameters kollector: r=0.7, s=0.5, k=116, n=25000000

| kmer size | Reads read  (bbmaker) | Hit cap | Reads in Abyss | Success rate (%) | time BBT | time Abyss |
|:---------:|---------------------:|:----:|--------------:|:------------:|:--------:|:----------:|
|     25    |       978489452       | yes|   225192806   |       -       | 10.5h    |	AdjList - std::bad_alloc   |
|     30    |       232415256       | yes|   360270678   |       -       | 8.16h    |	AdjList - std::bad_alloc   |
|     35    |       1357510607      | yes|   165205072   |        -      | 11.1h    |   abyss-scaffold (RepFil35_250-6.path) - std::bad_alloc         |
|     40    |       1351266439      | yes|   132713566   |      35.4        | 11.1h    |      151h      |


parameters kollector: r=0.9, s=0.7, k=116, n=25000000

| kmer size | Reads read  (bbmaker) | Hit cap |Reads in Abyss | Success rate (%) | time BBT | time Abyss|
|:---------:|---------------------:|:---:|--------------:|:------------:|:--------:|:---------:|
|     25    |       3086004679      | no |   41811908     |    33.2      | 13.7h    |  38.82h   | 	
|     30    |       780822074       | yes|   510510780    |    -          | 9h       |	AdjList - std::bad_alloc   |
|     35    |       3086004679      | no |  48426602     |    33.1      | 13.8h    |  32.4h	   |	
|     40    |       3086004679      | no |  30229938     |    30.7      | 20.1h    |  13.h     |

______________________

parameters kollector: r=0.9, s=0.7, k=116, n=50000000

| kmer size | Reads read  (bbmaker) | Hit cap | Reads in Abyss | Success rate (%) | time BBT | time Abyss|
|:---------:|---------------------:|:---------:|--------------:|:------------:|:--------:|:---------:|
|     25    |       3086004679    | no |    41811904     |    33.2      | 13.6h    |  38.5h    |
|     30    |       812373807     | yes | 680357468    |      -        | 9.32h    |   AdjList - std::bad_alloc        |
|     35    |       3086004679    | no | 48425234     |    33.1      | 13.7h    |  32.3h    |
|     40    |       3086004679    | no |  30229838     |    30.7      | 13.9h    |  32.1h    |




