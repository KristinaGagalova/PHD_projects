GO-termsTest
================

## Test go terms mapped to GO slim - Based on Pfam domains

``` r
q903 <- read.delim("/projects/btl/kgagalova/PHD_projects2/SpruceUp/SpruceUpPaper/data/GoTerms/MapToPfam/Pfam/my_Slimgo_terms_plantQ903.counts", header=FALSE)
ws <- read.delim("/projects/btl/kgagalova/PHD_projects2/SpruceUp/SpruceUpPaper/data/GoTerms/MapToPfam/Pfam/my_Slimgo_terms_plantWS77111.counts", header=FALSE)
pg29 <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/Pfam/PG29/my_Slimgo_terms_plant.counts", header=FALSE)

#######################################################
########Q903 vs WS77111

##bind the 2 classes
q903_wsSlimAll = merge(q903,ws,by="V2")
q903_wsSlimAll = q903_wsSlimAll[,c(1,2,5,3,4)]
colnames(q903_wsSlimAll) = c("Go_term","Q903","WS77111","Domain","Go_name")
dim(q903_wsSlimAll)
```

    ## [1] 76  5

``` r
#############All domains together - Q903 vs WS77111
colSums_q903_wsSlimAll = unname(apply(q903_wsSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_q903_wsSlimAll)
p1 =colSums_q903_wsSlimAll / GrandTotal

p1
```

    ## [1] 0.4503037 0.5496963

``` r
#filter for low counts
q903_wsSlimAllSubs = subset(q903_wsSlimAll,Q903 > 5 & WS77111 > 5)
dim(q903_wsSlimAllSubs)
```

    ## [1] 72  5

``` r
prob_slim = c()
for (i in 1:nrow(q903_wsSlimAllSubs)){
  #print(i)
  obs = q903_wsSlimAllSubs[i,c("Q903","WS77111")]
  #print(obs)
  prob_slim[i] = chisq.test(obs, p=p1)$p.value
}

#sum(p.adjust(prob_slim, method = "BH", n = length(prob_slim)) < 0.05)
#sum(p.adjust(prob_slim, method = "BH", n = length(prob_slim)) < 0.001)
#q903_wsSlimAllSubs$prob_slim = prob_slim

q903_wsSlimAllSubs$prob_slimBH = p.adjust(prob_slim, method = "bonferroni", n = length(prob_slim))
subset(q903_wsSlimAllSubs,q903_wsSlimAllSubs$prob_slimBH < 0.05)
```

    ##       Go_term Q903 WS77111                                Domain
    ## 2  GO:0000166 2952    3125                    nucleotide_binding
    ## 16 GO:0005515 2593    4320                       protein_binding
    ## 24 GO:0005654   25      69                           nucleoplasm
    ## 40 GO:0006464 2177    2220 cellular_protein_modification_process
    ## 45 GO:0007165  347     264                   signal_transduction
    ## 54 GO:0009579  204     404                             thylakoid
    ## 61 GO:0015979  251     405                        photosynthesis
    ## 62 GO:0016020 1919    2734                              membrane
    ## 65 GO:0016301 1936    1980                       kinase_activity
    ##               Go_name  prob_slimBH
    ## 2  molecular_function 1.982208e-06
    ## 16 molecular_function 2.240629e-34
    ## 24 cellular_component 2.359107e-02
    ## 40 biological_process 1.689305e-07
    ## 45 biological_process 3.678405e-07
    ## 54 cellular_component 9.230336e-07
    ## 61 biological_process 3.553357e-02
    ## 62 cellular_component 1.483991e-05
    ## 65 molecular_function 2.127137e-06

``` r
subset(q903_wsSlimAllSubs,q903_wsSlimAllSubs$prob_slimBH < 0.001)
```

    ##       Go_term Q903 WS77111                                Domain
    ## 2  GO:0000166 2952    3125                    nucleotide_binding
    ## 16 GO:0005515 2593    4320                       protein_binding
    ## 40 GO:0006464 2177    2220 cellular_protein_modification_process
    ## 45 GO:0007165  347     264                   signal_transduction
    ## 54 GO:0009579  204     404                             thylakoid
    ## 62 GO:0016020 1919    2734                              membrane
    ## 65 GO:0016301 1936    1980                       kinase_activity
    ##               Go_name  prob_slimBH
    ## 2  molecular_function 1.982208e-06
    ## 16 molecular_function 2.240629e-34
    ## 40 biological_process 1.689305e-07
    ## 45 biological_process 3.678405e-07
    ## 54 cellular_component 9.230336e-07
    ## 62 cellular_component 1.483991e-05
    ## 65 molecular_function 2.127137e-06

``` r
#######################################################
########PG - Q903

##bind the 2 classes
q903_pgSlimAll = merge(q903,pg29,by="V2")
q903_pgSlimAll = q903_pgSlimAll[,c(1,2,5,3,4)]
colnames(q903_pgSlimAll) = c("Go_term","Q903","PG29","Domain","Go_name")
dim(q903_pgSlimAll)
```

    ## [1] 76  5

``` r
#############All domains together - Q903 vs PG29
colSums_q903_pgSlimAll = unname(apply(q903_pgSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_q903_pgSlimAll)
p1 =colSums_q903_pgSlimAll / GrandTotal

p1
```

    ## [1] 0.4731404 0.5268596

``` r
#filter for low counts
q903_pgSlimAllSubs = subset(q903_pgSlimAll,Q903 > 5 & PG29 > 5)
dim(q903_pgSlimAllSubs)
```

    ## [1] 70  5

``` r
prob_slim = c()
for (i in 1:nrow(q903_pgSlimAllSubs)){
  #print(i)
  obs = q903_pgSlimAllSubs[i,c("Q903","PG29")]
  #print(obs)
  prob_slim[i] = chisq.test(obs, p=p1)$p.value
}

sum(p.adjust(prob_slim, method = "bonferroni", n = length(prob_slim)) < 0.05)
```

    ## [1] 10

``` r
sum(p.adjust(prob_slim, method = "bonferroni", n = length(prob_slim)) < 0.001)
```

    ## [1] 8

``` r
q903_pgSlimAllSubs$prob_slim = prob_slim

q903_pgSlimAllSubs$prob_slimBH = p.adjust(prob_slim, method = "bonferroni", n = length(prob_slim))
subset(q903_pgSlimAllSubs,prob_slimBH < 0.05) 
```

    ##       Go_term Q903 PG29                                Domain
    ## 3  GO:0003676  496  695                  nucleic_acid_binding
    ## 14 GO:0005488 2850 2850                               binding
    ## 15 GO:0005515 2593 4122                       protein_binding
    ## 23 GO:0005654   25   74                           nucleoplasm
    ## 39 GO:0006464 2177 2059 cellular_protein_modification_process
    ## 44 GO:0007165  347  259                   signal_transduction
    ## 48 GO:0008152 2901 2832                     metabolic_process
    ## 53 GO:0009579  204  349                             thylakoid
    ## 62 GO:0016020 1919 2492                              membrane
    ## 65 GO:0016301 1936 1867                       kinase_activity
    ##               Go_name    prob_slim  prob_slimBH
    ## 3  molecular_function 8.926370e-05 6.248459e-03
    ## 14 molecular_function 4.874339e-05 3.412037e-03
    ## 15 molecular_function 3.025699e-46 2.117990e-44
    ## 23 cellular_component 1.099928e-05 7.699496e-04
    ## 39 biological_process 1.054973e-07 7.384811e-06
    ## 44 biological_process 9.378735e-07 6.565114e-05
    ## 48 biological_process 6.166667e-07 4.316667e-05
    ## 53 cellular_component 9.114166e-07 6.379916e-05
    ## 62 cellular_component 4.040035e-07 2.828024e-05
    ## 65 molecular_function 9.076803e-06 6.353762e-04

``` r
subset(q903_pgSlimAllSubs,prob_slimBH < 0.001)
```

    ##       Go_term Q903 PG29                                Domain
    ## 15 GO:0005515 2593 4122                       protein_binding
    ## 23 GO:0005654   25   74                           nucleoplasm
    ## 39 GO:0006464 2177 2059 cellular_protein_modification_process
    ## 44 GO:0007165  347  259                   signal_transduction
    ## 48 GO:0008152 2901 2832                     metabolic_process
    ## 53 GO:0009579  204  349                             thylakoid
    ## 62 GO:0016020 1919 2492                              membrane
    ## 65 GO:0016301 1936 1867                       kinase_activity
    ##               Go_name    prob_slim  prob_slimBH
    ## 15 molecular_function 3.025699e-46 2.117990e-44
    ## 23 cellular_component 1.099928e-05 7.699496e-04
    ## 39 biological_process 1.054973e-07 7.384811e-06
    ## 44 biological_process 9.378735e-07 6.565114e-05
    ## 48 biological_process 6.166667e-07 4.316667e-05
    ## 53 cellular_component 9.114166e-07 6.379916e-05
    ## 62 cellular_component 4.040035e-07 2.828024e-05
    ## 65 molecular_function 9.076803e-06 6.353762e-04

``` r
#######################################################
########WS - PG

##bind the 2 classes
ws_pgSlimAll = merge(ws,pg29,by="V2")
ws_pgSlimAll = ws_pgSlimAll[,c(1,2,5,3,4)]
colnames(ws_pgSlimAll) = c("Go_term","WS77111","PG29","Domain","Go_name")
dim(ws_pgSlimAll)
```

    ## [1] 75  5

``` r
#############All domains together - WS vs PG29
colSums_ws_pgSlimAll = unname(apply(ws_pgSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_ws_pgSlimAll)
p1 =colSums_ws_pgSlimAll / GrandTotal

p1
```

    ## [1] 0.5229655 0.4770345

``` r
#filter for low counts
ws_pgSlimAllSubs = subset(ws_pgSlimAll,WS77111 > 5 & PG29 > 5)
dim(ws_pgSlimAllSubs)
```

    ## [1] 70  5

``` r
prob_slim = c()
for (i in 1:nrow(ws_pgSlimAllSubs)){
  #print(ws_pgSlimAllSubs[i,c("WS77111","PG29")])
  obs = ws_pgSlimAllSubs[i,c("WS77111","PG29")]
  print(chisq.test(obs, p=p1))
  #print(obs)
  prob_slim[i] = chisq.test(obs, p=p1)$p.value
}
```

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 10.788, df = 1, p-value = 0.001022
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 3.6056, df = 1, p-value = 0.05759
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 18.226, df = 1, p-value = 1.962e-05
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.96467, df = 1, p-value = 0.326
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.7551, df = 1, p-value = 0.3849
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.90394, df = 1, p-value = 0.3417
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.061063, df = 1, p-value = 0.8048
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.17447, df = 1, p-value = 0.6762
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.849, df = 1, p-value = 0.3568
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.96015, df = 1, p-value = 0.3271
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 2.9887, df = 1, p-value = 0.08385
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.93517, df = 1, p-value = 0.3335
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.012189, df = 1, p-value = 0.9121
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 15.441, df = 1, p-value = 8.513e-05
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 4.274, df = 1, p-value = 0.0387
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 5.4398, df = 1, p-value = 0.01968
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.3633, df = 1, p-value = 0.243
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 4.6193, df = 1, p-value = 0.03161
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.0494, df = 1, p-value = 0.3057
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.16698, df = 1, p-value = 0.6828
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.18625, df = 1, p-value = 0.6661
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.9378, df = 1, p-value = 0.3328
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.55896, df = 1, p-value = 0.4547
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.0192, df = 1, p-value = 0.3127
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.0028317, df = 1, p-value = 0.9576
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.90915, df = 1, p-value = 0.3403
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.013325, df = 1, p-value = 0.9081
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.36866, df = 1, p-value = 0.5437
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.93916, df = 1, p-value = 0.3325
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.67201, df = 1, p-value = 0.4124
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.2892, df = 1, p-value = 0.2562
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.074292, df = 1, p-value = 0.7852
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.2669, df = 1, p-value = 0.6054
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.39693, df = 1, p-value = 0.5287
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 3.3369, df = 1, p-value = 0.06774
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.1492, df = 1, p-value = 0.2837
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.29579, df = 1, p-value = 0.5865
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.41336, df = 1, p-value = 0.5203
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.004055, df = 1, p-value = 0.9492
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 2.0887, df = 1, p-value = 0.1484
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.1655, df = 1, p-value = 0.2803
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.69331, df = 1, p-value = 0.405
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.7865, df = 1, p-value = 0.3752
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.38238, df = 1, p-value = 0.5363
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.0406, df = 1, p-value = 0.3077
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 3.5665, df = 1, p-value = 0.05896
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 4.6641, df = 1, p-value = 0.0308
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.335, df = 1, p-value = 0.2479
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.29136, df = 1, p-value = 0.5894
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.8705, df = 1, p-value = 0.3508
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.55459, df = 1, p-value = 0.4564
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.063424, df = 1, p-value = 0.8012
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.011615, df = 1, p-value = 0.9142
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.096142, df = 1, p-value = 0.7565
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 3.6056, df = 1, p-value = 0.05759
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.25627, df = 1, p-value = 0.6127
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.2693, df = 1, p-value = 0.2599
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.54816, df = 1, p-value = 0.4591
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.00073998, df = 1, p-value = 0.9783
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.11043, df = 1, p-value = 0.7397
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.89846, df = 1, p-value = 0.3432
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 1.0569, df = 1, p-value = 0.3039
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 7.3594, df = 1, p-value = 0.006671
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.088416, df = 1, p-value = 0.7662
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.21943, df = 1, p-value = 0.6395
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.060788, df = 1, p-value = 0.8053
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.084564, df = 1, p-value = 0.7712
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.16187, df = 1, p-value = 0.6874
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.26822, df = 1, p-value = 0.6045
    ## 
    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 0.10233, df = 1, p-value = 0.7491

``` r
sum(p.adjust(prob_slim, method = "bonferroni", n = length(prob_slim)) < 0.05)
```

    ## [1] 2

``` r
sum(p.adjust(prob_slim, method = "bonferroni", n = length(prob_slim)) < 0.001)
```

    ## [1] 0

``` r
ws_pgSlimAllSubs$prob_slim = prob_slim

ws_pgSlimAllSubs$prob_slimBH = p.adjust(prob_slim, method = "bonferroni", n = length(prob_slim))
subset(ws_pgSlimAllSubs,prob_slimBH < 0.05)
```

    ##       Go_term WS77111 PG29               Domain            Go_name
    ## 3  GO:0003676     601  695 nucleic_acid_binding molecular_function
    ## 14 GO:0005488    3451 2850              binding molecular_function
    ##       prob_slim prob_slimBH
    ## 3  1.962335e-05 0.001373634
    ## 14 8.512888e-05 0.005959022

``` r
subset(ws_pgSlimAllSubs,prob_slimBH < 0.001)
```

    ## [1] Go_term     WS77111     PG29        Domain      Go_name     prob_slim  
    ## [7] prob_slimBH
    ## <0 rows> (or 0-length row.names)

## Test go terms mapped to GO slim - Based on BLAST with id \> 50%, eval 1e-10 and aa align \> 100

``` r
q903_50id <- read.delim("/projects/btl/kgagalova/PHD_projects2/SpruceUp/SpruceUpPaper/data/GoTerms/MapToPfam/BLAST/my_Slimgo_terms_plant.Q903IdLen.counts", header=FALSE)
ws_50id <- read.delim("/projects/btl/kgagalova/PHD_projects2/SpruceUp/SpruceUpPaper/data/GoTerms/MapToPfam/BLAST/my_Slimgo_terms_plant.WS77111IdLen.counts", header=FALSE)

##bind the 2 classes
q903_50id_wsSlimAll = merge(q903_50id,ws_50id,by="V2")
q903_50id_wsSlimAll = q903_50id_wsSlimAll[,c(1,2,5,3,4)]
colnames(q903_50id_wsSlimAll) = c("Go_term","Q903","WS77111","Domain","Go_name")
dim(q903_50id_wsSlimAll)
```

    ## [1] 98  5

``` r
#############All domains together
colSums_q903_50id_wsSlimAll = unname(apply(q903_50id_wsSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_q903_50id_wsSlimAll)
p1 =colSums_q903_50id_wsSlimAll / GrandTotal

p1
```

    ## [1] 0.4877626 0.5122374

``` r
#filter for low counts
q903_50id_wsSlimAllSubs = subset(q903_50id_wsSlimAll,q903_50id_wsSlimAll$Q903 > 5 & q903_50id_wsSlimAll$WS77111 > 5)
dim(q903_50id_wsSlimAllSubs)
```

    ## [1] 97  5

``` r
prob_slim50id = c()
for (i in 1:nrow(q903_50id_wsSlimAllSubs)){
  #print(i)
  obs = q903_50id_wsSlimAllSubs[i,c("Q903","WS77111")]
  #print(obs)
  prob_slim50id[i] = chisq.test(obs, p=p1)$p.value
}

sum(p.adjust(prob_slim50id, method = "BH", n = length(prob_slim50id)) < 0.05)
```

    ## [1] 45

``` r
sum(p.adjust(prob_slim50id, method = "BH", n = length(prob_slim50id)) < 0.001)
```

    ## [1] 31

``` r
q903_50id_wsSlimAllSubs$prob_slim50id = prob_slim50id
q903_50id_wsSlimAllSubs$prob_slim50BH = p.adjust(prob_slim50id, method = "BH", n = length(prob_slim50id))
subset(q903_50id_wsSlimAllSubs,q903_50id_wsSlimAllSubs$prob_slim50BH < 0.05) 
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166  80177   78945
    ## 3  GO:0003674   2582    3097
    ## 5  GO:0003677  66667   68491
    ## 8  GO:0003723  11593   13906
    ## 10 GO:0003824  94325  105825
    ## 11 GO:0004518    347     473
    ## 12 GO:0005102    584     514
    ## 13 GO:0005198   7788    8690
    ## 14 GO:0005215  29289   28962
    ## 15 GO:0005488  83329   98353
    ## 17 GO:0005575   5656    6345
    ## 21 GO:0005622  16403   14347
    ## 23 GO:0005634  40896   42213
    ## 24 GO:0005635    297     454
    ## 31 GO:0005773   8301    9240
    ## 32 GO:0005777   1343    1649
    ## 33 GO:0005783   6048    5935
    ## 34 GO:0005794   6429    7063
    ## 36 GO:0005840  10055   12519
    ## 37 GO:0005856   9703    9806
    ## 38 GO:0005886  20970   22951
    ## 39 GO:0005975  20462   23871
    ## 40 GO:0006091  24336   27473
    ## 41 GO:0006139  68248   61355
    ## 42 GO:0006259   1915    1818
    ## 43 GO:0006412   7593    9096
    ## 44 GO:0006464  33249   35679
    ## 46 GO:0006810  70695   63677
    ## 50 GO:0007165  21838   22116
    ## 55 GO:0008152  56512   65635
    ## 58 GO:0009056  17161   19099
    ## 59 GO:0009058  88960   85970
    ## 60 GO:0009507  71065   75570
    ## 61 GO:0009536  42966   46343
    ## 64 GO:0009606    319     408
    ## 68 GO:0009719  37490   38043
    ## 71 GO:0009835    557     688
    ## 76 GO:0009987  83223   88705
    ## 78 GO:0015979  36948   42802
    ## 79 GO:0016020 148148  157269
    ## 83 GO:0016740  51094   50388
    ## 88 GO:0019825     62      40
    ## 89 GO:0030154  22413   23013
    ## 90 GO:0030234   1742    1976
    ## 94 GO:0038023   2583    2347
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 5                                       DNA_binding molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 11                                nuclease_activity molecular_function
    ## 12                       signaling_receptor_binding molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 17                               cellular_component cellular_component
    ## 21                                    intracellular cellular_component
    ## 23                                          nucleus cellular_component
    ## 24                                 nuclear_envelope cellular_component
    ## 31                                          vacuole cellular_component
    ## 32                                       peroxisome cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 34                                  Golgi_apparatus cellular_component
    ## 36                                         ribosome cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 38                                  plasma_membrane cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 42                            DNA_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 44            cellular_protein_modification_process biological_process
    ## 46                                        transport biological_process
    ## 50                              signal_transduction biological_process
    ## 55                                metabolic_process biological_process
    ## 58                                catabolic_process biological_process
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 64                                          tropism biological_process
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 71                                   fruit_ripening biological_process
    ## 76                                 cellular_process biological_process
    ## 78                                   photosynthesis biological_process
    ## 79                                         membrane cellular_component
    ## 83                             transferase_activity molecular_function
    ## 88                                   oxygen_binding molecular_function
    ## 89                             cell_differentiation biological_process
    ## 90                        enzyme_regulator_activity molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ##    prob_slim50id prob_slim50BH
    ## 2   8.026084e-38  8.650335e-37
    ## 3   6.005545e-07  2.912690e-06
    ## 5   5.398189e-05  2.181768e-04
    ## 8   3.697916e-26  2.989149e-25
    ## 10  2.652534e-49  3.675655e-48
    ## 11  2.152868e-04  6.960941e-04
    ## 12  3.451485e-03  9.565545e-03
    ## 13  1.018404e-04  3.799431e-04
    ## 14  3.755602e-13  2.276834e-12
    ## 15 5.078348e-136 1.641999e-134
    ## 17  3.070123e-04  9.606514e-04
    ## 21  9.075394e-58  1.467189e-56
    ## 23  1.284223e-02  3.038283e-02
    ## 24  4.196725e-07  2.142539e-06
    ## 31  1.183467e-04  4.251715e-04
    ## 32  2.074119e-05  8.747373e-05
    ## 33  2.051686e-04  6.862537e-04
    ## 34  8.893152e-03  2.156589e-02
    ## 36  4.224762e-37  4.098019e-36
    ## 37  7.320950e-03  1.868769e-02
    ## 38  1.528300e-05  6.738414e-05
    ## 39  2.432344e-28  2.144885e-27
    ## 40  2.146529e-16  1.388089e-15
    ## 41 4.165226e-172 2.020134e-170
    ## 42  2.043084e-03  6.193099e-03
    ## 43  2.348024e-17  1.626845e-16
    ## 44  4.641921e-03  1.216936e-02
    ## 46 4.812357e-174 4.667986e-172
    ## 50  1.410482e-04  4.886313e-04
    ## 55  5.467074e-69  1.325765e-67
    ## 58  3.416683e-08  1.949519e-07
    ## 59  9.722157e-68  1.886099e-66
    ## 60  1.670354e-02  3.682372e-02
    ## 61  6.687700e-05  2.594827e-04
    ## 64  8.248947e-03  2.051661e-02
    ## 68  2.408234e-06  1.112375e-05
    ## 71  4.372690e-03  1.178197e-02
    ## 76  2.114279e-03  6.214700e-03
    ## 78  1.880862e-43  2.280546e-42
    ## 79  2.889377e-03  8.243224e-03
    ## 83  1.297077e-23  9.678193e-23
    ## 88  1.525641e-02  3.523503e-02
    ## 89  1.630624e-02  3.678384e-02
    ## 90  1.897835e-02  4.090889e-02
    ## 94  3.751600e-07  2.021695e-06

``` r
subset(q903_50id_wsSlimAllSubs,q903_50id_wsSlimAllSubs$prob_slim50BH < 0.001)
```

    ##       Go_term  Q903 WS77111
    ## 2  GO:0000166 80177   78945
    ## 3  GO:0003674  2582    3097
    ## 5  GO:0003677 66667   68491
    ## 8  GO:0003723 11593   13906
    ## 10 GO:0003824 94325  105825
    ## 11 GO:0004518   347     473
    ## 13 GO:0005198  7788    8690
    ## 14 GO:0005215 29289   28962
    ## 15 GO:0005488 83329   98353
    ## 17 GO:0005575  5656    6345
    ## 21 GO:0005622 16403   14347
    ## 24 GO:0005635   297     454
    ## 31 GO:0005773  8301    9240
    ## 32 GO:0005777  1343    1649
    ## 33 GO:0005783  6048    5935
    ## 36 GO:0005840 10055   12519
    ## 38 GO:0005886 20970   22951
    ## 39 GO:0005975 20462   23871
    ## 40 GO:0006091 24336   27473
    ## 41 GO:0006139 68248   61355
    ## 43 GO:0006412  7593    9096
    ## 46 GO:0006810 70695   63677
    ## 50 GO:0007165 21838   22116
    ## 55 GO:0008152 56512   65635
    ## 58 GO:0009056 17161   19099
    ## 59 GO:0009058 88960   85970
    ## 61 GO:0009536 42966   46343
    ## 68 GO:0009719 37490   38043
    ## 78 GO:0015979 36948   42802
    ## 83 GO:0016740 51094   50388
    ## 94 GO:0038023  2583    2347
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 5                                       DNA_binding molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 11                                nuclease_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 17                               cellular_component cellular_component
    ## 21                                    intracellular cellular_component
    ## 24                                 nuclear_envelope cellular_component
    ## 31                                          vacuole cellular_component
    ## 32                                       peroxisome cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 38                                  plasma_membrane cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 46                                        transport biological_process
    ## 50                              signal_transduction biological_process
    ## 55                                metabolic_process biological_process
    ## 58                                catabolic_process biological_process
    ## 59                             biosynthetic_process biological_process
    ## 61                                          plastid cellular_component
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 78                                   photosynthesis biological_process
    ## 83                             transferase_activity molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ##    prob_slim50id prob_slim50BH
    ## 2   8.026084e-38  8.650335e-37
    ## 3   6.005545e-07  2.912690e-06
    ## 5   5.398189e-05  2.181768e-04
    ## 8   3.697916e-26  2.989149e-25
    ## 10  2.652534e-49  3.675655e-48
    ## 11  2.152868e-04  6.960941e-04
    ## 13  1.018404e-04  3.799431e-04
    ## 14  3.755602e-13  2.276834e-12
    ## 15 5.078348e-136 1.641999e-134
    ## 17  3.070123e-04  9.606514e-04
    ## 21  9.075394e-58  1.467189e-56
    ## 24  4.196725e-07  2.142539e-06
    ## 31  1.183467e-04  4.251715e-04
    ## 32  2.074119e-05  8.747373e-05
    ## 33  2.051686e-04  6.862537e-04
    ## 36  4.224762e-37  4.098019e-36
    ## 38  1.528300e-05  6.738414e-05
    ## 39  2.432344e-28  2.144885e-27
    ## 40  2.146529e-16  1.388089e-15
    ## 41 4.165226e-172 2.020134e-170
    ## 43  2.348024e-17  1.626845e-16
    ## 46 4.812357e-174 4.667986e-172
    ## 50  1.410482e-04  4.886313e-04
    ## 55  5.467074e-69  1.325765e-67
    ## 58  3.416683e-08  1.949519e-07
    ## 59  9.722157e-68  1.886099e-66
    ## 61  6.687700e-05  2.594827e-04
    ## 68  2.408234e-06  1.112375e-05
    ## 78  1.880862e-43  2.280546e-42
    ## 83  1.297077e-23  9.678193e-23
    ## 94  3.751600e-07  2.021695e-06

``` r
q903_50id_wsSlimAllSubs$prob_slimBon = p.adjust(prob_slim50id, method = "bonferroni", n = length(prob_slim50id))
subset(q903_50id_wsSlimAllSubs,q903_50id_wsSlimAllSubs$prob_slimBon < 0.05)
```

    ##       Go_term  Q903 WS77111
    ## 2  GO:0000166 80177   78945
    ## 3  GO:0003674  2582    3097
    ## 5  GO:0003677 66667   68491
    ## 8  GO:0003723 11593   13906
    ## 10 GO:0003824 94325  105825
    ## 11 GO:0004518   347     473
    ## 13 GO:0005198  7788    8690
    ## 14 GO:0005215 29289   28962
    ## 15 GO:0005488 83329   98353
    ## 17 GO:0005575  5656    6345
    ## 21 GO:0005622 16403   14347
    ## 24 GO:0005635   297     454
    ## 31 GO:0005773  8301    9240
    ## 32 GO:0005777  1343    1649
    ## 33 GO:0005783  6048    5935
    ## 36 GO:0005840 10055   12519
    ## 38 GO:0005886 20970   22951
    ## 39 GO:0005975 20462   23871
    ## 40 GO:0006091 24336   27473
    ## 41 GO:0006139 68248   61355
    ## 43 GO:0006412  7593    9096
    ## 46 GO:0006810 70695   63677
    ## 50 GO:0007165 21838   22116
    ## 55 GO:0008152 56512   65635
    ## 58 GO:0009056 17161   19099
    ## 59 GO:0009058 88960   85970
    ## 61 GO:0009536 42966   46343
    ## 68 GO:0009719 37490   38043
    ## 78 GO:0015979 36948   42802
    ## 83 GO:0016740 51094   50388
    ## 94 GO:0038023  2583    2347
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 5                                       DNA_binding molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 11                                nuclease_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 17                               cellular_component cellular_component
    ## 21                                    intracellular cellular_component
    ## 24                                 nuclear_envelope cellular_component
    ## 31                                          vacuole cellular_component
    ## 32                                       peroxisome cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 38                                  plasma_membrane cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 46                                        transport biological_process
    ## 50                              signal_transduction biological_process
    ## 55                                metabolic_process biological_process
    ## 58                                catabolic_process biological_process
    ## 59                             biosynthetic_process biological_process
    ## 61                                          plastid cellular_component
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 78                                   photosynthesis biological_process
    ## 83                             transferase_activity molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ##    prob_slim50id prob_slim50BH  prob_slimBon
    ## 2   8.026084e-38  8.650335e-37  7.785301e-36
    ## 3   6.005545e-07  2.912690e-06  5.825379e-05
    ## 5   5.398189e-05  2.181768e-04  5.236244e-03
    ## 8   3.697916e-26  2.989149e-25  3.586978e-24
    ## 10  2.652534e-49  3.675655e-48  2.572958e-47
    ## 11  2.152868e-04  6.960941e-04  2.088282e-02
    ## 13  1.018404e-04  3.799431e-04  9.878521e-03
    ## 14  3.755602e-13  2.276834e-12  3.642934e-11
    ## 15 5.078348e-136 1.641999e-134 4.925998e-134
    ## 17  3.070123e-04  9.606514e-04  2.978019e-02
    ## 21  9.075394e-58  1.467189e-56  8.803132e-56
    ## 24  4.196725e-07  2.142539e-06  4.070824e-05
    ## 31  1.183467e-04  4.251715e-04  1.147963e-02
    ## 32  2.074119e-05  8.747373e-05  2.011896e-03
    ## 33  2.051686e-04  6.862537e-04  1.990136e-02
    ## 36  4.224762e-37  4.098019e-36  4.098019e-35
    ## 38  1.528300e-05  6.738414e-05  1.482451e-03
    ## 39  2.432344e-28  2.144885e-27  2.359374e-26
    ## 40  2.146529e-16  1.388089e-15  2.082134e-14
    ## 41 4.165226e-172 2.020134e-170 4.040269e-170
    ## 43  2.348024e-17  1.626845e-16  2.277583e-15
    ## 46 4.812357e-174 4.667986e-172 4.667986e-172
    ## 50  1.410482e-04  4.886313e-04  1.368168e-02
    ## 55  5.467074e-69  1.325765e-67  5.303061e-67
    ## 58  3.416683e-08  1.949519e-07  3.314183e-06
    ## 59  9.722157e-68  1.886099e-66  9.430493e-66
    ## 61  6.687700e-05  2.594827e-04  6.487069e-03
    ## 68  2.408234e-06  1.112375e-05  2.335987e-04
    ## 78  1.880862e-43  2.280546e-42  1.824436e-41
    ## 83  1.297077e-23  9.678193e-23  1.258165e-21
    ## 94  3.751600e-07  2.021695e-06  3.639052e-05

``` r
subset(q903_50id_wsSlimAllSubs,q903_50id_wsSlimAllSubs$prob_slimBon < 0.001)
```

    ##       Go_term  Q903 WS77111
    ## 2  GO:0000166 80177   78945
    ## 3  GO:0003674  2582    3097
    ## 8  GO:0003723 11593   13906
    ## 10 GO:0003824 94325  105825
    ## 14 GO:0005215 29289   28962
    ## 15 GO:0005488 83329   98353
    ## 21 GO:0005622 16403   14347
    ## 24 GO:0005635   297     454
    ## 36 GO:0005840 10055   12519
    ## 39 GO:0005975 20462   23871
    ## 40 GO:0006091 24336   27473
    ## 41 GO:0006139 68248   61355
    ## 43 GO:0006412  7593    9096
    ## 46 GO:0006810 70695   63677
    ## 55 GO:0008152 56512   65635
    ## 58 GO:0009056 17161   19099
    ## 59 GO:0009058 88960   85970
    ## 68 GO:0009719 37490   38043
    ## 78 GO:0015979 36948   42802
    ## 83 GO:0016740 51094   50388
    ## 94 GO:0038023  2583    2347
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 21                                    intracellular cellular_component
    ## 24                                 nuclear_envelope cellular_component
    ## 36                                         ribosome cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 46                                        transport biological_process
    ## 55                                metabolic_process biological_process
    ## 58                                catabolic_process biological_process
    ## 59                             biosynthetic_process biological_process
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 78                                   photosynthesis biological_process
    ## 83                             transferase_activity molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ##    prob_slim50id prob_slim50BH  prob_slimBon
    ## 2   8.026084e-38  8.650335e-37  7.785301e-36
    ## 3   6.005545e-07  2.912690e-06  5.825379e-05
    ## 8   3.697916e-26  2.989149e-25  3.586978e-24
    ## 10  2.652534e-49  3.675655e-48  2.572958e-47
    ## 14  3.755602e-13  2.276834e-12  3.642934e-11
    ## 15 5.078348e-136 1.641999e-134 4.925998e-134
    ## 21  9.075394e-58  1.467189e-56  8.803132e-56
    ## 24  4.196725e-07  2.142539e-06  4.070824e-05
    ## 36  4.224762e-37  4.098019e-36  4.098019e-35
    ## 39  2.432344e-28  2.144885e-27  2.359374e-26
    ## 40  2.146529e-16  1.388089e-15  2.082134e-14
    ## 41 4.165226e-172 2.020134e-170 4.040269e-170
    ## 43  2.348024e-17  1.626845e-16  2.277583e-15
    ## 46 4.812357e-174 4.667986e-172 4.667986e-172
    ## 55  5.467074e-69  1.325765e-67  5.303061e-67
    ## 58  3.416683e-08  1.949519e-07  3.314183e-06
    ## 59  9.722157e-68  1.886099e-66  9.430493e-66
    ## 68  2.408234e-06  1.112375e-05  2.335987e-04
    ## 78  1.880862e-43  2.280546e-42  1.824436e-41
    ## 83  1.297077e-23  9.678193e-23  1.258165e-21
    ## 94  3.751600e-07  2.021695e-06  3.639052e-05
