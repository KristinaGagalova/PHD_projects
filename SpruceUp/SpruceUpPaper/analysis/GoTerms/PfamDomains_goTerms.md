GO-termsTest
================

Test go terms mapped to GO slim - Based on Pfam domains
-------------------------------------------------------

``` r
q903 <- read.delim("/projects/btl/kgagalova/PHD_projects2/SpruceUp/SpruceUpPaper/data/GoTerms/MapToPfam/my_Slimgo_terms_plantQ903.counts", header=FALSE)
ws <- read.delim("/projects/btl/kgagalova/PHD_projects2/SpruceUp/SpruceUpPaper/data/GoTerms/MapToPfam/my_Slimgo_terms_plantWS77111.counts", header=FALSE)

##bind the 2 classes
q903_wsSlimAll = merge(q903,ws,by="V2")
q903_wsSlimAll = q903_wsSlimAll[,c(1,2,5,3,4)]
colnames(q903_wsSlimAll) = c("Go_term","Q903","WS77111","Domain","Go_name")
dim(q903_wsSlimAll)
```

    ## [1] 76  5

``` r
#############All domains together
colSums_q903_wsSlimAll = unname(apply(q903_wsSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_q903_wsSlimAll)
p1 =colSums_q903_wsSlimAll / GrandTotal

p1
```

    ## [1] 0.4503037 0.5496963

``` r
#filter for low counts
q903_wsSlimAllSubs = subset(q903_wsSlimAll,q903_wsSlimAll$Q903 > 5 & q903_wsSlimAll$WS77111 > 5)
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

sum(p.adjust(prob_slim, method = "BH", n = length(prob_slim)) < 0.05)
```

    ## [1] 11

``` r
sum(p.adjust(prob_slim, method = "BH", n = length(prob_slim)) < 0.001)
```

    ## [1] 7

``` r
q903_wsSlimAllSubs$prob_slim = prob_slim

q903_wsSlimAllSubs$prob_slimBH = p.adjust(prob_slim, method = "BH", n = length(prob_slim))

subset(q903_wsSlimAllSubs,q903_wsSlimAllSubs$prob_slimBH < 0.05) 
```

    ##       Go_term Q903 WS77111                                Domain
    ## 2  GO:0000166 2952    3125                    nucleotide_binding
    ## 16 GO:0005515 2593    4320                       protein_binding
    ## 24 GO:0005654   25      69                           nucleoplasm
    ## 40 GO:0006464 2177    2220 cellular_protein_modification_process
    ## 45 GO:0007165  347     264                   signal_transduction
    ## 49 GO:0008152 2901    3259                     metabolic_process
    ## 54 GO:0009579  204     404                             thylakoid
    ## 61 GO:0015979  251     405                        photosynthesis
    ## 62 GO:0016020 1919    2734                              membrane
    ## 65 GO:0016301 1936    1980                       kinase_activity
    ## 72 GO:0030246  166     141                  carbohydrate_binding
    ##               Go_name    prob_slim  prob_slimBH
    ## 2  molecular_function 2.753067e-08 3.545228e-07
    ## 16 molecular_function 3.111985e-36 2.240629e-34
    ## 24 cellular_component 3.276537e-04 2.948884e-03
    ## 40 biological_process 2.346257e-09 8.446526e-08
    ## 45 biological_process 5.108895e-09 1.226135e-07
    ## 49 biological_process 1.131212e-03 8.144728e-03
    ## 54 cellular_component 1.281991e-08 2.307584e-07
    ## 61 biological_process 4.935219e-04 3.948175e-03
    ## 62 cellular_component 2.061099e-07 2.119988e-06
    ## 65 molecular_function 2.954357e-08 3.545228e-07
    ## 72 molecular_function 1.452086e-03 9.504561e-03

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
    ##               Go_name    prob_slim  prob_slimBH
    ## 2  molecular_function 2.753067e-08 3.545228e-07
    ## 16 molecular_function 3.111985e-36 2.240629e-34
    ## 40 biological_process 2.346257e-09 8.446526e-08
    ## 45 biological_process 5.108895e-09 1.226135e-07
    ## 54 cellular_component 1.281991e-08 2.307584e-07
    ## 62 cellular_component 2.061099e-07 2.119988e-06
    ## 65 molecular_function 2.954357e-08 3.545228e-07

Test go terms mapped to GO slim - Based on BLAST and coverage 50% from subject and eval 1e-5
--------------------------------------------------------------------------------------------

``` r
q903_50 <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/Q903/ThresholdCovIdentity2/my_Slimgo_terms_plant.Q903cov50.counts", header=FALSE)
ws_50 <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/WS77111/ThresholdCovIdentity2/my_Slimgo_terms_plant.WS77111cov50.counts", header=FALSE)

##bind the 2 classes
q903_50_wsSlimAll = merge(q903_50,ws_50,by="V2")
q903_50_wsSlimAll = q903_50_wsSlimAll[,c(1,2,5,3,4)]
colnames(q903_50_wsSlimAll) = c("Go_term","Q903","WS77111","Domain","Go_name")
dim(q903_50_wsSlimAll)
```

    ## [1] 98  5

``` r
#############All domains together
colSums_q903_50_wsSlimAll = unname(apply(q903_50_wsSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_q903_50_wsSlimAll)
p1 =colSums_q903_50_wsSlimAll / GrandTotal

p1
```

    ## [1] 0.5241383 0.4758617

``` r
#filter for low counts
q903_50_wsSlimAllSubs = subset(q903_50_wsSlimAll,q903_50_wsSlimAll$Q903 > 5 & q903_50_wsSlimAll$WS77111 > 5)
dim(q903_50_wsSlimAllSubs)
```

    ## [1] 98  5

``` r
prob_slim50 = c()
for (i in 1:nrow(q903_50_wsSlimAllSubs)){
  #print(i)
  obs = q903_50_wsSlimAllSubs[i,c("Q903","WS77111")]
  #print(obs)
  prob_slim50[i] = chisq.test(obs, p=p1)$p.value
}

sum(p.adjust(prob_slim50, method = "BH", n = length(prob_slim50)) < 0.05)
```

    ## [1] 75

``` r
sum(p.adjust(prob_slim50, method = "BH", n = length(prob_slim50)) < 0.001)
```

    ## [1] 64

``` r
q903_50_wsSlimAllSubs$prob_slim50 = prob_slim50

q903_50_wsSlimAllSubs$prob_slim50BH = p.adjust(prob_slim50, method = "BH", n = length(prob_slim50))

subset(q903_50_wsSlimAllSubs,q903_50_wsSlimAllSubs$prob_slim50BH < 0.05) 
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166 497561  419196
    ## 3  GO:0003674  33330   31271
    ## 4  GO:0003676   6335    6254
    ## 5  GO:0003677  37207   36884
    ## 7  GO:0003700  13200   12897
    ## 8  GO:0003723  46711   47375
    ## 9  GO:0003774    496     372
    ## 10 GO:0003824 538560  516847
    ## 11 GO:0004518   2335    2325
    ## 12 GO:0005102  12488    9696
    ## 13 GO:0005198  20197   20604
    ## 14 GO:0005215  86105   82538
    ## 15 GO:0005488 482975  463095
    ## 16 GO:0005515 173113  155755
    ## 17 GO:0005575  35286   31247
    ## 18 GO:0005576  98916   91707
    ## 19 GO:0005615   3679    3125
    ## 20 GO:0005618  51135   47806
    ## 21 GO:0005622  57607   51155
    ## 22 GO:0005623   8575    7294
    ## 23 GO:0005634  83159   81913
    ## 25 GO:0005654   3327    3264
    ## 26 GO:0005730   3637    3569
    ## 27 GO:0005737  93411   90460
    ## 29 GO:0005764   1238     941
    ## 30 GO:0005768  16321   15514
    ## 31 GO:0005773  35745   33807
    ## 32 GO:0005777   6538    6239
    ## 33 GO:0005783  54017   51133
    ## 34 GO:0005794  23256   22065
    ## 35 GO:0005829  29841   28379
    ## 36 GO:0005840  28211   29793
    ## 38 GO:0005886 208432  178771
    ## 39 GO:0005975  53502   52508
    ## 40 GO:0006091  28831   28821
    ## 41 GO:0006139 206689  181357
    ## 43 GO:0006412  28066   28511
    ## 44 GO:0006464 282505  249784
    ## 45 GO:0006629 100279   94358
    ## 46 GO:0006810 150901  143646
    ## 47 GO:0006950 321704  276566
    ## 48 GO:0007049  28832   28226
    ## 51 GO:0007267    396     435
    ## 52 GO:0007275  68302   60575
    ## 53 GO:0008135   2327    2300
    ## 55 GO:0008152 393454  363668
    ## 56 GO:0008219  10037    7658
    ## 59 GO:0009058 178203  166054
    ## 60 GO:0009507 150086  144818
    ## 61 GO:0009536  93814   92219
    ## 62 GO:0009579  78643   75720
    ## 63 GO:0009605 102581   85325
    ## 64 GO:0009606   1574    1620
    ## 65 GO:0009607 102051   85252
    ## 66 GO:0009628  80933   76737
    ## 67 GO:0009653  36351   31954
    ## 68 GO:0009719  86313   80132
    ## 69 GO:0009790  15476   13181
    ## 71 GO:0009835   4857    4606
    ## 72 GO:0009838   3415    2820
    ## 75 GO:0009908  22575   19378
    ## 76 GO:0009987 463805  417436
    ## 78 GO:0015979  40959   43255
    ## 80 GO:0016043  87456   82604
    ## 82 GO:0016301 650502  556870
    ## 83 GO:0016740 552764  471897
    ## 84 GO:0016787 240853  225597
    ## 85 GO:0019538  24081   22963
    ## 86 GO:0019725   3313    3545
    ## 87 GO:0019748  38030   36457
    ## 90 GO:0030234  11813   11531
    ## 91 GO:0030246  27493   23692
    ## 94 GO:0038023  68075   55439
    ## 95 GO:0040007  15459   13045
    ## 98 GO:0140110    556     630
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 4                              nucleic_acid_binding molecular_function
    ## 5                                       DNA_binding molecular_function
    ## 7         DNA-binding_transcription_factor_activity molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 9                                    motor_activity molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 11                                nuclease_activity molecular_function
    ## 12                       signaling_receptor_binding molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 16                                  protein_binding molecular_function
    ## 17                               cellular_component cellular_component
    ## 18                             extracellular_region cellular_component
    ## 19                              extracellular_space cellular_component
    ## 20                                        cell_wall cellular_component
    ## 21                                    intracellular cellular_component
    ## 22                                             cell cellular_component
    ## 23                                          nucleus cellular_component
    ## 25                                      nucleoplasm cellular_component
    ## 26                                        nucleolus cellular_component
    ## 27                                        cytoplasm cellular_component
    ## 29                                         lysosome cellular_component
    ## 30                                         endosome cellular_component
    ## 31                                          vacuole cellular_component
    ## 32                                       peroxisome cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 34                                  Golgi_apparatus cellular_component
    ## 35                                          cytosol cellular_component
    ## 36                                         ribosome cellular_component
    ## 38                                  plasma_membrane cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 44            cellular_protein_modification_process biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 47                               response_to_stress biological_process
    ## 48                                       cell_cycle biological_process
    ## 51                              cell-cell_signaling biological_process
    ## 52               multicellular_organism_development biological_process
    ## 53         translation_factor_activity,_RNA_binding molecular_function
    ## 55                                metabolic_process biological_process
    ## 56                                       cell_death biological_process
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 63                    response_to_external_stimulus biological_process
    ## 64                                          tropism biological_process
    ## 65                      response_to_biotic_stimulus biological_process
    ## 66                     response_to_abiotic_stimulus biological_process
    ## 67               anatomical_structure_morphogenesis biological_process
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 69                               embryo_development biological_process
    ## 71                                   fruit_ripening biological_process
    ## 72                                       abscission biological_process
    ## 75                               flower_development biological_process
    ## 76                                 cellular_process biological_process
    ## 78                                   photosynthesis biological_process
    ## 80                  cellular_component_organization biological_process
    ## 82                                  kinase_activity molecular_function
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ## 85                        protein_metabolic_process biological_process
    ## 86                             cellular_homeostasis biological_process
    ## 87                      secondary_metabolic_process biological_process
    ## 90                        enzyme_regulator_activity molecular_function
    ## 91                             carbohydrate_binding molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ## 95                                           growth biological_process
    ## 98                 transcription_regulator_activity molecular_function
    ##      prob_slim50 prob_slim50BH
    ## 2  1.451071e-278 1.422049e-276
    ## 3   2.990124e-05  5.232718e-05
    ## 4   2.598797e-06  5.197594e-06
    ## 5   5.225994e-33  2.438797e-32
    ## 7   3.026615e-09  7.234348e-09
    ## 8   9.312424e-65  6.518697e-64
    ## 9   5.274523e-03  7.280328e-03
    ## 10 1.398657e-178 3.426710e-177
    ## 11  1.617431e-03  2.365794e-03
    ## 12  5.953258e-31  2.536606e-30
    ## 13  4.939786e-32  2.200450e-31
    ## 14  6.974583e-29  2.531515e-28
    ## 15 2.644299e-155 5.182826e-154
    ## 16  9.703375e-03  1.285042e-02
    ## 17  1.327426e-03  1.971027e-03
    ## 18  4.841928e-06  9.125173e-06
    ## 19  6.194619e-03  8.316064e-03
    ## 20  4.079489e-06  7.839019e-06
    ## 21  2.653108e-04  4.127058e-04
    ## 22  4.273401e-05  7.098191e-05
    ## 23  1.210118e-61  7.906107e-61
    ## 25  1.649590e-03  2.377351e-03
    ## 26  9.637489e-04  1.453037e-03
    ## 27  1.562362e-43  9.006556e-43
    ## 29  3.892658e-05  6.692640e-05
    ## 30  4.212441e-05  7.098191e-05
    ## 31  7.061020e-08  1.472298e-07
    ## 32  4.876965e-03  6.827750e-03
    ## 33  1.300526e-11  3.267988e-11
    ## 34  2.753257e-06  5.396384e-06
    ## 35  2.194308e-08  4.887322e-08
    ## 36  3.794354e-74  3.718467e-73
    ## 38  1.072246e-69  8.756679e-69
    ## 39  7.601677e-37  3.724822e-36
    ## 40  6.313427e-31  2.577983e-30
    ## 41  2.825675e-26  9.548831e-26
    ## 43  9.119932e-41  4.965296e-40
    ## 44  5.494027e-22  1.794716e-21
    ## 45  3.101287e-15  8.939004e-15
    ## 46  8.837211e-38  4.558141e-37
    ## 47  2.779994e-98  3.891991e-97
    ## 48  2.151502e-19  6.801524e-19
    ## 51  6.000210e-03  8.166952e-03
    ## 52  2.694025e-05  4.800263e-05
    ## 53  3.848633e-03  5.466175e-03
    ## 55  7.023507e-15  1.966582e-14
    ## 56  1.747617e-30  6.587170e-30
    ## 59  2.379677e-14  6.478009e-14
    ## 60  2.049377e-61  1.255243e-60
    ## 61  6.919265e-66  5.216062e-65
    ## 62  8.185165e-31  3.208585e-30
    ## 63  1.077113e-79  1.172856e-78
    ## 64  3.904692e-04  5.979060e-04
    ## 65  5.386996e-72  4.799323e-71
    ## 66  7.160336e-18  2.192853e-17
    ## 67  2.533805e-05  4.598386e-05
    ## 68  5.348479e-06  9.889640e-06
    ## 69  7.007913e-08  1.472298e-07
    ## 71  3.413468e-02  4.460264e-02
    ## 72  1.933099e-04  3.055544e-04
    ## 75  1.022413e-08  2.330151e-08
    ## 76  4.500916e-05  7.351495e-05
    ## 78 9.220117e-107 1.505952e-105
    ## 80  3.573927e-16  1.061348e-15
    ## 82 1.569032e-227 7.688258e-226
    ## 83 9.476026e-212 3.095502e-210
    ## 84  1.815824e-26  6.355386e-26
    ## 85  1.022553e-07  2.087712e-07
    ## 86  9.940943e-12  2.563717e-11
    ## 87  1.163243e-13  3.081023e-13
    ## 90  3.080445e-08  6.708526e-08
    ## 91  3.971000e-09  9.265667e-09
    ## 94  1.410929e-80  1.728388e-79
    ## 95  7.513134e-10  1.840718e-09
    ## 98  1.357531e-04  2.180951e-04

``` r
subset(q903_50_wsSlimAllSubs,q903_50_wsSlimAllSubs$prob_slim50BH < 0.001)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166 497561  419196
    ## 3  GO:0003674  33330   31271
    ## 4  GO:0003676   6335    6254
    ## 5  GO:0003677  37207   36884
    ## 7  GO:0003700  13200   12897
    ## 8  GO:0003723  46711   47375
    ## 10 GO:0003824 538560  516847
    ## 12 GO:0005102  12488    9696
    ## 13 GO:0005198  20197   20604
    ## 14 GO:0005215  86105   82538
    ## 15 GO:0005488 482975  463095
    ## 18 GO:0005576  98916   91707
    ## 20 GO:0005618  51135   47806
    ## 21 GO:0005622  57607   51155
    ## 22 GO:0005623   8575    7294
    ## 23 GO:0005634  83159   81913
    ## 27 GO:0005737  93411   90460
    ## 29 GO:0005764   1238     941
    ## 30 GO:0005768  16321   15514
    ## 31 GO:0005773  35745   33807
    ## 33 GO:0005783  54017   51133
    ## 34 GO:0005794  23256   22065
    ## 35 GO:0005829  29841   28379
    ## 36 GO:0005840  28211   29793
    ## 38 GO:0005886 208432  178771
    ## 39 GO:0005975  53502   52508
    ## 40 GO:0006091  28831   28821
    ## 41 GO:0006139 206689  181357
    ## 43 GO:0006412  28066   28511
    ## 44 GO:0006464 282505  249784
    ## 45 GO:0006629 100279   94358
    ## 46 GO:0006810 150901  143646
    ## 47 GO:0006950 321704  276566
    ## 48 GO:0007049  28832   28226
    ## 52 GO:0007275  68302   60575
    ## 55 GO:0008152 393454  363668
    ## 56 GO:0008219  10037    7658
    ## 59 GO:0009058 178203  166054
    ## 60 GO:0009507 150086  144818
    ## 61 GO:0009536  93814   92219
    ## 62 GO:0009579  78643   75720
    ## 63 GO:0009605 102581   85325
    ## 64 GO:0009606   1574    1620
    ## 65 GO:0009607 102051   85252
    ## 66 GO:0009628  80933   76737
    ## 67 GO:0009653  36351   31954
    ## 68 GO:0009719  86313   80132
    ## 69 GO:0009790  15476   13181
    ## 72 GO:0009838   3415    2820
    ## 75 GO:0009908  22575   19378
    ## 76 GO:0009987 463805  417436
    ## 78 GO:0015979  40959   43255
    ## 80 GO:0016043  87456   82604
    ## 82 GO:0016301 650502  556870
    ## 83 GO:0016740 552764  471897
    ## 84 GO:0016787 240853  225597
    ## 85 GO:0019538  24081   22963
    ## 86 GO:0019725   3313    3545
    ## 87 GO:0019748  38030   36457
    ## 90 GO:0030234  11813   11531
    ## 91 GO:0030246  27493   23692
    ## 94 GO:0038023  68075   55439
    ## 95 GO:0040007  15459   13045
    ## 98 GO:0140110    556     630
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 4                              nucleic_acid_binding molecular_function
    ## 5                                       DNA_binding molecular_function
    ## 7         DNA-binding_transcription_factor_activity molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 12                       signaling_receptor_binding molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 18                             extracellular_region cellular_component
    ## 20                                        cell_wall cellular_component
    ## 21                                    intracellular cellular_component
    ## 22                                             cell cellular_component
    ## 23                                          nucleus cellular_component
    ## 27                                        cytoplasm cellular_component
    ## 29                                         lysosome cellular_component
    ## 30                                         endosome cellular_component
    ## 31                                          vacuole cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 34                                  Golgi_apparatus cellular_component
    ## 35                                          cytosol cellular_component
    ## 36                                         ribosome cellular_component
    ## 38                                  plasma_membrane cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 44            cellular_protein_modification_process biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 47                               response_to_stress biological_process
    ## 48                                       cell_cycle biological_process
    ## 52               multicellular_organism_development biological_process
    ## 55                                metabolic_process biological_process
    ## 56                                       cell_death biological_process
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 63                    response_to_external_stimulus biological_process
    ## 64                                          tropism biological_process
    ## 65                      response_to_biotic_stimulus biological_process
    ## 66                     response_to_abiotic_stimulus biological_process
    ## 67               anatomical_structure_morphogenesis biological_process
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 69                               embryo_development biological_process
    ## 72                                       abscission biological_process
    ## 75                               flower_development biological_process
    ## 76                                 cellular_process biological_process
    ## 78                                   photosynthesis biological_process
    ## 80                  cellular_component_organization biological_process
    ## 82                                  kinase_activity molecular_function
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ## 85                        protein_metabolic_process biological_process
    ## 86                             cellular_homeostasis biological_process
    ## 87                      secondary_metabolic_process biological_process
    ## 90                        enzyme_regulator_activity molecular_function
    ## 91                             carbohydrate_binding molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ## 95                                           growth biological_process
    ## 98                 transcription_regulator_activity molecular_function
    ##      prob_slim50 prob_slim50BH
    ## 2  1.451071e-278 1.422049e-276
    ## 3   2.990124e-05  5.232718e-05
    ## 4   2.598797e-06  5.197594e-06
    ## 5   5.225994e-33  2.438797e-32
    ## 7   3.026615e-09  7.234348e-09
    ## 8   9.312424e-65  6.518697e-64
    ## 10 1.398657e-178 3.426710e-177
    ## 12  5.953258e-31  2.536606e-30
    ## 13  4.939786e-32  2.200450e-31
    ## 14  6.974583e-29  2.531515e-28
    ## 15 2.644299e-155 5.182826e-154
    ## 18  4.841928e-06  9.125173e-06
    ## 20  4.079489e-06  7.839019e-06
    ## 21  2.653108e-04  4.127058e-04
    ## 22  4.273401e-05  7.098191e-05
    ## 23  1.210118e-61  7.906107e-61
    ## 27  1.562362e-43  9.006556e-43
    ## 29  3.892658e-05  6.692640e-05
    ## 30  4.212441e-05  7.098191e-05
    ## 31  7.061020e-08  1.472298e-07
    ## 33  1.300526e-11  3.267988e-11
    ## 34  2.753257e-06  5.396384e-06
    ## 35  2.194308e-08  4.887322e-08
    ## 36  3.794354e-74  3.718467e-73
    ## 38  1.072246e-69  8.756679e-69
    ## 39  7.601677e-37  3.724822e-36
    ## 40  6.313427e-31  2.577983e-30
    ## 41  2.825675e-26  9.548831e-26
    ## 43  9.119932e-41  4.965296e-40
    ## 44  5.494027e-22  1.794716e-21
    ## 45  3.101287e-15  8.939004e-15
    ## 46  8.837211e-38  4.558141e-37
    ## 47  2.779994e-98  3.891991e-97
    ## 48  2.151502e-19  6.801524e-19
    ## 52  2.694025e-05  4.800263e-05
    ## 55  7.023507e-15  1.966582e-14
    ## 56  1.747617e-30  6.587170e-30
    ## 59  2.379677e-14  6.478009e-14
    ## 60  2.049377e-61  1.255243e-60
    ## 61  6.919265e-66  5.216062e-65
    ## 62  8.185165e-31  3.208585e-30
    ## 63  1.077113e-79  1.172856e-78
    ## 64  3.904692e-04  5.979060e-04
    ## 65  5.386996e-72  4.799323e-71
    ## 66  7.160336e-18  2.192853e-17
    ## 67  2.533805e-05  4.598386e-05
    ## 68  5.348479e-06  9.889640e-06
    ## 69  7.007913e-08  1.472298e-07
    ## 72  1.933099e-04  3.055544e-04
    ## 75  1.022413e-08  2.330151e-08
    ## 76  4.500916e-05  7.351495e-05
    ## 78 9.220117e-107 1.505952e-105
    ## 80  3.573927e-16  1.061348e-15
    ## 82 1.569032e-227 7.688258e-226
    ## 83 9.476026e-212 3.095502e-210
    ## 84  1.815824e-26  6.355386e-26
    ## 85  1.022553e-07  2.087712e-07
    ## 86  9.940943e-12  2.563717e-11
    ## 87  1.163243e-13  3.081023e-13
    ## 90  3.080445e-08  6.708526e-08
    ## 91  3.971000e-09  9.265667e-09
    ## 94  1.410929e-80  1.728388e-79
    ## 95  7.513134e-10  1.840718e-09
    ## 98  1.357531e-04  2.180951e-04

``` r
q903_50_wsSlimAllSubs$prob_slim50Bon = p.adjust(prob_slim50, method = "bonferroni", n = length(prob_slim50))
subset(q903_50_wsSlimAllSubs,q903_50_wsSlimAllSubs$prob_slim50Bon < 0.05)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166 497561  419196
    ## 3  GO:0003674  33330   31271
    ## 4  GO:0003676   6335    6254
    ## 5  GO:0003677  37207   36884
    ## 7  GO:0003700  13200   12897
    ## 8  GO:0003723  46711   47375
    ## 10 GO:0003824 538560  516847
    ## 12 GO:0005102  12488    9696
    ## 13 GO:0005198  20197   20604
    ## 14 GO:0005215  86105   82538
    ## 15 GO:0005488 482975  463095
    ## 18 GO:0005576  98916   91707
    ## 20 GO:0005618  51135   47806
    ## 21 GO:0005622  57607   51155
    ## 22 GO:0005623   8575    7294
    ## 23 GO:0005634  83159   81913
    ## 27 GO:0005737  93411   90460
    ## 29 GO:0005764   1238     941
    ## 30 GO:0005768  16321   15514
    ## 31 GO:0005773  35745   33807
    ## 33 GO:0005783  54017   51133
    ## 34 GO:0005794  23256   22065
    ## 35 GO:0005829  29841   28379
    ## 36 GO:0005840  28211   29793
    ## 38 GO:0005886 208432  178771
    ## 39 GO:0005975  53502   52508
    ## 40 GO:0006091  28831   28821
    ## 41 GO:0006139 206689  181357
    ## 43 GO:0006412  28066   28511
    ## 44 GO:0006464 282505  249784
    ## 45 GO:0006629 100279   94358
    ## 46 GO:0006810 150901  143646
    ## 47 GO:0006950 321704  276566
    ## 48 GO:0007049  28832   28226
    ## 52 GO:0007275  68302   60575
    ## 55 GO:0008152 393454  363668
    ## 56 GO:0008219  10037    7658
    ## 59 GO:0009058 178203  166054
    ## 60 GO:0009507 150086  144818
    ## 61 GO:0009536  93814   92219
    ## 62 GO:0009579  78643   75720
    ## 63 GO:0009605 102581   85325
    ## 64 GO:0009606   1574    1620
    ## 65 GO:0009607 102051   85252
    ## 66 GO:0009628  80933   76737
    ## 67 GO:0009653  36351   31954
    ## 68 GO:0009719  86313   80132
    ## 69 GO:0009790  15476   13181
    ## 72 GO:0009838   3415    2820
    ## 75 GO:0009908  22575   19378
    ## 76 GO:0009987 463805  417436
    ## 78 GO:0015979  40959   43255
    ## 80 GO:0016043  87456   82604
    ## 82 GO:0016301 650502  556870
    ## 83 GO:0016740 552764  471897
    ## 84 GO:0016787 240853  225597
    ## 85 GO:0019538  24081   22963
    ## 86 GO:0019725   3313    3545
    ## 87 GO:0019748  38030   36457
    ## 90 GO:0030234  11813   11531
    ## 91 GO:0030246  27493   23692
    ## 94 GO:0038023  68075   55439
    ## 95 GO:0040007  15459   13045
    ## 98 GO:0140110    556     630
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 4                              nucleic_acid_binding molecular_function
    ## 5                                       DNA_binding molecular_function
    ## 7         DNA-binding_transcription_factor_activity molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 12                       signaling_receptor_binding molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 18                             extracellular_region cellular_component
    ## 20                                        cell_wall cellular_component
    ## 21                                    intracellular cellular_component
    ## 22                                             cell cellular_component
    ## 23                                          nucleus cellular_component
    ## 27                                        cytoplasm cellular_component
    ## 29                                         lysosome cellular_component
    ## 30                                         endosome cellular_component
    ## 31                                          vacuole cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 34                                  Golgi_apparatus cellular_component
    ## 35                                          cytosol cellular_component
    ## 36                                         ribosome cellular_component
    ## 38                                  plasma_membrane cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 44            cellular_protein_modification_process biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 47                               response_to_stress biological_process
    ## 48                                       cell_cycle biological_process
    ## 52               multicellular_organism_development biological_process
    ## 55                                metabolic_process biological_process
    ## 56                                       cell_death biological_process
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 63                    response_to_external_stimulus biological_process
    ## 64                                          tropism biological_process
    ## 65                      response_to_biotic_stimulus biological_process
    ## 66                     response_to_abiotic_stimulus biological_process
    ## 67               anatomical_structure_morphogenesis biological_process
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 69                               embryo_development biological_process
    ## 72                                       abscission biological_process
    ## 75                               flower_development biological_process
    ## 76                                 cellular_process biological_process
    ## 78                                   photosynthesis biological_process
    ## 80                  cellular_component_organization biological_process
    ## 82                                  kinase_activity molecular_function
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ## 85                        protein_metabolic_process biological_process
    ## 86                             cellular_homeostasis biological_process
    ## 87                      secondary_metabolic_process biological_process
    ## 90                        enzyme_regulator_activity molecular_function
    ## 91                             carbohydrate_binding molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ## 95                                           growth biological_process
    ## 98                 transcription_regulator_activity molecular_function
    ##      prob_slim50 prob_slim50BH prob_slim50Bon
    ## 2  1.451071e-278 1.422049e-276  1.422049e-276
    ## 3   2.990124e-05  5.232718e-05   2.930322e-03
    ## 4   2.598797e-06  5.197594e-06   2.546821e-04
    ## 5   5.225994e-33  2.438797e-32   5.121474e-31
    ## 7   3.026615e-09  7.234348e-09   2.966082e-07
    ## 8   9.312424e-65  6.518697e-64   9.126176e-63
    ## 10 1.398657e-178 3.426710e-177  1.370684e-176
    ## 12  5.953258e-31  2.536606e-30   5.834193e-29
    ## 13  4.939786e-32  2.200450e-31   4.840990e-30
    ## 14  6.974583e-29  2.531515e-28   6.835092e-27
    ## 15 2.644299e-155 5.182826e-154  2.591413e-153
    ## 18  4.841928e-06  9.125173e-06   4.745090e-04
    ## 20  4.079489e-06  7.839019e-06   3.997899e-04
    ## 21  2.653108e-04  4.127058e-04   2.600046e-02
    ## 22  4.273401e-05  7.098191e-05   4.187933e-03
    ## 23  1.210118e-61  7.906107e-61   1.185916e-59
    ## 27  1.562362e-43  9.006556e-43   1.531115e-41
    ## 29  3.892658e-05  6.692640e-05   3.814805e-03
    ## 30  4.212441e-05  7.098191e-05   4.128192e-03
    ## 31  7.061020e-08  1.472298e-07   6.919800e-06
    ## 33  1.300526e-11  3.267988e-11   1.274515e-09
    ## 34  2.753257e-06  5.396384e-06   2.698192e-04
    ## 35  2.194308e-08  4.887322e-08   2.150422e-06
    ## 36  3.794354e-74  3.718467e-73   3.718467e-72
    ## 38  1.072246e-69  8.756679e-69   1.050802e-67
    ## 39  7.601677e-37  3.724822e-36   7.449644e-35
    ## 40  6.313427e-31  2.577983e-30   6.187158e-29
    ## 41  2.825675e-26  9.548831e-26   2.769161e-24
    ## 43  9.119932e-41  4.965296e-40   8.937533e-39
    ## 44  5.494027e-22  1.794716e-21   5.384147e-20
    ## 45  3.101287e-15  8.939004e-15   3.039261e-13
    ## 46  8.837211e-38  4.558141e-37   8.660467e-36
    ## 47  2.779994e-98  3.891991e-97   2.724394e-96
    ## 48  2.151502e-19  6.801524e-19   2.108472e-17
    ## 52  2.694025e-05  4.800263e-05   2.640145e-03
    ## 55  7.023507e-15  1.966582e-14   6.883037e-13
    ## 56  1.747617e-30  6.587170e-30   1.712664e-28
    ## 59  2.379677e-14  6.478009e-14   2.332083e-12
    ## 60  2.049377e-61  1.255243e-60   2.008389e-59
    ## 61  6.919265e-66  5.216062e-65   6.780880e-64
    ## 62  8.185165e-31  3.208585e-30   8.021462e-29
    ## 63  1.077113e-79  1.172856e-78   1.055570e-77
    ## 64  3.904692e-04  5.979060e-04   3.826598e-02
    ## 65  5.386996e-72  4.799323e-71   5.279256e-70
    ## 66  7.160336e-18  2.192853e-17   7.017129e-16
    ## 67  2.533805e-05  4.598386e-05   2.483128e-03
    ## 68  5.348479e-06  9.889640e-06   5.241509e-04
    ## 69  7.007913e-08  1.472298e-07   6.867755e-06
    ## 72  1.933099e-04  3.055544e-04   1.894437e-02
    ## 75  1.022413e-08  2.330151e-08   1.001965e-06
    ## 76  4.500916e-05  7.351495e-05   4.410897e-03
    ## 78 9.220117e-107 1.505952e-105  9.035714e-105
    ## 80  3.573927e-16  1.061348e-15   3.502449e-14
    ## 82 1.569032e-227 7.688258e-226  1.537652e-225
    ## 83 9.476026e-212 3.095502e-210  9.286505e-210
    ## 84  1.815824e-26  6.355386e-26   1.779508e-24
    ## 85  1.022553e-07  2.087712e-07   1.002102e-05
    ## 86  9.940943e-12  2.563717e-11   9.742124e-10
    ## 87  1.163243e-13  3.081023e-13   1.139979e-11
    ## 90  3.080445e-08  6.708526e-08   3.018837e-06
    ## 91  3.971000e-09  9.265667e-09   3.891580e-07
    ## 94  1.410929e-80  1.728388e-79   1.382710e-78
    ## 95  7.513134e-10  1.840718e-09   7.362871e-08
    ## 98  1.357531e-04  2.180951e-04   1.330380e-02

``` r
subset(q903_50_wsSlimAllSubs,q903_50_wsSlimAllSubs$prob_slim50Bon < 0.001)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166 497561  419196
    ## 4  GO:0003676   6335    6254
    ## 5  GO:0003677  37207   36884
    ## 7  GO:0003700  13200   12897
    ## 8  GO:0003723  46711   47375
    ## 10 GO:0003824 538560  516847
    ## 12 GO:0005102  12488    9696
    ## 13 GO:0005198  20197   20604
    ## 14 GO:0005215  86105   82538
    ## 15 GO:0005488 482975  463095
    ## 18 GO:0005576  98916   91707
    ## 20 GO:0005618  51135   47806
    ## 23 GO:0005634  83159   81913
    ## 27 GO:0005737  93411   90460
    ## 31 GO:0005773  35745   33807
    ## 33 GO:0005783  54017   51133
    ## 34 GO:0005794  23256   22065
    ## 35 GO:0005829  29841   28379
    ## 36 GO:0005840  28211   29793
    ## 38 GO:0005886 208432  178771
    ## 39 GO:0005975  53502   52508
    ## 40 GO:0006091  28831   28821
    ## 41 GO:0006139 206689  181357
    ## 43 GO:0006412  28066   28511
    ## 44 GO:0006464 282505  249784
    ## 45 GO:0006629 100279   94358
    ## 46 GO:0006810 150901  143646
    ## 47 GO:0006950 321704  276566
    ## 48 GO:0007049  28832   28226
    ## 55 GO:0008152 393454  363668
    ## 56 GO:0008219  10037    7658
    ## 59 GO:0009058 178203  166054
    ## 60 GO:0009507 150086  144818
    ## 61 GO:0009536  93814   92219
    ## 62 GO:0009579  78643   75720
    ## 63 GO:0009605 102581   85325
    ## 65 GO:0009607 102051   85252
    ## 66 GO:0009628  80933   76737
    ## 68 GO:0009719  86313   80132
    ## 69 GO:0009790  15476   13181
    ## 75 GO:0009908  22575   19378
    ## 78 GO:0015979  40959   43255
    ## 80 GO:0016043  87456   82604
    ## 82 GO:0016301 650502  556870
    ## 83 GO:0016740 552764  471897
    ## 84 GO:0016787 240853  225597
    ## 85 GO:0019538  24081   22963
    ## 86 GO:0019725   3313    3545
    ## 87 GO:0019748  38030   36457
    ## 90 GO:0030234  11813   11531
    ## 91 GO:0030246  27493   23692
    ## 94 GO:0038023  68075   55439
    ## 95 GO:0040007  15459   13045
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 4                              nucleic_acid_binding molecular_function
    ## 5                                       DNA_binding molecular_function
    ## 7         DNA-binding_transcription_factor_activity molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 12                       signaling_receptor_binding molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 18                             extracellular_region cellular_component
    ## 20                                        cell_wall cellular_component
    ## 23                                          nucleus cellular_component
    ## 27                                        cytoplasm cellular_component
    ## 31                                          vacuole cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 34                                  Golgi_apparatus cellular_component
    ## 35                                          cytosol cellular_component
    ## 36                                         ribosome cellular_component
    ## 38                                  plasma_membrane cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 44            cellular_protein_modification_process biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 47                               response_to_stress biological_process
    ## 48                                       cell_cycle biological_process
    ## 55                                metabolic_process biological_process
    ## 56                                       cell_death biological_process
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 63                    response_to_external_stimulus biological_process
    ## 65                      response_to_biotic_stimulus biological_process
    ## 66                     response_to_abiotic_stimulus biological_process
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 69                               embryo_development biological_process
    ## 75                               flower_development biological_process
    ## 78                                   photosynthesis biological_process
    ## 80                  cellular_component_organization biological_process
    ## 82                                  kinase_activity molecular_function
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ## 85                        protein_metabolic_process biological_process
    ## 86                             cellular_homeostasis biological_process
    ## 87                      secondary_metabolic_process biological_process
    ## 90                        enzyme_regulator_activity molecular_function
    ## 91                             carbohydrate_binding molecular_function
    ## 94                      signaling_receptor_activity molecular_function
    ## 95                                           growth biological_process
    ##      prob_slim50 prob_slim50BH prob_slim50Bon
    ## 2  1.451071e-278 1.422049e-276  1.422049e-276
    ## 4   2.598797e-06  5.197594e-06   2.546821e-04
    ## 5   5.225994e-33  2.438797e-32   5.121474e-31
    ## 7   3.026615e-09  7.234348e-09   2.966082e-07
    ## 8   9.312424e-65  6.518697e-64   9.126176e-63
    ## 10 1.398657e-178 3.426710e-177  1.370684e-176
    ## 12  5.953258e-31  2.536606e-30   5.834193e-29
    ## 13  4.939786e-32  2.200450e-31   4.840990e-30
    ## 14  6.974583e-29  2.531515e-28   6.835092e-27
    ## 15 2.644299e-155 5.182826e-154  2.591413e-153
    ## 18  4.841928e-06  9.125173e-06   4.745090e-04
    ## 20  4.079489e-06  7.839019e-06   3.997899e-04
    ## 23  1.210118e-61  7.906107e-61   1.185916e-59
    ## 27  1.562362e-43  9.006556e-43   1.531115e-41
    ## 31  7.061020e-08  1.472298e-07   6.919800e-06
    ## 33  1.300526e-11  3.267988e-11   1.274515e-09
    ## 34  2.753257e-06  5.396384e-06   2.698192e-04
    ## 35  2.194308e-08  4.887322e-08   2.150422e-06
    ## 36  3.794354e-74  3.718467e-73   3.718467e-72
    ## 38  1.072246e-69  8.756679e-69   1.050802e-67
    ## 39  7.601677e-37  3.724822e-36   7.449644e-35
    ## 40  6.313427e-31  2.577983e-30   6.187158e-29
    ## 41  2.825675e-26  9.548831e-26   2.769161e-24
    ## 43  9.119932e-41  4.965296e-40   8.937533e-39
    ## 44  5.494027e-22  1.794716e-21   5.384147e-20
    ## 45  3.101287e-15  8.939004e-15   3.039261e-13
    ## 46  8.837211e-38  4.558141e-37   8.660467e-36
    ## 47  2.779994e-98  3.891991e-97   2.724394e-96
    ## 48  2.151502e-19  6.801524e-19   2.108472e-17
    ## 55  7.023507e-15  1.966582e-14   6.883037e-13
    ## 56  1.747617e-30  6.587170e-30   1.712664e-28
    ## 59  2.379677e-14  6.478009e-14   2.332083e-12
    ## 60  2.049377e-61  1.255243e-60   2.008389e-59
    ## 61  6.919265e-66  5.216062e-65   6.780880e-64
    ## 62  8.185165e-31  3.208585e-30   8.021462e-29
    ## 63  1.077113e-79  1.172856e-78   1.055570e-77
    ## 65  5.386996e-72  4.799323e-71   5.279256e-70
    ## 66  7.160336e-18  2.192853e-17   7.017129e-16
    ## 68  5.348479e-06  9.889640e-06   5.241509e-04
    ## 69  7.007913e-08  1.472298e-07   6.867755e-06
    ## 75  1.022413e-08  2.330151e-08   1.001965e-06
    ## 78 9.220117e-107 1.505952e-105  9.035714e-105
    ## 80  3.573927e-16  1.061348e-15   3.502449e-14
    ## 82 1.569032e-227 7.688258e-226  1.537652e-225
    ## 83 9.476026e-212 3.095502e-210  9.286505e-210
    ## 84  1.815824e-26  6.355386e-26   1.779508e-24
    ## 85  1.022553e-07  2.087712e-07   1.002102e-05
    ## 86  9.940943e-12  2.563717e-11   9.742124e-10
    ## 87  1.163243e-13  3.081023e-13   1.139979e-11
    ## 90  3.080445e-08  6.708526e-08   3.018837e-06
    ## 91  3.971000e-09  9.265667e-09   3.891580e-07
    ## 94  1.410929e-80  1.728388e-79   1.382710e-78
    ## 95  7.513134e-10  1.840718e-09   7.362871e-08

Test go terms mapped to GO slim - Based on BLAST and coverage 95% from subject eval 1e-5
----------------------------------------------------------------------------------------

``` r
q903_95 <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/Q903/ThresholdCovIdentity3/my_Slimgo_terms_plant.Q903cov95.counts", header=FALSE)
ws_95 <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/WS77111/ThresholdCovIdentity3/my_Slimgo_terms_plant.WS77111cov95.counts", header=FALSE)

##bind the 2 classes
q903_95_wsSlimAll = merge(q903_95,ws_95,by="V2")
q903_95_wsSlimAll = q903_95_wsSlimAll[,c(1,2,5,3,4)]
colnames(q903_95_wsSlimAll) = c("Go_term","Q903","WS77111","Domain","Go_name")
dim(q903_95_wsSlimAll)
```

    ## [1] 98  5

``` r
#############All domains together
colSums_q903_95_wsSlimAll = unname(apply(q903_95_wsSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_q903_95_wsSlimAll)
p1 =colSums_q903_95_wsSlimAll / GrandTotal

p1
```

    ## [1] 0.5415536 0.4584464

``` r
#filter for low counts
q903_95_wsSlimAllSubs = subset(q903_95_wsSlimAll,q903_95_wsSlimAll$Q903 > 5 & q903_95_wsSlimAll$WS77111 > 5)
dim(q903_95_wsSlimAllSubs)
```

    ## [1] 97  5

``` r
prob_slim95 = c()
for (i in 1:nrow(q903_95_wsSlimAllSubs)){
  #print(i)
  obs = q903_95_wsSlimAllSubs[i,c("Q903","WS77111")]
  #print(obs)
  prob_slim95[i] = chisq.test(obs, p=p1)$p.value
}

sum(p.adjust(prob_slim95, method = "BH", n = length(prob_slim95)) < 0.05)
```

    ## [1] 40

``` r
sum(p.adjust(prob_slim95, method = "BH", n = length(prob_slim95)) < 0.001)
```

    ## [1] 27

``` r
q903_95_wsSlimAllSubs$prob_slim95 = prob_slim95

q903_95_wsSlimAllSubs$prob_slim95BH = p.adjust(prob_slim95, method = "BH", n = length(prob_slim95))

subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95BH < 0.05) 
```

    ##       Go_term   Q903 WS77111
    ## 1  GO:0000003   2370    1858
    ## 2  GO:0000166  25636   20233
    ## 3  GO:0003674   1480    1510
    ## 8  GO:0003723   3327    3844
    ## 10 GO:0003824  79250   72020
    ## 13 GO:0005198   2661    2620
    ## 14 GO:0005215  16279   14177
    ## 15 GO:0005488  65295   58154
    ## 21 GO:0005622  16919   12843
    ## 23 GO:0005634   5343    4743
    ## 28 GO:0005739   4283    3263
    ## 29 GO:0005764    202     119
    ## 32 GO:0005777    628     661
    ## 33 GO:0005783   7463    7101
    ## 36 GO:0005840   2682    4049
    ## 37 GO:0005856   3873    2220
    ## 40 GO:0006091   5326    5020
    ## 41 GO:0006139  14360   10681
    ## 43 GO:0006412   3416    3610
    ## 45 GO:0006629   9502    8712
    ## 46 GO:0006810  22449   19969
    ## 50 GO:0007165   4492    4061
    ## 52 GO:0007275   4802    3785
    ## 55 GO:0008152  39175   34679
    ## 57 GO:0008289    621     921
    ## 59 GO:0009058  25114   20348
    ## 60 GO:0009507  13065   14634
    ## 61 GO:0009536   7446    8686
    ## 62 GO:0009579  11393   12628
    ## 63 GO:0009605   6365    5134
    ## 67 GO:0009653   2883    2211
    ## 68 GO:0009719   5073    4076
    ## 70 GO:0009791   3073    2443
    ## 71 GO:0009835    267     295
    ## 75 GO:0009908   1230     932
    ## 78 GO:0015979   7206    7082
    ## 79 GO:0016020  64212   61390
    ## 83 GO:0016740 149762  102548
    ## 84 GO:0016787  28060   22432
    ## 86 GO:0019725    405     417
    ##                                              Domain            Go_name
    ## 1                                      reproduction biological_process
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 21                                    intracellular cellular_component
    ## 23                                          nucleus cellular_component
    ## 28                                    mitochondrion cellular_component
    ## 29                                         lysosome cellular_component
    ## 32                                       peroxisome cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 50                              signal_transduction biological_process
    ## 52               multicellular_organism_development biological_process
    ## 55                                metabolic_process biological_process
    ## 57                                    lipid_binding molecular_function
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 63                    response_to_external_stimulus biological_process
    ## 67               anatomical_structure_morphogenesis biological_process
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 70                       post-embryonic_development biological_process
    ## 71                                   fruit_ripening biological_process
    ## 75                               flower_development biological_process
    ## 78                                   photosynthesis biological_process
    ## 79                                         membrane cellular_component
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ## 86                             cellular_homeostasis biological_process
    ##      prob_slim95 prob_slim95BH
    ## 1   1.318190e-02  3.455795e-02
    ## 2   9.038790e-14  5.479767e-13
    ## 3   3.209772e-07  1.353687e-06
    ## 8   1.021985e-39  1.101472e-38
    ## 10  3.285246e-43  3.983361e-42
    ## 13  3.923563e-08  1.902928e-07
    ## 14  1.360961e-02  3.474033e-02
    ## 15  5.266688e-19  3.405791e-18
    ## 21  1.146308e-20  9.265993e-20
    ## 23  1.730090e-02  4.303045e-02
    ## 28  5.669397e-06  2.115121e-05
    ## 29  1.607553e-03  5.030084e-03
    ## 32  8.985508e-05  3.228127e-04
    ## 33  1.735162e-12  9.900629e-12
    ## 36 9.478089e-123 4.596873e-121
    ## 37  3.542167e-49  4.908432e-48
    ## 40  4.661058e-08  2.152965e-07
    ## 41  3.949209e-24  3.482485e-23
    ## 43  1.244879e-20  9.288714e-20
    ## 45  7.403862e-08  3.264430e-07
    ## 46  3.530313e-07  1.426835e-06
    ## 50  2.396475e-03  7.264315e-03
    ## 52  1.019691e-03  3.410689e-03
    ## 55  1.341788e-09  6.850179e-09
    ## 57  7.332310e-28  7.112341e-27
    ## 59  3.338711e-06  1.295420e-05
    ## 60 1.754950e-120 5.674337e-119
    ## 61  2.093704e-92  3.384821e-91
    ## 62  3.424592e-97  6.643709e-96
    ## 63  9.975459e-03  2.845940e-02
    ## 67  4.723484e-04  1.636350e-03
    ## 68  1.303818e-02  3.455795e-02
    ## 70  2.043550e-02  4.955609e-02
    ## 71  1.565645e-03  5.030084e-03
    ## 75  1.066327e-02  2.955249e-02
    ## 78  4.356859e-19  3.018681e-18
    ## 79 3.793034e-103 9.198107e-102
    ## 83  0.000000e+00  0.000000e+00
    ## 84  1.617813e-10  8.718216e-10
    ## 86  4.938759e-03  1.451696e-02

``` r
subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95BH < 0.001)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166  25636   20233
    ## 3  GO:0003674   1480    1510
    ## 8  GO:0003723   3327    3844
    ## 10 GO:0003824  79250   72020
    ## 13 GO:0005198   2661    2620
    ## 15 GO:0005488  65295   58154
    ## 21 GO:0005622  16919   12843
    ## 28 GO:0005739   4283    3263
    ## 32 GO:0005777    628     661
    ## 33 GO:0005783   7463    7101
    ## 36 GO:0005840   2682    4049
    ## 37 GO:0005856   3873    2220
    ## 40 GO:0006091   5326    5020
    ## 41 GO:0006139  14360   10681
    ## 43 GO:0006412   3416    3610
    ## 45 GO:0006629   9502    8712
    ## 46 GO:0006810  22449   19969
    ## 55 GO:0008152  39175   34679
    ## 57 GO:0008289    621     921
    ## 59 GO:0009058  25114   20348
    ## 60 GO:0009507  13065   14634
    ## 61 GO:0009536   7446    8686
    ## 62 GO:0009579  11393   12628
    ## 78 GO:0015979   7206    7082
    ## 79 GO:0016020  64212   61390
    ## 83 GO:0016740 149762  102548
    ## 84 GO:0016787  28060   22432
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 15                                          binding molecular_function
    ## 21                                    intracellular cellular_component
    ## 28                                    mitochondrion cellular_component
    ## 32                                       peroxisome cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 55                                metabolic_process biological_process
    ## 57                                    lipid_binding molecular_function
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 78                                   photosynthesis biological_process
    ## 79                                         membrane cellular_component
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ##      prob_slim95 prob_slim95BH
    ## 2   9.038790e-14  5.479767e-13
    ## 3   3.209772e-07  1.353687e-06
    ## 8   1.021985e-39  1.101472e-38
    ## 10  3.285246e-43  3.983361e-42
    ## 13  3.923563e-08  1.902928e-07
    ## 15  5.266688e-19  3.405791e-18
    ## 21  1.146308e-20  9.265993e-20
    ## 28  5.669397e-06  2.115121e-05
    ## 32  8.985508e-05  3.228127e-04
    ## 33  1.735162e-12  9.900629e-12
    ## 36 9.478089e-123 4.596873e-121
    ## 37  3.542167e-49  4.908432e-48
    ## 40  4.661058e-08  2.152965e-07
    ## 41  3.949209e-24  3.482485e-23
    ## 43  1.244879e-20  9.288714e-20
    ## 45  7.403862e-08  3.264430e-07
    ## 46  3.530313e-07  1.426835e-06
    ## 55  1.341788e-09  6.850179e-09
    ## 57  7.332310e-28  7.112341e-27
    ## 59  3.338711e-06  1.295420e-05
    ## 60 1.754950e-120 5.674337e-119
    ## 61  2.093704e-92  3.384821e-91
    ## 62  3.424592e-97  6.643709e-96
    ## 78  4.356859e-19  3.018681e-18
    ## 79 3.793034e-103 9.198107e-102
    ## 83  0.000000e+00  0.000000e+00
    ## 84  1.617813e-10  8.718216e-10

``` r
q903_95_wsSlimAllSubs$prob_slim95Bon = p.adjust(prob_slim95, method = "bonferroni", n = length(prob_slim95))
subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95Bon < 0.05)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166  25636   20233
    ## 3  GO:0003674   1480    1510
    ## 8  GO:0003723   3327    3844
    ## 10 GO:0003824  79250   72020
    ## 13 GO:0005198   2661    2620
    ## 15 GO:0005488  65295   58154
    ## 21 GO:0005622  16919   12843
    ## 28 GO:0005739   4283    3263
    ## 32 GO:0005777    628     661
    ## 33 GO:0005783   7463    7101
    ## 36 GO:0005840   2682    4049
    ## 37 GO:0005856   3873    2220
    ## 40 GO:0006091   5326    5020
    ## 41 GO:0006139  14360   10681
    ## 43 GO:0006412   3416    3610
    ## 45 GO:0006629   9502    8712
    ## 46 GO:0006810  22449   19969
    ## 55 GO:0008152  39175   34679
    ## 57 GO:0008289    621     921
    ## 59 GO:0009058  25114   20348
    ## 60 GO:0009507  13065   14634
    ## 61 GO:0009536   7446    8686
    ## 62 GO:0009579  11393   12628
    ## 67 GO:0009653   2883    2211
    ## 78 GO:0015979   7206    7082
    ## 79 GO:0016020  64212   61390
    ## 83 GO:0016740 149762  102548
    ## 84 GO:0016787  28060   22432
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 15                                          binding molecular_function
    ## 21                                    intracellular cellular_component
    ## 28                                    mitochondrion cellular_component
    ## 32                                       peroxisome cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 55                                metabolic_process biological_process
    ## 57                                    lipid_binding molecular_function
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 67               anatomical_structure_morphogenesis biological_process
    ## 78                                   photosynthesis biological_process
    ## 79                                         membrane cellular_component
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ##      prob_slim95 prob_slim95BH prob_slim95Bon
    ## 2   9.038790e-14  5.479767e-13   8.767626e-12
    ## 3   3.209772e-07  1.353687e-06   3.113479e-05
    ## 8   1.021985e-39  1.101472e-38   9.913252e-38
    ## 10  3.285246e-43  3.983361e-42   3.186689e-41
    ## 13  3.923563e-08  1.902928e-07   3.805856e-06
    ## 15  5.266688e-19  3.405791e-18   5.108687e-17
    ## 21  1.146308e-20  9.265993e-20   1.111919e-18
    ## 28  5.669397e-06  2.115121e-05   5.499315e-04
    ## 32  8.985508e-05  3.228127e-04   8.715943e-03
    ## 33  1.735162e-12  9.900629e-12   1.683107e-10
    ## 36 9.478089e-123 4.596873e-121  9.193746e-121
    ## 37  3.542167e-49  4.908432e-48   3.435902e-47
    ## 40  4.661058e-08  2.152965e-07   4.521226e-06
    ## 41  3.949209e-24  3.482485e-23   3.830733e-22
    ## 43  1.244879e-20  9.288714e-20   1.207533e-18
    ## 45  7.403862e-08  3.264430e-07   7.181746e-06
    ## 46  3.530313e-07  1.426835e-06   3.424404e-05
    ## 55  1.341788e-09  6.850179e-09   1.301534e-07
    ## 57  7.332310e-28  7.112341e-27   7.112341e-26
    ## 59  3.338711e-06  1.295420e-05   3.238549e-04
    ## 60 1.754950e-120 5.674337e-119  1.702301e-118
    ## 61  2.093704e-92  3.384821e-91   2.030893e-90
    ## 62  3.424592e-97  6.643709e-96   3.321854e-95
    ## 67  4.723484e-04  1.636350e-03   4.581780e-02
    ## 78  4.356859e-19  3.018681e-18   4.226153e-17
    ## 79 3.793034e-103 9.198107e-102  3.679243e-101
    ## 83  0.000000e+00  0.000000e+00   0.000000e+00
    ## 84  1.617813e-10  8.718216e-10   1.569279e-08

``` r
subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95Bon < 0.001)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166  25636   20233
    ## 3  GO:0003674   1480    1510
    ## 8  GO:0003723   3327    3844
    ## 10 GO:0003824  79250   72020
    ## 13 GO:0005198   2661    2620
    ## 15 GO:0005488  65295   58154
    ## 21 GO:0005622  16919   12843
    ## 28 GO:0005739   4283    3263
    ## 33 GO:0005783   7463    7101
    ## 36 GO:0005840   2682    4049
    ## 37 GO:0005856   3873    2220
    ## 40 GO:0006091   5326    5020
    ## 41 GO:0006139  14360   10681
    ## 43 GO:0006412   3416    3610
    ## 45 GO:0006629   9502    8712
    ## 46 GO:0006810  22449   19969
    ## 55 GO:0008152  39175   34679
    ## 57 GO:0008289    621     921
    ## 59 GO:0009058  25114   20348
    ## 60 GO:0009507  13065   14634
    ## 61 GO:0009536   7446    8686
    ## 62 GO:0009579  11393   12628
    ## 78 GO:0015979   7206    7082
    ## 79 GO:0016020  64212   61390
    ## 83 GO:0016740 149762  102548
    ## 84 GO:0016787  28060   22432
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 15                                          binding molecular_function
    ## 21                                    intracellular cellular_component
    ## 28                                    mitochondrion cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 55                                metabolic_process biological_process
    ## 57                                    lipid_binding molecular_function
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 78                                   photosynthesis biological_process
    ## 79                                         membrane cellular_component
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ##      prob_slim95 prob_slim95BH prob_slim95Bon
    ## 2   9.038790e-14  5.479767e-13   8.767626e-12
    ## 3   3.209772e-07  1.353687e-06   3.113479e-05
    ## 8   1.021985e-39  1.101472e-38   9.913252e-38
    ## 10  3.285246e-43  3.983361e-42   3.186689e-41
    ## 13  3.923563e-08  1.902928e-07   3.805856e-06
    ## 15  5.266688e-19  3.405791e-18   5.108687e-17
    ## 21  1.146308e-20  9.265993e-20   1.111919e-18
    ## 28  5.669397e-06  2.115121e-05   5.499315e-04
    ## 33  1.735162e-12  9.900629e-12   1.683107e-10
    ## 36 9.478089e-123 4.596873e-121  9.193746e-121
    ## 37  3.542167e-49  4.908432e-48   3.435902e-47
    ## 40  4.661058e-08  2.152965e-07   4.521226e-06
    ## 41  3.949209e-24  3.482485e-23   3.830733e-22
    ## 43  1.244879e-20  9.288714e-20   1.207533e-18
    ## 45  7.403862e-08  3.264430e-07   7.181746e-06
    ## 46  3.530313e-07  1.426835e-06   3.424404e-05
    ## 55  1.341788e-09  6.850179e-09   1.301534e-07
    ## 57  7.332310e-28  7.112341e-27   7.112341e-26
    ## 59  3.338711e-06  1.295420e-05   3.238549e-04
    ## 60 1.754950e-120 5.674337e-119  1.702301e-118
    ## 61  2.093704e-92  3.384821e-91   2.030893e-90
    ## 62  3.424592e-97  6.643709e-96   3.321854e-95
    ## 78  4.356859e-19  3.018681e-18   4.226153e-17
    ## 79 3.793034e-103 9.198107e-102  3.679243e-101
    ## 83  0.000000e+00  0.000000e+00   0.000000e+00
    ## 84  1.617813e-10  8.718216e-10   1.569279e-08

Test go terms mapped to GO slim - Based on BLAST with id &gt; 50%, cov subjet &gt; 50% and eval 1e-10
-----------------------------------------------------------------------------------------------------

``` r
q903_50id <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/Q903/ThresholdCovIdentity4/my_Slimgo_terms_plant.Q903IdLen.counts", header=FALSE)
ws_50id <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/WS77111/ThresholdCovIdentity4/my_Slimgo_terms_plant.WS77111IdLen.counts", header=FALSE)

##bind the 2 classes
q903_50id_wsSlimAll = merge(q903_50id,ws_50id,by="V2")
q903_50id_wsSlimAll = q903_50id_wsSlimAll[,c(1,2,5,3,4)]
colnames(q903_50id_wsSlimAll) = c("Go_term","Q903","WS77111","Domain","Go_name")
dim(q903_50id_wsSlimAll)
```

    ## [1] 97  5

``` r
#############All domains together
colSums_q903_50id_wsSlimAll = unname(apply(q903_50id_wsSlimAll[,c(2,3)],2,sum))
GrandTotal = sum(colSums_q903_50id_wsSlimAll)
p1 =colSums_q903_50id_wsSlimAll / GrandTotal

p1
```

    ## [1] 0.4894671 0.5105329

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

sum(p.adjust(prob_slim50id, method = "bonferroni", n = length(prob_slim50id)) < 0.05)
```

    ## [1] 32

``` r
sum(p.adjust(prob_slim50id, method = "bonferroni", n = length(prob_slim50id)) < 0.001)
```

    ## [1] 26

``` r
#q903_95_wsSlimAllSubs$prob_slim95 = prob_slim95
#q903_95_wsSlimAllSubs$prob_slim95BH = p.adjust(prob_slim95, method = "BH", n = length(prob_slim95))
#subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95BH < 0.05) 
#subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95BH < 0.001)

q903_50id_wsSlimAllSubs$prob_slimBon = p.adjust(prob_slim50id, method = "bonferroni", n = length(prob_slim50id))
subset(q903_50id_wsSlimAllSubs,q903_50id_wsSlimAllSubs$prob_slimBon < 0.05)
```

    ##       Go_term  Q903 WS77111
    ## 2  GO:0000166 38253   35364
    ## 8  GO:0003723 14651   16404
    ## 10 GO:0003824 61144   68526
    ## 14 GO:0005215 20032   19683
    ## 15 GO:0005488 50755   60523
    ## 18 GO:0005576 11872   11798
    ## 21 GO:0005622 12933   12609
    ## 27 GO:0005737 17618   17661
    ## 33 GO:0005783  3920    3615
    ## 37 GO:0005856  5658    5162
    ## 39 GO:0005975 11995   13801
    ## 40 GO:0006091 16101   18599
    ## 41 GO:0006139 26492   24316
    ## 42 GO:0006259   942     806
    ## 45 GO:0006629  7535    7149
    ## 46 GO:0006810 42215   39116
    ## 47 GO:0006950 22372   22296
    ## 48 GO:0007049  2485    2132
    ## 55 GO:0008152 31378   38260
    ## 57 GO:0008289  1151     955
    ## 59 GO:0009058 41953   39287
    ## 60 GO:0009507 49671   55272
    ## 61 GO:0009536 30240   34665
    ## 62 GO:0009579 45605   49848
    ## 68 GO:0009719  9978    9829
    ## 76 GO:0009987 48804   49108
    ## 78 GO:0015979 24516   30268
    ## 82 GO:0016301 15032   16371
    ## 83 GO:0016740 29956   25347
    ## 84 GO:0016787 36583   34809
    ## 85 GO:0019538  4753    4549
    ## 87 GO:0019748  4875    4504
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 14                             transporter_activity molecular_function
    ## 15                                          binding molecular_function
    ## 18                             extracellular_region cellular_component
    ## 21                                    intracellular cellular_component
    ## 27                                        cytoplasm cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 39                   carbohydrate_metabolic_process biological_process
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 42                            DNA_metabolic_process biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 47                               response_to_stress biological_process
    ## 48                                       cell_cycle biological_process
    ## 55                                metabolic_process biological_process
    ## 57                                    lipid_binding molecular_function
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 68                  response_to_endogenous_stimulus biological_process
    ## 76                                 cellular_process biological_process
    ## 78                                   photosynthesis biological_process
    ## 82                                  kinase_activity molecular_function
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ## 85                        protein_metabolic_process biological_process
    ## 87                      secondary_metabolic_process biological_process
    ##     prob_slimBon
    ## 2   3.187149e-58
    ## 8   4.337249e-08
    ## 10  3.495518e-36
    ## 14  2.589735e-07
    ## 15 8.807416e-108
    ## 18  1.911078e-02
    ## 21  6.638839e-06
    ## 27  1.867445e-02
    ## 33  8.846129e-06
    ## 37  3.273811e-10
    ## 39  3.641875e-13
    ## 40  2.285879e-19
    ## 41  4.651492e-45
    ## 42  3.449908e-03
    ## 45  9.216586e-07
    ## 46  6.332657e-62
    ## 47  1.442958e-04
    ## 48  3.301114e-09
    ## 55  1.261569e-91
    ## 57  1.567126e-05
    ## 59  2.887999e-51
    ## 60  1.177805e-23
    ## 61  3.252148e-31
    ## 62  4.803500e-11
    ## 68  5.541700e-03
    ## 76  1.837132e-06
    ## 78  5.764738e-84
    ## 82  1.274471e-02
    ## 83 3.415956e-131
    ## 84  1.261051e-32
    ## 85  3.255973e-03
    ## 87  4.169287e-07

``` r
subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95Bon < 0.001)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166  25636   20233
    ## 3  GO:0003674   1480    1510
    ## 8  GO:0003723   3327    3844
    ## 10 GO:0003824  79250   72020
    ## 13 GO:0005198   2661    2620
    ## 15 GO:0005488  65295   58154
    ## 21 GO:0005622  16919   12843
    ## 28 GO:0005739   4283    3263
    ## 33 GO:0005783   7463    7101
    ## 36 GO:0005840   2682    4049
    ## 37 GO:0005856   3873    2220
    ## 40 GO:0006091   5326    5020
    ## 41 GO:0006139  14360   10681
    ## 43 GO:0006412   3416    3610
    ## 45 GO:0006629   9502    8712
    ## 46 GO:0006810  22449   19969
    ## 55 GO:0008152  39175   34679
    ## 57 GO:0008289    621     921
    ## 59 GO:0009058  25114   20348
    ## 60 GO:0009507  13065   14634
    ## 61 GO:0009536   7446    8686
    ## 62 GO:0009579  11393   12628
    ## 78 GO:0015979   7206    7082
    ## 79 GO:0016020  64212   61390
    ## 83 GO:0016740 149762  102548
    ## 84 GO:0016787  28060   22432
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 15                                          binding molecular_function
    ## 21                                    intracellular cellular_component
    ## 28                                    mitochondrion cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 55                                metabolic_process biological_process
    ## 57                                    lipid_binding molecular_function
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 78                                   photosynthesis biological_process
    ## 79                                         membrane cellular_component
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ##      prob_slim95 prob_slim95BH prob_slim95Bon
    ## 2   9.038790e-14  5.479767e-13   8.767626e-12
    ## 3   3.209772e-07  1.353687e-06   3.113479e-05
    ## 8   1.021985e-39  1.101472e-38   9.913252e-38
    ## 10  3.285246e-43  3.983361e-42   3.186689e-41
    ## 13  3.923563e-08  1.902928e-07   3.805856e-06
    ## 15  5.266688e-19  3.405791e-18   5.108687e-17
    ## 21  1.146308e-20  9.265993e-20   1.111919e-18
    ## 28  5.669397e-06  2.115121e-05   5.499315e-04
    ## 33  1.735162e-12  9.900629e-12   1.683107e-10
    ## 36 9.478089e-123 4.596873e-121  9.193746e-121
    ## 37  3.542167e-49  4.908432e-48   3.435902e-47
    ## 40  4.661058e-08  2.152965e-07   4.521226e-06
    ## 41  3.949209e-24  3.482485e-23   3.830733e-22
    ## 43  1.244879e-20  9.288714e-20   1.207533e-18
    ## 45  7.403862e-08  3.264430e-07   7.181746e-06
    ## 46  3.530313e-07  1.426835e-06   3.424404e-05
    ## 55  1.341788e-09  6.850179e-09   1.301534e-07
    ## 57  7.332310e-28  7.112341e-27   7.112341e-26
    ## 59  3.338711e-06  1.295420e-05   3.238549e-04
    ## 60 1.754950e-120 5.674337e-119  1.702301e-118
    ## 61  2.093704e-92  3.384821e-91   2.030893e-90
    ## 62  3.424592e-97  6.643709e-96   3.321854e-95
    ## 78  4.356859e-19  3.018681e-18   4.226153e-17
    ## 79 3.793034e-103 9.198107e-102  3.679243e-101
    ## 83  0.000000e+00  0.000000e+00   0.000000e+00
    ## 84  1.617813e-10  8.718216e-10   1.569279e-08

Test go terms mapped to GO slim - Based on BLAST with id &gt; 50%, eval 1e-10 and aa align &gt; 100
---------------------------------------------------------------------------------------------------

``` r
q903_50id <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/Q903/ThresholdCovIdentity5/my_Slimgo_terms_plant.Q903IdLen.counts", header=FALSE)
ws_50id <- read.delim("/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/BLAST/WS77111/ThresholdCovIdentity5/my_Slimgo_terms_plant.WS77111IdLen.counts", header=FALSE)

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

sum(p.adjust(prob_slim50id, method = "bonferroni", n = length(prob_slim50id)) < 0.05)
```

    ## [1] 31

``` r
sum(p.adjust(prob_slim50id, method = "bonferroni", n = length(prob_slim50id)) < 0.001)
```

    ## [1] 21

``` r
#q903_95_wsSlimAllSubs$prob_slim95 = prob_slim95
#q903_95_wsSlimAllSubs$prob_slim95BH = p.adjust(prob_slim95, method = "BH", n = length(prob_slim95))
#subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95BH < 0.05) 
#subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95BH < 0.001)

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
    ##     prob_slimBon
    ## 2   7.785301e-36
    ## 3   5.825379e-05
    ## 5   5.236244e-03
    ## 8   3.586978e-24
    ## 10  2.572958e-47
    ## 11  2.088282e-02
    ## 13  9.878521e-03
    ## 14  3.642934e-11
    ## 15 4.925998e-134
    ## 17  2.978019e-02
    ## 21  8.803132e-56
    ## 24  4.070824e-05
    ## 31  1.147963e-02
    ## 32  2.011896e-03
    ## 33  1.990136e-02
    ## 36  4.098019e-35
    ## 38  1.482451e-03
    ## 39  2.359374e-26
    ## 40  2.082134e-14
    ## 41 4.040269e-170
    ## 43  2.277583e-15
    ## 46 4.667986e-172
    ## 50  1.368168e-02
    ## 55  5.303061e-67
    ## 58  3.314183e-06
    ## 59  9.430493e-66
    ## 61  6.487069e-03
    ## 68  2.335987e-04
    ## 78  1.824436e-41
    ## 83  1.258165e-21
    ## 94  3.639052e-05

``` r
subset(q903_95_wsSlimAllSubs,q903_95_wsSlimAllSubs$prob_slim95Bon < 0.001)
```

    ##       Go_term   Q903 WS77111
    ## 2  GO:0000166  25636   20233
    ## 3  GO:0003674   1480    1510
    ## 8  GO:0003723   3327    3844
    ## 10 GO:0003824  79250   72020
    ## 13 GO:0005198   2661    2620
    ## 15 GO:0005488  65295   58154
    ## 21 GO:0005622  16919   12843
    ## 28 GO:0005739   4283    3263
    ## 33 GO:0005783   7463    7101
    ## 36 GO:0005840   2682    4049
    ## 37 GO:0005856   3873    2220
    ## 40 GO:0006091   5326    5020
    ## 41 GO:0006139  14360   10681
    ## 43 GO:0006412   3416    3610
    ## 45 GO:0006629   9502    8712
    ## 46 GO:0006810  22449   19969
    ## 55 GO:0008152  39175   34679
    ## 57 GO:0008289    621     921
    ## 59 GO:0009058  25114   20348
    ## 60 GO:0009507  13065   14634
    ## 61 GO:0009536   7446    8686
    ## 62 GO:0009579  11393   12628
    ## 78 GO:0015979   7206    7082
    ## 79 GO:0016020  64212   61390
    ## 83 GO:0016740 149762  102548
    ## 84 GO:0016787  28060   22432
    ##                                              Domain            Go_name
    ## 2                                nucleotide_binding molecular_function
    ## 3                                molecular_function molecular_function
    ## 8                                       RNA_binding molecular_function
    ## 10                               catalytic_activity molecular_function
    ## 13                     structural_molecule_activity molecular_function
    ## 15                                          binding molecular_function
    ## 21                                    intracellular cellular_component
    ## 28                                    mitochondrion cellular_component
    ## 33                            endoplasmic_reticulum cellular_component
    ## 36                                         ribosome cellular_component
    ## 37                                     cytoskeleton cellular_component
    ## 40   generation_of_precursor_metabolites_and_energy biological_process
    ## 41 nucleobase-containing_compound_metabolic_process biological_process
    ## 43                                      translation biological_process
    ## 45                          lipid_metabolic_process biological_process
    ## 46                                        transport biological_process
    ## 55                                metabolic_process biological_process
    ## 57                                    lipid_binding molecular_function
    ## 59                             biosynthetic_process biological_process
    ## 60                                      chloroplast cellular_component
    ## 61                                          plastid cellular_component
    ## 62                                        thylakoid cellular_component
    ## 78                                   photosynthesis biological_process
    ## 79                                         membrane cellular_component
    ## 83                             transferase_activity molecular_function
    ## 84                               hydrolase_activity molecular_function
    ##      prob_slim95 prob_slim95BH prob_slim95Bon
    ## 2   9.038790e-14  5.479767e-13   8.767626e-12
    ## 3   3.209772e-07  1.353687e-06   3.113479e-05
    ## 8   1.021985e-39  1.101472e-38   9.913252e-38
    ## 10  3.285246e-43  3.983361e-42   3.186689e-41
    ## 13  3.923563e-08  1.902928e-07   3.805856e-06
    ## 15  5.266688e-19  3.405791e-18   5.108687e-17
    ## 21  1.146308e-20  9.265993e-20   1.111919e-18
    ## 28  5.669397e-06  2.115121e-05   5.499315e-04
    ## 33  1.735162e-12  9.900629e-12   1.683107e-10
    ## 36 9.478089e-123 4.596873e-121  9.193746e-121
    ## 37  3.542167e-49  4.908432e-48   3.435902e-47
    ## 40  4.661058e-08  2.152965e-07   4.521226e-06
    ## 41  3.949209e-24  3.482485e-23   3.830733e-22
    ## 43  1.244879e-20  9.288714e-20   1.207533e-18
    ## 45  7.403862e-08  3.264430e-07   7.181746e-06
    ## 46  3.530313e-07  1.426835e-06   3.424404e-05
    ## 55  1.341788e-09  6.850179e-09   1.301534e-07
    ## 57  7.332310e-28  7.112341e-27   7.112341e-26
    ## 59  3.338711e-06  1.295420e-05   3.238549e-04
    ## 60 1.754950e-120 5.674337e-119  1.702301e-118
    ## 61  2.093704e-92  3.384821e-91   2.030893e-90
    ## 62  3.424592e-97  6.643709e-96   3.321854e-95
    ## 78  4.356859e-19  3.018681e-18   4.226153e-17
    ## 79 3.793034e-103 9.198107e-102  3.679243e-101
    ## 83  0.000000e+00  0.000000e+00   0.000000e+00
    ## 84  1.617813e-10  8.718216e-10   1.569279e-08
