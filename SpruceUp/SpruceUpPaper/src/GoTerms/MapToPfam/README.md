## Commands to use the scripts in the directory

### Select the Longest isoform only from Interproscan
```
grep -F -f  list_TransFinal.in /projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.maker.proteinsLoweAEDPolished.rename.tsv > Q903.all.maker.proteinsLoweAEDPolishedFinal.tsv
```

### Transcript ID and GO terms
```
awk -F $'\t' '{print $1,$14}' Q903.all.maker.proteinsLoweAEDPolishedFinal.tsv | tr ' ' '\t' | awk  -F $'\t' '$2!=""' > Q903.all.maker.proteinsLoweAEDPolishedFinal_go.tsv
```

### Reshape
```
bash ../Reshape_goAllEntities.sh Q903.all.maker.proteinsLoweAEDPolishedFinal_go.tsv > Q903.all.maker.proteinsLoweAEDPolishedFinal_goAll.tsv
bash ../Transform2GAF.sh Q903.all.maker.proteinsLoweAEDPolishedFinal_goAll.tsv E0M31 > Q903.all.maker.proteinsLoweAEDPolishedFinal_goAll.gaf
```

### Map to Go slim
```
bash ../MapToGoSplim.sh Q903.all.maker.proteinsLoweAEDPolishedFinal_goAll.gaf
```

### Assign GO terms
```
join -t $'\t' -1 1 -2 3 <(sort -k1,1 /projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/Dependencies/goslim_plant.nams) <(awk '{print $2,$3,$4}' my_Slimgo_terms_plant.gaf | tr ' ' '\t' | sort -k3,3) | tr ' ' '\t' > my_Slimgo_terms_plant.nams
```

### Count the occurrences
```
awk -F $'\t' '{print $1,$2,$3}' my_Slimgo_terms_plant.nams | sort | uniq -c | sed -e 's/^[ \t]*//' | tr ' ' '\t' > my_Slimgo_terms_plant.counts
````

