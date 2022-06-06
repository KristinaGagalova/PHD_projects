## Extract, download metadata for the IDS

Download enzyme names
```
while read id; do echo -n $id && echo -n "," && efetch -db nucleotide -id $id -format gp | grep "product=" | cut -d"=" -f2 |sed 's/ //g;s/"//g'; done < GenBankIds.in > IdsEnzymeNames.txt
```

Download classification (some manual curation needed)
```
while read id; do echo -n $id && echo -n "," && efetch -db nucleotide -id $id -format gp | grep "note=" | cut -d"=" -f2 |sed 's/ //g;s/"//g'; done < GenBankIds.in > IdsEnzymeType.txt
```

```
cut -d"," -f2 IdsEnzymeType.txt | paste -d, IdsEnzymeNames.txt - > MetadataEnzymes.txt
```
