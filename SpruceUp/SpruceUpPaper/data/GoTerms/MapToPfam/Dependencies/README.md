## Downloaded on the 22nd January 2019

```
wget http://snapshot.geneontology.org/ontology/go-basic.obo

```

# go slims 
```
http://current.geneontology.org/ontology/subsets/goslim_plant.obo
http://current.geneontology.org/ontology/subsets/goslim_generic.obo
```

### The other files are reshaped based on the following blog
```
http://avrilomics.blogspot.com/2015/06/creating-go-slim-and-mapping-go-terms.html
```

```
grep "id: GO:" goslim_generic.obo | grep -v alt | cut -d" " -f2 > goslim_terms.txt
```

### Names and GO terms
```
wget http://www.geneontology.org/doc/GO.terms_alt_ids
```

### Add names to GO terms - manually remove the last rows because not contaning GO terms - 98 terms
```
paste -d '\t' <(grep "^id: " goslim_plant.obo | sed 's/id: //') <(grep "name: " goslim_plant.obo | sed 's/name: //')  <(grep "namespace: " goslim_plant.obo | sed 's/namespace: //')  > goslim_plant.nams 
44 biological_process
28 cellular_component
26 molecular_function
```
