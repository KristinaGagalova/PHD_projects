## Contigs selection

Select the contigs that have at least 1000/500bp after masking and gaps


```
awk '{print ($2 - $3 - $4),$1 }' WS77111-v2_all.stats | awk '$1 >= 1000 {print $2}' > WS77111-v2_sel1000.in
awk '{print ($2 - $3 - $4),$1 }' WS77111-v2_all.stats | awk '$1 >= 500 {print $2}' > WS77111-v2_sel500.in

```
