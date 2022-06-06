## Samtools

Index multiple files at the same time with xargs

```
#P5 - run 5 jobs at the same time
find `pwd` -name "*bam" | xargs -n1 -P5 samtools index
```
