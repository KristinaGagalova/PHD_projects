# Pipeline for multeplicity evaluation

**input**: fa file from cdhit
**out**: table counts for different identity and coverage thresholds

## Pipeline steps
- Run BLAST all-vs-all
- Filter for self-hits and reverse complement
- Count the hits per query and output in tab format
