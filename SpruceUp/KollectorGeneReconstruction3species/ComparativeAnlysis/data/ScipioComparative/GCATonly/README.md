## Output files from the comparative analysis - general and introns

- **ScipioParsedAllRaw.txt** - Scipio parsed output from yaml
- **ScipioParsedScore09Raw.txt** - Scipio parsed output from yaml, filtered for score >= 0.9
- **ScipioParsedScore09StatusRaw.txt** - Scipio parsed output from yaml, filtered for score >= 0.9 and kept status = auto (complete) and incomplete = mismatches
- **ScipioParsedScore09StatusMinIter.txt** - Scipio parsed output from yaml, filtered for score >= 0.9 and kept status = auto (complete) and incomplete = mismatches, if multiple hits per target keep only the ones that are assembled in earlier iteration.
- **CommonReconstructionTop.txt** - From the previous filtering steps, list of names that are commonly reconstructed in the species
- **NumbExons.txt** - From the previous filtering steps, number of exons per species including id and iteration
- **NumbDiffExons.txt** - From the previous	filtering steps, differential number of exons per species including id and iteration
