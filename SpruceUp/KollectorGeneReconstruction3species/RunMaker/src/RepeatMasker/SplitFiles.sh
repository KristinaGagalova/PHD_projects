#!/usr/bin/env bash

#UCSC faSplit needed:http://hgdownload.cse.ucsc.edu/admin/jksrc.zip

fasta=/projects/spruceup_scratch/pglauca/WS77111/assemblies/comparative_genomics/kollector/kollector_runs/kollector_cdhit4AllTargets/FastaFiles/FilteredFasta/WS77111_allIterations_allbinsShortNams.fasta

faSplit sequence $fasta 21 WS77111batch
