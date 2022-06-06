#!/usr/bin/bash


genomeWS=/projects/spruceup/pglauca/WS77111/assemblies/releases/version2/WS77111v2_release/HardMasked/WS77111-v2_1000plus_hardMasked.fa
genomeQ903=/projects/spruceup/psitchensis/Q903/assembly/releases/version1/HardMasked/Q903_v1_1000plus_hardMasked.fa

#Extract og and gene spruce
#find /projects/spruceup_scratch/dev/SprucePaper2018/GeneFamilies/OthofinderFinalRuns/AllConifers_noDF1/Results_Out_6Species_pine_spruce/Orthogroup_Sequences/ -type f | xargs grep "ABT39\|DB47\|E0M31" | tr ':>' '\t' | cut -d"/" -f11 | sed 's/.fa//' > OG_spruce.in


####################
#WS77111
####################
grep "DB47" ../OG_spruce.in | cut -d"-" -f1 > OG_spruceDB47.in

#get genes only
grep -F -f <(awk '{print $2}' OG_spruceDB47.in) <(awk '$3 == "gene" {print $0}' WS77111.all.polished.genesLoweAED.function_domainAnnotatedCompletenoFrag.gff) > WS77111.all.genesOG.gff


paste -d'\t' WS77111.all.genesOG.gff <(awk '{print $9}' WS77111.all.genesOG.gff | cut -d";" -f1 | sed "s/ID=//") > WS77111.all.genesOG2.gff

#sorting works only with the following command
LOCALE=C
join -t $'\t' -1 10 -2 2 <(sort -t $'\t' -k10,10 WS77111.all.genesOG2.gff) <(sort -k2,2 OG_spruceDB47.in) > WS77111.all.genesOG2sort.gff

awk -F $'\t' -v OFS='\t' '{print $2,$5 - 1,$6,$11,$1,$8}' WS77111.all.genesOG2sort.gff > WS77111.all.genesOG2sort.bed

bedtools getfasta -s -name -fi $genomeWS -bed WS77111.all.genesOG2sort.bed | cut -d"(" -f1 > WS77111.all.genesOGExtract.fa

####################
#Q903
####################

grep "E0M31" ../OG_spruce.in | cut -d"-" -f1 > OG_spruceE0M31.in

grep -F -f <(awk '{print $2}' OG_spruceE0M31.in) <(awk '$3 == "gene" {print $0}' Q903.all.polished.genesLoweAED.function_domainAnnotatedCompletenoFragNams.gff) > Q903.all.genesOG.gff

paste -d'\t' Q903.all.genesOG.gff <(awk '{print $9}' Q903.all.genesOG.gff | cut -d";" -f1 | sed "s/ID=//") > Q903.all.genesOG2.gff

#sorting works only with the following command
LOCALE=C
join -t $'\t' -1 10 -2 2 <(sort -t $'\t' -k10,10 Q903.all.genesOG2.gff) <(sort -k2,2 OG_spruceE0M31.in) >Q903.all.genesOG2sort.gff

awk -F $'\t' -v OFS='\t' '{print $2,$5 - 1,$6,$11,$1,$8}' Q903.all.genesOG2sort.gff > Q903.all.genesOG2sort.bed 

bedtools getfasta -s -name -fi $genomeQ903 -bed Q903.all.genesOG2sort.bed | cut -d"(" -f1 > Q903.all.genesOGExtract.fa

cat Q903.all.genesOGExtract.fa WS77111.all.genesOGExtract.fa > WS_Q903.all.genesOGExtract.fa

cat WS_Q903.all.genesOGExtract.fa | paste - - | sort -t $'\t' -k1 | datamash -g 1 collapse 2 | sed 's/,/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g'| tr "\t" "\n" > WS_Q903.all.genesOGExtractConc.fa

rm *genesOG2sort.gff *genesOG2.gff

tar -czvf IntermediateFiles.tar.gz *gff *bed *.all.genesOGExtract.fa OG*in

rm *gff *bed *.all.genesOGExtract.fa OG*in
