#############################
####Kristina Gagalova########
#############################
######3 May 2017#############
#############################

#Description: given a list of contigs, check all the sam files in the directory and extract coverage, name target, contig, alignment, CIGAR and concatenate to ND tag

LIST_TARGEN=/extscratch/btl/kgagalova/EasterRun/ExtractResults/ExtractSuccessful/Out.results.PG29.top

while read line
do
target=$(cut -f2 -d' ' <<< $line)
gene=$(cut -f3 -d' ' <<< $line)
grep $target *sam | grep $gene | sed 's/.sam:/ /' | sort | uniq | awk '$2 >= 0.9 {print $0}' | grep -o 'NM:[^[:space:]]*\|MD:[^[:space:]]*' | tr '\n' '\t' | sed '$s/\t$/\n/' >> md_nm.tmp
grep $target *sam | grep $gene | sed 's/.sam:/ /' | sort | uniq | awk '$2 >= 0.9 {print $1,$2,$3,$4,$5,$8}' >> Out.results.tmp
done < $LIST_TARGEN

paste -d$'\t' Out.results.tmp md_nm.tmp > Out.results.PG29.top_sam.out
rm *tmp

