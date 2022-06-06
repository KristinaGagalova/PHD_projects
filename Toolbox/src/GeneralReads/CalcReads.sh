################################
######Kristina Gagalova#########
################################
#########Apr 26 2017############
################################

#Description: calculate the number of reads for each fastq file in list. Output is a table of file and number of reads
#git repo: https://github.com/bcgsc/PHD_projects/blob/master/Toolbox/src/GeneralReads/CalcReads.sh

#one per line list of fastq files
READS_FILE=$1
nam=$2

for file in $(cat $READS_FILE); do
	echo $file >>  "$nam".files
	echo $file && gunzip -c $file | wc -l | awk '{d=$1; print d/4;}' >> "$nam".reads
done

paste -d'\t' "$nam".files "$nam".reads > "$nam".out
rm "$nam".files "$nam".reads

