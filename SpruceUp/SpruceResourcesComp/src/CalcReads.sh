READS_FILE=/projects/spruceup/pglauca/WS77111/data/reads/all-pet-filer-1.in

for file in $(cat $READS_FILE); do
	echo $file >>  Gsc_PET.files
	echo $file && gunzip -c $file | wc -l | awk '{d=$1; print d/4;}' >> Gsc_PET.reads
done

paste -d'\t' Gsc_PET.files Gsc_PET.reads > Gsc_PET.out
rm Gsc_PET.reads Gsc_PET.files
