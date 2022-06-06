# The script performs ParseBLAST with the selected thresholds for the files in specified directories

ident=80
cov=25

base_dir=/projects/spruceup/pglauca/WS77111/assemblies/kollector/target-sequences/evaluation/alignments/usr/kgagalova/SameSeq/ParseBlast

cd $base_dir
mkdir -p output

for f in $base_dir/*.tsv
do
	python ../ParseBLAST.py -f $f -o $base_dir/output -i $ident -c $cov   
done
