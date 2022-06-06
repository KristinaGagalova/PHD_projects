# The script performs ParseBLAST with the selected thresholds for the files in specified directories

ident=(80 85 90)
cov=(25 50 75)

base_dir=/projects/spruceup/pglauca/WS77111/assemblies/kollector/target-sequences/evaluation/alignments/usr/kgagalova/SameSeq/ParseBlast

cd $base_dir
mkdir -p output_v2

for f in $base_dir/*.tsv

do
	echo $f
	for ((i=0; i<=${#ident[@]} - 1; i++)); do echo "-i: ${ident[$i]}, -c: ${cov[$i]}"; python ../ParseBLASTv2.py -f $f -o $base_dir/output_v2 -i ${ident[$i]} -c ${cov[$i]} -t count; done   
done
