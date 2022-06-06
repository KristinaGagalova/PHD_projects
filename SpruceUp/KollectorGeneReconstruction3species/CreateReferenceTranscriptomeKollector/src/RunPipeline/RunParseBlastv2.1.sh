# The script performs ParseBLAST with the selected thresholds for the files in specified directories

in=$1
out=$2

ident=(80 85 90)
cov=(25 50 75)

for f in *$in
do
	echo $f
	for ((i=0; i<=${#ident[@]} - 1; i++)); do echo "-i: ${ident[$i]}, -c: ${cov[$i]}"; python ParseBLASTv2.1.py -f $f -o $out -i ${ident[$i]} -c ${cov[$i]} -t count; done   
done

