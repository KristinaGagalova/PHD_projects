N_TIMES=$1
NUM_SEQ=18243
SEQ=$2

SEQ_NEW=${SEQ%.*}

for ((i=1; i<=N_TIMES; i++)); 
do 	
	echo $SEQ_NEW
	Rscript GenRandom.R $SEQ "${SEQ_NEW}_${i}ran.fa" $NUM_SEQ  
done

