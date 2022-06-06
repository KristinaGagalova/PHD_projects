#!/bin/bash

#SBATCH --job-name=q64_sealer
#SBATCH --partition=all
#SBATCH --ntasks=48
#SBATCH --mem=304G
#SBATCH --exclusive

export PATH=/gsc/btl/linuxbrew/Cellar/abyss/2.0.1-k256/bin/:$PATH

echo "PATH for job was:"
echo $PATH

k=80

cd /projects/spruceup/interior_spruce/PG29/bloom/hiseq_miseq/k80

list_reads=/projects/spruceup/interior_spruce/PG29/data/reads/fastq/PG29_reads.in

reads=$(cat $list_reads | tr "\t" " ")

echo "abyss-bloom build command:"
echo "abyss-bloom build -v -v -k$k -j48 -b300G -l2 -q15 - $reads |gzip -c>k$k.bloom.hiseq_miseq.z"
abyss-bloom build -v -v -k$k -j48 -b300G -l2 -q15 - $reads |gzip -c>k$k.bloom.hiseq_miseq.z

echo "Job ended at $(date)"

exit

