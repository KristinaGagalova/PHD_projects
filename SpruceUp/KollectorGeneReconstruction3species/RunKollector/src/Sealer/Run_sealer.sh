#!/bin/bash
#SBATCH --job-name=ws_sealer
#SBATCH --partition=all
#SBATCH --ntasks=48
#SBATCH --mem=350G
#SBATCH --exclusive

export PATH=/gsc/btl/linuxbrew/Cellar/abyss/2.0.1-k256/bin/:$PATH

export paths_dir=/projects/spruceup/pglauca/WS77111/data/reads

export reads=$(cat \
        $paths_dir/chromium-reads-split.in \
)

abyss-sealer -v -S $1 \
        -t $1-sealed-trace.txt \
        -o $1-sealed \
        -L 260 \
        -j 48 \
        -P 10  \
        -k135 --input-bloom=<(zcat /projects/spruceup/pglauca/WS77111/bloom/chromium-only/k135.bloom.all.z) \
        -k125 --input-bloom=<(zcat /projects/spruceup/pglauca/WS77111/bloom/chromium-only/k125.bloom.all.z) \
        -k115 --input-bloom=<(zcat /projects/spruceup/pglauca/WS77111/bloom/chromium-only/k115.bloom.all.z) \
        -k105 --input-bloom=<(zcat /projects/spruceup/pglauca/WS77111/bloom/chromium-only/k105.bloom.all.z) \
        -k95 --input-bloom=<(zcat /projects/spruceup/pglauca/WS77111/bloom/chromium-only/k95.bloom.all.z) \
        -k85 --input-bloom=<(zcat /projects/spruceup/pglauca/WS77111/bloom/chromium-only/k85.bloom.all.z) \
        -k75 --input-bloom=<(zcat /projects/spruceup/pglauca/WS77111/bloom/chromium-only/k75.bloom.all.z) \

### EOF ###

