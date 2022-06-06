#!/bin/bash -l

#SBATCH --job-name=makerWS1
#SBATCH --partition=all
#SBATCH --ntasks=24
#SBATCH --mem=350G

export PATH=/gsc/btl/linuxbrew/bin:$PATH

nam=_RePlaCeNam_

maker -fix_nucleotides -base $nam maker_opts.ctl maker_bopts.ctl maker_exe.ctl

