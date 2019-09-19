#!/bin/bash
#SBATCH -J md_iter
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-8:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=38G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../_slurm/medip_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../../_slurm/medip_%j.err  # File to which STDERR will be written, %j inserts jobid

export WINDOWSIZE
export iter
export ntop

if [ $iter -eq 0 ]; then
  R CMD BATCH --quiet --no-restore --no-save 20190712_MEDIPS_RCC.R log/20190712_MEDIPS_RCC\_$WINDOWSIZE.Rout
else
  R CMD BATCH --quiet --no-restore --no-save 20190712_DMRs_RCC.R log/20190712_DMRs_RCC\_$WINDOWSIZE\_$iter.Rout
fi
