#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p commons,shared,serial_requeue   # Partition to submit to
#SBATCH --mem=20G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../_slurm/downloaddropbox_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../../_slurm/downloaddropbox_%j.err  # File to which STDERR will be written, %j inserts jobid

# BLCA samples, healthy control

module load rclone
cd /n/irizarryfs01/kkorthauer/cfMeDIPseq/data

# cf MeDIP data from Sandor
rclone copy dropbox:cfDNA .
