#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p commons,shared,serial_requeue   # Partition to submit to
#SBATCH --mem=20G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../_slurm/downloadnovogene_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../../_slurm/downloadnovogene_%j.err  # File to which STDERR will be written, %j inserts jobid

# RCC Samples

cd /n/irizarryfs01/kkorthauer/cfMeDIPseq/data

wget -r --user=P202SC18123014-01_20190326_YTSG0A --password=IYQ8dx ftp://128.120.88.242/