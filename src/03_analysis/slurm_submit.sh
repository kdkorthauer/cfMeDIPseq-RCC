# Main slurm controller workflow

# initial MEDIPS run
sbatch slurm_run.sh


# iterate DMR detection
ntop=300
export ntop
WINDOWSIZE=300
export WINDOWSIZE
iter=0
export iter
sbatch slurm_one_iteration.sh

for iter in {1..100}; do
  export iter
  sbatch slurm_one_iteration.sh
  sleep 0.5
done

# summarize and make figures
summary=1
echo summary
sbatch slurm_run.sh


