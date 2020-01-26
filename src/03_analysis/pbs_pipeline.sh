# Main slurm controller workflow

# initial MEDIPS run
WINDOWSIZE=300
export WINDOWSIZE
ntop=300
export ntop
summary=0
export summary
echo summary
qsub pbs_run.sh


# iterate DMR detection
iter=0
export iter
qsub pbs_one_iteration.sh

for iter in {1..75}; do
  export iter
  qsub pbs_one_iteration.sh
  sleep 0.5
done

# summarize and make figures
summary=1
echo summary
export summary
qsub pbs_run.sh


