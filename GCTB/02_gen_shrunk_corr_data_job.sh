#!/bin/bash

# Make sure job is submitted directly
if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch."
   exit 1
fi

# allows user to specify a job to wait for completion before running any of these scripts
if [ ! -z "$1" ]; then
  previous_job=$1
else
  previous_job=1 # run first job immediately
fi

# Create logging directory
log_dir=logs/02_gen_shrunk_corr/
mkdir -p $log_dir


gen_scorr_job=$(sbatch --dependency afterany:$previous_job \
                  --parsable \
                  --account INOUYE-COVID19-SL2-CPU \
                  --job-name "gen_scorr" \
                  --array 1-22 \
                  --time 12:0:0 \
                  --mem 60000 \
                  --output $log_dir/gen_scorr_%A_%a.o \
                  --error $log_dir/gen_scorr_%A_%a.e \
                  --partition skylake-himem \
                  02_gen_shrunk_corr_data/gctb_gen_shrunk_corr_matrix.sh)

echo "Submitted jobs $gen_scorr_job"

