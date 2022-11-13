#!/bin/bash

# Make sure job is submitted directly
if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch."
   exit 1
fi


# Create logging directory
log_dir=logs/03_shrunk_sparse_corr/
mkdir -p $log_dir


sparse_job=$(sbatch  \
                  --parsable \
                  --account INOUYE-COVID19-SL2-CPU \
                  --job-name "sparse" \
                  --array 1-22 \
                  --time 2:0:0 \
                  --mem 60000 \
                  --output $log_dir/shrunk_sparse_%A_%a.o \
                  --error $log_dir/shrunk_sparse_%A_%a.e \
                  --partition skylake,skylake-himem \
                  03_convert_to_sparse/gctb_convert_to_sparse_shrunk_corr_matrix.sh)

echo "Submitted jobs $sparse_job"

