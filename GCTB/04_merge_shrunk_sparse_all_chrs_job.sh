#!/bin/bash

# Make sure job is submitted directly
if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch."
   exit 1
fi


# Create logging directory
log_dir=logs/04_merge_all_chrs/
mkdir -p $log_dir


merge_all_job=$(sbatch  \
                  --parsable \
                  --account INOUYE-COVID19-SL2-CPU \
                  --job-name "all_merge" \
                  --time 2:0:0 \
                  --mem 60000 \
                  --output $log_dir/M_all_chrs_%A_%a.o \
                  --error $log_dir/M_all_chrs_%A_%a.e \
                  --partition skylake-himem \
                  04_merge_shrunk_sparse_all_chrs/gctb_merge_shrunk_sparse_bins_all_chrs.sh)

echo "Submitted jobs $merge_all_job"

