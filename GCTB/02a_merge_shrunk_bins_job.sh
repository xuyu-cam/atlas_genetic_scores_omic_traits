#!/bin/bash

# Make sure job is submitted directly
if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch."
   exit 1
fi


# Create logging directory
log_dir=logs/02_gen_shrunk_corr/
mkdir -p $log_dir


merge_job=$(sbatch  \
                  --parsable \
                  --account INOUYE-COVID19-SL2-CPU \
                  --job-name "merge_bin" \
                  --array 1-6 \
                  --time 2:0:0 \
                  --mem 60000 \
                  --output $log_dir/merge_shrunk_%A_%a.o \
                  --error $log_dir/merge_shrunk_%A_%a.e \
                  --partition skylake,skylake-himem \
                  02a_gen_shrunk_corr_data_mcpu/gctb_merge_shrunk_bins_by_chr.sh)

echo "Submitted jobs $merge_job"

