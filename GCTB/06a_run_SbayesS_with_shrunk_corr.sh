#!/bin/bash

# Make sure job is submitted directly
if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch."
   exit 1
fi


# Create logging directory
log_dir=logs/06_run_SbayesS/
mkdir -p $log_dir


run_SBayesS=$(sbatch  \
                  --parsable \
                  --account INOUYE-SL2-CPU \
                  --job-name "SS${1}" \
                  --time 12:0:0 \
                  --mem 240000 \
                  --output $log_dir/run_SBayesS_shrunk_${1}_%A_%a.o \
                  --error $log_dir/run_SBayesS_shrunk_${1}_%A_%a.e \
                  --partition skylake-himem \
                  06a_run_SbayesS_with_shrunk_corr/run_SbayesS_with_shrunk.sh ${1})

echo "Submitted jobs $run_SBayesS"

