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
                  --job-name "SbayesS" \
                  --time 8:0:0 \
                  --mem 60000 \
                  --output $log_dir/run_SBayesS_${1}_%A_%a.o \
                  --error $log_dir/run_SBayesS_${1}_%A_%a.e \
                  --partition skylake \
                  06_run_SbayesS/run_SbayesS.sh ${1})

echo "Submitted jobs $run_SBayesS"

