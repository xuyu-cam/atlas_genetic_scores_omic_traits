#!/bin/bash

# Make sure job is submitted directly
if [ ! -z ${SLURM_JOB_ID+x} ]; then
   echo "This script should be executed directly, not with sbatch."
   exit 1
fi


# Create logging directory
log_dir=logs/02_gen_shrunk_corr/
mkdir -p $log_dir


gen_scorr_job=$(sbatch  \
                  --parsable \
                  --account INOUYE-COVID19-SL2-CPU \
                  --job-name "gen_scorr" \
                  --array 1-${2} \
                  --time 2:0:0 \
                  --mem 11000 \
                  --output $log_dir/gen_scorr_chr${1}_%A_%a.o \
                  --error $log_dir/gen_scorr_chr${1}_%A_%a.e \
                  --partition skylake,skylake-himem \
                  02a_gen_shrunk_corr_data_mcpu/gctb_gen_shrunk_corr_matrix_mcpu.sh ${1})

echo "Submitted jobs $gen_scorr_job"

