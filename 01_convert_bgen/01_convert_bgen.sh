#!/bin/bash

ref_dir=$HOME/rds/rds-jmmh2-post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
out_dir=geno_files/genotype_data/
mkdir -p $out_dir

chr=$SLURM_ARRAY_TASK_ID

plink2 --bgen $ref_dir/impute_${chr}_interval.bgen \
       --sample $ref_dir/interval.samples \
       --threads $SLURM_CPUS_ON_NODE \
       --memory $SLURM_MEM_PER_NODE \
       --silent \
       --out $out_dir/impute_${chr}_interval

