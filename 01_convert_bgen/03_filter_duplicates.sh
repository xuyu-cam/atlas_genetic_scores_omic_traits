#!/bin/bash

out_dir=geno_files/genotype_data/
chr=$SLURM_ARRAY_TASK_ID

# Identify variants to remove
grep 'remove' $out_dir/impute_${chr}_interval.pvar | cut -f 3 > $out_dir/chr${chr}_duplicates.txt

# Exclude these and create new pgen files
plink2 --pfile $out_dir/impute_${chr}_interval \
       --exclude $out_dir/chr${chr}_duplicates.txt \
       --threads $SLURM_CPUS_ON_NODE \
       --memory $SLURM_MEM_PER_NODE \
       --silent \
       --make-pgen \
       --out $out_dir/impute_${chr}_interval_dedup

# Remove old pgen files.
rm $out_dir/impute_${chr}_interval.pgen
rm $out_dir/impute_${chr}_interval.psam
rm $out_dir/impute_${chr}_interval.pvar

# Remove temporary exclusion file
rm $out_dir/chr${chr}_duplicates.txt

