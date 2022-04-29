#!/bin/bash

ref_dir=geno_files/genotype_data
out_dir=$ref_dir/ldthinned
chr=$SLURM_ARRAY_TASK_ID

# Get LD thinned set of variants with MAF > 0.5%
plink2 --pfile $ref_dir/impute_${chr}_interval_dedup \
       --extract $out_dir/chr${chr}_keep.txt \
       --maf 0.005 \
       --indep-pairwise 1000kb 0.8 \
       --threads $SLURM_CPUS_ON_NODE \
       --memory $SLURM_MEM_PER_NODE \
       --silent \
       --out $out_dir/chr${chr}_ldthinned

# Extract those variants
plink2 --pfile $ref_dir/impute_${chr}_interval_dedup \
       --extract $out_dir/chr${chr}_ldthinned.prune.in \
       --threads $SLURM_CPUS_ON_NODE \
       --memory $SLURM_MEM_PER_NODE \
       --silent \
       --make-pgen \
       --out $out_dir/impute_${chr}_interval_dedup_unambig_SNPs_maf0.005_ldthin0.8

# Remove temporary file
rm $out_dir/chr${chr}_keep.txt
rm $out_dir/chr${chr}_ldthinned.prune.in
rm $out_dir/chr${chr}_ldthinned.prune.out

