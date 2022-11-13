#!/bin/sh

plink=/home/yx322/plink_2.0/plink2

bed_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/interval_genetics/interval_impute_chr${SLURM_ARRAY_TASK_ID}

output_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_INTERVAL_genetics/filtered_interval_chr${SLURM_ARRAY_TASK_ID}

variant_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/OmicsPred_LDpred_benchmarks/data/HapMap3_transformed/HapMap1kg_variants_matched2INTERVAL_rsid.txt

sample_id_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/OmicsPred_LDpred_benchmarks/data/INTERVAL_genotypes/all_externalIDs.txt


$plink  --bfile $bed_file  \
        --extract $variant_file \
        --geno 0.05  \
        --hwe 1e-6 \
        --mach-r2-filter 0.3 \
        --maf 0.01 \
        --remove $sample_id_file \
        --make-bed \
        --out $output_file
        
