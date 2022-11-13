#!/bin/bash


gctb=/home/yx322/GCTB/gctb_2.02_Linux/gctb

plink_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_INTERVAL_genetics/filtered_interval_chr${SLURM_ARRAY_TASK_ID}

genetic_map_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_genetic_map/genetic_map_chr${SLURM_ARRAY_TASK_ID}.txt

output_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix/interval_shrunk_chr${SLURM_ARRAY_TASK_ID}

$gctb --bfile $plink_file \
      --make-shrunk-ldm  \
      --gen-map $genetic_map_file \
      --out $output_file