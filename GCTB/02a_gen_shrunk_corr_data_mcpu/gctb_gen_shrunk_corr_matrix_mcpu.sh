#!/bin/bash

k=5000
i=${SLURM_ARRAY_TASK_ID}

gctb=/home/yx322/GCTB/gctb_2.02_Linux/gctb

plink_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_INTERVAL_genetics/filtered_interval_chr${1}

genetic_map_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_genetic_map/genetic_map_chr${1}.txt

output_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix/interval_shrunk_chr${1}

$gctb --bfile $plink_file \
      --make-shrunk-ldm  \
      --gen-map $genetic_map_file \
      --snp $((k*(i-1)+1))-$((k*i)) \
      --out ${output_file}