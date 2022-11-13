#!/bin/bash


gctb=/home/yx322/GCTB/gctb_2.02_Linux/gctb

shrunk_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix/interval_shrunk_chr${SLURM_ARRAY_TASK_ID}.ldm.shrunk

output_file=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix/interval_shrunk_chr${SLURM_ARRAY_TASK_ID}

$gctb --ldm $shrunk_file \
      --make-sparse-ldm  \
      --chisq 0  \
      --out $output_file


