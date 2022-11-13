#!/bin/bash

gctb=/home/yx322/GCTB/gctb_2.02_Linux/gctb

PWD=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix

out=interval_shrunk_chr${SLURM_ARRAY_TASK_ID}



$gctb --mldm ${PWD}/${out}.mldmlist --make-shrunk-ldm --out ${PWD}/${out}