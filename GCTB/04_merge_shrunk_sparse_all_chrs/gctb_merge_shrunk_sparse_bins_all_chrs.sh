#!/bin/bash

gctb=/home/yx322/GCTB/gctb_2.02_Linux/gctb

PWD=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix


$gctb --mldm ${PWD}/interval_shrunk_chr_all.mldmlist --make-sparse-ldm --chisq 0 --out ${PWD}/interval_shrunk_chr_all