#!/bin/bash

gctb=/home/yx322/GCTB/gctb_2.03beta_Linux/gctb

PWD=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix


$gctb --sbayes S  \
      --ldm ${PWD}/interval_shrunk_chr_all.ldm.shrunk  \
      --gwas-summary /rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_GWAS_summ_stats/Metabolon/gwas_filteredCleaned_${1}.txt \
      --pi 0.01 \
      --hsq 0.5 \
    --num-chains 4 \
    --chain-length 25000 \
    --burn-in 2000 \
    --seed 12345 \
    --thread 18 \
    --no-mcmc-bin \
    --out-freq 10 \
    --out /rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/results/GCTB/${1}