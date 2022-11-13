#!/bin/bash


PWD=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix

out=interval_shrunk_chr


for i in $( seq 1 22 )
do

echo "${PWD}/${out}${i}.ldm.shrunk" >> "${PWD}/interval_shrunk_only_chr_all.mldmlist"

done