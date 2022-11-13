#!/bin/bash


k=5000
PWD=/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/methods_benchmark/data/GCTB_corr_matrix


out=interval_shrunk_chr${1}


for i in $( seq 1 ${2} )
do

echo "${PWD}/${out}.snp$((k*(i-1)+1))-$((k*i)).ldm.shrunk" >> "${PWD}/${out}.mldmlist"

done