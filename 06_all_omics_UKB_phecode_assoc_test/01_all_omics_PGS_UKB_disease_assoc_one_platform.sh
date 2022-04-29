#!/bin/bash

module load miniconda3
source activate ml

python /rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/scripts_gene_expressions/08_all_omics_UKB_phecode_assoc_test/01_all_omics_PGS_UKB_disease_assoc.py ${1} ${2}
