#!/bin/bash

if [[ $1 = "" ]]; then 
  out_dir=geno_files/ml_inputs
else
  out_dir=$1
fi

# Get array of phenotypes
phenotypes=($(ls $out_dir/*_variant_effects.txt | sed "s#$out_dir/##" | sed 's/_variant_effects.txt//'))
nphen=${#phenotypes[@]}

# Determine which phenotypes we're working with for this task:
task=$SLURM_ARRAY_TASK_ID
ntasks=$SLURM_ARRAY_TASK_MAX

nphen_per_task=$(echo "a=$nphen; b=$ntasks; if ( a%b )  a/b+1 else a/b" | bc)
task_start=$(echo "a=$task; b=$nphen_per_task; (a-1)*b" | bc)
task_end=$(echo "a=$task_start; b=$nphen_per_task; c=$nphen; if ( (a+b)>c ) c - 1 else a+b-1" | bc)

if [ $task_start -gt $task_end ]; then
  echo "No remaining phenotypes for task $task of $ntasks."
  exit 0
else 
  echo "Task $task of $ntasks extracting dosages of phenotypes $task_start - $task_end of $nphen."
fi

# Iterate through each protein for this task to extract the genotype dosages
for phenIdx in $(seq $task_start $task_end); do
  # What phenotype are we working with?
  phen=${phenotypes[$phenIdx]}

  # Determine the chromosomes we need to extract for this phenotype
  chrs=( $(tail -n +2 $out_dir/${phen}_variant_effects.txt | cut -f 2 | sort | uniq) )

  # Determine the samples we need to keep for this phenotype
  grep -P "^${phen}\t" $out_dir/phenotypes.txt | cut -f 3 > $out_dir/${phen}_IID.txt
  echo "#FID"$'\t'"IID" > $out_dir/${phen}.samples
  paste $out_dir/${phen}_IID.txt $out_dir/${phen}_IID.txt >> $out_dir/${phen}.samples
  rm $out_dir/${phen}_IID.txt
  
  # For each chromosome, extract the dosages of the pQTLs
  for chr in ${chrs[@]}; do
    # Get the list of pQTLs on this chromosome to extract
    cut -f 1,2 $out_dir/${phen}_variant_effects.txt | grep -w "${chr}"'$' | cut -f 1 > $out_dir/${phen}_chr${chr}_variant_ids.txt

    # Get their effect alleles
    grep -f $out_dir/${phen}_chr${chr}_variant_ids.txt -w $out_dir/${phen}_variant_effects.txt | cut -f 1,4 > $out_dir/${phen}_chr${chr}_variant_effect_alleles.txt 

    # Extract the dosages of the effect alleles
    plink2 --pfile geno_files/genotype_data/ldthinned/impute_${chr}_interval_dedup_unambig_SNPs_maf0.005_ldthin0.8 \
           --out $out_dir/${phen}_chr${chr}_dosages \
           --keep $out_dir/${phen}.samples \
           --extract $out_dir/${phen}_chr${chr}_variant_ids.txt \
           --export A --export-allele $out_dir/${phen}_chr${chr}_variant_effect_alleles.txt \
           --memory $SLURM_MEM_PER_NODE \
           --threads $SLURM_CPUS_ON_NODE \
           --silent

    # remove temporary files
    rm $out_dir/${phen}_chr${chr}_variant_ids.txt 
    rm $out_dir/${phen}_chr${chr}_variant_effect_alleles.txt
    rm $out_dir/${phen}_chr${chr}_dosages.log
  done

  # Combine chromosome data
  for chr in ${chrs[@]}; do
    paste $out_dir/${phen}_dosages.txt $out_dir/${phen}_chr${chr}_dosages.raw > $out_dir/${phen}_dosage_tmpfile
    mv $out_dir/${phen}_dosage_tmpfile $out_dir/${phen}_dosages.txt
    rm $out_dir/${phen}_chr${chr}_dosages.raw
  done

  # get the varID column and remove extra sample identifier information from each chromosome
  Rscript scripts/04_extract_QTLs/02_helpers/reformat_dosages.R $out_dir $phen
  rm $out_dir/${phen}_dosages.txt

  # sample file no longer needed
  rm $out_dir/${phen}.samples
done

