#!/bin/sh
#SBATCH -J LD-prune
#SBATCH --partition=highmem
#SBATCH -N 1
#SBATCH -n 16

for i in {1..22}; do

plink2 --pfile /data/lfahey/ukb_qc/qcd_files2/${input}_qcd --rm-dup force-first --keep bmi4plink.txt --pheno bmi4plink.txt --make-bed --out /data/aforde/${input}_qcd_bmi
plink2 --bfile /data/aforde/${i}_qcd_bmi --indep-pairwise 50 5 0.5 --out ${i}_bmi;

done

