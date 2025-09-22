#! /bin/sh -l
#SBATCH -J gwas-bmi
#SBATCH --partition=highmem
#SBATCH -N 1
#SBATCH -n 16

for i in {1..22}; do

plink2 --pfile /data/lfahey/ukb_qc/qcd_files2/${i}_qcd --rm-dup force-first --keep /data/aforde/bmi1A.txt --pheno /data/aforde/bmi1A.txt --make-bed --out ${i}_qcd_bmi
plink2 --bfile ${i}_qcd_bmi --extract /data/aforde/pruned/${i}_bmi.prune.in --glm omit-ref hide-covar --covar covariates.txt --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 AGE --out bmi_res_${i}

rm ${i}_qcd_bmi*;

done


Rscript summary-stats-BMI.R bmi_res_1.PHENO1.glm.linear bmi_res_2.PHENO1.glm.linear bmi_res_3.PHENO1.glm.linear bmi_res_4.PHENO1.glm.linear bmi_res_5.PHENO1.glm.linear bmi_res_6.PHENO1.glm.linear bmi_res_7.PHENO1.glm.linear bmi_res_8.PHENO1.glm.linear bmi_res_9.PHENO1.glm.linear bmi_res_10.PHENO1.glm.linear bmi_res_11.PHENO1.glm.linear bmi_res_12.PHENO1.glm.linear bmi_res_13.PHENO1.glm.linear bmi_res_14.PHENO1.glm.linear bmi_res_15.PHENO1.glm.linear bmi_res_16.PHENO1.glm.linear bmi_res_17.PHENO1.glm.linear bmi_res_18.PHENO1.glm.linear bmi_res_19.PHENO1.glm.linear bmi_res_20.PHENO1.glm.linear bmi_res_21.PHENO1.glm.linear bmi_res_22.PHENO1.glm.linear summary-stats-BMI1A.txt
