## Simulated sample splitting approach to address selection bias in Mendelian Randomization studies

This repository contains the scripts used to provide results and figures for the manuscript, ***"Simulated sample splitting approach to address selection bias in Mendelian Randomization studies"***. In our work, we propose MR Simulated Sample Splitting (MR-SimSS), a novel method that corrects Winner's Curse bias in Mendelian Randomization (MR) studies using only GWAS summary data. MR-SimSS is also designed to remain robust under varying degrees of sample overlap between exposure and outcome datasets. The performance of MR-SimSS is evaluated first using a **simulation study** with a wide variety of simulated data sets and then, by means of two **empirical analyses** using real data sets of two different traits, body mass index (BMI) and type 2 diabetes (T2D).

Scripts related to simulations are contained in `simulations/scripts` while scripts used for the empirical analyses are contained in `real_data/scripts`.

&nbsp;

### Simulation study

**00-sims-funs.R:** This script provides all functions required for simulations.

**01-sims.R:** This script runs the main simulation study which evaluates
and compares four variants of MR-SimSS together with other classical MR approaches, across varying parameters such as exposure-outcome correlation, degree of sample overlap and sample size.

- *Output:* sims-res-200.csv, sims-res-500.csv, sims-res-50.csv

**02-sims-threshold.R:** This scripts run a small simulation study which assesses the
performance of SimSS-3-RAPS in low-power settings with different significance
thresholds applied in the selection step.

- *Output:* sims-res-thresh.csv

&nbsp;

### Empirical analyses

#### Quality control

The genotypic data used was collected, processed and imputed by UK Biobank (UKBB) ([http://www.ukbiobank.ac.uk/](http://www.ukbiobank.ac.uk/)). The required quality control steps were implemented using the code available at [https://github.com/coggene/UK-Biobank-QC](https://github.com/coggene/UK-Biobank-QC).

These steps included the removal of both related individuals and those that were of non-European ancestry. These non-European samples were identified by principal component analysis (PCA) using 1000 Genomes Project (1KGP) data. Furthermore, samples which had been identified as outliers with respect to heterozygosity and missingness, together with samples with discordant sex information and those suffering from chromosomal aneuploidy, were also discarded. The total number of samples remaining after the execution of these quality control steps were 332,618 and 333,642 for BMI and T2D, respectively. With respect to variants, only those with an information score greater than 0.8, a minor allele frequency greater than 0.01, a genotyping rate of at least 98% and those that passed the Hardy-Weinberg test at the specified significance threshold of 1 <span>&#215;</span> 10<sup>-8</sup> were included. This resulted in a total of 7,915,560 SNPs that were considered suitable for our analyses. 

This QC process provided us with the following files: bmi4plink.txt, T2D4plink.txt, covariates.txt and \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome.

&nbsp;

#### LD-pruning

The total set of 7,915,560 SNPs were **pruned** to provide a *sub set of approximately independent variants*. Pruning occurred by first calculating LD between each pair of SNPs in a window of 50 SNPs. If an LD value greater than 0.5 was observed, then one SNP out of this pair was removed. The window was shifted 5 SNPs forward and the process repeated. 1,589,295 genetic variants remained after this procedure, a data set about 20% of the size of the original.

The following code was run: `sbatch LD-prune.sh`

**LD-prune.sh:** This shell script performs pruning using PLINK 2.0 and `--indep-pairwise 50 5 0.5` to obtain a full list of pruned SNPs.

- *Input:* \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, bmi4plink.txt

- *Output:* \*_bmi.prune.in for each chromosome

&nbsp;

#### **Same-trait empirical analysis**

##### 1. Obtaining 10 pairs of independent samples

The following code was run: `Rscript 00-split-BMI.R bmi4plink.txt`

**00-split-BMI.R:** This script takes the full UK Biobank BMI data set and randomly splits it into two independent samples 10 times, providing 10 pairs of independent sample sets. 

- *Input:* bmi4plink.txt

- *Output:* bmi$\dagger$A.txt, bmi$\dagger$B.txt for $\dagger$ $\in$ {1,2,...,10}

&nbsp;

##### 2. Conducting GWASs for each set of samples

With respect to each sample set contained in bmi$\dagger$A.txt or bmi$\dagger$B.txt, the following code was run: `sbatch 01-GWAS-BMI.sh`

**01-GWAS-BMI.sh:** This script uses PLINK 2.0 to produce \*_qcd.bim, \*_qcd.bed, \*_qcd.fam files, with the phenotype data being specified in a given sample set, e.g. bmi1A.txt. A linear model is then fitted for every variant in which 8 principal components are included as well as age. Results from the files bmi_res_PHENO1.glm.linear for chromosomes 1 to 22 are then combined into a summarized text file, e.g. summary-stats-BMI1A.txt, which contains information regarding chromosome, position, rsID, effect size and corresponding standard error for each genetic variant.

- *Input:* \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, bmi$\dagger$A.txt or bmi$\dagger$B.txt for $\dagger$ $\in$ {1,2,...,10}, summary-stats-BMI.R

- *Output:* summary-stats-BMI$\dagger$A.txt or summary-stats-BMI$\dagger$B.txt for $\dagger$ $\in$ {1,2,...,10}


&nbsp;

##### 3. Performing MR analyses

**02-BMI-BMI.R:** This script is used to conduct same-trait MR analyses by applying four variants of MR-SimSS, together with other summary-level MR methods, to BMI GWAS summary statistics.

- *Input:* summary-stats-BMI$\dagger$A.txt, summary-stats-BMI$\dagger$B.txt for $\dagger$ $\in$ {1,2,...,10}
- *Output:* res-BMI-BMI.csv

&nbsp;


#### **Effect of BMI on T2D**

##### 1. Obtaining sample sets with varying overlap

The following code was run: `Rscript 00-sample-T2D-BMI.R T2D4plink.txt bmi4plink.txt`

**00-sample-T2D-BMI.R:** This script randomly produces UK Biobank sample sets, one for T2D and five for BMI, in which the BMI sample sets have varying degrees of sample overlap with the T2D sample set.  

- *Input:* T2D4plink.txt bmi4plink.txt

- *Output:* T2D.txt, bmi$\ddagger$.txt for $\ddagger$ $\in$ {0,25,50,75,100}

&nbsp;

##### 2. Conducting GWASs

The following code was run: `sbatch 01-GWAS-T2D-BMI.sh`

**01-GWAS-T2D-BMI.sh:** This script uses PLINK 2.0 to produce \*_qcd.bim, \*_qcd.bed, \*_qcd.fam files, with the phenotype data being specified in a given sample set, e.g. T2D.txt. For T2D, a logistic model is then fitted for every variant in which 8 principal components are included as well as age. Results from the files T2D_res_PHENO1.glm.logistic.hybrid for chromosomes 1 to 22 are then combined into a summarized text file, e.g. summary-stats-T2D.txt, which contains information regarding chromosome, position, rsID, effect size and corresponding standard error for each genetic variant. A similar process is repeated for each of the BMI sample sets.

- *Input:* \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, T2D.txt, bmi$\ddagger$.txt for $\ddagger$ $\in$ {0,25,50,75,100}, summary-stats.R

- *Output:* summary-stats-T2D.txt, summary-stats-BMI$\ddagger$.txt for $\ddagger$ $\in$ {0,25,50,75,1}

&nbsp;

##### 3. Performing MR analyses

**02-BMI-T2D.R:** This script applies four variants of MR-SimSS, together with other summary-level MR methods, to GWAS summary statistics in order to estimate the effect of BMI on T2D, under varying degrees of sample overlap.

- *Input:* summary-stats-T2D.txt, summary-stats-BMI$\ddagger$.txt for $\ddagger$ $\in$ {0,25,50,75,1}
- *Output:* res-BMI-T2D.csv

&nbsp;

### Figures

The script below was used to produce all figures in the manuscript, as well as supplementary figures, and can be found in the `figures` folder.

**figs.R** This script  was used to generate the following figures: Figs 1-4, Supplementary Figs 1-10
