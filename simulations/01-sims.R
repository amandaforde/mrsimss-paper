## MAIN SIMULATION STUDY SCRIPT:

## This script allows us to run the entire main simulation study which evaluates
## and compares four variants of MR-SimSS together with other classical MR
## approaches.


## Load all required packages for simulations.
library(devtools)
devtools::install_github("amandaforde/mr.simss")
library(mr.simss)
library(tidyverse)
library(parallel)
library(mr.raps)
library(mr.divw)
library(stats)


## Run 'sims_funs.R' in order to define additional functions required for
## simulations below.
source("simulations/00-sims-funs.R")


################################################################################

## Total number of simulations:
tot_sim <- 100
n_snps <- 10^6
beta_xy <- 0.3

################################################################################
## PART A) Sample Size = 200,000

## Set of parameters:
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  h2 = c(0.3,0.7),
  prop_effect = c(0.01,0.001),
  cor_xy = c(-0.1,0.1,0.3,0.5),
  frac_overlap = c(0,0.25,0.5,0.75,1),
  n_x = c(200000)
)

set.seed(1998)
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

## Organise results:
total_res <- data.frame(sim=c(),h2=c(),prop_effect=c(),cor_xy=c(),frac_overlap=c(),n_x=c(),method=c(),nsnp=c(),b=c(),se=c(),pval=c())

for(i in 1:length(res)){
  for(j in 2:length(res[[1]])){
    if(is.null(res[[i]][j][[1]]) == FALSE){
      res1 <- cbind(data.frame(res[[i]][1]$params),res[[i]][j][[1]])
      if(j == 2 | j == 3 ){res1$method <- paste("SimSS2", res1$method, sep = "-")
      res1 <- res1[,c(1:8,10:12)]}
      if(j == 4 | j == 5 ){res1$method <- paste("SimSS3", res1$method, sep = "-")
      res1 <- res1[,c(1:8,10:12)]}
      total_res <- rbind(total_res,res1)
    }
  }
}

write.table(total_res,"simulations/results/sims-res-200.csv",quote=FALSE,row.names=FALSE)


################################################################################
## PART B) Sample Size = 500,000

## Set of parameters:
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  h2 = c(0.3,0.7),
  prop_effect = c(0.01,0.001),
  cor_xy = c(-0.1,0.1,0.3,0.5),
  frac_overlap = c(0,0.25,0.5,0.75,1),
  n_x = c(500000)
)

set.seed(1998)
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

## Organise results:
total_res <- data.frame(sim=c(),h2=c(),prop_effect=c(),cor_xy=c(),frac_overlap=c(),n_x=c(),method=c(),nsnp=c(),b=c(),se=c(),pval=c())

for(i in 1:length(res)){
  for(j in 2:length(res[[1]])){
    if(is.null(res[[i]][j][[1]]) == FALSE){
      res1 <- cbind(data.frame(res[[i]][1]$params),res[[i]][j][[1]])
      if(j == 2 | j == 3 ){res1$method <- paste("SimSS2", res1$method, sep = "-")
      res1 <- res1[,c(1:8,10:12)]}
      if(j == 4 | j == 5 ){res1$method <- paste("SimSS3", res1$method, sep = "-")
      res1 <- res1[,c(1:8,10:12)]}
      total_res <- rbind(total_res,res1)
    }
  }
}

write.table(total_res,"simulations/results/sims-res-500.csv",quote=FALSE,row.names=FALSE)

################################################################################
## PART C) Sample Size = 50,000

## Set of parameters:
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  h2 = c(0.3,0.7),
  prop_effect = c(0.01,0.001),
  cor_xy = c(-0.1,0.1,0.3,0.5),
  frac_overlap = c(0,0.25,0.5,0.75,1),
  n_x = c(500000)
)

set.seed(1998)
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

## Organise results:
total_res <- data.frame(sim=c(),h2=c(),prop_effect=c(),cor_xy=c(),frac_overlap=c(),n_x=c(),method=c(),nsnp=c(),b=c(),se=c(),pval=c())

for(i in 1:length(res)){
  for(j in 2:length(res[[1]])){
    if(is.null(res[[i]][j][[1]]) == FALSE){
      res1 <- cbind(data.frame(res[[i]][1]$params),res[[i]][j][[1]])
      if(j == 2 | j == 3 ){res1$method <- paste("SimSS2", res1$method, sep = "-")
      res1 <- res1[,c(1:8,10:12)]}
      if(j == 4 | j == 5 ){res1$method <- paste("SimSS3", res1$method, sep = "-")
      res1 <- res1[,c(1:8,10:12)]}
      total_res <- rbind(total_res,res1)
    }
  }
}

write.table(total_res,"simulations/results/sims-res-50.csv",quote=FALSE,row.names=FALSE)
