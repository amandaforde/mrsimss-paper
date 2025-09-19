## SIMULATION STUDY SCRIPT - varying significance thresholds:

## This script allows us to run a small simulation study which assesses the
## performance of SimSS-3-RAPS in low-power settings with different significance
## thresholds applied in the selection step.

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

## Set of parameters:
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  h2 = c(0.3),
  prop_effect = c(0.01),
  cor_xy = c(0.5),
  frac_overlap = c(0,0.25,0.5,0.75,1),
  n_x = c(50000)
)

set.seed(1353)

run_sim <- function(h2,prop_effect,cor_xy,frac_overlap,n_x,sim){
  data <- sim_mr_stats(n_snps, prop_effect, h2, frac_overlap, n_x, n_x, cor_xy, beta_xy)

  ## SimSS-3-RAPS-5e-8
  res.sim.3.raps.1 <- mr.simss::mr_simss(data,mr_method="mr_raps",threshold=5e-8,n.iter=1000,splits=3,est.lambda=TRUE)$summary

  ## SimSS-3-RAPS-5e-6
  res.sim.3.raps.2 <- mr.simss::mr_simss(data,mr_method="mr_raps",threshold=5e-6,n.iter=1000,splits=3,est.lambda=TRUE)$summary

  ## SimSS-3-RAPS-5e-5
  res.sim.3.raps.3 <- mr.simss::mr_simss(data,mr_method="mr_raps",threshold=5e-5,n.iter=1000,splits=3,est.lambda=TRUE)$summary

  ## SimSS-3-RAPS-5e-4
  res.sim.3.raps.4 <- mr.simss::mr_simss(data,mr_method="mr_raps",threshold=5e-4,n.iter=1000,splits=3,est.lambda=TRUE)$summary

  ## compile results
  results <- list(params = data.frame(sim=sim,h2=h2,prop_effect=prop_effect,cor_xy=cor_xy,frac_overlap=frac_overlap,n_x=n_x),res.sim.3.raps.1,res.sim.3.raps.2,res.sim.3.raps.3,res.sim.3.raps.4)
  return(results)
}

res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

## Organise results:
total_res <- data.frame(sim=c(),h2=c(),prop_effect=c(),cor_xy=c(),frac_overlap=c(),n_x=c(),method=c(),nsnp=c(),b=c(),se=c(),pval=c())

for(i in 1:length(res)){
  for(j in 2:length(res[[1]])){
    if(is.null(res[[i]][j][[1]]) == FALSE){
      res1 <- cbind(data.frame(res[[i]][1]$params),res[[i]][j][[1]])
      if(j == 2 ){res1$method <- paste("5e-8", res1$method, sep = "-")}
      if(j == 3 ){res1$method <- paste("5e-6", res1$method, sep = "-")}
      if(j == 4 ){res1$method <- paste("5e-5", res1$method, sep = "-")}
      if(j == 5 ){res1$method <- paste("5e-4", res1$method, sep = "-")}
      res1 <- res1[,c(1:8,10:12)]}
    total_res <- rbind(total_res,res1)
  }
}

write.table(total_res,"simulations/results/sims-res-thresh.csv",quote=FALSE,row.names=FALSE)

