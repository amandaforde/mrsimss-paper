## SIMULATION STUDY SCRIPT - FUNCTIONS:

## This script provides all functions required for simulations.

################################################################################
# 1. sim_mr_stats() - Given values for the total number of variants (n_snps),
# proportion of effect variants relative to exposure (prop_effect), exposure
# heritability (h2), fraction of sample overlap (frac_overlap), sample sizes for
# exposure and outcome GWASs (n_x, n_y), exposure-outcome correlation (cor_xy)
# and the causal effect of exposure on outcome (beta_xy), GWAS summary
# statistics, i.e. variant-exposure and variant-outcome association estimates and
# corresponding standard errors, are simulated for each variant.


sim_mr_stats <- function(n_snps, prop_effect, h2, frac_overlap, n_x, n_y, cor_xy, beta_xy){
  n_overlap <- frac_overlap*min(n_x, n_y)
  maf <- runif(n_snps, 0.01, 0.5)
  effect_snps <- n_snps*prop_effect
  index <- sample(1:n_snps, ceiling(effect_snps), replace=FALSE) # random sampling
  beta_gx <- rep(0,n_snps)
  beta_gx[index] <- rnorm(length(index),0,1)
  var_x <- sum(2*maf*(1-maf)*beta_gx^2)/h2
  if(var_x != 0){beta_gx <- beta_gx/sqrt(var_x)} # scaling to represent an exposure with variance 1
  beta_gy <- beta_gx * beta_xy

  var_gx <- 1/(n_x*2*maf*(1-maf)) # var(X)=1
  var_gy <- 1/(n_y*2*maf*(1-maf)) # var(Y)=1
  cov_gx_gy <- ((n_overlap*cor_xy)/(n_x*n_y))*(1/(2*maf*(1-maf)))
  # create covariance matrix for each SNP
  cov_array <- array(dim=c(2, 2, n_snps))
  cov_array[1,1,] <- var_gx
  cov_array[2,1,] <- cov_gx_gy
  cov_array[1,2,] <- cov_array[2,1,]
  cov_array[2,2,] <- var_gy

  # summary_stats <- apply(cov_array, 3, function(x){MASS::mvrnorm(n=1, mu=c(0,0), Sigma=x)})
  # mvrnorm replaced by the following:
  tr_vec <- apply(cov_array,3,function(x){x[1,1]+x[2,2]})
  s_vec <- apply(cov_array,3,function(x){sqrt(sum(x[1,1]*x[2,2]-x[1,2]*x[2,1]))})
  t_vec <- sqrt(tr_vec+2*s_vec)
  s_vec <- rep(s_vec,times=rep(4,n_snps))
  I_vec <- rep(c(1,0,0,1),times=rep(n_snps))
  t_vec <- rep(t_vec,times=rep(4,n_snps))
  sqrt_array <- array((as.vector(cov_array)+s_vec*I_vec)/t_vec,dim=c(2,2,n_snps))

  Z_array <- array(rnorm(2*n_snps),dim=c(1,2,n_snps))
  Z_array <- abind::abind(Z_array,Z_array,along=1)
  sqrt_array_normal <- sqrt_array*Z_array
  # Rearrange array so to use matrix multiplication:
  sqrt_array_normal <- aperm(a=sqrt_array_normal,perm=c(3,1,2))
  dim1 <- sqrt_array_normal[,1,] %*% matrix(c(1,1),nrow=2)
  dim2 <- sqrt_array_normal[,2,] %*% matrix(c(1,1),nrow=2)
  summary_stats <- cbind(dim1,dim2)

  summary_stats <- (summary_stats + cbind(beta_gx, beta_gy))

  data <- tibble(
    SNP = 1:n_snps,
    beta.exposure = summary_stats[,1],
    beta.outcome = summary_stats[,2],
    se.exposure = sqrt(var_gx),
    se.outcome = sqrt(var_gy)
  )
  return(data)
}


################################################################################
# 2. mr_ivw() - Given a set of GWAS summary statistics and a significance threshold
# (threshold=5e-8), this function applies the classical inverse variance weighted
# (IVW) method.

mr_ivw <- function(data, threshold=5e-8){
  ## subset data set
  data_sig <- data %>% dplyr::filter(2*(stats::pnorm(abs(data$beta.exposure/data$se.exposure), lower.tail=FALSE)) < threshold)
  if(nrow(data_sig) < 3){return(NULL)}else{
    ivw.res <- summary(stats::lm(data_sig$beta.outcome ~ -1 + data_sig$beta.exposure, weights = 1/data_sig$se.outcome^2))
    b <- ivw.res$coef[1,1]
    se <- ivw.res$coef[1,2]/min(1,ivw.res$sigma)
    pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)
    return(data.frame(method="mr_ivw", nsnp=nrow(data_sig), b=b, se=se, pval=pval))
  }
}


################################################################################
# 3. run_sim() - This is the main simulation function. Given certain parameters
# which we wish to alter over our set of simulations, run_sim() applies four
# variants of MR-SimSS, together with classical IVW, MR-RAPS and dIVW to
# simulated data set.

run_sim <- function(h2,prop_effect,cor_xy,frac_overlap,n_x,sim){

  data <- sim_mr_stats(n_snps, prop_effect, h2, frac_overlap, n_x, n_x, cor_xy, beta_xy)

  ## SimSS-2-IVW
  res.sim.2.ivw <- mr.simss::mr_simss(data,mr_method="mr_ivw",threshold=5e-8,n.iter=1000,splits=2,est.lambda=TRUE)$summary

  ## SimSS-2-RAPS
  res.sim.2.raps <- mr.simss::mr_simss(data,mr_method="mr_raps",threshold=5e-8,n.iter=1000,splits=2,est.lambda=TRUE)$summary

  ## SimSS-3-IVW
  res.sim.3.ivw <- mr.simss::mr_simss(data,mr_method="mr_ivw",threshold=5e-8,n.iter=1000,splits=3,est.lambda=TRUE)$summary

  ## SimSS-3-RAPS
  res.sim.3.raps <- mr.simss::mr_simss(data,mr_method="mr_raps",threshold=5e-8,n.iter=1000,splits=3,est.lambda=TRUE)$summary

  ## IVW
  res.ivw <- mr_ivw(data)

  ## MR-RAPS
  data_sig <- data %>% dplyr::filter(2*(stats::pnorm(abs(data$beta.exposure/data$se.exposure), lower.tail=FALSE)) < 5e-8)
  res.raps <- mr.raps::mr.raps(data_sig$beta.exposure,data_sig$beta.outcome,data_sig$se.exposure,data_sig$se.outcome)
  res.raps <- data.frame(method="mr_raps", nsnp=nrow(data_sig), b=res.raps$beta.hat, se=res.raps$beta.se, pval=res.raps$beta.p.value)

  ## dIVW
  res.divw <- mr.divw::mr.divw(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome, diagnostics=FALSE)
  res.divw <- data.frame(method="mr_divw", nsnp=res.divw$n.IV, b=res.divw$beta.hat, se=res.divw$beta.se, pval=2*(stats::pnorm(abs(res.divw$beta.hat/res.divw$beta.se), lower.tail=FALSE)))

  ## compile results
  results <- list(params = data.frame(sim=sim,h2=h2,prop_effect=prop_effect,cor_xy=cor_xy,frac_overlap=frac_overlap,n_x=n_x), res.sim.2.ivw, res.sim.2.raps, res.sim.3.ivw, res.sim.3.raps, res.ivw, res.raps, res.divw)
  return(results)
}


