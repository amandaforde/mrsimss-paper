## SAME-TRAIT EMPIRICAL MR ANALYSES


fun <- function(data){
  ## SimSS-2-IVW
  res.sim.2.ivw <- mr.simss::mr_simss(data,splits=2,mr_method = "mr_ivw")$summary
  ## SimSS-2-RAPS
  res.sim.2.raps <- mr.simss::mr_simss(data,splits=2,mr_method = "mr_raps")$summary
  ## SimSS-3-IVW
  res.sim.3.ivw <- mr.simss::mr_simss(data,splits=3,mr_method = "mr_ivw")$summary
  ## SimSS-3-RAPS
  res.sim.3.raps <- mr.simss::mr_simss(data,splits=3,mr_method = "mr_raps")$summary
  ## IVW
  res.ivw <- mr_ivw(data)
  ## MR-RAPS
  data_sig <- data %>% dplyr::filter(2*(stats::pnorm(abs(data$beta.exposure/data$se.exposure), lower.tail=FALSE)) < 5e-8)
  res.raps <- mr.raps::mr.raps(data_sig$beta.exposure,data_sig$beta.outcome,data_sig$se.exposure,data_sig$se.outcome)
  res.raps <- data.frame(method="mr_raps", nsnp=nrow(data_sig), b=res.raps$beta.hat, se=res.raps$beta.se, pval=res.raps$beta.p.value)
  ## MR-Egger
  res.egger <- mr.egger(data)
  ## Weighted median
  res.median <- mr.weighted.median(data)
  ## dIVW
  res.divw <- mr.divw::mr.divw(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome, diagnostics=FALSE)
  res.divw <- data.frame(method="mr_divw", nsnp=res.divw$n.IV, b=res.divw$beta.hat, se=res.divw$beta.se, pval=2*(stats::pnorm(abs(res.divw$beta.hat/res.divw$beta.se), lower.tail=FALSE)))
  ## compile results
  res.bmi <- list(res.sim.2.ivw, res.sim.2.raps, res.sim.3.ivw, res.sim.3.raps, res.ivw, res.raps, res.egger, res.median,res.divw)
  total_res <- data.frame(method=c(rep("0",9)))
  total_res$method <- c("SimSS-2-IVW","SimSS-2-RAPS","SimSS-3-IVW","SimSS-3-RAPS","IVW","RAPS","Egger","Weighted median","dIVW")
  total_res$method <- factor(total_res$method, levels=c("SimSS-2-IVW","SimSS-3-IVW","SimSS-2-RAPS","SimSS-3-RAPS","IVW","RAPS","Egger","Weighted median","dIVW"))

  snp <- function(i){res.bmi[[i]]$nsnp}
  b <- function(i){res.bmi[[i]]$b}
  se <- function(i){res.bmi[[i]]$se}
  total_res$b <- c(b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),b(9))
  total_res$se <- c(se(1),se(2),se(3),se(4),se(5),se(6),se(7),se(8),se(9))
  total_res$CI.lower <- total_res$b - 1.96*total_res$se
  total_res$CI.upper <- total_res$b + 1.96*total_res$se
  total_res$n.IV <- c(snp(1),snp(2),snp(3),snp(4),snp(5),snp(6),snp(7),snp(8),snp(9))
  return(total_res)
}


set.seed(1996)

###############################################################################
## BMI-1A-BMI-1B
bmi1A <-  read.table('summary-stats-BMI1A.txt',header=TRUE)
bmi1B <-  read.table('summary-stats-BMI1B.txt',header=TRUE)
data <- data.frame(SNP=bmi1A$rsid,beta.exposure=bmi1A$beta,beta.outcome=bmi1B$beta,se.exposure=bmi1A$se,se.outcome=bmi1B$se)
res.bmi1AB <- fun(data)

## BMI-1B-BMI-1A
data <- data.frame(SNP=bmi1B$rsid,beta.exposure=bmi1B$beta,beta.outcome=bmi1A$beta,se.exposure=bmi1B$se,se.outcome=bmi1A$se)
res.bmi1BA  <- fun(data)


###############################################################################
## BMI-2A-BMI-2B
bmi2A <-  read.table('summary-stats-BMI2A.txt',header=TRUE)
bmi2B <-  read.table('summary-stats-BMI2B.txt',header=TRUE)
data <- data.frame(SNP=bmi2A$rsid,beta.exposure=bmi2A$beta,beta.outcome=bmi2B$beta,se.exposure=bmi2A$se,se.outcome=bmi2B$se)
res.bmi2AB <- fun(data)

## BMI-2B-BMI-2A
data <- data.frame(SNP=bmi2B$rsid,beta.exposure=bmi2B$beta,beta.outcome=bmi2A$beta,se.exposure=bmi2B$se,se.outcome=bmi2A$se)
res.bmi2BA  <- fun(data)

###############################################################################
## BMI-3A-BMI-3B
bmi3A <-  read.table('summary-stats-BMI3A.txt',header=TRUE)
bmi3B <-  read.table('summary-stats-BMI3B.txt',header=TRUE)
data <- data.frame(SNP=bmi3A$rsid,beta.exposure=bmi3A$beta,beta.outcome=bmi3B$beta,se.exposure=bmi3A$se,se.outcome=bmi3B$se)
res.bmi3AB <- fun(data)

## BMI-3B-BMI-3A
data <- data.frame(SNP=bmi3B$rsid,beta.exposure=bmi3B$beta,beta.outcome=bmi3A$beta,se.exposure=bmi3B$se,se.outcome=bmi3A$se)
res.bmi3BA  <- fun(data)

###############################################################################
## BMI-4A-BMI-4B
bmi4A <-  read.table('summary-stats-BMI4A.txt',header=TRUE)
bmi4B <-  read.table('summary-stats-BMI4B.txt',header=TRUE)
data <- data.frame(SNP=bmi4A$rsid,beta.exposure=bmi4A$beta,beta.outcome=bmi4B$beta,se.exposure=bmi4A$se,se.outcome=bmi4B$se)
res.bmi4AB <- fun(data)

## BMI-4B-BMI-4A
data <- data.frame(SNP=bmi4B$rsid,beta.exposure=bmi4B$beta,beta.outcome=bmi4A$beta,se.exposure=bmi4B$se,se.outcome=bmi4A$se)
res.bmi4BA  <- fun(data)

###############################################################################
## BMI-5A-BMI-5B
bmi5A <-  read.table('summary-stats-BMI5A.txt',header=TRUE)
bmi5B <-  read.table('summary-stats-BMI5B.txt',header=TRUE)
data <- data.frame(SNP=bmi5A$rsid,beta.exposure=bmi5A$beta,beta.outcome=bmi5B$beta,se.exposure=bmi5A$se,se.outcome=bmi5B$se)
res.bmi5AB <- fun(data)

## BMI-5B-BMI-5A
data <- data.frame(SNP=bmi5B$rsid,beta.exposure=bmi5B$beta,beta.outcome=bmi5A$beta,se.exposure=bmi5B$se,se.outcome=bmi5A$se)
res.bmi5BA  <- fun(data)

###############################################################################
## BMI-6A-BMI-6B
bmi6A <-  read.table('summary-stats-BMI6A.txt',header=TRUE)
bmi6B <-  read.table('summary-stats-BMI6B.txt',header=TRUE)
data <- data.frame(SNP=bmi6A$rsid,beta.exposure=bmi6A$beta,beta.outcome=bmi6B$beta,se.exposure=bmi6A$se,se.outcome=bmi6B$se)
res.bmi6AB <- fun(data)

## BMI-6B-BMI-6A
data <- data.frame(SNP=bmi6B$rsid,beta.exposure=bmi6B$beta,beta.outcome=bmi6A$beta,se.exposure=bmi6B$se,se.outcome=bmi6A$se)
res.bmi6BA  <- fun(data)

###############################################################################
## BMI-7A-BMI-7B
bmi7A <-  read.table('summary-stats-BMI7A.txt',header=TRUE)
bmi7B <-  read.table('summary-stats-BMI7B.txt',header=TRUE)
data <- data.frame(SNP=bmi7A$rsid,beta.exposure=bmi7A$beta,beta.outcome=bmi7B$beta,se.exposure=bmi7A$se,se.outcome=bmi7B$se)
res.bmi7AB <- fun(data)

## BMI-7B-BMI-7A
data <- data.frame(SNP=bmi7B$rsid,beta.exposure=bmi7B$beta,beta.outcome=bmi7A$beta,se.exposure=bmi7B$se,se.outcome=bmi7A$se)
res.bmi7BA  <- fun(data)

###############################################################################
## BMI-8A-BMI-8B
bmi8A <-  read.table('summary-stats-BMI8A.txt',header=TRUE)
bmi8B <-  read.table('summary-stats-BMI8B.txt',header=TRUE)
data <- data.frame(SNP=bmi8A$rsid,beta.exposure=bmi8A$beta,beta.outcome=bmi8B$beta,se.exposure=bmi8A$se,se.outcome=bmi8B$se)
res.bmi8AB <- fun(data)

## BMI-8B-BMI-8A
data <- data.frame(SNP=bmi8B$rsid,beta.exposure=bmi8B$beta,beta.outcome=bmi8A$beta,se.exposure=bmi8B$se,se.outcome=bmi8A$se)
res.bmi8BA  <- fun(data)

###############################################################################
## BMI-9A-BMI-9B
bmi9A <-  read.table('summary-stats-BMI9A.txt',header=TRUE)
bmi9B <-  read.table('summary-stats-BMI9B.txt',header=TRUE)
data <- data.frame(SNP=bmi9A$rsid,beta.exposure=bmi9A$beta,beta.outcome=bmi9B$beta,se.exposure=bmi9A$se,se.outcome=bmi9B$se)
res.bmi9AB <- fun(data)

## BMI-9B-BMI-9A
data <- data.frame(SNP=bmi9B$rsid,beta.exposure=bmi9B$beta,beta.outcome=bmi9A$beta,se.exposure=bmi9B$se,se.outcome=bmi9A$se)
res.bmi9BA  <- fun(data)

###############################################################################
## BMI-10A-BMI-10B
bmi10A <-  read.table('summary-stats-BMI10A.txt',header=TRUE)
bmi10B <-  read.table('summary-stats-BMI10B.txt',header=TRUE)
data <- data.frame(SNP=bmi10A$rsid,beta.exposure=bmi10A$beta,beta.outcome=bmi10B$beta,se.exposure=bmi10A$se,se.outcome=bmi10B$se)
res.bmi10AB <- fun(data)

## BMI-10B-BMI-10A
data <- data.frame(SNP=bmi10B$rsid,beta.exposure=bmi10B$beta,beta.outcome=bmi10A$beta,se.exposure=bmi10B$se,se.outcome=bmi10A$se)
res.bmi10BA  <- fun(data)


res.all <- rbind(res.bmi1AB,res.bmi1BA,res.bmi2AB,res.bmi2BA,res.bmi3AB,res.bmi3BA,res.bmi4AB,res.bmi4BA,res.bmi5AB,res.bmi5BA,res.bmi6AB,res.bmi6BA,res.bmi7AB,res.bmi7BA,res.bmi8AB,res.bmi8BA,res.bmi9AB,res.bmi9BA,res.bmi10AB,res.bmi10BA)
res.all <- res.all %>%
  mutate(method = recode(method,
                         "SimSS2-mr_ivw"  = "SimSS-2-IVW",
                         "SimSS2-mr_raps" = "SimSS-2-RAPS",
                         "SimSS3-mr_ivw"  = "SimSS-3-IVW",
                         "SimSS3-mr_raps" = "SimSS-3-RAPS",
                         "mr_ivw" = "IVW",
                         "mr_raps" = "RAPS",
                         "Egger" = "Egger",
                         "Weighted median" = "Weighted median",
                         "mr_divw" = "dIVW"
  ))


write.csv(res.all,"real_data/results/res-BMI-BMI.csv",row.names=FALSE)
