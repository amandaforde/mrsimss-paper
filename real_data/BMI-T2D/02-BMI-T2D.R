## BMI-T2D MR ANALYSES with varying sample overlap

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
## T2D-BMI-zero-overlap
T2D <-  read.table('summary-stats-T2D.txt',header=TRUE)
T2D$beta <- log(T2D$beta)
bmi0 <-  read.table('summary-stats-BMI-0.txt',header=TRUE)
data <- data.frame(SNP=bmi0$rsid,beta.exposure=bmi0$beta,beta.outcome=T2D$beta,se.exposure=bmi0$se,se.outcome=T2D$se)
res.0 <- fun(data)

###############################################################################
## T2D-BMI-25%-overlap
bmi25 <-  read.table('summary-stats-BMI-25.txt',header=TRUE)
data <- data.frame(SNP=bmi25$rsid,beta.exposure=bmi25$beta,beta.outcome=T2D$beta,se.exposure=bmi25$se,se.outcome=T2D$se)
res.25 <- fun(data)

###############################################################################
## T2D-BMI-50%-overlap
bmi50 <-  read.table('summary-stats-BMI-50.txt',header=TRUE)
data <- data.frame(SNP=bmi50$rsid,beta.exposure=bmi50$beta,beta.outcome=T2D$beta,se.exposure=bmi50$se,se.outcome=T2D$se)
res.50 <- fun(data)

###############################################################################
## T2D-BMI-75%-overlap
bmi75 <-  read.table('summary-stats-BMI-75.txt',header=TRUE)
data <- data.frame(SNP=bmi75$rsid,beta.exposure=bmi75$beta,beta.outcome=T2D$beta,se.exposure=bmi75$se,se.outcome=T2D$se)
res.75 <- fun(data)

###############################################################################
## T2D-BMI-100%-overlap
bmi100 <-  read.table('summary-stats-BMI-100.txt',header=TRUE)
data <- data.frame(SNP=bmi100$rsid,beta.exposure=bmi100$beta,beta.outcome=T2D$beta,se.exposure=bmi100$se,se.outcome=T2D$se)
res.100 <- fun(data)

###############################################################################

res.0$overlap <- c(rep("0%",11))
res.25$overlap <- c(rep("25%",11))
res.50$overlap <- c(rep("50%",11))
res.75$overlap <- c(rep("75%",11))
res.100$overlap <- c(rep("100%",11))

res.all <- rbind(res.0,res.25,res.50,res.75,res.100)
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


write.csv(res.all,"real_data/results/res-BMI-T2D.csv",row.names=FALSE)


