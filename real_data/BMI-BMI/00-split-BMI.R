#! /usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

bmi_full <- read.table(args[1],header=TRUE,fill=TRUE)
bmi_full <- bmi_full[is.na(bmi_full$f) == FALSE,]
sample_size <- floor(0.5*nrow(bmi_full))

set.seed(123)

## PAIR 1
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi1A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi1B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 2
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi2A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi2B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 3
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi3A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi3B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 4
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi4A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi4B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 5
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi5A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi5B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 6
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi6A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi6B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 7
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi7A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi7B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 8
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi8A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi8B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 9
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi9A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi9B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

## PAIR 10
picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]
write.table(gwas_A,"bmi10A.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,"bmi10B.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))

