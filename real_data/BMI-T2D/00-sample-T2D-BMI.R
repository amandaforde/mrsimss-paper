
#! /usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)


## read in both full BMI and T2D
T2D <- read.table(args[1], header=TRUE,fill=TRUE)
BMI <- read.table(args[2], header=TRUE, fill=TRUE)

set.seed(2345)

## Pick T2D "base" dataset, 249,749 individuals at random
pick <- sample(nrow(T2D),249749)
T2D_base <- T2D[pick,]

## Full overlap
BMI_100 <- BMI[pick,]

## Zero overlap
BMI_0 <- BMI[-pick,]

## 25% overlap
BMI_same <- BMI_100[sample(nrow(BMI_100),62437),]
BMI_diff <- BMI_0[sample(nrow(BMI_0),187312),]

BMI_25 <- rbind(BMI_same,BMI_diff)

## 50% overlap
BMI_same <- BMI_100[sample(nrow(BMI_100),124874),]
BMI_diff <- BMI_0[sample(nrow(BMI_0),124875),]
BMI_50 <- rbind(BMI_same,BMI_diff)

## 75% overlap
BMI_same <- BMI_100[sample(nrow(BMI_100),187312),]
BMI_diff <- BMI_0[sample(nrow(BMI_0),62437),]
BMI_75 <- rbind(BMI_same,BMI_diff)


write.table(T2D_base,"T2D.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','T2D'))
write.table(BMI_0,"bmi0.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(BMI_25,"bmi25.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(BMI_50,"bmi50.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(BMI_75,"bmi75.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(BMI_100,"bmi100.txt",quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
