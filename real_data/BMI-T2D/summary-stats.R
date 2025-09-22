#! usr/bin/env Rscript


args=commandArgs(trailingOnly=TRUE)


## function to create a subsetted version of results for each chromosome, suitable for use with package

sub <- function(x){
  x <- x[,c(1,2,3,4,5,6,9,10)]
  names(x)[1] <- 'chr'
  names(x)[2] <- 'pos'
  names(x)[3] <- 'rsid'
  names(x)[4] <- 'ref_allele'
  names(x)[5] <- 'alt_allele'
  names(x)[6] <- 'A1_allele'
  names(x)[7] <- 'beta'
  names(x)[8] <- 'se'
  return(x)
}

chr1 <- read.table(args[1],fill=TRUE)
chr1 <- sub(chr1)

chr2 <- read.table(args[2],fill=TRUE)
chr2 <- sub(chr2)

chr3 <- read.table(args[3],fill=TRUE)
chr3 <- sub(chr3)

chr4 <- read.table(args[4],fill=TRUE)
chr4 <- sub(chr4)

chr5 <- read.table(args[5],fill=TRUE)
chr5 <- sub(chr5)

chr6 <- read.table(args[6],fill=TRUE)
chr6 <- sub(chr6)

chr7 <- read.table(args[7],fill=TRUE)
chr7 <- sub(chr7)

chr8 <- read.table(args[8],fill=TRUE)
chr8 <- sub(chr8)

chr9 <- read.table(args[9],fill=TRUE)
chr9 <- sub(chr9)

chr10 <- read.table(args[10],fill=TRUE)
chr10 <- sub(chr10)

chr11 <- read.table(args[11],fill=TRUE)
chr11 <- sub(chr11)

chr12 <- read.table(args[12],fill=TRUE)
chr12 <- sub(chr12)

chr13 <- read.table(args[13],fill=TRUE)
chr13 <- sub(chr13)

chr14 <- read.table(args[14],fill=TRUE)
chr14 <- sub(chr14)

chr15 <- read.table(args[15],fill=TRUE)
chr15 <- sub(chr15)

chr16 <- read.table(args[16],fill=TRUE)
chr16 <- sub(chr16)

chr17 <- read.table(args[17],fill=TRUE)
chr17 <- sub(chr17)

chr18 <- read.table(args[18],fill=TRUE)
chr18 <- sub(chr18)

chr19 <- read.table(args[19],fill=TRUE)
chr19 <- sub(chr19)

chr20 <- read.table(args[20],fill=TRUE)
chr20 <- sub(chr20)

chr21 <- read.table(args[21],fill=TRUE)
chr21 <- sub(chr21)

chr22 <- read.table(args[22],fill=TRUE)
chr22 <- sub(chr22)

chr.all <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)

write.table(chr.all, args[23], quote=FALSE, row.names=FALSE)
