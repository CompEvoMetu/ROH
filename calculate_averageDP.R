#!/usr/bin/env Rscript
# @Kanat Gürün

# calculate average DP using vcf files
library(vcfR)

for (i in files){
  x <- read.table(paste(i, ".vcf",sep=""))
  w <- matrix(i, ncol = 1, nrow = nrow(w))
  w <- data.frame(w)
  y <- read.vcfR(paste(i, ".vcf",sep=""))
  noofsnps <- length(x$V1)
  averageDP <- mean(extract.info(y, element="DP", as.numeric=TRUE))
  w$No_of_SNPs <- noofsnps
  w$Average_DP <- averageDP
  write.table(w, file = "_ROH_SNPnoDP.csv", col.names = F, row.names = T, append = T, quote = F, sep = ",")
}
