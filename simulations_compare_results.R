# Dec 12, 2020
# Kanat Gürün
# The output hom files of the mapping and ROH calling pipeline will be used as input.

files <- list.files(, pattern="*genotype.all.hs37d5.cons.90perc.trimBAM.hom")
files <- gsub("_genotype.all.hs37d5.cons.90perc.trimBAM.hom", "", files)
head(files)

# to view all the detected ROH. NROH, SROH and FROH can be calculated using the resulting Excel file.
for (i in files){
  x <- read.table(paste(i, "_genotype.all.hs37d5.cons.90perc.trimBAM.hom",sep=""), head = TRUE)
  x <- data.frame(POS1 = x$POS1, POS2 = x$POS2, KB = x$KB)
  j <- matrix(i, ncol = 1, nrow = nrow(x))
  x$ids <- j
  write.table(x, file = "all.csv", col.names = F, row.names = T, append = T, quote = F, sep = ",")
}


library(GenomicRanges)

# to find the numbers of true and false positives
for (i in files){
  x <- read.table(paste(i, "_genotype.all.hs37d5.cons.90perc.trimBAM.hom",sep=""), head = TRUE) # observed ROH
  x <- data.frame(chr = x$CHR, start = x$POS1, end = x$POS2, KB = x$KB)
  x2 <- makeGRangesFromDataFrame(x)
  y <- read.table(paste(i, "_ROH_positions.txt",sep="")) # inserted/expected ROH
  y <- data.frame(chr = 1, start = y$V1, end = y$V2)
  y2 <- makeGRangesFromDataFrame(y)
  w <- overlapsAny(x2,y2)==TRUE & x$KB>=1000 # or x$KB>=10000 - ROH longer than a threshold value observed in expected positions
  TruePos <- length(w[w ==TRUE]) # no of true positives
  w <- data.frame(w)
  v <- x$KB>=1000 # or x$KB>=10000
  AllPos <- length(v[v ==TRUE]) # no of all positives
  FalsePos <- AllPos - TruePos  # no of false positives
  v <- data.frame(v)
  j <- matrix(i, ncol = 1, nrow = nrow(w)) # sample IDs
  k <- matrix(TruePos, ncol = 1, nrow = nrow(w))
  l <- matrix(FalsePos, ncol = 1, nrow = nrow(w))
  w$v <- v$v
  w$ids <- j
  w$true_positives <- k
  w$false_positives <- l
  write.table(w, file = "results.csv", col.names = F, row.names = T, append = T, quote = F, sep = ",")
}
