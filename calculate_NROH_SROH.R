#!/usr/bin/env Rscript
# @Kanat Gürün

files <- list.files(, pattern="*.hom")
files <- gsub(".hom", "", files)
head(files)
for (i in files){
  x <- read.table(paste(i, ".hom",sep=""), head = TRUE)
  w <- x$KB>=1000
  y <- x$KB>=1500
  NROH1000 <- length(w[w ==TRUE])
  NROH1500 <- length(y[y ==TRUE])
  SROH1500 <- sum(x$KB[x$KB>=1500])
  SROH1000 <- sum(x$KB[x$KB>=1000])
  SROH1_2 <- sum(x$KB[x$KB>=1000&x$KB<2000])
  SROH2_4 <- sum(x$KB[x$KB>=2000&x$KB<4000])
  SROH4_8 <- sum(x$KB[x$KB>=4000&x$KB<8000])
  SROH8 <- sum(x$KB[x$KB>=8000])
  w <- data.frame(w)
  w$ID <- matrix(i, ncol = 1, nrow = nrow(w))
  w$NROH1000 <- matrix(NROH1000, ncol = 1, nrow = nrow(w))
  w$NROH1500 <- matrix(NROH1500, ncol = 1, nrow = nrow(w))
  w$SROH1500 <- matrix(SROH1500, ncol = 1, nrow = nrow(w))
  w$SROH1000 <- matrix(SROH1000, ncol = 1, nrow = nrow(w))
  w$SROH1_2 <- matrix(SROH1_2, ncol = 1, nrow = nrow(w))
  w$SROH2_4 <- matrix(SROH2_4, ncol = 1, nrow = nrow(w))
  w$SROH4_8 <- matrix(SROH4_8, ncol = 1, nrow = nrow(w))
  w$SROH8 <- matrix(SROH8, ncol = 1, nrow = nrow(w))
  w <- w[!duplicated(w[3:4]),]
  write.table(w, file = "NROH_SROH_allSROHcategories.csv", col.names = F, row.names = T, append = T, quote = F, sep = ",")
}
