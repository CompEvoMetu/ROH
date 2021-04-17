#!/usr/bin/env Rscript
# Nov 16, 2020 @Mehmet Somel

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (bed file).\n", call.=FALSE)
}

# input 1: genotype file, "geno", with the 2 alleles per 1240K SNP on chr1. the order does not matter.
geno <- read.table(args[1])  # from Kanat Gurun, received Sep 12, 2020
head(geno)
nsnps <- nrow(geno)
nsnps

# input 2: list of simulated ROH lengths.

# the list of 4 ROH scenarios:
# 1) 2 x 1 Mb
# 2) 10 x 1 Mb

# simple estimation of the size of ROHs in SNP number, that should correspond or 1 and 10 Mb
# this also allows some variance among ROH sizes
shortROHlen <- 1e6
longROHlen <- 1e7

scenarios <- list(
  rep(shortROHlen, 2),
  rep(shortROHlen, 10))

scenarios
# we plan to run: 4 scenarios x 50 indv x 1 repeat (x 4 coverages)
# so we create 50 ROH lists per scenario

# we divide chr1 into 10 segments. each segment will have 1 ROH max
seg_starts <- round(seq(1, nsnps, length.out = 11))[-11]  # start positions of segments
seg_ends <- round(seq(1, nsnps, length.out = 11))[-1] - 1  # end positions of segments

cbind(seg_starts, seg_ends)

# centromere start and end
centromere <- c(121311798, 143979325)

###########

# scenario_1
i = 1
list_scenario_1 <- list()
x = 1; while(x < 51) {
  ROHlengths <- scenarios[[i]]  # for that scenario
  # define empty list and matrices
  ROH_pos_list <- list()
  start_end_snp_mat <- matrix(,0,2)
  start_end_bp_mat <- matrix(,0,2)
  numROHs <- length(ROHlengths)  # number of ROHs to be simulated in that scenario
  segs_chosen <- sample( 1:10, numROHs) # choose the segments, as many as the number of ROHs in that scenario
  for (j in 1:numROHs) {  # in each segment
    ROHseg <- segs_chosen[j]    # the segment
    ROHlength <- ROHlengths[j]  # the length of that ROH

    # the farthest 3' position we can go in that segment, given ROHlength (1 Mb or 10 Mb)
    end_pos_bp <- geno[ seg_ends[ROHseg], 2] - ROHlength
    # the snp closest to the farthest 3' position
    end_pos_bp_snp <- which.min(abs(geno[, 2] - end_pos_bp))
    # upper and lower boundaries for sampling the ROH
    up_bound_bp <- geno[ seg_starts[ROHseg], 2]
    down_bound_bp <- geno[ end_pos_bp_snp, 2]

    # randomly choose the new ROH start position in bp in that segment
    start_bp <- sample( up_bound_bp:down_bound_bp, 1)
    # define the end position in bp
    end_bp <- start_bp + ROHlength
    # the closest SNPs to the start and end positions of the new ROH
    start <- which.min(abs(geno[, 2] - start_bp))
    end <- which.min(abs(geno[, 2] - end_bp))
    newROH <- list( start:end ) # list the positions of the SNPs inside the new ROH

    ROH_pos_list <- c(ROH_pos_list, newROH)
    start_end_bp <- c(geno[start,2], geno[end,2])  # the start + end positions in bp
    start_end_bp_mat <- rbind(start_end_bp_mat, start_end_bp)
  }
  # check that there is >1.5 Mb overlap between conseq ROHs, all ROHs are size <15 Mb, and ROHs don't overlap the centromere
  start_end_bp_mat <- apply(start_end_bp_mat, 2, sort)
  # distance between conseq ROHs
  ROHdistance <- start_end_bp_mat[-1,1] - start_end_bp_mat[-numROHs,2]
  # size of ROHs
  ROHsize <- start_end_bp_mat[,2] - start_end_bp_mat[,1]
  # do the ROHs overlap with the centromere? a boolean vector
  CENTROMERE <- (start_end_bp_mat[,1] <= centromere[2] & start_end_bp_mat[,2] >= centromere[1])
  # if distance between conseq ROHs <1.5Mb or any ROH size >15 Mb or if they overlap with the centromere, skip that turn
  if ((sum(ROHdistance <= 1.5e6) > 0) | (sum(ROHsize > 15e6) > 0) | sum(CENTROMERE)) {
    print ("proximal ROHs or large ROHs or centromeric ROHs, repeating")
  } else {
    # if not, add the new list and matrix
    list_scenario_1 <- c(list_scenario_1,
                         list(list(sort(unlist(ROH_pos_list)), start_end_bp_mat)))
    x <- x + 1
    print (x)
  }
}
length(list_scenario_1)

# plot of sizes
matplot(sapply(list_scenario_1, function(x) {
  start_end_bp_mat <- x[[2]]
  size <- start_end_bp_mat[,2] - start_end_bp_mat[,1]
  sort(size)
}), ylab="ROH size")

# create the 20 individuals with ROHs under scenario_1
# include realistic proportions of heterozygosity at these SNPs
geno_indv_list <- list()
for (i in 1:20) {
  heterozygosity <- 0.10  # expected heterozygosity at those SNPs
  geno_indv <- t(apply(as.matrix(geno[,4:5]), 1, function(x) {
    x2 <- sample(x)
    first_allele <- x2[1]
    second_allele <- sample(x2, prob = c(1 - heterozygosity, heterozygosity))[1]
    return( c(first_allele, second_allele) )
  }))
  # summary(geno_indv[,1]==geno_indv[,2])
  # head(geno_indv)
  # now create ROHs at the simulated ROH positions
  ROHpos <- list_scenario_1[[i]][[1]]
  geno_indv[ROHpos, 2] <- geno_indv[ROHpos, 1]  # copy 1 allele on the other. both are random, so should be ok
  # summary(geno_indv[,1]==geno_indv[,2])
  geno_indv_list <- c(geno_indv_list, list(geno_indv))
  # write the genotype matrix into a file
  write.table(geno_indv_list[[i]], file=paste("individual_", i, "_ROH_scenario_1_genotype.txt", sep=""), quote=F, col.names = F, row.names = F)
  # write the ROH positions into a file (to be later checked with ROH estimates)
  write.table(list_scenario_1[[i]][[2]], file=paste("individual_", i, "_ROH_scenario_1_ROH_positions.txt", sep=""), quote=F, col.names = F, row.names = F)
}


# scenario_2
i = 2
list_scenario_2 <- list()
x = 1; while(x < 51) {
  ROHlengths <- scenarios[[i]]  # for that scenario
  # define empty list and matrices
  ROH_pos_list <- list()
  start_end_snp_mat <- matrix(,0,2)
  start_end_bp_mat <- matrix(,0,2)
  numROHs <- length(ROHlengths)  # number of ROHs to be simulated in that scenario
  segs_chosen <- sample( 1:10, numROHs) # choose the segments, as many as the number of ROHs in that scenario
  for (j in 1:numROHs) {  # in each segment
    ROHseg <- segs_chosen[j]    # the segment
    ROHlength <- ROHlengths[j]  # the length of that ROH

    # the farthest 3' position we can go in that segment, given ROHlength (1 Mb or 10 Mb)
    end_pos_bp <- geno[ seg_ends[ROHseg], 2] - ROHlength
    # the snp closest to the farthest 3' position
    end_pos_bp_snp <- which.min(abs(geno[, 2] - end_pos_bp))
    # upper and lower boundaries for sampling the ROH
    up_bound_bp <- geno[ seg_starts[ROHseg], 2]
    down_bound_bp <- geno[ end_pos_bp_snp, 2]

    # randomly choose the new ROH start position in bp in that segment
    start_bp <- sample( up_bound_bp:down_bound_bp, 1)
    # define the end position in bp
    end_bp <- start_bp + ROHlength
    # the closest SNPs to the start and end positions of the new ROH
    start <- which.min(abs(geno[, 2] - start_bp))
    end <- which.min(abs(geno[, 2] - end_bp))
    newROH <- list( start:end ) # list the positions of the SNPs inside the new ROH

    ROH_pos_list <- c(ROH_pos_list, newROH)
    start_end_bp <- c(geno[start,2], geno[end,2])  # the start + end positions in bp
    start_end_bp_mat <- rbind(start_end_bp_mat, start_end_bp)
  }
  # check that there is >1.5 Mb overlap between conseq ROHs, all ROHs are size <15 Mb, and ROHs don't overlap the centromere
  start_end_bp_mat <- apply(start_end_bp_mat, 2, sort)
  # distance between conseq ROHs
  ROHdistance <- start_end_bp_mat[-1,1] - start_end_bp_mat[-numROHs,2]
  # size of ROHs
  ROHsize <- start_end_bp_mat[,2] - start_end_bp_mat[,1]
  # do the ROHs overlap with the centromere? a boolean vector
  CENTROMERE <- (start_end_bp_mat[,1] <= centromere[2] & start_end_bp_mat[,2] >= centromere[1])
  # if distance between conseq ROHs <1.5Mb or any ROH size >15 Mb or if they overlap with the centromere, skip that turn
  if ((sum(ROHdistance <= 1.5e6) > 0) | (sum(ROHsize > 15e6) > 0) | sum(CENTROMERE)) {
    print ("proximal ROHs or large ROHs or centromeric ROHs, repeating")
  } else {
    # if not, add the new list and matrix
    list_scenario_2 <- c(list_scenario_2,
                         list(list(sort(unlist(ROH_pos_list)), start_end_bp_mat)))
    x <- x + 1
    print (x)
  }
}
length(list_scenario_2)

# plot of sizes
matplot(sapply(list_scenario_2, function(x) {
  start_end_bp_mat <- x[[2]]
  size <- start_end_bp_mat[,2] - start_end_bp_mat[,1]
  sort(size)
}), ylab="ROH size")

# create the 50 individuals with ROHs under scenario_2
# include realistic proportions of heterozygosity at these SNPs
geno_indv_list <- list()
for (i in 1:20) {
  heterozygosity <- 0.10  # expected heterozygosity at those SNPs
  geno_indv <- t(apply(as.matrix(geno[,4:5]), 1, function(x) {
    x2 <- sample(x)
    first_allele <- x2[1]
    second_allele <- sample(x2, prob = c(1 - heterozygosity, heterozygosity))[1]
    return( c(first_allele, second_allele) )
  }))
  # summary(geno_indv[,1]==geno_indv[,2])
  # head(geno_indv)
  # now create ROHs at the simulated ROH positions
  ROHpos <- list_scenario_2[[i]][[1]]
  geno_indv[ROHpos, 2] <- geno_indv[ROHpos, 1]  # copy 1 allele on the other. both are random, so should be ok
  # summary(geno_indv[,1]==geno_indv[,2])
  geno_indv_list <- c(geno_indv_list, list(geno_indv))
  # write the genotype matrix into a file
  write.table(geno_indv_list[[i]], file=paste("individual_", i, "_ROH_scenario_2_genotype.txt", sep=""), quote=F, col.names = F, row.names = F)
  # write the ROH positions into a file (to be later checked with ROH estimates)
  write.table(list_scenario_2[[i]][[2]], file=paste("individual_", i, "_ROH_scenario_2_ROH_positions.txt", sep=""), quote=F, col.names = F, row.names = F)
}
