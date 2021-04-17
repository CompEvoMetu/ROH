# Runs of Homozygosity analysis for ancient genomes

**Pipeline for ROH analysis for ancient genomes (Ceballos, Gürün et al. 2021, Current Biology)**

Here we describe a pipeline to call ROH, starting from bam files as the first type of input files.

### dependencies

**for real data**
- samtools
- bcftools
- PLINK
- [trimBam](https://genome.sph.umich.edu/wiki/BamUtil:_trimBam)

**for simulations**
- [gargamel](https://github.com/grenaud/gargammel)

## 1. Genotyping and ROH calling

Before calling ROH, we clipped 10 bases of the ends of each read for each sample prepared by shotgun sequencing and 2 bases for each sample prepared by Uracil-DNA-glycosylase (UDG) protocol using trimBam (Jun et al. 2015).

```bash
sh ROH_calling.sh in.bam 1240K.bed hs37d5.fa
```

> To calculate NROH and SROH from hom files, we used ```calculate_NROH_SROH.R```

> To calculate average DP from vcf files, we used ```calculate_averageDP.R```


## 2. Simulation

> **NOTE**: These scripts work just with chromosome 1

### 2.1. Creating simulated genotypes with ROHs of different number and size on ***chromosome 1***, using the 1240K SNP list


```bash
Rscript --vanilla simulations_ROH_to_bed.R in.bed
```

The output will be 20 individual for each scenario

### 2.2. Generating simulated reads of size 70 bp with ancient DNA damage by using gargamel

```bash
ls individual_{1..20}_ROH_scenario_1_genotype.txt > simulated_individual_list.txt
ls individual_{1..20}_ROH_scenario_2_genotype.txt >> simulated_individual_list.txt

# Text output files of Simulations_1_ROH_to_bed.R will be used as input
for line in $(cat path/simulated_individual_list.txt)
do
paste -d"\t" 1240K_only_positions.bed $line.out > $line.bed
done

for line in $(cat path/simulated_individual_list.txt)
do
path/modifyFasta.sh $line.bed # this creates two fasta files for each simulated individual with inserted ROH
done

# gargammel to generate simulated reads of size 70 bp with ancient DNA damage

for line in $(cat path/simulated_individual_list.txt)
do
cp path/$line.1.fa path/gargammel/data/endo
cp path/$line.2.fa path/gargammel/data/endo
samtools faidx path/gargammel/data/endo/$line.1.fa
samtools faidx path/gargammel/data/endo/$line.2.fa
./gargammel.pl -c 3 --comp 0,0,1 -l 70 -damage 0.03,0.4,0.01,0.3 -o data/simulation data/
rm path/gargammel/data/endo/$line.1.fa
rm path/gargammel/data/endo/$line.2.fa
rm /gargammel/data/endo/$line.1.fa.fai
rm path/gargammel/data/endo/$line.2.fa.fai
rm simulation_a.fa.gz
rm simulation.b.fa.gz
rm simulation.c.fa.gz
rm simulation_d.fa.gz
rm simulation.e.fa.gz
mv simulation_s1.fq.gz $line.s1.fq.gz
mv simulation_s2.fq.gz $line.s2.fq.gz
done
```

### 2.3. Comparing results

These simulated fastq files were mapped to human reference genome (hs37d5.fa) and clipped bam files by using trimBam. Then we called ROH by using plink (see step 1).
Finally we got hom file for each simulated indivual.

The compare results following R code were used.

```r
#!/usr/bin/env Rscript
# Dec 12, 2020
# Kanat Gürün
# The output hom files of the mapping and ROH calling pipeline will be used as input.

files <- list.files(, pattern="*.hom")
files <- gsub("*.hom", "", files)
head(files)

# to view all the detected ROH. NROH, SROH and FROH can be calculated using the resulting Excel file.
for (i in files){
  x <- read.table(paste(i, ".hom",sep=""), head = TRUE)
  x <- data.frame(POS1 = x$POS1, POS2 = x$POS2, KB = x$KB)
  j <- matrix(i, ncol = 1, nrow = nrow(x))
  x$ids <- j
  write.table(x, file = "all.csv", col.names = F, row.names = T, append = T, quote = F, sep = ",")
}

library(GenomicRanges)

# to find the numbers of true and false positives
for (i in files){
  x <- read.table(paste(i, ".hom",sep=""), head = TRUE) # observed ROH
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
```
