### G1 - ANALYSIS 1 - CT~FT + HD~FT
#- FUNCTION ------------------------------------#
#gc               Genotype Characterization  
#forch            Format and Check
#grmx             Genetic Relationship MatriX
#snpa             Single SNP Analysis
#emmax            EMMAX analysis
#plotting         Plots building


# installing packages
# install.packages("data.table")
# install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# install.packages("ggplot2")
# install.packages("qqman")
# install.packages("tidyverse")
# install.packages("plotly")
# install.packages('R.utils')

#upload libraries for the analysis
library(data.table)  
library(SNPRelate)
library(ggplot2)
library(qqman)
library(gridExtra)

setwd("//ad.helsinki.fi/home/g/gfriso/Desktop/RStudioNEW/tribolium")

source("./R/functions.R")

# # #CT and HD fitness table with histograms
fitness.G1 <- read.table("./data/offspr_G1.txt")                               # from Counts.cpm-fit_CT244.txt.gz and Counts.cpm-fit_HD242.txt.gz
FTCT <- fitness.G1[grepl("CT$", rownames(fitness.G1)), , drop = F]
FTHD <- fitness.G1[grepl("HD$", rownames(fitness.G1)), , drop = F]
# hist(FTCT[,OFFSPR])
# hist(FTHD[,OFFSPR])

#matrix from CT and HD vcf file with genotype data (GT: 1|1, 1|0, 0|1, 0|0)
gtct <- fread("./data/Tribolium_castaneum_CT_G1_Tcas3.30_imputed.vcf.gz")
gthd <- fread("./data/Tribolium_castaneum_HD_G1_Tcas3.30_imputed.vcf.gz")

#FUNCTION genotype characterization ############################################
gc <- function(m) {                                                           
  mfil <- as.matrix(m[,10:ncol(m)])                               
  MAT <- sub(":.*","", mfil)
  MAT[MAT %in% c("1|0","0|1")] <- 1     
  MAT[MAT == "1|1"] <- 2
  MAT[MAT == "0|0"] <- 0
  MAT <- matrix(as.numeric(MAT), nrow = nrow(mfil), ncol = ncol(mfil))
  colnames(MAT) <- colnames(mfil)
  rownames(MAT) <- unlist(m[,3])
  return(MAT) 
}
################################################################################

cat("## Creating Genotype Matrix .. .. ..\n\n")
GTCT <- gc(gtct)
GTHD <- gc(gthd)                           

FTCT <- FTCT[-47, , drop = F]          # removing a missing sample

# FUNCTION format and check ##############no sense function######################
forch <- function(gt0, ft0, gt1, ft1) {
  # table(colnames(gt0) %in% unlist(ft0[,1]))
  # table(unlist(ft0[,1]) %in% colnames(gt0))
  GTCT <<- gt0[, unlist(ft0[,1])]
  
  
  FTCT$OFFSPR <<- as.numeric(ft0$OFFSPR)
  
  # table(colnames(gt1) %in% unlist(ft1[,1]))
  # table(unlist(ft1[,1]) %in% colnames(gt1))
  GTHD <<- gt1[, unlist(ft1[,1])]

  
  FTHD$OFFSPR <<- as.numeric(ft1$OFFSPR)

  en0 <- as.data.frame(ft0)
  en0$OFFSPR <- as.numeric(0) 
  colnames(en0)[colnames(en0) == "OFFSPR"] <- "ENVIRONMENT"
  en1 <- as.data.frame(ft1)
  en1$OFFSPR <- as.numeric(1)
  colnames(en1)[colnames(en1) == "OFFSPR"] <- "ENVIRONMENT"
  ENV <<- rbind(en0, en1) 
  
  dgt0 <- as.data.frame(t(gt0))
  dgt1 <- as.data.frame(t(gt1))
  GTGT <<- rbind(dgt0, dgt1)
  
  df0 <- as.data.frame(ft0)
  df1 <- as.data.frame(ft1)
  FTFT <<- rbind(df0, df1)
}
################################################################################

cat("## Formatting and checking genotype|fitness matrixes\n\n") 
forch(GTCT, FTCT, GTHD, FTHD)

#FUNCTION genetic relationship matrix ##########################################
grmx <- function(m) {
  mname <- as.character(deparse(substitute(m)))
  name <- paste0(mname,".gds")
  snpgdsCreateGeno(name, genmat = m, sample.id = colnames(m), snp.id = rownames(m), snpfirstdim = TRUE)
  o <- snpgdsOpen(name)
  GRM <- snpgdsGRM(o, method="GCTA",verbose = FALSE)$grm #remove.monosnp = TRUE,
  snpgdsClose(o)
  system(paste0("rm ", name))
  return(GRM)
}
################################################################################

cat("## Creating Genetic Relationship Matrix\n\n")
GRMCT <- grmx(GTCT)           
GRMHD <- grmx(GTHD)           
GTGT <- as.matrix(t(GTGT))
GRMGT <- grmx(GTGT)

# ------------------------------------------------------------------------
# library(rrBLUP)
# GRMCT1 <- A.mat(GTCT)
# plot(GRMCT1,GRMCT)
# ------------------------------------------------------------------------
cat("## Performing PCA\n\n")
# PCAa <- prcomp(GRMCT,1)
# plot(PCAa$x[,1:2])
# screeplot(PCAa)
# PCAb <- prcomp(GRMHD,1)
# plot(PCAb$x[,1:2])
# screeplot(PCAb)

#### FUNCTION single SNP analyses #######################################################
snpa <- function(gt, ft, snp) {
  gt = t(gt)
  snp = as.integer(snp)
  plot(gt[,snp],ft$OFFSPR, xlim = c(0, 2))
  abline(lm(ft$OFFSPR~gt[,snp]))
  singleSNP <- cor.test(gt[,snp],ft$OFFSPR)
  # assign(paste0(deparse(substitute(gt)), "--", deparse(substitute(snp))), singleSNP, envir = .GlobalEnv)
}
################################################################################

cat("## Performing single SNP outlier analyes\n\n") 
# snpa(GTCT, FTCT, 9)             
# snpa(GTHD, FTHD, 9)             

GTCT = t(GTCT)
GTHD = t(GTHD)
GTGT = t(GTGT)

### emmax analyses including GRM as a random effect
cat("## Performing EMMAX analyes\n\n") 
emx_ct <- emmax(X = GTCT, Y = FTCT$OFFSPR, K = GRMCT)     
emx_hd <- emmax(X = GTHD, Y = FTHD$OFFSPR, K = GRMHD)     

## FUNCTION plotting ###########################################################
plotting <- function(emx, gt, ft, lm = "no", beta = "no", r2 = "no") {
  par(mfrow = c(2, 2))
  plot(-log10(emx$pval), ylab = paste0("-log10_", deparse(substitute(emx))), main = "Emmax P-values")
  
  pj <- p.adjust(emx$pval,"fdr")
  plot(-log10(pj), ylab = paste0("adj-log10_", deparse(substitute(emx))), main = "Adjusted pval 'fdr'")
  
  qqplot(-log10(ppoints(emx$pval)),-log10(emx$pval), ylab = paste0(deparse(substitute(emx))), main = "Against Uniform Distribution")
  abline(0,1)
  
  if (any(pj <= 0.05))  {
  assign(paste0("adj_", deparse(substitute(emx))), pj, envir = .GlobalEnv)
  qqplot(-log10(ppoints(pj)),-log10(pj), ylab = paste0(deparse(substitute(emx))), main = "Against Uniform Distribution (adj_pval)")
  abline(0,1)
  }
  #pjss <- as.matrix(pj[-log10(pj) > 1.3])
  #assign(paste0("adj_up_05_", deparse(substitute(emx))), pjss, envir = .GlobalEnv)
  #pjmax <- as.matrix(sort(-log10(pj), decreasing = TRUE)[1:10])
  #assign(paste0("adj_max10_", deparse(substitute(emx))), pjmax, envir = .GlobalEnv)
  
  if (lm == "yes") {
  p_lm <- apply(gt,2,function(x) cor.test(x,ft[,OFFSPR])$p.value)   #-> -> -> for ENV dataset substitute "OFFSPR" with "2"
  plot(-log10(p_lm),-log10(emx$pval), ylab = deparse(substitute(emx)), main = "Linear model vs Emmax model")
  # cor.test(p_lm, emx$pval)  
  }
  
  if (beta == "yes") {
  beta_lm <- apply(gt,2,function(x) lm(ft[,OFFSPR]~x)$coefficients[2])
  plot(abs(beta_lm),-log10(p_lm), ylab = paste0("lm_", deparse(substitute(gt)), "|", deparse(substitute(ft))), main = "Beta LM")
  assign(paste0("Beta_", deparse(substitute(ft))), beta_lm, envir = .GlobalEnv)
  }
  
  if (r2 == "yes")  {
  r2_lm <- apply(gt,2,function(x) cor(x,ft[,OFFSPR])^2)
  plot(r2_lm,-log10(p_lm), ylab = paste0("lm_", deparse(substitute(gt)), "|", deparse(substitute(ft))), main = "r2 LM")
  }
  
  # qqplot(-log10(ppoints(emx$pval)),-log10(p_lm))
  # abline(0,1)
}
# ##############################################################################

cat("## Building plots\n\n")
plotting(emx_ct, GTCT, FTCT, lm = "yes", beta = "no", r2 = "no")
plotting(emx_hd, GTHD, FTHD, lm = "yes")

# SAVE DATASET and FUNCTION for next anaysis
save(GTGT, GRMGT, ENV, FTFT, file = "G1_A1_MERGED_DATASETS.RData")
save(GTCT, FTCT, GRMCT, emx_ct, GTHD, FTHD, emx_hd, GRMHD, file = "G1_A1_DATASETS.RData")
save(gc, forch, grmx, snpa, plotting, file = "G1_A1_FUNCTIONS.R")

