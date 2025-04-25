#######################
  ### G21 -- ANALYSIS
#######################
  # - - BASEMENT - - 
      ## 1 ## Building GENOTYPE MATRIX GTG21
      ## 2 ## Building FITNESS DATASET
      ## 3 ## Building corrected GENOTYPE MATRIXES
      ## 4 ## Building ENVIRONMENT DATASET (0 = CT, 1 = HD)
      ## 5 ## Building GRM - Genetic relationship matrix
  # - - ANALYSIS  - -  
      ## 6 ## Building Emmax Linear Model
      ## 7 ## RESULTS - adjpvgn - plotting - mantx - annot -genet
              # -> -> # ..ct      = CT(fitness)strict, CT_CT            HIDDEN
                      #  ..hd     = HD(fitness)strict, HD_HD            HIDDEN
                      #  ..ctA    = CT(fitness)all, CT_xyz
                      #  ..hdA    = HD(fitness)all, HD_xyz
                      #  ENV..    = ENV(HD~CT)all, HD_xyz vs CT_xyz
                      #  FT..     = FT(HD~CT)all, HD_xyz vs CT_xyz
                      #  ENV..z   = ENV(HD~CT)strict, HD_HD vs CT_CT 
                      #  FT..z    = FT(HD~CT)strict, HD_HD vs CT_CT

# setwd("C:/Users/gfriso/RStudioNEW/tribolium")
setwd("//ad.helsinki.fi/home/g/gfriso/Desktop/RStudioNEW/tribolium")

library(data.table)  
library(SNPRelate)
library(ggplot2)
library(qqman)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(cowplot)
library(ggrepel) 
library(gridExtra)

source("./R/functions.R")
load("./OUTPUTs/G1_A1_FUNCTIONS.R")
load("./OUTPUTs/G1_A2_FUNCTIONS.R")
load("./OUTPUTs/plotting2.R")
# load("./OUTPUTs/G1_A1_DATASETS.RData")
# load("./OUTPUTs/G1_A1_MERGED_DATASETS.RData")
load("./OUTPUTs/G21_CTHD.RData")    #raw vcf file of all CT and HD samples G21 (+ID : last column) previously shortened (G1 data deleted) to # - SAVE SPACE

G21FTR <- read.csv("./data/G21-fitness+count.csv.gz") # raw fitness data G21

SNP_dataset <- readRDS("./data/SNP_dataset.rds")
GENE_dataset <- readRDS("./data/gene_dataset.rds")

# - SAVE SPACE
# env <- new.env()
# load("./OUTPUTs/G1_A1_MERGED_DATASETS.RData", envir = env)
# GTGT <- env$GTGT
# G21_vcf <- fread("./data/Tribolium_castaneum_ALL_Tcas3.30_imputed.vcf.gz")
# G21_vcf_filtered <- G21_vcf[, !colnames(G21_vcf) %in% rownames(GTGT), with = FALSE]
# G21 <- G21_vcf_filtered[, grepl("HD|CT", colnames(G21_vcf_filtered)) & grepl("L|M", colnames(G21_vcf_filtered)), with = FALSE]
# G21 <- G21 %>%
#   mutate(ID = G21_vcf$ID)
# save(G21, file = "G21_CTHD.RData")

# - - BASEMENT - - - --- - - - --- - - - --- - - - --- - - - --- - - - --- - - -

## 1 ## Building GENOTYPE MATRIX GTG21
G21T <- as.matrix(G21[, grepl(".*HD-.*|.*CT-.*", colnames(G21)), with = FALSE])
rownames(G21T) <- G21$ID
GTG21r <- sub(":.*", "", G21T)
GTG21r[GTG21r %in% c("1|0", "0|1")] <- 1
GTG21r[GTG21r == "0|0"] <- 0
GTG21r[GTG21r == "1|1"] <- 2
GTG21 <- matrix(as.numeric(GTG21r), nrow = nrow(GTG21r), ncol = ncol(GTG21r))   
colnames(GTG21) <- colnames(GTG21r)
rownames(GTG21) <- rownames(GTG21r)        #158 of 159(ft.dataset) but not full correspondence (true = 154 of 159)

## 2 ## Building FITNESS DATASETS
G21FT <- data.frame(ID = colnames(G21FTR)[-1], OFFSPR = as.numeric(t(G21FTR[1,-1])))
G21FT <- G21FT %>%
  mutate(ID = sub("\\.", "-", ID)) %>%
  mutate(ID = sub("\\.", "-", ID)) %>%
  filter(grepl(".*HD-|.*CT-", ID)) %>%
  filter(ID %in% colnames(GTG21))                   #FT-ALL = 154 of 159

G21FTct <- G21FT %>%
  filter(grepl(".*CT-", ID ) & grepl("-.CT", ID))   #CT + CT = 25 of 26

G21FThd <- G21FT %>%
  filter(grepl(".*HD-", ID ) & grepl("-.HD", ID))   #HD + HD = 28 of 28
  
G21FTctA <- G21FT %>%
  filter(grepl(".*CT-", ID ))                       #CT + CT-H-D-HD = 87 of 91

G21FThdA <- G21FT %>%
  filter(grepl(".*HD-", ID ))                       #HD + CT-H-D-HD = 67 of 68

G21FTz <- as.data.frame(rbind(G21FTct, G21FThd))    #FT-STRICT = 53 of 54

## 3 ## Building corrected GENOTYPE MATRIXES              
GTG21 <- GTG21[, colnames(GTG21) %in% G21FT$ID]                                         
                 
GTG21ct <- GTG21[, colnames(GTG21) %in% G21FTct$ID]

GTG21hd <- GTG21[, colnames(GTG21) %in% G21FThd$ID]

GTG21ctA <- GTG21[, colnames(GTG21) %in% G21FTctA$ID]  

GTG21hdA <- GTG21[, colnames(GTG21) %in% G21FThdA$ID]   

GTG21z <- GTG21[, colnames(GTG21) %in% G21FTz$ID]       

## 4 ## Building ENVIRONMENT DATASET (0 = CT, 1 = HD)
ENVG21 <- data.frame(ID = colnames(GTG21))
ENVG21 <- ENVG21 %>%
  mutate(ENV = ID) %>%
  mutate(ENV = sub(".*CT-.*", "0", ENV)) %>%
  mutate(ENV = sub(".*HD-.*", "1", ENV)) %>%
  mutate(ENV = as.numeric(ENV))

ENVG21z <- data.frame(ID = colnames(GTG21z))
ENVG21z <- ENVG21z %>%
  mutate(ENV = ID) %>%
  mutate(ENV = sub(".*CT-.*", "0", ENV)) %>%
  mutate(ENV = sub(".*HD-.*", "1", ENV)) %>%
  mutate(ENV = as.numeric(ENV))

## 5 ## Building GRM - Genetic relationship matrix
GRMGTG21 <- grmx(GTG21)
GRMGTG21ct <- grmx(GTG21ct)
GRMGTG21hd <- grmx(GTG21hd)
GRMGTG21z <- grmx(GTG21z)
GRMGTG21ctA <- grmx(GTG21ctA)
GRMGTG21hdA <- grmx(GTG21hdA)
 
# - - ANALYSIS - - - - - ------------------ - - - - ----- - - ---------- - - -
GTG21 <- t(GTG21)
GTG21ct <- t(GTG21ct)
GTG21hd <- t(GTG21hd)
GTG21z <- t(GTG21z)
GTG21ctA <- t(GTG21ctA)
GTG21hdA <- t(GTG21hdA)

# # 6 ## Building Emmax Linear Model
emx_G21_ct <- emmax(X = GTG21ct, Y = G21FTct$OFFS, K = GRMGTG21ct )
emx_G21_hd <- emmax(X = GTG21hd, Y = G21FThd$OFFS, K = GRMGTG21hd )
emx_G21_env <- emmax(X = GTG21, Y = ENVG21$ENV, K = GRMGTG21 )
emx_G21_ft <- emmax(X = GTG21, Y = G21FT$OFFS, K = GRMGTG21 )
emx_G21_ctA <- emmax(X = GTG21ctA, Y = G21FTctA$OFFS, K = GRMGTG21ctA )
emx_G21_hdA <- emmax(X = GTG21hdA, Y = G21FThdA$OFFS, K = GRMGTG21hdA )
emx_G21_envz <- emmax(X = GTG21z, Y = ENVG21z$ENV, K = GRMGTG21z )
emx_G21_ftz <- emmax(X = GTG21z, Y = G21FTz$OFFS, K = GRMGTG21z )

## 7 ## RESULTS - CT(fitness)all, HD(fitness)all, ENV(HD~CT)all, FT(HD~CT)all, CT(fitness)strict, HD(fitness)strict, ENV(HD~CT)strict, FT(HD~CT)strict

# ct
adjpvgn(emx_G21_ct, SNP_dataset)
# plotting(emx_G21_ct, GTG21ct, G21FTct)
# mantx(adj_emx_G21_ct, GTG21ct)
# annot(dtfr_adj_emx_G21_ct, SNP_dataset, adj_spv_emx_G21_ct)
# genet(dtfr_adj_emx_G21_ct, GENE_dataset, adj_spv_emx_G21_ct)

# hd
adjpvgn(emx_G21_hd, SNP_dataset)
# plotting(emx_G21_hd, GTG21hd, G21FThd)
# mantx(adj_emx_G21_hd, GTG21hd)
# annot(dtfr_adj_emx_G21_hd, SNP_dataset, adj_spv_emx_G21_hd)
# genet(dtfr_adj_emx_G21_hd, GENE_dataset, adj_spv_emx_G21_hd)

# ENV
adjpvgn(emx_G21_env, SNP_dataset)
plotting2(emx_G21_env, GTG21, ENVG21, lm = "no", beta = "no")
mantx(adj_emx_G21_env, GTG21)
annot(dtfr_adj_emx_G21_env, SNP_dataset, adj_spv_emx_G21_env)
genet(dtfr_adj_emx_G21_env, GENE_dataset, adj_spv_emx_G21_env)

# FT
adjpvgn(emx_G21_ft, SNP_dataset)
plotting2(emx_G21_ft, GTG21, G21FT, lm = "no", beta = "no") #any significant adjusted p-value - stop analysis
# mantx(adj_emx_G21_ft, GTG21)
# annot(dtfr_adj_emx_G21_ft, SNP_dataset, adj_spv_emx_G21_ft)
# genet(dtfr_adj_emx_G21_ft, GENE_dataset, adj_spv_emx_G21_ft)

# ctA
adjpvgn(emx_G21_ctA, SNP_dataset)
plotting2(emx_G21_ctA, GTG21ctA, G21FTctA, lm = "no", beta = "no") #any significant adjusted p-value - stop analysis
# mantx(adj_emx_G21_ctA, GTG21ctA)
# annot(dtfr_adj_emx_G21_ctA, SNP_dataset, adj_spv_emx_G21_ctA)
# genet(dtfr_adj_emx_G21_ctA, GENE_dataset, adj_spv_emx_G21_ctA)

# hdA
adjpvgn(emx_G21_hdA, SNP_dataset)
plotting(emx_G21_hdA, GTG21hdA, G21FThdA, lm = "no", beta = "no") #any significant adjusted p-value - stop analysis
# mantx(adj_emx_G21_hdA, GTG21hdA)
# annot(dtfr_adj_emx_G21_hdA, SNP_dataset, adj_spv_emx_G21_hdA)
# genet(dtfr_adj_emx_G21_hdA, GENE_dataset, adj_spv_emx_G21_hdA)

# ENVz
adjpvgn(emx_G21_envz, SNP_dataset)
plotting2(emx_G21_envz, GTG21z, ENVG21z, lm = "no", beta = "no") 
mantx(adj_emx_G21_envz, GTG21z)
annot(dtfr_adj_emx_G21_envz, SNP_dataset, adj_spv_emx_G21_envz)
genet(dtfr_adj_emx_G21_envz, GENE_dataset, adj_spv_emx_G21_envz)

# FTz
adjpvgn(emx_G21_ftz, SNP_dataset)
plotting2(emx_G21_ftz, GTG21z, G21FTz, lm = "no", beta = "no") #any significant adjusted p-value - stop analysis
# mantx(adj_emx_G21_ftz, GTG21z)
# annot(dtfr_adj_emx_G21_ftz, SNP_datasetz, adj_spv_emx_G21_ftz)
# genet(dtfr_adj_emx_G21_ftz, GENE_dataset, adj_spv_emx_G21_ftz)

##########
  # save(emx_G21_env, ENVG21, adj_spv_emx_G21_env, dtfr_adj_emx_G21_env, SNPsoi_dtfr_adj_emx_G21_env, GoI_dtfr_adj_emx_G21_env, G21FT, emx_G21_ft, GTG21, GRMGTG21, file = "G21_ENV_FT_Results.RData")
 # save(emx_G21_envz, ENVG21z, SNPsoi_dtfr_adj_emx_G21_envz, GoI_dtfr_adj_emx_G21_envz, G21FTz, emx_G21_ftz, GTG21z, GRMGTG21z, file = "G21_ENV_FT_zeta_Results.RData")
# save(emx_G21_ctA, GTG21ctA, G21FTctA, emx_G21_hdA, GTG21hdA, GRMGTG21ctA, GRMGTG21hdA, file = "G21_CTa_HDa_Results.RData")

adjpvgnB <- function(emx, ann_dataset) {
  adj <- p.adjust(emx$pval, "bonferroni")
  n <- length(setdiff(names(adj), rownames(ann_dataset)))
  cat("\nThe number of SNPs not annotated from ", deparse(substitute(emx)), " are:\n", n)
  adj1 <- adj[names(adj) %in% rownames(ann_dataset)]
  pvgn <- ann_dataset[rownames(ann_dataset) %in% names(adj1), c("gene", "ID")]
  adj1_ordered <- adj1[match(rownames(pvgn), names(adj1))]
  pvgn <- pvgn %>%
    mutate(pval = adj1_ordered) %>%
    filter(-log10(pval) > 1.3)
  nn <- length(adj[-log10(adj) > 1.3])
  nnn <- nrow(pvgn)
  cat("\nThe number of total significant SNPs [a] and of annotated significant SNPs [b] from ", deparse(substitute(emx)), " are:\n", "[a] ", nn, "\n[b] ", nnn, "\n")
  assign(paste0("adj_", deparse(substitute(emx))), adj1, envir = .GlobalEnv)
  assign(paste0("adj_spv_", deparse(substitute(emx))), pvgn, envir = .GlobalEnv)
}

adjpvgnH <- function(emx, ann_dataset) {
  adj <- p.adjust(emx$pval, "holm")
  n <- length(setdiff(names(adj), rownames(ann_dataset)))
  cat("\nThe number of SNPs not annotated from ", deparse(substitute(emx)), " are:\n", n)
  adj1 <- adj[names(adj) %in% rownames(ann_dataset)]
  pvgn <- ann_dataset[rownames(ann_dataset) %in% names(adj1), c("gene", "ID")]
  adj1_ordered <- adj1[match(rownames(pvgn), names(adj1))]
  pvgn <- pvgn %>%
    mutate(pval = adj1_ordered) %>%
    filter(-log10(pval) > 1.3)
  nn <- length(adj[-log10(adj) > 1.3])
  nnn <- nrow(pvgn)
  cat("\nThe number of total significant SNPs [a] and of annotated significant SNPs [b] from ", deparse(substitute(emx)), " are:\n", "[a] ", nn, "\n[b] ", nnn, "\n")
  assign(paste0("adj_", deparse(substitute(emx))), adj1, envir = .GlobalEnv)
  assign(paste0("adj_spv_", deparse(substitute(emx))), pvgn, envir = .GlobalEnv)
}

adjpvgnBH <- function(emx, ann_dataset) {
  adj <- p.adjust(emx$pval, "BH")
  n <- length(setdiff(names(adj), rownames(ann_dataset)))
  cat("\nThe number of SNPs not annotated from ", deparse(substitute(emx)), " are:\n", n)
  adj1 <- adj[names(adj) %in% rownames(ann_dataset)]
  pvgn <- ann_dataset[rownames(ann_dataset) %in% names(adj1), c("gene", "ID")]
  adj1_ordered <- adj1[match(rownames(pvgn), names(adj1))]
  pvgn <- pvgn %>%
    mutate(pval = adj1_ordered) %>%
    filter(-log10(pval) > 1.3)
  nn <- length(adj[-log10(adj) > 1.3])
  nnn <- nrow(pvgn)
  cat("\nThe number of total significant SNPs [a] and of annotated significant SNPs [b] from ", deparse(substitute(emx)), " are:\n", "[a] ", nn, "\n[b] ", nnn, "\n")
  assign(paste0("adj_", deparse(substitute(emx))), adj1, envir = .GlobalEnv)
  assign(paste0("adj_spv_", deparse(substitute(emx))), pvgn, envir = .GlobalEnv)
}

adjpvgn(emx_G21_env, SNP_dataset)
adjpvgnB(emx_G21_env, SNP_dataset)
adjpvgnH(emx_G21_env, SNP_dataset)
adjpvgnBH(emx_G21_env, SNP_dataset)

adjpvgn(emx_G21_ft, SNP_dataset)
adjpvgnB(emx_G21_ft, SNP_dataset)
adjpvgnH(emx_G21_ft, SNP_dataset)
adjpvgnBH(emx_G21_ft, SNP_dataset)


test <- intersect(adj_spv_emx_G21_env$gene, adj_spv_emx_G21_envz$gene)
write.table(test, file = "intersect_18.txt", row.names = F, col.names = F, quote = F)

# 
# sets <- list(sg_EN22 = sg_EN22, PC_G1 = sig.fit2.G1, EC_hd = sig.fit2.EC_hd)
# venn.plot <- venn.diagram(
#   x = sets,
#   filename = NULL, 
#   fill = c("blue", "darkred", "red"),
#   alpha = 0.5,
#   cex = 1.5,
#   cat.cex = 1.2,
#   main = "GROUPS"
# )

