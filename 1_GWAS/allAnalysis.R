####
#####
###### GWAS ANALYSIS ON GENERATION 1 AND 22, ON AND BETWEEN CONTROL AND HOT-DRY CONDITIONS

# 1- Library source
# 2- Datasets import
# 3- Function loading

# A- G1 ANALYSIS ---------------------------------------------------------------
## A1- Building fitness datasets
## A2 Building genotype matrixes 
## A3- Formatting fitness and genotype datasets
## A4- Building total environment, total genotype, and total fitness datasets
## A5- Building a genetic relationship matrix
## A6- Building emmax linear regression model
## A7- Building adjusted emmax p-values vector and significant p-values dataframe
## A8- Building significant genes table for G1-Environment and G1-Fitness
## A9- Performing a Fisher´s exact test on the gene category specified for G1-Environment and G1-Fitness
## A10- Performing a permutation test for the category specified on significant genes of G1-Environment and G1-Fitness



# B- G22 ANALYSIS --------------------------------------------------------------
## B1- Filtering genotype data
## B2- Building total genotype matrix
## B3- Building fitness dataset
## B4- Formatting total genotype matrix and building CT and HD genotype matrixes  
## B5- Building environment dataset
## B6- Building genetic relationship matrixes
## B7- Building emmax linear regression model
## B8- Building adjusted emmax p-values vector and significant p-values dataframe
## B9- Building significant genes table for G22-Environment
## B10- Performing a Fisher´s exact test on the gene category specified of significant genes for G22-Environment
## B11- Performing a permutation test for the category specified on significant genes of G22-Environment





# setwd("//ad.helsinki.fi/home/g/gfriso/Desktop/RStudioNEW/tribolium")


# 1- Library source
library(data.table)  

# 2- Datasets import
ftct <- fread("./data/Counts.cpm-fit_CT244.txt.gz") # fitness data G1-CT
fthd <- fread("./data/Counts.cpm-fit_HD242.txt.gz") # fitness data G1-HD

gtct <- fread("./data/Tribolium_castaneum_CT_G1_Tcas3.30_imputed.vcf.gz")   # genotype data G1-CT
gthd <- fread("./data/Tribolium_castaneum_HD_G1_Tcas3.30_imputed.vcf.gz")   # genotype data G1-HD

G21FTR <- read.csv("./data/G21-fitness+count.csv.gz") # fitness data G22-ALL

G21_vcf <- fread("./data/Tribolium_castaneum_ALL_Tcas3.30_imputed.vcf.gz")  # genotype data G22-ALL

SNP_dataset <- readRDS("./data/SNP_dataset.rds")    # SNPs annotated dataset T. castaneum
GENE_dataset <- readRDS("./data/gene_dataset.rds")  # Genes annotated dataset T. castaneum

# 3- Function loading [ emmax(), gc(), grmx(), adjpvgn(), ftest() ]
source("./allFunctions.R")  

# A- G1 ANALYSIS ------------------------------------------------------------------------------------------------------------------------------------------------------

## A1- Building fitness datasets
FTCT <- ftct[,1:2] # G1-fitness-CT
FTCT <- FTCT[-47,] # removing the missing sample of the genotype data
  
FTHD <- fthd[,1:2] # G1-fitness-HD

## A2- Building genotype matrixes 
GTCT <- gc(gtct)   # G1-genotype-CT

GTHD <- gc(gthd)   # G1-genotype-HD

## A3- Formatting fitness and genotype datasets
GTCT <- GTCT[, unlist(FTCT[,1])]        # G1-genotype-CT
FTCT$OFFSPR <- as.numeric(FTCT$OFFSPR)  # G1-fitness-CT

GTHD <- GTHD[, unlist(FTHD[,1])]        # G1-genotype-HD
FTHD$OFFSPR <- as.numeric(FTHD$OFFSPR)  # G1-fitness-HD

## A4- Building total environment, total genotype, and total fitness datasets
enCT <- as.data.frame(FTCT)
enCT$OFFSPR <- as.numeric(0) 
colnames(enCT)[colnames(enCT) == "OFFSPR"] <- "ENVIRONMENT"
enHD <- as.data.frame(FTHD)
enHD$OFFSPR <- as.numeric(1)
colnames(enHD)[colnames(enHD) == "OFFSPR"] <- "ENVIRONMENT"
ENV <- rbind(enCT, enHD)    # G1-Environment

dgtCT <- as.data.frame(t(GTCT))
dgtHD <- as.data.frame(t(GTHD))
GTGT <- rbind(dgtCT, dgtHD) # G1-Genotype-ALL

dfCT <- as.data.frame(FTCT)
dfHD <- as.data.frame(FTHD)
FTFT <- rbind(dfCT, dfHD)   # G1-Fitness-ALL

## A5- Building genetic relationship matrixes
GRMCT <- grmx(GTCT) # G1-grm-CT

GRMHD <- grmx(GTHD) # G1-grm-HD      

GTGT <- as.matrix(t(GTGT))
GRMGT <- grmx(GTGT) # G1-grm-ALL

## A6- Building emmax linear regression model
GTCT = t(GTCT)
GTHD = t(GTHD)
GTGT = t(GTGT)

emx_ct <- emmax(X = GTCT, Y = FTCT$OFFSPR, K = GRMCT)     # G1-emmax-CT      testing correlation between genotype and fitness

emx_hd <- emmax(X = GTHD, Y = FTHD$OFFSPR, K = GRMHD)     # G1-emmax-HD      testing correlation between genotype and fitness

emx_EN <- emmax(X = GTGT, Y = ENV$ENVIRONMENT, K = GRMGT) # G1-emmax-Environment       testing correlation between genotype and environment

emx_FT <- emmax(X = GTGT, Y = FTFT$OFFSPR, K = GRMGT)     # G1-emmax-Fitness              testing correlation between genotype and fitness

## A7- Building adjusted emmax p-values vector and significant p-values dataframe
adjpvgn(emx_ct, SNP_dataset) # G1-adj-CT  

adjpvgn(emx_hd, SNP_dataset) # G1-adj-HD

adjpvgn(emx_EN, SNP_dataset) # G1-adj-Environment  ---256  -> 100

adjpvgn(emx_FT, SNP_dataset) # G1-adj-Fitness      ---19  -> 13

## A8- Building significant genes table for G1-Environment and G1-Fitness
sg_genes_EN <- GENE_dataset %>%
  filter(GENE_dataset$gene %in% adj_spv_emx_EN$gene)  # G1-sg_genes-Environment

sg_genes_FT <- GENE_dataset %>%
  filter(GENE_dataset$gene %in% adj_spv_emx_FT$gene)  # G1-sg_genes-Fitness

print(rownames(sg_genes_FT) %in% rownames(sg_genes_EN))

## A9- Performing a Fisher´s exact test on the gene category specified for G1-Environment and G1-Fitness
ftest(GENE_dataset, sg_genes_EN, "hub_gene")          # G1-sg_genes-Environment
ftest(GENE_dataset, sg_genes_EN, "highest_sg_pool")   # G1-sg_genes-Environment
ftest(GENE_dataset, sg_genes_EN, "highest_de_pool")   # G1-sg_genes-Environment  ++++ 6-94, 100-11116, 0.00034
ftest(GENE_dataset, sg_genes_EN, "eQTL_carrier")      # G1-sg_genes-Environment

ftest(GENE_dataset, sg_genes_FT, "hub_gene")          # G1-sg_genes-Fitness
ftest(GENE_dataset, sg_genes_FT, "highest_sg_pool")   # G1-sg_genes-Fitness
ftest(GENE_dataset, sg_genes_FT, "highest_de_pool")   # G1-sg_genes-Fitness  + 2-11, 100-11116, 0.0059
ftest(GENE_dataset, sg_genes_FT, "eQTL_carrier")      # G1-sg_genes-Fitness

## A10- Performing a permutation test for the category specified on significant genes of G1-Environment and G1-Fitness
permtest(GENE_dataset, sg_genes_EN, "logFC")          # G1-sg_genes-Environment ++++ 0
permtest(GENE_dataset, sg_genes_EN, "net_selection")  # G1-sg_genes-Environment ++ 0.0376

permtest(GENE_dataset, sg_genes_FT, "logFC")          # G1-sg_genes-Fitness  ++++ 0
permtest(GENE_dataset, sg_genes_FT, "net_selection")  # G1-sg_genes-Fitness  +++ 0.0108

# B- G22 ANALYSIS ------------------------------------------------------------------------------------------------------------------------------------------------------

## B1- Filtering genotype data
G21_vcf_filtered <- G21_vcf[, !colnames(G21_vcf) %in% rownames(GTGT), with = FALSE]
G21 <- G21_vcf_filtered[, grepl("HD|CT", colnames(G21_vcf_filtered)) & grepl("L|M", colnames(G21_vcf_filtered)), with = FALSE]
G21 <- G21 %>%
  mutate(ID = G21_vcf$ID) 

## B2- Building total genotype matrix
G21T <- as.matrix(G21[, grepl(".*HD-.*|.*CT-.*", colnames(G21)), with = FALSE])
rownames(G21T) <- G21$ID
GTG21r <- sub(":.*", "", G21T)
GTG21r[GTG21r %in% c("1|0", "0|1")] <- 1
GTG21r[GTG21r == "0|0"] <- 0
GTG21r[GTG21r == "1|1"] <- 2
GTG21 <- matrix(as.numeric(GTG21r), nrow = nrow(GTG21r), ncol = ncol(GTG21r))   
colnames(GTG21) <- colnames(GTG21r)
rownames(GTG21) <- rownames(GTG21r)    # G22-Genotype

## B3- Building fitness dataset
G21FT <- data.frame(ID = colnames(G21FTR)[-1], OFFSPR = as.numeric(t(G21FTR[1,-1])))
G21FT <- G21FT %>%
  mutate(ID = sub("\\.", "-", ID)) %>%
  mutate(ID = sub("\\.", "-", ID)) %>%
  filter(grepl(".*HD-|.*CT-", ID)) %>%
  filter(ID %in% colnames(GTG21))      # G22-Fitness

G21FTctA <- G21FT %>%
  filter(grepl(".*CT-", ID ))          # G22-fitness-CT

G21FThdA <- G21FT %>%
  filter(grepl(".*HD-", ID ))          # G22-fitness-HD

## B4- Formatting total genotype matrix and building CT and HD genotype matrixes           
GTG21 <- GTG21[, colnames(GTG21) %in% G21FT$ID]          # G22-Genotype                               

GTG21ctA <- GTG21[, colnames(GTG21) %in% G21FTctA$ID]    # G22-Genotype-CT

GTG21hdA <- GTG21[, colnames(GTG21) %in% G21FThdA$ID]    # G22-Genotype-HD

## B5- Building environment dataset
ENVG21 <- data.frame(ID = colnames(GTG21))
ENVG21 <- ENVG21 %>%
  mutate(ENV = ID) %>%
  mutate(ENV = sub(".*CT-.*", "0", ENV)) %>%
  mutate(ENV = sub(".*HD-.*", "1", ENV)) %>%
  mutate(ENV = as.numeric(ENV))                # G22-Environment

## B6- Building Building genetic relationship matrixes
GRMGTG21 <- grmx(GTG21)        # G22-grm-ALL

GRMGTG21ctA <- grmx(GTG21ctA)  # G22-grm-CT

GRMGTG21hdA <- grmx(GTG21hdA)  # G22-grm-HD

## B7- Building emmax linear regression model
GTG21 <- t(GTG21)
GTG21ctA <- t(GTG21ctA)
GTG21hdA <- t(GTG21hdA)

emx_G21_ctA <- emmax(X = GTG21ctA, Y = G21FTctA$OFFS, K = GRMGTG21ctA )   # G22-emmax-CT      testing correlation between genotype and fitness

emx_G21_hdA <- emmax(X = GTG21hdA, Y = G21FThdA$OFFS, K = GRMGTG21hdA )   # G22-emmax-HD      testing correlation between genotype and fitness

emx_G21_env <- emmax(X = GTG21, Y = ENVG21$ENV, K = GRMGTG21 )            # G22-emmax-Environment        testing correlation between genotype and environment

emx_G21_ft <- emmax(X = GTG21, Y = G21FT$OFFS, K = GRMGTG21 )             # G22-emmax-Fitness            testing correlation between genotype and fitness

## B8- Building adjusted emmax p-values vector and significant p-values dataframe
adjpvgn(emx_G21_ctA, SNP_dataset)  # G22-adj-CT  

adjpvgn(emx_G21_hdA, SNP_dataset)  # G22-adj-HD

adjpvgn(emx_G21_env, SNP_dataset)  # G22-adj-Environment  ---1699 -> 902

adjpvgn(emx_G21_ft, SNP_dataset)   # G22-adj-Fitness

## B9- Building significant genes table for G22-Environment
sg_genes_EN22 <- GENE_dataset %>%
  filter(GENE_dataset$gene %in% adj_spv_emx_G21_env$gene)  #G22-sg_genes-Environment

## B10- Performing a Fisher´s exact test on the gene category specified for G22-Environment
ftest(GENE_dataset, sg_genes_EN22, "hub_gene")        # ++++ 40-862, 259-10957, 0.00020
ftest(GENE_dataset, sg_genes_EN22, "highest_sg_pool")
ftest(GENE_dataset, sg_genes_EN22, "highest_de_pool")
ftest(GENE_dataset, sg_genes_EN22, "eQTL_carrier")    # +++++ 74-828, 564-10652, 7.604e-05

## B11- Performing a permutation test for the category specified on significant genes of G22-Environment
permtest(GENE_dataset, sg_genes_EN22, "logFC")          
permtest(GENE_dataset, sg_genes_EN22, "net_selection")  # +++ 0.0016


