## DIFFERENTIAL EXPRESSION ANALYSIS on G1 and G21 with Limma

# PLASTIC GENES G1  =  DE (HD vs CT)
# MOST PLASTIC GENES G1  =  PLASTIC GENES G1 with logFC > 1
# EVOLVED GENES G21  =  DE (CT.HD vs HD.HD)
# PLASTIC GENES G21  =  DE (HD.CT vs HD.HD)
# MOST PLASTIC GENES G21  =  MOST PLASTIC GENES G21 with logFC > 1
# GWAS GENES  =  GENES of SNPs from GWAS (CT vs HD - G21)
# LOST PLASTICITY GENES  =  PLASTIC in G1 - NO PLASTIC in G21
# KEEP PALSTICITY GENES  =  PLASTIC in G1 - PLASTIC in G21
# NEW PLASTIC GENES  =  NO PLASTIC in G1 - PLASTIC in G21
# PLASTIC GENES EVOLVED  =  PLASTIC in G1 - EVOLVED GENES G21
# NEW PLASTIC GENES EVOLVED  =  NO PLASTIC in G1 - PLASTIC in G21 - EVOLVED GENES G21

# # # DE Limma G1  + SELECTION - - - - - - - - - - - - - - - - - - - - - - - - - 
# # # DE Limma G21 + SELECTION - - - - - - - - - - - - - - - - - - - - - - - - - 
# # # DE TABLE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # # RESULT NUMBERS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# # # TESTs - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - 
# # # PLOT - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


setwd("//ad.helsinki.fi/home/g/gfriso/Desktop/RStudioNEW/tribolium")

library(data.table)  
library(dplyr)
library(tidyr)
library(edgeR)
library(limma)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(VennDiagram)
library(grid)
library(cowplot)
library(gghalves)


source("./allFunctions.R")  

# GENE DATASET EVA
GENE_dataset <- readRDS("./data/gene_dataset.rds")

# GWAS GENES
sg_genes_EN22 <- read.table("./data/sg_genes_EN22.txt")

# # # DE Limma G1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#####
offspr.G1 <- read.table("./data/offspr_G1.txt")
offspr.G1.CT <- offspr.G1[grepl("CT$", rownames(offspr.G1)), , drop = F]
offspr.G1.HD <- offspr.G1[grepl("HD$", rownames(offspr.G1)), , drop = F]

offspr.G1.CT <- offspr.G1.CT[-47, , drop = F]

count.G1 <- read.table("./data/raw_count_CT_HD_G1.txt",  header = TRUE, check.names = FALSE)

batch <- read.table("./data/BatchEffect.txt", T)
batch <- batch[-52, ]

order <- c()
for (i in 1:dim(batch)[1]) {
  order[i] <- which(batch$Sample==colnames(count.G1)[i])
}

batch <- batch[order,]

# DGEList object  
DGE.count.G1 <- DGEList(counts = count.G1, group = batch$Condition)

# Filtering + Normalization for Limma
keep <- rowSums(cpm(DGE.count.G1) > 1) >= 2
DGE.count.G1 <- DGE.count.G1[keep, , keep.lib.sizes = FALSE]
DGE.count.G1 <- DGEList(DGE.count.G1, group =  batch$Condition)
DGE.count.G1 <- calcNormFactors(DGE.count.G1)                                   # TMM normalization

# Design + Contrast Matrix
design <- model.matrix(~0 + Sequencing + Condition, data = batch)
contrast.matrix <- makeContrasts(ConditionHD, levels = design)

# DE Limma Analisys
G1.norm <- voom(DGE.count.G1, design)
fit.G1 <- lmFit(G1.norm, design)
fit2.G1 <- contrasts.fit(fit.G1, contrast.matrix)
fit2.G1 <- eBayes(fit2.G1)

test.DE_G1 <- decideTests(fit2.G1)
table(test.DE_G1@.Data)

table.fit2.G1 <- topTable(fit2.G1, adjust.method = "BH", number = Inf)
sig.table.fit2.G1 <- table.fit2.G1[table.fit2.G1$adj.P.Val < 0.05 & abs(table.fit2.G1$logFC) > 0.136,] # ~10% up/down

sig.fit2.G1 <- row.names(sig.table.fit2.G1)

# cpm tranformation for Covariance
logCPM <- cpm(DGE.count.G1, log=TRUE, prior.count=5, normalized.lib.sizes = T)
logCPMc <- removeBatchEffect(logCPM, batch$Sequencing)

logCPMc.CT <- logCPMc[, colnames(logCPMc) %in% rownames(offspr.G1.CT)]
logCPMc.HD <- logCPMc[, colnames(logCPMc) %in% rownames(offspr.G1.HD)]

# Genetic Selection Coefficient G1   --  cov(exp.lv, rel.ft) G1
compute_covariance <- function(gene_expression, fitness) {
  relative_fitness <- fitness/mean(fitness, na.rm = TRUE)
  return(cov(gene_expression, relative_fitness, use = "complete.obs")) 
}
# CT
cov_results.CT1 <- apply(logCPMc.CT, 1, compute_covariance, fitness = offspr.G1.CT$OFFSPR)
cov_table.CT1 <- data.frame(Gene = rownames(logCPMc.CT), CovarianceCT = cov_results.CT1)
# HD
cov_results.HD1 <- apply(logCPMc.HD, 1, compute_covariance, fitness = offspr.G1.HD$OFFSPR)
cov_table.HD1 <- data.frame(Gene = rownames(logCPMc.HD), CovarianceHD = cov_results.HD1)

# Genetic Selection Coefficient G1 standardized (effect size)  --  cor(exp.lv, rel.ft) 
compute_correlation <- function(gene_expression, fitness) {
  relative_fitness <- fitness/mean(fitness, na.rm = TRUE)
  return(cor(gene_expression, relative_fitness, use = "complete.obs"))
}
# CT
cor_results.CT1 <- apply(logCPMc.CT, 1, compute_correlation, fitness = offspr.G1.CT$OFFSPR)
cor_table.CT1 <- data.frame(Gene = rownames(logCPMc.CT), Correlation1 = cor_results.CT1)
# HD
cor_results.HD1 <- apply(logCPMc.HD, 1, compute_correlation, fitness = offspr.G1.HD$OFFSPR)
cor_table.HD1 <- data.frame(Gene = rownames(logCPMc.HD), Correlation1 = cor_results.HD1)
#####

# # # DE Limma G21 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#####
offspr.G21 <- as.data.frame(t(read.table("./data/offspring_data_G21.txt")))

data2 <- read.table("./data/Counts_Transpl_2018.txt") 

sample.data0 <- read.table("./data/Samples_Transpl_2018.txt", row.names = 1, header = T) 
sample.data0$Sample2 <- gsub("-", ".", sample.data0$Sample2)
rownames(sample.data0) <- sample.data0$Sample2

data.count <- data2[, colnames(data2) != "Mx1CT.1.4CT"] # remove outlier
sample.data0 <- sample.data0[rownames(sample.data0) != "Mx1CT.1.4CT",]

sample.data <- sample.data0 %>%
  filter(Selection %in% c("CT","HD")) %>%
  filter(Treatment %in% c("CT","HD")) %>%
  mutate(line = ifelse(grepl("L",Sample2), substr(sample.data0$Sample2,1,2), substr(sample.data0$Sample2,1,3))) %>%
  mutate(line = paste(line, Selection, sep=""))

data.count <- data.count[, match(sample.data$Sample2, colnames(data.count))]

identical(colnames(data.count), sample.data$Sample2)
table(sample.data$Sample2 == colnames(data.count))

# Factors
sample.data$group <- factor(sample.data$group, levels=c("CT.CT", "HD.CT", "CT.HD", "HD.HD"))
sample.data$line <- factor(sample.data$line)
sample.data$batch <- factor(sample.data$batch)

# DEGList object 
dge.count <- DGEList(data.count, remove.zeros = TRUE, group = sample.data$group)

# design matrix
design = model.matrix(~0 + batch + group, data = sample.data)

# Filtering + Normalization + logCPM transformation + Limma linear model
keep <- rowSums(cpm(dge.count) > 1) >= 2
dge.count <- dge.count[keep, , keep.lib.sizes = FALSE]
dge.count2 <- DGEList(dge.count, group = sample.data$group)
dge.count2 <- calcNormFactors(dge.count2)                                       # determine the appropriate scale values = normalization TMM

V1     = voom(dge.count2, design)                                               # transorm the normalized data in logCPM
corfit = duplicateCorrelation(V1, design, block = sample.data$line)                 
V2     = voom(dge.count2, design, block = sample.data$line, correlation = corfit$consensus) 
corfit <- duplicateCorrelation(V2, design, block = sample.data$line)

fit <- lmFit(V2, design, block = sample.data$line, correlation = corfit$consensus)

# PLASTIC GENES G21  -  HD.CT_HD.HD 
contrast_PC21 <- makeContrasts(Plasticity = groupHD.CT-groupHD.HD, levels = design)

fit2.PC21 <- contrasts.fit(fit, contrast_PC21)                      
fit2.PC21 <- eBayes(fit2.PC21)

test.DE_PC21 <- decideTests(fit2.PC21, adjust.method = "bonferroni")
table(test.DE_PC21@.Data)

topTable(fit2.PC21, adjust.method = "bonferroni")                                     

table.fit2.PC21 <- topTable(fit2.PC21, adjust.method = "BH", number = Inf)
sig.table.fit2.PC21 <- table.fit2.PC21[table.fit2.PC21$adj.P.Val < 0.05 & abs(table.fit2.PC21$logFC) > 0.136,] # ~10% up/down

sig.fit2.PC21 <- row.names(sig.table.fit2.PC21)

# EVOLVED GENES G21  -  CT.HD-HD.HD 
contrast_EC21 <- makeContrasts(groupHD.HD-groupCT.HD, levels = design)

fit2.EC21 <- contrasts.fit(fit, contrast_EC21)                      
fit2.EC21 <- eBayes(fit2.EC21)

test.DE_EC21 <- decideTests(fit2.EC21)
table(test.DE_EC21@.Data)

topTable(fit2.EC21, adjust.method = "bonferroni")                                     

table.fit2.EC21 <- topTable(fit2.EC21, adjust.method = "BH", number = Inf)
sig.table.fit2.EC21 <- table.fit2.EC21[table.fit2.EC21$adj.P.Val < 0.05 & abs(table.fit2.PC21$logFC) > 0.136,] # ~10% up/down

sig.fit2.EC21 <- row.names(sig.table.fit2.EC21)

# Offspring CT.HD and HD.HD
offspr.G21.CT.HD <- na.omit(offspr.G21[match(rownames(sample.data)[grepl("CT.HD", sample.data$group)], rownames(offspr.G21)), , drop = F ])
offspr.G21.HD.HD <- na.omit(offspr.G21[match(rownames(sample.data)[grepl("HD.HD", sample.data$group)], rownames(offspr.G21)), , drop = F ])

# cpm tranformation for Covariance
logCPM <- cpm(dge.count2, log=TRUE, prior.count=5, normalized.lib.sizes = T)
logCPMc <- removeBatchEffect(logCPM, sample.data$batch)

logCPMc.CT.HD <- logCPMc[, colnames(logCPMc) %in% rownames(offspr.G21.CT.HD)]
logCPMc.HD.HD <- logCPMc[, colnames(logCPMc) %in% rownames(offspr.G21.HD.HD)]

# Genetic Selection Coefficient G21   --  cov(exp.lv, rel.ft) G21    
# CT
cov_results.CT.HD21 <- apply(logCPMc.CT.HD, 1, compute_covariance, fitness = offspr.G21.CT.HD$offspring.num)
cov_table.CT.HD21 <- data.frame(Gene = rownames(logCPMc.CT.HD), CovarianceCT = cov_results.CT.HD21)
# HD
cov_results.HD.HD21 <- apply(logCPMc.HD.HD, 1, compute_covariance, fitness = offspr.G21.HD.HD$offspring.num)
cov_table.HD.HD21 <- data.frame(Gene = rownames(logCPMc.HD.HD), CovarianceHD = cov_results.HD.HD21)

# Genetic Selection Coefficient G21 standardized (effect size)  --  cor(exp.lv, rel.ft)  
# CT
cor_results.CT.HD21 <- apply(logCPMc.CT.HD, 1, compute_correlation, fitness = offspr.G21.CT.HD$offspring.num)
cor_table.CT.HD21 <- data.frame(Gene = rownames(logCPMc.CT.HD), Correlation1 = cor_results.CT.HD21)
# HD
cor_results.HD.HD21 <- apply(logCPMc.HD.HD, 1, compute_correlation, fitness = offspr.G21.HD.HD$offspring.num)
cor_table.HD.HD21 <- data.frame(Gene = rownames(logCPMc.HD.HD), Correlation1 = cor_results.HD.HD21)

#####

# # # DE TABLE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#####
# MOST PLASTIC GENES G1
mostP_G1 <- filter(sig.table.fit2.G1, abs(sig.table.fit2.G1$logFC) > 1) # = 100% 
mostP_G1.50 <- filter(sig.table.fit2.G1, abs(sig.table.fit2.G1$logFC) > 0.585) # = 50%

# MOST PLASTIC GENES G21
mostPC21_G21 <- filter(sig.table.fit2.PC21, abs(sig.table.fit2.PC21$logFC) > 1) # = 100%
mostPC21_G21.50 <- filter(sig.table.fit2.PC21, abs(sig.table.fit2.PC21$logFC) > 0.585) # = 50%

# DE_table
DE_row <- unique(c(rownames(table.fit2.G1), rownames(table.fit2.PC21)))

DE_table <- data.frame(
  row.names = DE_row,
  PC1 = as.integer(DE_row %in% rownames(sig.table.fit2.G1)),    # plastic genes G1
  PC21 = as.integer(DE_row %in% rownames(sig.table.fit2.PC21)), # plastic genes G21
  EC21 = as.integer(DE_row %in% rownames(sig.table.fit2.EC21))  # evolved genes G21
)

DE_table <- DE_table %>%
  mutate(genes = rownames(DE_table)) %>%
  # Plastic genes G1: lost, keep, new
  mutate(LSPC = ifelse(PC1 == 1 & PC21 == 0, 1, 0)) %>%
  mutate(PRPC = ifelse(PC1 + PC21 == 2, 1, 0)) %>%
  mutate(NWPC = ifelse(PC1 == 0 & PC21 == 1, 1, 0)) %>%
  mutate(G1.logFC_PC1 = table.fit2.G1$logFC[match(DE_row, rownames(table.fit2.G1))]) %>%
  mutate(G1.AveExpr_PC1 = table.fit2.G1$AveExpr[match(DE_row, rownames(table.fit2.G1))]) %>%
  mutate(G1.up_down_PC1 = sign(G1.logFC_PC1)) %>%
  
  # Most plastic genes G1
  mutate(MPC1 = ifelse(DE_row %in% rownames(mostP_G1), 1, 0)) %>%
  mutate(MPC1.50 = ifelse(DE_row %in% rownames(mostP_G1.50), 1, 0)) %>%
  
  # G1-Plastic evolved genes
  mutate(PCEC = ifelse(PC1 + EC21 == 2, 1, 0)) %>%
  mutate(PCECmost = ifelse(MPC1 + EC21 == 2, 1, 0)) %>%
  mutate(PCECmost.50 = ifelse(MPC1.50 + EC21 == 2, 1, 0)) %>%
  
  # Plastic genes G21
  mutate(logFC_PC21 = -table.fit2.PC21$logFC[match(DE_row, rownames(table.fit2.PC21))]) %>%
  mutate(AveExpr_PC21 = table.fit2.PC21$AveExpr[match(DE_row, rownames(table.fit2.PC21))]) %>%
  mutate(up_down_PC21 = sign(logFC_PC21)) %>%
  
  # Most plastic genes G21
  mutate(PC21most = ifelse(DE_row %in% rownames(mostPC21_G21), 1, 0)) %>%
  mutate(PC21most.50 = ifelse(DE_row %in% rownames(mostPC21_G21.50), 1, 0)) %>%
  
  # New Plastic Evolved genes G21
  mutate(NWPCEC = ifelse(NWPC == 1 & EC21 == 1, 1, 0)) %>%
  
  # Evolved genes G21
  mutate(logFC_EC21 = table.fit2.EC21$logFC[match(DE_row, rownames(table.fit2.EC21))]) %>%
  mutate(up_down_EC21 = sign(logFC_EC21)) %>%

  # Up_Down : PC1 vs EC_G21, 
  mutate(Contrast = G1.up_down_PC1 + up_down_PC21) %>%
  
  # comparing the direction of PC1 responses to the direction of EC21
  mutate(PC_resp = sign(G1.logFC_PC1) - sign(logFC_PC21)) %>%
  
  # GENE_dataset & sg_genes_EN22(GWAS)
  mutate(module_gene = GENE_dataset$module_gene[match(DE_row, rownames(GENE_dataset))]) %>%
  mutate(hub_gene = GENE_dataset$hub_gene[match(DE_row, rownames(GENE_dataset))]) %>%
  mutate(direct_selection = GENE_dataset$direct_selection[match(DE_row, rownames(GENE_dataset))]) %>%
  mutate(diff_selection = GENE_dataset$diff_selection[match(DE_row, rownames(GENE_dataset))]) %>%
  mutate(eQTL_carrier = GENE_dataset$eQTL_carrier[match(DE_row, rownames(GENE_dataset))]) %>%
  mutate(Pleiotropy = GENE_dataset$Pleiotropy[match(DE_row, rownames(GENE_dataset))]) %>%
  mutate(SNPsEN22 = ifelse(DE_row %in% rownames(sg_genes_EN22), 1, 0)) %>%
  
  # Selection G1 & G21 + Standardized SelG1 & SelG21 (cor)
  mutate(sel.cov.CT.1 = cov_table.CT1$CovarianceCT[match(DE_row, rownames(cov_table.CT1))]) %>%
  mutate(sel.cov.HD.1 = cov_table.HD1$CovarianceHD[match(DE_row, rownames(cov_table.HD1))]) %>%
  mutate(sel.cov.CT.HD.21 = cov_table.CT.HD21$CovarianceCT[match(DE_row, rownames(cov_table.CT.HD21))]) %>%
  mutate(sel.cov.HD.HD.21 = cov_table.HD.HD21$CovarianceHD[match(DE_row, rownames(cov_table.HD.HD21))]) %>%
  mutate(sign_S1 = sign(sel.cov.HD.1)) %>%
  mutate(sign_S21 = sign(sel.cov.HD.HD.21)) %>%
  mutate(contrastSEL = sign(sel.cov.HD.1) - sign(sel.cov.HD.HD.21)) %>%
  mutate(contrastPS1 = G1.up_down_PC1 + sign_S1) %>%
  mutate(maladaptive = ifelse(contrastPS1 == 0, 1, 0)) %>%
  mutate(adaptive = ifelse(contrastPS1 != 0, 1, 0)) %>%
  mutate(contrastPS21 = G1.up_down_PC1 + sign_S21) %>%
  mutate(contrastP21S1 = up_down_PC21 + sign_S1) %>%
  mutate(contrastP21S21 = up_down_PC21 + sign_S21) %>%
  mutate(sel.cor.CT.1 = cor_table.CT1$Correlation1[match(DE_row, rownames(cor_table.CT1))]) %>%
  mutate(sel.cor.HD.1 = cor_table.HD1$Correlation1[match(DE_row, rownames(cor_table.HD1))]) %>%
  mutate(sel.cor.CT.HD.21 = cor_table.CT.HD21$Correlation1[match(DE_row, rownames(cor_table.CT.HD21))]) %>%
  mutate(sel.cor.HD.HD.21 = cor_table.HD.HD21$Correlation1[match(DE_row, rownames(cor_table.HD.HD21))]) 
#####

# # # RESULT NUMBERS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#####

table(DE_table$PC1)            #- -5037   PLASTIC GENES G1               
table(DE_table$MPC1)           #- -235    MOST PLASTIC GENES G1             
table(DE_table$MPC1.50)        #- -766    MOST.50 PLASTIC GENES G1              
table(DE_table$LSPC)           #- -3398   LOST PLASTICITY GENES            
table(DE_table$NWPC)           #- -643    NEW PLASTIC GENES                 
table(DE_table$PRPC)           #- -1639   KEEP PLASTICITY GENES            
table(DE_table$PCEC)           #- -336    PLASTIC GENES EVOLVED             
table(DE_table$PCECmost)       #- -14     MOST PLASTIC GENES EVOLVED           
table(DE_table$PCECmost.50)    #- -48     MOST.50 PLASTIC GENES EVOLVED        
table(DE_table$EC21)           #- 731     EVOLVED GENES                    
table(DE_table$PC21)           #- 2282    PLASTIC GENES G21                
table(DE_table$PC21most)       #- 139     MOST PLASTIC GENES G21                
table(DE_table$PC21most.50)    #- 398     MOST.50 PLASTIC GENES G21  
table(DE_table$NWPCEC)         #- 49      NEW PLASTIC EVOLVED GENES                    
table(DE_table$SNPsEN22)       #- 894     GWAS GENES                    
#

table(DE_table$Contrast[DE_table$PRPC == 1])
table(DE_table$PC_resp[DE_table$PC1 == 1])
table(DE_table$PC_resp[DE_table$EC21 == 1])
table(DE_table$PC_resp[DE_table$EC21 == 1 & DE_table$PC1 == 1])
table(DE_table$PC_resp[DE_table$MPC1.50 == 1])
table(DE_table$PC_resp[DE_table$EC21 == 1 & DE_table$MPC1.50 == 1])
table(DE_table$contrastSEL)
table(DE_table$contrastSEL[DE_table$SNPsEN22 == 1])
table(DE_table$contrastSEL[DE_table$PRPC == 1])
table(DE_table$contrastSEL[DE_table$LSPC == 1])
table(DE_table$contrastSEL[DE_table$NWPC == 1])
table(DE_table$contrastSEL[DE_table$EC21 == 1])
#####

# # # ADAPTIVE - MALADAPTIVE RESPONSE - - - - - - - - - - - - - - - - - 
#####
# Selection intensity of UP - DOWN regulated genes vs Non-Responding genes

sel.resp <- function(resp, noresp, sign = c("up", "down")) {
  D_mean_RESP_NORSP <- mean(resp) - mean(noresp)
  
  n_perm = 10000
  D_perm_RESP <- numeric(n_perm)
  
  rep.RESP_NORSP <- c(rep("sign", length(resp)),
                      rep("norsp", length(noresp)))
  comb.RESP_NORSP <- c(resp, noresp)
  
  for (i in 1:n_perm) {
    perm_labels <- sample(rep.RESP_NORSP)
    
    perm.RESP <- comb.RESP_NORSP[perm_labels == "sign"]
    perm.NORSP <- comb.RESP_NORSP[perm_labels == "norsp"]
    
    D_perm_RESP[i] <- mean(perm.RESP) - mean(perm.NORSP)
  }
  if(sign[1] == "up") {
    pval_RESP <- sum(D_perm_RESP >= D_mean_RESP_NORSP) / n_perm
  } else {
    pval_RESP <- sum(D_perm_RESP <= D_mean_RESP_NORSP) / n_perm
  }
  
  if(pval_RESP < 0.05) {
    cat(deparse(substitute(resp)), " genes are under significantly more ", 
        if(sign[1] == "up") "positive" else "negative", 
        " selection compared to ",  deparse(substitute(noresp)), " genes\n\n")
  } else {
    cat(deparse(substitute(resp)), " genes are NOT under significantly more ", 
        if(sign[1] == "up") "positive" else "negative", 
        " selection compared to ",  deparse(substitute(noresp)), " genes\n\n")
  }
  cat("p-value is", pval_RESP, "\n\n")
  cat("Observed difference in means:", D_mean_RESP_NORSP, "\n\n")
  
  return(pval_RESP)
}

# PC1
obs.sel_PC1 <- DE_table$sel.cov.HD.1[!is.na(DE_table$G1.logFC_PC1) & DE_table$PC1 == 1] #2604
obs.sel_UP_G1 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 > 0 & DE_table$PC1 == 1] 
obs.sel_DOWN_G1 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 < 0 & DE_table$PC1 == 1] 
obs.sel_NORSP_G1<- DE_table$sel.cov.HD.1[!is.na(DE_table$G1.logFC_PC1) & DE_table$PC1 == 0] # 7006

sel.resp(obs.sel_PC1, obs.sel_NORSP_G1, sign = "up")     # significant positive (0)     --------------

sel.resp(obs.sel_UP_G1, obs.sel_NORSP_G1, sign = "up")  # significant positive (0.0089) --------------
sel.resp(obs.sel_DOWN_G1, obs.sel_NORSP_G1, sign = "down") # not significant negative (1)
sel.resp(obs.sel_DOWN_G1, obs.sel_NORSP_G1, sign = "up") # significant negative (0)     --------------
# .M50
obs.sel_UP_G1.M50 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 > 0 & DE_table$MPC1.50 == 1] #2604
obs.sel_DOWN_G1.M50 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 < 0 & DE_table$MPC1.50 == 1] #2433
obs.sel_NORSP_G1.M50<- DE_table$sel.cov.HD.1[!is.na(DE_table$G1.logFC_PC1) & DE_table$MPC1.50 == 0] # 

sel.resp(obs.sel_UP_G1.M50, obs.sel_NORSP_G1.M50, sign = "up")  # significant positive (0)  --------------
sel.resp(obs.sel_DOWN_G1.M50, obs.sel_NORSP_G1.M50, sign = "down") # not significant negative (1)
sel.resp(obs.sel_DOWN_G1.M50, obs.sel_NORSP_G1.M50, sign = "up") # significant negative (0)   ----------------
# LOST
obs.sel_LOST_G1 <- DE_table$sel.cov.HD.1[DE_table$LSPC == 1] 
obs.sel_UPLOST_G1 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 > 0 & DE_table$LSPC == 1] 
obs.sel_DOWNLOST_G1 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 < 0 & DE_table$LSPC == 1] 
obs.sel_NOLOST_G1 <- DE_table$sel.cov.HD.1[!is.na(DE_table$G1.logFC_PC1) & DE_table$LSPC == 0] 

sel.resp(obs.sel_LOST_G1, obs.sel_NOLOST_G1, sign = "up")  # significant positive (0) -------------  
sel.resp(obs.sel_UPLOST_G1, obs.sel_NOLOST_G1, sign = "up") # not significant positive (0.0524)
sel.resp(obs.sel_DOWNLOST_G1, obs.sel_NOLOST_G1, sign = "down") # not significant negative (1)
sel.resp(obs.sel_DOWNLOST_G1, obs.sel_NOLOST_G1, sign = "up") # significant positive (0) -----------

# NEW
obs.sel_NEW_G1 <- na.omit(DE_table$sel.cov.HD.1[DE_table$NWPC == 1])
obs.sel_UPNEW_G1 <- na.omit(DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 > 0 & DE_table$NWPC == 1])
obs.sel_DOWNNEW_G1 <- na.omit(DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 < 0 & DE_table$NWPC == 1])
obs.sel_NONEW_G1 <- na.omit(DE_table$sel.cov.HD.1[!is.na(DE_table$G1.logFC_PC1) & DE_table$NWPC == 0])

sel.resp(obs.sel_NEW_G1, obs.sel_NONEW_G1, sign = "up")  #  not significant positive (0.5) 
sel.resp(obs.sel_UPNEW_G1, obs.sel_NONEW_G1, sign = "up")  # not significant positive (0.8) 
sel.resp(obs.sel_DOWNNEW_G1, obs.sel_NONEW_G1, sign = "down") # not significant negative (0.8)
sel.resp(obs.sel_DOWNNEW_G1, obs.sel_NONEW_G1, sign = "up") # not significant negative (0.17)
# EC
obs.sel_EC_G1 <- DE_table$sel.cov.HD.1[DE_table$EC21 == 1] 
obs.sel_UPEC_G1 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 > 0 & DE_table$EC21 == 1] 
obs.sel_DOWNEC_G1 <- DE_table$sel.cov.HD.1[DE_table$G1.logFC_PC1 < 0 & DE_table$EC21 == 1] 
obs.sel_NOEC_G1 <- DE_table$sel.cov.HD.1[!is.na(DE_table$G1.logFC_PC1) & DE_table$EC21 == 0] 

sel.resp(obs.sel_EC_G1, obs.sel_NOEC_G1, sign = "up")     # significant positive (0.01) ------------
sel.resp(obs.sel_EC_G1, obs.sel_NOEC_G1, sign = "down")   # not significant negative (0.98)

sel.resp(obs.sel_UPEC_G1, obs.sel_NOEC_G1, sign = "up")  # not significant positive (0.7) 
sel.resp(obs.sel_DOWNEC_G1, obs.sel_NOEC_G1, sign = "down") # significant negative (0.008) -----------
### G21 PC
obs.sel_PC21 <- DE_table$sel.cov.HD.HD.21[!is.na(DE_table$logFC_PC21) & DE_table$PC21 == 1] #
obs.sel_UP_PCG21 <- DE_table$sel.cov.HD.HD.21[DE_table$logFC_PC21 > 0 & DE_table$PC21 == 1] 
obs.sel_DOWN_PCG21 <- DE_table$sel.cov.HD.HD.21[DE_table$logFC_PC21 < 0 & DE_table$PC21 == 1] 
obs.sel_NORSP_PCG21<- DE_table$sel.cov.HD.HD.21[!is.na(DE_table$logFC_PC21) & DE_table$PC21 == 0] # 

sel.resp(obs.sel_PC21, obs.sel_NORSP_PCG21, sign = "up")     # NOT significant positive (0.9) 
sel.resp(obs.sel_PC21, obs.sel_NORSP_PCG21, sign = "down")   # not significant negative (0.08)

sel.resp(obs.sel_UP_PCG21, obs.sel_NORSP_PCG21, sign = "down")  # significant negative (0) ------------
sel.resp(obs.sel_DOWN_PCG21, obs.sel_NORSP_PCG21, sign = "up") # significant positive (0) ------------
# EC
obs.sel_EC_G21 <- DE_table$sel.cov.HD.HD.21[DE_table$EC21 == 1] 
obs.sel_UPEC_G21 <- DE_table$sel.cov.HD.HD.21[DE_table$logFC_EC21 > 0 & DE_table$EC21 == 1] 
obs.sel_DOWNEC_G21 <- DE_table$sel.cov.HD.HD.21[DE_table$logFC_EC21 < 0 & DE_table$EC21 == 1] 
obs.sel_NOEC_G21 <- DE_table$sel.cov.HD.HD.21[!is.na(DE_table$logFC_EC21) & DE_table$EC21 == 0] 

sel.resp(obs.sel_EC_G21, obs.sel_NOEC_G21, sign = "up")     # significant positive (0) ------------
sel.resp(obs.sel_EC_G21, obs.sel_NOEC_G21, sign = "down")   # not significant negative (1)

sel.resp(obs.sel_UPEC_G21, obs.sel_NOEC_G21, sign = "up")  # significant positive (0) ----------
sel.resp(obs.sel_DOWNEC_G21, obs.sel_NOEC_G21, sign = "down") #  significant negative (0.0036) -----------
# NEW
obs.sel_NEW_G21 <- na.omit(DE_table$sel.cov.HD.HD.21[DE_table$NWPC == 1])
obs.sel_UPNEW_G21 <- na.omit(DE_table$sel.cov.HD.HD.21[DE_table$logFC_PC21 > 0 & DE_table$NWPC == 1])
obs.sel_DOWNNEW_G21 <- na.omit(DE_table$sel.cov.HD.HD.21[DE_table$logFC_PC21 < 0 & DE_table$NWPC == 1])
obs.sel_NONEW_G21 <- na.omit(DE_table$sel.cov.HD.HD.21[!is.na(DE_table$logFC_PC21) & DE_table$NWPC == 0])

sel.resp(obs.sel_NEW_G21, obs.sel_NONEW_G21, sign = "up")     # not significant positive (0.1) 
sel.resp(obs.sel_NEW_G21, obs.sel_NONEW_G21, sign = "down")   # not significant negative (0.8)

sel.resp(obs.sel_UPNEW_G21, obs.sel_NONEW_G21, sign = "down")  # significant negative (0)   --------------
sel.resp(obs.sel_DOWNNEW_G21, obs.sel_NONEW_G21, sign = "up") # significant positive (0.0008) -----------
#####

# # # TEST (Permutation + Fisher)- - - - - - - - - - - - - - - - - - - - - - - - 
#####

## SELECTION INTENSITY

# 5037 # are the plastic genes (G1) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$PC1 == 1, ], "sel.cov.HD.1")                                    #yes 0.00009
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$PC1 == 1, ], "sel.cov.HD.HD.21")                                #yes 0.00009

# 3398 # are the genes with lost plasticity under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$LSPC == 1, ], "sel.cov.HD.1")                                   #yes 0.00009
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$LSPC == 1, ], "sel.cov.HD.HD.21")                               #yes 0.00009

# 235 # are the most plastic genes (G1) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$MPC1 == 1, ], "sel.cov.HD.1")                                   #yes 0.00009
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$MPC1 == 1, ], "sel.cov.HD.HD.21")                               #yes 0.00009

# 766 # are the most.50 plastic genes (G1) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$MPC1.50 == 1, ], "sel.cov.HD.1")                                #yes 0.00009
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$MPC1.50 == 1, ], "sel.cov.HD.HD.21")                            #yes 0.00009

# 1639 # are the genes with preserved plasticity under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$PRPC == 1, ], "sel.cov.HD.1")                                   #yes 0.00009          
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$PRPC == 1, ], "sel.cov.HD.HD.21")                               #yes 0.00009

# 643 are the new plastic genes under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$NWPC == 1, ], "sel.cov.HD.1")                                   #no 1
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$NWPC == 1, ], "sel.cov.HD.HD.21")                               #no 1

# 731 # are the evolved genes (G21) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$EC21 == 1, ], "sel.cov.HD.1")                                   #no 0.08
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$EC21 == 1, ], "sel.cov.HD.HD.21")                               #yes 0.00459 

# 336 # are the plastic genes evolved (G21) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$PCEC == 1,], "sel.cov.HD.1")                                    #yes 0.0006
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$PCEC == 1, ], "sel.cov.HD.HD.21")                               #yes 0.00009

# 49 # are the new plastic genes evolved (G21) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$NWPCEC == 1,], "sel.cov.HD.1")                                  #no 0.9
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$NWPCEC == 1, ], "sel.cov.HD.HD.21")                             #no 0.88

# 2282 # are the plastic genes (G21) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$PC21 == 1,], "sel.cov.HD.1")                                    #yes 0.0008 
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$PC21 == 1, ], "sel.cov.HD.HD.21")                               #yes 0.00009  

# 893 # are the genetically different genes (GWAS) under stronger selection in G1?
permtest(DE_table, DE_table[DE_table$SNPsEN22 == 1,], "sel.cov.HD.1")                                #no 1
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$SNPsEN22 == 1, ], "sel.cov.HD.HD.21")                           #no 1   

# 259 # are the hub genes under stronger selection in G1?     
permtest(DE_table, DE_table[DE_table$hub_gene == 1,], "sel.cov.HD.1")                                #yes 0.001
# under stronger selection in G21?
permtest(DE_table, DE_table[DE_table$hub_gene == 1, ], "sel.cov.HD.HD.21")                           #yes 0.00009

## ENRICHMENT IN THE HUB GENES
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "PC1")                 #no 0.11
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "MPC1")                #yes 0.03
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "MPC1.50")             #yes 0.02
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "PC21")                #no 0.19
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "NWPC")                #no 0.8
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "LSPC")                #no 0.5
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "PRPC")                #no 0.053
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "PCEC")                #yes 0.01
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "EC21")                #yes 0.00000008
ftest(DE_table[DE_table$hub_gene == 0, ], DE_table[DE_table$hub_gene == 1, ], "SNPsEN22")            #yes 0.00001

## ENRICHMENT IN THE EVOLVED (evolved/new plastic/plastic evolved) GENES
ftest(DE_table[DE_table$EC21 == 0,], DE_table[DE_table$EC21 == 1,], "SNPsEN22")                      #yes 0.0000001
ftest(DE_table[DE_table$NWPC == 0,], DE_table[DE_table$NWPC == 1,], "SNPsEN22")                      #yes 0.01 
ftest(DE_table[DE_table$PCEC == 0,], DE_table[DE_table$PCEC == 1, ], "SNPsEN22")                     #yes 0.02   

ftest(DE_table[DE_table$EC21 == 0,], DE_table[DE_table$EC21 == 1,], "PC1")                           #yes 0.009
ftest(DE_table[DE_table$EC21 == 0,], DE_table[DE_table$EC21 == 1,], "MPC1")                          #no 0.5
ftest(DE_table[DE_table$EC21 == 0,], DE_table[DE_table$EC21 == 1,], "MPC1.50")                       #no 0.4

# # # TEST (t.test/wilcox + cor.test)- - - - - - - - - - - - - - - - - - - - - - 

# SELECTION G1 VS SELECTION G21
t.test(DE_table$sel.cov.HD.1, DE_table$sel.cov.HD.HD.21, paired = T) #sign different distribution 1.072e-13
wilcox.test(DE_table$sel.cov.HD.1, DE_table$sel.cov.HD.HD.21, paired = T) #sign different distribution 3.286e-13

# # # PLOT  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#####

## SEL G1 vs SEL G21    < < < < < < < < < < < < < < < 
# box
boxplot(DE_table$sel.cor.HD.1, DE_table$sel.cor.HD.HD.21, 
        names = c("cor_selection G1", "cor_selection G21"), 
        col = c("blue", "red"),
        main = "cor.selG1 vs cor.selG21")
abline(h = 0, col = "black", lwd = 2, lty = 2)  

# distribution
df <- data.frame(
  valore = c(DE_table$sel.cor.HD.1, DE_table$sel.cor.HD.HD.21),
  gruppo = rep(c("cor.selG1", "cor.selG21"), each = length(DE_table$sel.cor.HD.1))
)

ggplot(df, aes(x = valore, fill = gruppo)) +
  geom_density(alpha = 0.5) +  
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Selection(cor) G1 vs G21 - not tranformed", x = "Value", y = "Density") +
  theme_minimal()



## Intersection venn.plot ( GWAS + PG1 + PG21 + EVOLVED )    < < < < < < < < < < < < < < < 
sets <- list(sg_EN22 = rownames(sg_genes_EN22), PC1 = sig.fit2.G1, EC21 = sig.fit2.EC21, PC21 = sig.fit2.PC21)
venn.plot <- venn.diagram(
  x = sets,
  filename = NULL, 
  fill = c("blue", "darkred", "yellow", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "GROUPS"
)

grid.newpage()
grid.draw(venn.plot)



## Permutations distribution | Selection intensity  < < < < < < < < < < < < < < < 

# EC21 sel.cov.HD.1
perm_stats1 <- perm_stats_DE_table_DE_table_EC21____1____sel.cov.HD.1
observed_stat1 <- observed_DE_table_DE_table_EC21____1____sel.cov.HD.1

S1.EC21 <- ggplot(data.frame(stat = perm_stats1), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat1, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution EC21 sel.cov.HD.1",
       x = "mean value", y = "Frequency") +
  theme_minimal()

# EC21 sel.cov.HD.HD.21
perm_stats2 <- perm_stats_DE_table_DE_table_EC21____1____sel.cov.HD.HD.21
observed_stat2 <- observed_DE_table_DE_table_EC21____1____sel.cov.HD.HD.21

S21.EC21 <- ggplot(data.frame(stat = perm_stats2), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat2, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution EC21 sel.cov.HD.HD.21",
       x = "mean value", y = "Frequency") +
  theme_minimal()

grid.arrange(S1.EC21, S21.EC21, ncol = 2)




# PC1 sel.cov.HD.1
perm_stats7 <- perm_stats_DE_table_DE_table_PC1____1____sel.cov.HD.1
observed_stat7 <- observed_DE_table_DE_table_PC1____1____sel.cov.HD.1

S1.PC1 <- ggplot(data.frame(stat = perm_stats7), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat7, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution PC1 sel.cov.HD.1",
       x = "mean value", y = "Frequency") +
  theme_minimal()

# PC1 sel.cov.HD.HD.21
perm_stats8 <- perm_stats_DE_table_DE_table_PC1____1____sel.cov.HD.HD.21
observed_stat8 <- observed_DE_table_DE_table_PC1____1____sel.cov.HD.HD.21

S21.PC1 <- ggplot(data.frame(stat = perm_stats8), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat8, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution PC1 sel.cov.HD.HD.21",
       x = "mean value", y = "Frequency") +
  theme_minimal()

grid.arrange(S1.PC1, S21.PC1, ncol = 2)




# PC21 sel.cov.HD.1
perm_stats9 <- perm_stats_DE_table_DE_table_PC21____1____sel.cov.HD.1
observed_stat9 <- observed_DE_table_DE_table_PC21____1____sel.cov.HD.1

S1.PC21 <- ggplot(data.frame(stat = perm_stats9), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat9, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution PC21 sel.cov.HD.1",
       x = "mean value", y = "Frequency") +
  theme_minimal()

# PC21 sel.cov.HD.HD.21
perm_stats10 <- perm_stats_DE_table_DE_table_PC21____1____sel.cov.HD.HD.21
observed_stat10 <- observed_DE_table_DE_table_PC21____1____sel.cov.HD.HD.21

S21.PC21 <- ggplot(data.frame(stat = perm_stats10), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat10, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution PC21 sel.cov.HD.HD.21",
       x = "mean value", y = "Frequency") +
  theme_minimal()

grid.arrange(S1.PC21, S21.PC21, ncol = 2)




# SNPsEN22 sel.cov.HD.1
perm_stats9 <- perm_stats_DE_table_DE_table_SNPsEN22____1____sel.cov.HD.1
observed_stat9 <- observed_DE_table_DE_table_SNPsEN22____1____sel.cov.HD.1

S1.ENV22 <- ggplot(data.frame(stat = perm_stats9), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat9, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution SNPsEN22 sel.cov.HD.1",
       x = "mean value", y = "Frequency") +
  theme_minimal()
# SNPsEN22 sel.cov.HD.HD.21
perm_stats10 <-  perm_stats_DE_table_DE_table_SNPsEN22____1____sel.cov.HD.HD.21
observed_stat10 <- observed_DE_table_DE_table_SNPsEN22____1____sel.cov.HD.HD.21

S21.ENV22 <- ggplot(data.frame(stat = perm_stats10), aes(x = stat)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = observed_stat10, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Permutations distribution SNPsEN22 sel.cov.HD.HD.21",
       x = "mean value", y = "Frequency") +
  theme_minimal()

grid.arrange(S1.ENV22, S21.ENV22, ncol = 2)  


## Selection distribution | Adaptive/Maladptive response   < < < < < < < < < < < < < < < 


plot_half_violin_vertical <- function(DE_table, sel_col, logFC_col = "G1.logFC_PC1", CAT_col = "PC1") {
  
  plot_data <- DE_table %>%
    filter(!is.na(.data[[sel_col]]), !is.na(.data[[logFC_col]]), !is.na(.data[[CAT_col]])) %>%
    mutate(Regulation = case_when(
      .data[[CAT_col]] == 0 ~ "Non-Responding",
      .data[[logFC_col]] > 0 & .data[[CAT_col]] == 1 ~ "Up-Regulated",
      .data[[logFC_col]] < 0 & .data[[CAT_col]] == 1 ~ "Down-Regulated"
    )) %>%
    filter(Regulation %in% c("Up-Regulated", "Down-Regulated", "Non-Responding"))
  
  df_up <- plot_data %>% filter(Regulation == "Up-Regulated")
  df_down <- plot_data %>% filter(Regulation == "Down-Regulated")
  df_nrs <- plot_data %>% filter(Regulation == "Non-Responding")
  
  stats_up <- df_up %>% summarise(Q1 = quantile(.data[[sel_col]], 0.25),
                                  Median = median(.data[[sel_col]]),
                                  Q3 = quantile(.data[[sel_col]], 0.75))
  stats_down <- df_down %>% summarise(Q1 = quantile(.data[[sel_col]], 0.25),
                                      Median = median(.data[[sel_col]]),
                                      Q3 = quantile(.data[[sel_col]], 0.75))
  stats_nrs <- df_nrs %>% summarise(Q1 = quantile(.data[[sel_col]], 0.25),
                                    Median = median(.data[[sel_col]]),
                                    Q3 = quantile(.data[[sel_col]], 0.75))
  
  box_width <- 0.06  
  
  ggplot() +

    gghalves::geom_half_violin(
      data = df_down, aes(x = "Selection", y = .data[[sel_col]], fill = "Down-Regulated"),
      side = "l", trim = FALSE, adjust = 1.2, alpha = 0.8, color = "black"
    ) +
    gghalves::geom_half_violin(
      data = df_up, aes(x = "Selection", y = .data[[sel_col]], fill = "Up-Regulated"),
      side = "r", trim = FALSE, adjust = 1.2, alpha = 0.8, color = "black"
    ) +
  
    geom_rect(aes(xmin = 1 - box_width, xmax = 1 - 0.01,
                  ymin = stats_down$Q1, ymax = stats_down$Q3),
              fill = "#0072B2", alpha = 0.8) +
    geom_rect(aes(xmin = 1 + 0.01, xmax = 1 + box_width,
                  ymin = stats_up$Q1, ymax = stats_up$Q3),
              fill = "#D55E00", alpha = 0.8) +
    geom_rect(aes(xmin = 1 - box_width/2, xmax = 1 + box_width/2,
                  ymin = stats_nrs$Q1, ymax = stats_nrs$Q3),
              fill = "grey40", alpha = 0.8) +
   
    annotate("point", x = 1 - 0.03, y = stats_down$Median, shape = 21, fill = "white", size = 3, stroke = 0.5, color = "black") +
    annotate("point", x = 1 + 0.03, y = stats_up$Median, shape = 21, fill = "white", size = 3, stroke = 0.5, color = "black") +
    annotate("point", x = 1, y = stats_nrs$Median, shape = 21, fill = "white", size = 3, stroke = 0.5, color = "black") +

    scale_fill_manual(values = c("Up-Regulated" = "#D55E00", "Down-Regulated" = "#0072B2")) +
    labs(x = NULL, y = sel_col, fill = "Regulation",
         title = paste("Selection Coefficient Distribution [", sel_col, "] of UP-DOWN genes in ", CAT_col, " respect to ", logFC_col)) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "top")
}




plot_half_violin_vertical(DE_table, sel_col = "sel.cov.HD.1", logFC_col = "G1.logFC_PC1", CAT_col = "PC1")
plot_half_violin_vertical(DE_table, sel_col = "sel.cov.HD.1", logFC_col = "G1.logFC_PC1", CAT_col = "MPC1")
plot_half_violin_vertical(DE_table, sel_col = "sel.cov.HD.HD.21", logFC_col = "logFC_PC21", CAT_col = "PC21")
plot_half_violin_vertical(DE_table, sel_col = "sel.cov.HD.HD.21", logFC_col = "logFC_EC21", CAT_col = "EC21")
plot_half_violin_vertical(DE_table, sel_col = "sel.cov.HD.HD.21", logFC_col = "logFC_PC21", CAT_col = "NWPC")





