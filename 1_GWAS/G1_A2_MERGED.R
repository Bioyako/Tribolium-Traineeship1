### GWAS MERGED DATASETS ( ENVIRONMENT || FITNESS ) WITH UNIQUE SNPS
    #   1- Creating Emmax Linear Regression model
    #   2- FUNCTION adjpvgn - adjusted significant p-values
    #   3- Creating Plots (EmmaxLM_p-value, adjusted_p-value, EmmaxLM_vs_LM, Beta_LM, r2_LM, qqplot_EmmaxLM)
    #   4- Beta and Emmax_p-values correlation ENV-FTFT
    #   5- FUNCTION  mantx - Manhattan Plot 
    #   6- FUNCTION annot - SNPs of Interest
    #   7- FUNCTION genet -Gene of Interest

# setwd("\\\\ad.helsinki.fi/home/g/gfriso/Desktop/RStudioNEW/tribolium")
# setwd("C:/Users/gfriso/RStudioNEW/tribolium")

#install.packages("cowplot")
#install.packages("ggrepel")  
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
load("./OUTPUTs/plotting2.R")
load("./OUTPUTs/G1_A1_FUNCTIONS.R")
SNP_dataset <- readRDS("./data/SNP_dataset.rds")
GENE_dataset <- readRDS("./data/gene_dataset.rds")
load("./OUTPUTs/G1_A1_MERGED_DATASETS.RData")


#-----------------------------------------------------------------------------##

## 1 ## --- Emmax Linear Regression model
emx_EN <- emmax(X = GTGT, Y = ENV$ENVIRONMENT, K = GRMGT)
emx_FT <- emmax(X = GTGT, Y = FTFT$OFFSPR, K = GRMGT)   

## 2 ## FUNCTION Adjusted p-value & Significant p-value + gene dataframe #######
adjpvgn <- function(emx, ann_dataset) {
  adj <- p.adjust(emx$pval, "fdr")
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
################################################################################
adjpvgn(emx_EN, SNP_dataset)
adjpvgn(emx_FT, SNP_dataset)

## 3 ## Creating Plots (EmmaxLM_p-value, Emmax_adjusted_p-value, EmmaxLM_vs_LM, Beta_LM, r2_LM, qqplot_EmmaxLM)
plotting2(emx_EN, GTGT, ENV, lm = "yes", beta = "yes")
plotting2(emx_FT, GTGT, FTFT, lm = "yes", beta = "yes")

## 4 ## Correlation p-value and beta(LM) ENV-FTFT
par(mfrow = c(2, 2))
plot(-log10(emx_EN$pval), -log10(emx_FT$pval), main = "P-val_FTvsENV")
plot(-log10(adj_emx_EN), -log10(adj_emx_FT), main = "adjP-val_FTvsENV")
plot(abs(Beta_ENV),abs(Beta_FTFT), main = "BETA_FTvsENV")
cor.test(Beta_ENV, Beta_FTFT) 

#SNPS p-value above 0.05 in common EN-FT = n.23
wowSNP_05 <- cbind(adj_spv_emx_EN[rownames(adj_spv_emx_EN) %in% rownames(adj_spv_emx_FT), ], adj_spv_emx_FT[rownames(adj_spv_emx_FT) %in% rownames(adj_spv_emx_EN), ])
nw <- nrow(wowSNP_05)
cat("\nThe number of significant SNPs common to the two datasets are:\n", nw)

## 5 ## FUNCTION Manhattan Plot ################################################
mantx <- function(pv_vct_r, gt_r, wow_snp = NULL) {  
  pv_vct <- pv_vct_r[names(pv_vct_r) %in% rownames(SNP_dataset)]
  gt <- colnames(gt_r)[colnames(gt_r) %in% rownames(SNP_dataset)]
  pv_vector <- pv_vct
  ch_vector <- gsub(".*G", "", gt)
  ch_vector <- as.numeric(as.factor(gsub("-.*", "", ch_vector)))
  bp_vector <- as.numeric(gsub(".*-", "", gt))
  gn_vector <- SNP_dataset[,1][rownames(SNP_dataset) %in% gt]
  dataframe <- data.frame(CHR = ch_vector, BP = bp_vector, P = pv_vector, SNP = gt, GENE = gn_vector) 

  if (!is.null(wow_snp)) {
    significant_snps <- rownames(wow_snp) }
  else { wow_snp <- pv_vct_r[names(pv_vct_r) %in% rownames(SNP_dataset)]
  significant_snps <- names(wow_snp)[wow_snp < 0.05]
  }
  
  assign(paste0("pv_vector_", deparse(substitute(pv_vct_r))), pv_vector, envir = .GlobalEnv)
  assign(paste0("dtfr_", deparse(substitute(pv_vct_r))), dataframe, envir = .GlobalEnv)
  
  # manhattan(dataframe, chr ="CHR", bp ="BP", snp ="SNP", p ="P", highlight = significant_snps, annotatePval = 0.05)
  
  # dataframe <- dataframe %>%
  #   filter(-log10(P) > 0.6)   ##otherwise too heavy
  
  cps <- dataframe %>%
    group_by(CHR) %>%    ##COMPUTE CHROMOSOME SIZE
    summarise(chr_len = max(BP)) %>%  # 1 row for each combination - chr lenght 
    mutate(tot = cumsum(chr_len)-chr_len) %>% # create new columns - start position ##CALCULATE COMULATIVE POSITION OF EACH CHROMOSOME
    #select(-chr_len) %>% # keep or drop columns - just the start
    left_join(dataframe, ., by = c("CHR" = "CHR")) %>%    # add columns (keep observation in x) -  ## ADD INFO (start) TO THE INITIAL DATASET
    arrange(CHR, BP) %>% # order rows using column values
    mutate(BPcum = BP+tot) %>% ## ADD COMULATIVE POSITION FOR EACH SNP
    mutate( is_highlight=ifelse(SNP %in% significant_snps, "yes", "no"))
  axisdf = cps %>%
    group_by(CHR) %>%
    #filter(-log10(P) > 1.3) %>%
    summarise(center = ((max(BPcum)+min(BPcum))/2)) 
  
  cps$text <- paste("SNP: ", cps$SNP, "\nPosition: ", cps$BP, "\nChromosome: ", cps$CHR, "\nLOD score:", -log10(cps$P) %>% round(2), "\nGENE: ", cps$GENE)   # Prepare text description for each SNP
  cps$text1 <- paste("SNP: ", cps$SNP, "\nGENE: ", cps$GENE)  #for normal plot
  
  cps_annot <- cps %>%   # for normal plot
    arrange(desc(-log10(P))) %>%  
    slice(1:10) # choose the number of SNPs you want annotated
  
  PLOT1 <- ggplot(cps, aes(x = BPcum, y = -log10(P))) + # for normal plot
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) + 
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) + 
    scale_y_continuous(expand = c(0, 0) ) +     
    ylim(0,12) +
    scale_y_continuous(breaks = seq(0, 11.7, by = 1.3)) +
    geom_point(data=subset(cps, is_highlight=="yes"), color="orange", alpha=0.8, size=1.8) +  
    ggtitle(paste0("Manhattan Plot ", deparse(substitute(pv_vct_r)))) +
    theme_bw() +  
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold") 
      )  +  #:for normal plot
    geom_text_repel(data = cps_annot, aes(label = text1), color = "darkblue", size = 1.5,
                    nudge_y = 0.3, nudge_x = 0.3, max.overlaps = Inf, box.padding = 0.5, 
                    point.padding = 0.3, segment.size = 0.1, force = 1, force_pull = 0.5 )  +
    labs(x = "Chromosome")

  print(PLOT1) #normal plot
  
  PLOT <- ggplot(cps, aes(x = BPcum, y = -log10(P), text = text)) + #interactive plot
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) + 
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) + 
    scale_y_continuous(expand = c(0, 0) ) +     
    ylim(0,12) +
    scale_y_continuous(breaks = seq(0, 11.7, by = 1.3)) +
    geom_point(data=subset(cps, is_highlight=="yes"), color="orange", alpha=0.8, size=1.8) +  
    ggtitle(paste0("Manhattan Plot ", deparse(substitute(pv_vct_r)))) +
    theme_bw() +  
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold") 
      )  +
    labs(x = "Chromosome")
  p <- ggplotly(PLOT, tooltip="text")
  #print(p)    #interactive plot 
  file_name <- paste0("HtmlWidget/interactiveManhattanPlot_", deparse(substitute(pv_vct_r)), ".html")
  saveWidget(p, file = file_name) #saving interactive plot html
}

################################################################################

cat("## Building Manhattan plots\n\n")
mantx(adj_emx_EN, GTGT)
mantx(adj_emx_FT, GTGT)   #, wowSNP_05)

# How much differ the observed p-value distribution from an expected one by chance?
# qq(emx_EN$pval)
# qq(emx_FT$pval)

#ggplot correlation p-value FT-EN (creating single dataframe)
tot_dtf <- copy(dtfr_adj_emx_EN)
setnames(tot_dtf, "P", "P_EN")
P_FT <- dtfr_adj_emx_FT[3]
setnames(P_FT, "P", "P_FT")
tot_dtf <- cbind(tot_dtf, P_FT)

ggplot(tot_dtf, aes(-log10(P_EN), -log10(P_FT))) + geom_point() + labs(x = "P_EN", y = "P_FT") + ggtitle("Correlation p-value FT-EN")


## 6 ## FUNCTION SNPs annotation Insights ######################################
annot <- function(dtfr1, ann_dataset, pv_05) {
  #(creating single dataframe)                                                                                                        
  dtfr <- cbind(dtfr1, ann_dataset[ann_dataset$ID %in% dtfr1$SNP, c("Parallelism", "abs_logFC", "hub_gene", "is_eQTL")])
  ## dtfr$text <- paste("SNP: ", dtfr$SNP, "\nGENE: ", dtfr$GENE)  
  
  # Evaluating SNP Parralelism
  ann_para <- dtfr %>%
    filter(Parallelism > 2) %>%
    filter(-log10(P) > 1.3) %>%
    arrange(desc(Parallelism), P) %>% 
    slice_head(n = 10) %>% 
    mutate(text = paste("SNP: ", SNP, "\nGENE: ", GENE))
  plot1 <- ggplot(dtfr %>% filter(-log10(P) > 1.3), aes(x = Parallelism, y = -log10(P))) +
    geom_point(size = 3, alpha = 0.5, aes(color=as.factor(CHR))) +
    geom_text_repel(data = ann_para, aes(label = text), color = "darkblue", size = 1.5,
                    nudge_y = 0.3, nudge_x = 0.3, max.overlaps = Inf, box.padding = 0.5, 
                    point.padding = 0.3, segment.size = 0.1, force = 1, force_pull = 0.5 ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Parallelism", y = paste0("-log10(p_value)_", deparse(substitute(dtfr1))), title = "Parallelism of SNPs in 6 lines", color = "CHR") +
    ylim(0,12) +
    scale_y_continuous(breaks = seq(0, 11.7, by = 1.3)) 
  #print(plot1)
  paral <- dtfr[dtfr$SNP %in% pv_05$ID,] %>%
    filter(Parallelism >= 4) 
  # Evaluating SNP abs_logFC
  ann_logFC <- dtfr %>%
    filter(abs_logFC > 0.4) %>%
    filter(-log10(P) > 1.3) %>%
    arrange(desc(abs_logFC)) %>% 
    slice_max(abs_logFC, n = 10) %>% 
    mutate(text = paste("SNP: ", SNP, "\nGENE: ", GENE))
  plot2 <- ggplot(dtfr %>% filter(-log10(P) > 1.3), aes(x = abs_logFC, y = -log10(P))) +
    geom_point(size = 3, alpha = 0.5, aes(color=as.factor(CHR))) +
    geom_text_repel(data = ann_logFC, aes(label = text), color = "darkblue", size = 1.5,
                    nudge_y = 0.3,nudge_x = 0.3, max.overlaps = Inf, box.padding = 0.5, 
                    point.padding = 0.3, segment.size = 0.1, force = 1, force_pull = 0.5 ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "abs_logFC", y = paste0("-log10(p_value)_", deparse(substitute(dtfr1))), title = "Fold Change values of SNPs", color = "CHR") +
    ylim(0,12) +
    scale_y_continuous(breaks = seq(0, 11.7, by = 1.3)) 
  #print(plot2)  
  logfc <- dtfr %>%
    filter(abs_logFC > 1.4 & rownames(dtfr) %in% rownames(pv_05)) 
  # Evaluating SNP outlier eQTL
  qtl <- dtfr[rownames(dtfr) %in% rownames(pv_05), c("P", "is_eQTL")]
  cat("\n## Displaying significant SNPs from ", deparse(substitute(dtfr1)), " that are eQTL: \n")
  print(dtfr[rownames(dtfr) %in% rownames(qtl) & dtfr$is_eQTL == 1, c("P", "GENE")])
  eqtl <- dtfr %>%
    filter(is_eQTL == 1 & rownames(dtfr) %in% rownames(pv_05)) 
  # Evaluating SNP hubgenes
  hub <- dtfr[rownames(dtfr) %in% rownames(pv_05), c("P", "hub_gene")]
  cat("\n## Displaying significant SNPs from ", deparse(substitute(dtfr1)), " that are within hubgene: \n")
  print(dtfr[rownames(dtfr) %in% rownames(hub) & dtfr$hub_gene == 1, c("P", "GENE")])
  hubg <- dtfr %>%
    filter(hub_gene == 1 & rownames(dtfr) %in% rownames(pv_05)) 
  # SNPS of interest
  soi <- rbind(paral, logfc, eqtl, hubg)
  n <- nrow(soi)
  cat("\nThe number of SNPs of interest are:\n", n, "\n")

  grid.arrange(plot1, plot2, ncol = 2)
  
  cps <<- dtfr %>%
    group_by(CHR) %>%   
    summarise(chr_len = max(BP)) %>%  
    mutate(tot = cumsum(chr_len)-chr_len) %>%
    left_join(dtfr, ., by = c("CHR" = "CHR")) %>% 
    arrange(CHR, BP) %>% 
    mutate(BPcum = BP+tot) 
  axisdf = cps %>%
    group_by(CHR) %>%
    #filter(-log10(P) > 1.3) %>%
    summarise(center = ((max(BPcum)+min(BPcum))/2)) 
  
  # indices <- match(soi$SNP, cps$SNP)
  # soi$BPcum <- cps$BPcum[indices] 
  soi2 <- left_join(soi, cps[, c("SNP", "BPcum")], by = "SNP") %>%
    mutate(
      abs_logFC = round(abs_logFC, 4),
      text = paste(
        "SNP: ", SNP,
        "\nGENE: ", GENE,
        "\nParallelism: ", Parallelism,
        "\nabs_logFC: ", abs_logFC,
        "\nhub_gene: ", hub_gene,
        "\neQTL: ", is_eQTL
      )
    )  
  
  assign(paste0("SNPsoi_", deparse(substitute(dtfr1))), soi, envir = .GlobalEnv)
  
  
  SoI_plot <- ggplot(cps, aes(x = BPcum, y = -log10(P))) + 
    #geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) + 
    scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) + 
    scale_y_continuous(expand = c(0, 0) ) +     
    ylim(0, 12) +
    #scale_y_continuous(breaks = seq(0, 11.7, by = 1.3)) +
    #geom_point(data = subset(cps, is_highlight == "yes"), color = "orange", size = 2) +  
    ggtitle(paste0("Manhattan SoI ", deparse(substitute(dtfr1)))) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_point(data = soi2, aes(x = BPcum, y = -log10(P), color = as.factor(CHR)), alpha = 0.8, size = 2) +
    geom_text_repel(data = soi2, aes(label = text), color = "darkblue", size = 1.5,
                    nudge_y = 0.3, nudge_x = 0.3, max.overlaps = Inf, box.padding = 0.5, 
                    point.padding = 0.3, segment.size = 0.1, force = 1, force_pull = 0.5 )  +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold") 
    )  +
    labs(x = "Chromosome")
  print(SoI_plot)
  
}
################################################################################
annot(dtfr_adj_emx_EN, SNP_dataset, adj_spv_emx_EN)
annot(dtfr_adj_emx_FT, SNP_dataset, adj_spv_emx_FT)

## do the snps variables losted from CT-G1 come back in the G21?
cmm_soi <- SNPsoi_dtfr_adj_emx_EN[SNPsoi_dtfr_adj_emx_EN$SNP %in% SNPsoi_dtfr_adj_emx_FT$SNP,]
cat("\nThe number of SNPs of Interest common are:\n")
print(cmm_soi)

## 7 ## FUNCTION GENE LEVEL INSIGHT ###########################################
genet <- function(dtfr, genes_dataset, pv_05) {
  sgoi <- genes_dataset[rownames(genes_dataset) %in% pv_05$gene,]
  n <- nrow(sgoi)
  sprint <- sgoi[, c("gene", "CHROM", "Parallelism", "abs_logFC", "highest_sg_pool", "highest_de_pool", "hub_gene", "net_selection", "eQTL_carrier", "Pleiotropy")]
  cat("\n## The number of significant Genes from ", substitute(dtfr), " are: \n", n, "\n\n")
  # Evaluating highest_sg_pool
  cat("\n## Displaying significant Genes from ", substitute(dtfr), " that are in highest_sg_pool: \n")
  print(sprint[sprint$highest_sg_pool == 1,])
  hgsg <- sgoi %>%
    filter(highest_sg_pool == 1) 
  # Evaluating highest_de_pool
  cat("\n## Displaying significant Genes from ", substitute(dtfr), " that are in highest_de_pool: \n")
  print(sprint[sprint$highest_de_pool == 1,])
  hgde <- sgoi %>%
    filter(highest_de_pool == 1)
  # Evaluating hub_gene
  cat("\n## Displaying significant Genes from ", substitute(dtfr), " that are hub_gene: \n")
  print(sprint[sprint$hub_gene == 1,])
  hub <- sgoi %>%
    filter(hub_gene == 1)
  # Evaluating net_selection
  cat("\n## Displaying significant Genes from ", substitute(dtfr), " that have higher net_selection (> 9): \n")
  print(sprint[abs(sprint$net_selection) > 9,])
  nsel <- sgoi %>%
    filter(abs(net_selection) > 9)
  # Evaluating eQTL carrier
  cat("\n## Displaying significant Genes from ", substitute(dtfr), " that are eQTL carrier: \n")
  print(sprint[sprint$eQTL_carrier == 1,])
  eqtl <- sgoi %>%
    filter(eQTL_carrier == 1)
  # Evaluating Pleiotropy
  cat("\n## Displaying significant Genes from ", substitute(dtfr), " that have pleiotropy (> 0): \n")
  print(sprint[sprint$Pleiotropy > 0,])
  plpy <- sgoi %>%
    filter(Pleiotropy > 0)
  # Gene of interest
  goi <- rbind(hgsg, hgde, hub, nsel, eqtl, plpy)
  nn <- nrow(goi)
  cat("\nThe number of GENE of Interest are:\n", nn, "\n\n")
  assign(paste0("GoI_", deparse(substitute(dtfr))), goi, envir = .GlobalEnv)
}
################################################################################
genet(dtfr_adj_emx_EN, GENE_dataset, adj_spv_emx_EN)  
genet(dtfr_adj_emx_FT, GENE_dataset, adj_spv_emx_FT)  

cmm_goi <- GoI_dtfr_adj_emx_EN[GoI_dtfr_adj_emx_EN$gene %in% GoI_dtfr_adj_emx_FT$gene, 
                               c("gene", "CHROM", "Parallelism", "abs_logFC", "highest_sg_pool", "highest_de_pool", "hub_gene", "net_selection", "eQTL_carrier", "Pleiotropy")]
cat("\nThe number of GENE of Interest common are:\n")
print(cmm_goi)


# save(emx_EN, ENV, dtfr_adj_emx_EN, adj_spv_emx_EN, SNPsoi_dtfr_adj_emx_EN, GoI_dtfr_adj_emx_EN, FTFT, emx_FT, dtfr_adj_emx_FT, adj_spv_emx_FT, SNPsoi_dtfr_adj_emx_FT, GoI_dtfr_adj_emx_FT, GTGT, GRMGT, file = "G1_A2_DATASETS_Results.RData")
# save(adjpvgn, mantx, annot, genet, file = "G1_A2_FUNCTIONS.R")




























