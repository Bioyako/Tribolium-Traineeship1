# [`DE_analysis`](./2_DE_analysis) + Genetic selection coefficient + Statistical tests

Different expression analysis between categories within generations 1 and 21 using *limma* (https://doi.org/10.1093/nar/gkv007).  
Genetic selection coefficient for generations 1 and 21 is calculated as covariance between expression level and relative fitness.  

## SUMMARY RESULTS TABLE
|   *category* |   *genes* |   *n.samples* | *filename* | *GO* |
|--------------|-----------|---------------|-----------------------|----------------------------------------------|
|    PC1       |       5037 |          485 | `PC1_genes.txt`       |    https://biit.cs.ut.ee/gplink/l/a1Vms8NN1Qt     |
|    PC21      |      2282 |          56   | `PC21_genes.txt`       |  https://biit.cs.ut.ee/gplink/l/aE4GjfbM7SM       |
|    EC21      |      731 |            54  | `EC21_genes.txt`       |  https://biit.cs.ut.ee/gplink/l/anbQlinpYS8    |
|    MPC1.50   |      766 |            485 | `MPC1.50_genes.txt`     |  https://biit.cs.ut.ee/gplink/l/a1wc-JQTTRB  |
|    MPC1      |       325 |           485 | `MPC1_genes.txt`        |   https://biit.cs.ut.ee/gplink/l/axVfsBjvQSP |
|    NWPC      |      643 |              - | `NWPC_genes.txt`        | https://biit.cs.ut.ee/gplink/l/avGYmoEgFQn     |
|    LSPC      |       3398 |            - | `LSPC_genes.txt`        | https://biit.cs.ut.ee/gplink/l/adB14YkdzT1      |
|    PRPC      |       1639 |            - | `PRPC_genes.txt`        | https://biit.cs.ut.ee/gplink/l/a_GkOdhnGSu  |
|    PCEC      |      336|              -  | `PCEC_genes.txt`        |  https://biit.cs.ut.ee/gplink/l/arz3Bqn4US4 |
|   GWAS.7     |       894 |           184 | `21ENV_ALL_902.txt`     | https://biit.cs.ut.ee/gplink/l/arLetGrJnRp   |
|  HUB.GENES   |       259 |             - | `HUB_genes.txt`          |             -                                |

## Methods information
 - P-value correction method used on DE is FDR (Benjamini-Hochberg)
 - Plastic and evolved genes (PC1, PC21, EC21) are definied as such if the difference in the expression level is higher than 10%
 - raw counts data are normalized with voom(limma) for DE analysis and with TMM(edgeR) for the selection coefficient

## DE Categories
1. Generation 1
   - [**PC1**] Plastic genes: CT vs HD
2. Generation 21
   - [**PC21**] Plastic genes: HD.HD vs HD.CT
   - [**EC21**] Evolved genes: CT.HD vs HD.HD

## Sub-categories from DE Categories
   - [**MPC1.50**] - Most.50 plastic genes: Plastic genes G1 with abs(logFC > 0.585) = expression level 50% higher or lower
   - [**MPC1**] - Most plastic genes: Plastic genes G1 with abs(logFC > 1) = expression level 100% higher or lower
   - [**NWPC**] - New plastic genes: Plastic genes G21 - Plastic genes G1
   - [**LSPC**] - Lost plasticity genes: Plastic genes G1 - Plastic genes G21
   - [**PRPC**] - Preserved plasticity genes: Plastic genes G1 ∩ Plastic genes G21
   - [**PCEC**] - Evolved plastic genes: Evolved genes G21 ∩ Plastic genes G1

## Category from GWAS.7 (21_ENV)
   - [**GN21ENV**] - Genes of the significant SNPs from [GWAS.7](../1_GWAS/README.md#gwas-information)

## Folder Content
- [`data`](./1a_data)
- [`Results-Output`](./2a_Results-Output)
- [`Plots`](./3a_Plots)
- [`gProfiler`](./4a_gProfiler)
