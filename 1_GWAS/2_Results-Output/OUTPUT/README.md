# CONTENT INFORMATION

## output datasets/functions

### from `G1_A1.R` script
`G1_A1_DATASETS_Results.RData` - [GWAS.1](../../README.md#gwas-information) and [GWAS.2](../../README.md#gwas-information)
 - GTCT > genotype CT G1 dataset
 - FTCT > fitness CT G1 dataset
 - GRMCT > genetic relationship matrix CT G1 
 - emx_ct > emmax results CT G1 (list F, pval, Rsq)
 - GTHD > genotype HD G1 dataset
 - FTHD > fitness HD G1 dataset
 - GRMHD > genetic relationship matrix HD G1 
 - emx_hd > emmax results HD G1 (list F, pval, Rsq)

 `G1_A1_MERGED_DATASETS.RData` - [GWAS.5](../../README.md#gwas-information) and [GWAS.6](../../README.md#gwas-information)
 - GTGT > total genotype CT+HD G1 dataset
 - GRMGT > genetic relationship matrix CT+HD G1  
 - ENV > total environment CT+HD G1 dataset (CT = 0; HD = 1)
 - FTFT > total fitness CT+HD G1 dataset

`G1_A1_FUNCTIONS.R`
 - gc(): Genotype Characterization (1|0, 0|1, 1|1, 0|0, --> 1, 1, 2, 0 )
 - forch(): Format and Check
 - grmx(): Build Genetic Relationship MatriX
 - snpa(): Single SNP Analysis
 - plotting(): Plots building

### from `G1_A2_MERGED.R` script
`G1_A2_DATASETS_Results.RData` - [GWAS.5](../../README.md#gwas-information) and [GWAS.6](../../README.md#gwas-information)
 - emx_EN > emmax results ENV G1 (list F, pval, Rsq)
 - ENV > total environment CT+HD G1 dataset (CT = 0; HD = 1)
 - dtfr_adj_emx_EN > emmax adjusted pvalues ENV G1 
 - adj_spv_emx_EN > emmax significant adjusted pvalues ENV G1
 - SNPsoi_dtfr_adj_emx_EN > SNPs of interest from adj_spv_emx_EN
 - GoI_dtfr_adj_emx_EN > genes of interest from adj_spv_emx_EN
 - FTFT > total fitness CT+HD G1 dataset
 - emx_FT > emmax results FT G1 (list F, pval, Rsq)
 - dtfr_adj_emx_FT > emmax adjusted pvalues FT G1
 - adj_spv_emx_FT > emmax significant adjusted pvalues FT G1
 - SNPsoi_dtfr_adj_emx_FT > SNPs of interest from adj_spv_emx_FT 
 - GoI_dtfr_adj_emx_FT > genes of interest from adj_spv_emx_FT
 - GTGT > total genotype CT+HD G1 dataset
 - GRMGT > genetic relationship matrix CT+HD G1

`G1_A2_FUNCTIONS.R` 
 - adjpvgn(): Adjusted p-value & Significant p-value + gene dataframe
 - mantx(): Build Manhattan Plot
 - annot(): SNPs annotation Insights; generate SNPsoi_dtfr
 - genet(): genes annotation Insights; generate GoI_dtfr

### from `G21_A1.R` script
`G21_ENV_FT_Results.RData`; *samples CT.XX, HD.XX* - [GWAS.7](../../README.md#gwas-information) and [GWAS.8](../../README.md#gwas-information)
 - emx_G21_env > emmax results ENV G21 (list F, pval, Rsq)
 - ENVG21 > total environment CT.XX+HD.XX G21 dataset (CT = 0; HD = 1)
 - adj_spv_emx_G21_env > emmax significant adjusted pvalues ENV G21
 - dtfr_adj_emx_G21_env > emmax adjusted pvalues ENV G1 
 - SNPsoi_dtfr_adj_emx_G21_env > SNPs of interest from adj_spv_emx_G21_env
 - GoI_dtfr_adj_emx_G21_env > genes of interest from adj_spv_emx_G21_env
 - G21FT > total fitness CT.XX+HD.XX G21 dataset
 - emx_G21_ft > emmax results FT G21 (list F, pval, Rsq)
 - GTG21 > total genotype CT.XX+HD.XX G21 dataset
 - GRMGTG21 > genetic relationship matrix CT.XX+HD.XX G21

 `G21_ENV_FT_zeta_Results.RData` - [GWAS.7z](../../README.md#gwas-information) and [GWAS.8z](../../README.md#gwas-information)
 Same information above; *samples CT.CT, HD.HD*
 - emx_G21_envz 
 - ENVG21z
 - SNPsoi_dtfr_adj_emx_G21_envz
 - GoI_dtfr_adj_emx_G21_envz
 - G21FTz
 - emx_G21_ftz
 - GTG21z
 - GRMGTG21z

 `G21_CTa_HDa_Results.RData`; *samples CT.XX, HD.XX* - [GWAS.7](../../README.md#gwas-information) and [GWAS.8](../../README.md#gwas-information)
 - emx_G21_ctA > emmax results CT G21 (list F, pval, Rsq)
 - GTG21ctA > genotype CT G21 dataset
 - G21FTctA > fitness CT G1 dataset
 - emx_G21_hdA > emmax results HD G21 (list F, pval, Rsq)
 - GTG21hdA > genotype HD G21 dataset
 - GRMGTG21ctA > genetic relationship matrix CT G21
 - GRMGTG21hdA > genetic relationship matrix HD G21

### other
`plotting2.R`: plot function plotting2()  
`G21_CTHD.RData`: raw vcf file of all CT and HD samples G21 (+ID : last column)
