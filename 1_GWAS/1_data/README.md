# CONTENT INFORMATION

## folder R
`function.R`: emmax() function + other functioin (source file: Petri)

## fitness files
`offspr_G1.txt`: offspring number of generation 1 in CT and HD; 486 samples  
`G21-fitness+count.csv.gz`: raw fitness+count data G21
[fread("./data/Counts.cpm-fit_CT244.txt.gz") # fitness data G1-CT] > used in `allAnalysis.R`  
[fread("./data/Counts.cpm-fit_HD242.txt.gz") # fitness data G1-HD] > used in `allAnalysis.R`

## vcf files
`Tribolium_castaneum_CT_G1_Tcas3.30_imputed.vcf.gz`: genotype data of generation 1 CT  
`Tribolium_castaneum_CT_G1_Tcas3.30_imputed.vcf.gz`: genotype data of generation 1 HD  
`Tribolium_castaneum_ALL_Tcas3.30_imputed.vcf.gz`: genotype data of generation 21 CT HD H D  (not able to upload - use the file below)

see Results-Output/OUTPUT for [`G21_CTHD.RData`](../2_Results-Output/OUTPUT/README.md#g21cthd): raw vcf file of all CT and HD samples G21 (+ID : last column))  

## Eva datasets
`SNP_dataset.rds`: SNP dataset  
`gene_dataset.rds`: gene dataset



