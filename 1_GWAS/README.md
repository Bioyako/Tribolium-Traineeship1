# Genome Wide Association Studies 

Assessment of genetic variants associates with Fitness or Environment.  
Here ´Selection Lines´, CT and HD, are defined as Environment.

Testing correlation between genotype and Fitness: [GWAS.1](#gwas1---no-significant), [GWAS.2](#gwas2---no-significant), [GWAS.3](#gwas3---no-significant), [GWAS.4](#gwas4---no-significant), [GWAS.6](#gwas6---significant), [GWAS.8](#gwas8---no-significant)  
Testing correlation between genotype and Environment: [GWAS.5](#gwas5---significant), [GWAS.7](#gwas7---significant)  

Use of Emmax(https://doi.org/10.1038/ng.548) to account for the genetic structure of the groups

For [`SCRIPTS INFORMATION`](#scripts-information) see below

## SUMMARY RESULTS TABLE

|  *datasets* | *SNPs* | *genes* | *n.samples* | *GWAS*  |  *filename*  | *GO*  |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
|  FTCT.G1 |  0   |   0   | 244   |  1  |  -  |  -  |
|  FTHD.G1 |  0   |   0   | 242   |  2  |  -  |  -  |
|  FT.G1 |  19   |   13   | 485   |  6  |  `1FT_13.txt`  |  no  |
| ENV.G1 |  256    |  100   | 485   |  5  | `1ENV_100.txt`  |  https://biit.cs.ut.ee/gplink/l/aWx28voXqRD  |
|  FTCT.G21 |  0   |   0   | 87   |  3  |  -  |  -  |
|  FTHD.G21 |  0   |   0   | 67   |  4  |  -  |  -  |
| FTCT.G21z |  boh  |  4   |  53   |  boh  |  `21FT_CT_CT_4.txt`  | no  |
| **ENV.G21** |  1699  |  902   |  154  |  7  |  `21ENV_ALL_902.txt`  | https://biit.cs.ut.ee/gplink/l/arLetGrJnRp  |
| ENV.G21z | 79    |  19    |  53  |  7z  |  `21ENV_CTCT_HDHD_19.txt`  | https://biit.cs.ut.ee/gplink/l/amhpEHh10S1  |
|  FT.G21 |  0   |   0   | 184   |  8  |  -  |  -  |
|  FT.G21z |  0   |   0   | 53   |  8z  |  -  |  -  |

## DIRECTORY CONTENT

- [`data`](./1_data)
- [`Results-Output`](./2_Results-Output)
    - [`GENE_LIST`]
    - [`OUTPUT`]
- [`Plots`](./3_Plots)
    - [`CT-HD.G1`](./3_Plots/CT-HD.G1)
    - [`CT-HD.G21`](./3_Plots/CT-HD.G21)
    - [`FT-ENV.G1`](./3_Plots/FT-ENV.G1)
      - [`ENV.G1`]
      - [`FT.G1`]
    - [`FT-ENV.G21`](./3_Plots/FT-ENV.G21)
      - [`ENV.G21`]
      - [`FT.G21`]
- [`Manhattan.html`](./4_Manhattan.html)
- [`gProfiler`](./5_gProfiler)

## <h2 id="gwas-information">GWAS INFORMATION</h2>

### `GWAS.1` - NO SIGNIFICANT  
SNPs associated with Fitness in CT lines in Generation_1  
Number of samples: 243

### `GWAS.2` - NO SIGNIFICANT  
SNPs associated with Fitness in HD lines in Generation_1  
Number of samples: 242

### `GWAS.3` - NO SIGNIFICANT  
SNPs associated with Fitness in CT (CT.XX) lines in Generation_21  
Number of samples: 87

### `GWAS.4` - NO SIGNIFICANT  
SNPs associated with Fitness in HD (HD.XX) lines in Generation_21  
Number of samples: 67

### `GWAS.5` - SIGNIFICANT  
SNPs associated with different environments (CT or HD) in Generation_1  
Number of samples: 485  
**Signicant SNPs: 256**

### `GWAS.6` - SIGNIFICANT  
SNPs associated with Fitness indipendently by the environment in Generation_1  
Number of samples: 485  
**Signicant SNPs: 19**

### `GWAS.7` - SIGNIFICANT  
SNPs associated with different environments (CT.XX or HD.XX) in Generation_21  
Number of samples: 184  
**Signicant SNPs: 1699**

### `GWAS.7z` - SIGNIFICANT  
SNPs associated with different environments (CT.CT or HD.HD) in Generation_21  
Number of samples: 53  
**Signicant SNPs: 79**

### `GWAS.8` - NO SIGNIFICANT  
SNPs associated with Fitness indipendently by the environment (CT.XX or HD.XX) in Generation_21  
Number of samples: 184

### `GWAS.8z` - NO SIGNIFICANT  
SNPs associated with Fitness indipendently by the environment (CT.CT or HD.HD) in Generation_21  
Number of samples: 53

## SCRIPTS INFORMATION

### Basic scripts

#### `allAnalysis.R`
- 1- Library source
- 2- Datasets import
- 3- Function loading
- 
- A- G1 ANALYSIS ---------------------------------------------------------------
- A1- Building fitness datasets
- A2 Building genotype matrixes 
- A3- Formatting fitness and genotype datasets
- A4- Building total environment, total genotype, and total fitness datasets
- A5- Building a genetic relationship matrix
- A6- Building emmax linear regression model
- A7- Building adjusted emmax p-values vector and significant p-values dataframe
- A8- Building significant genes table for G1-Environment and G1-Fitness
- A9- Performing a Fisher´s exact test on the gene category specified for G1-Environment and G1-Fitness
- A10- Performing a permutation test for the category specified on significant genes of G1-Environment and G1-Fitness
- 
- B- G22 ANALYSIS --------------------------------------------------------------
- B1- Filtering genotype data
- B2- Building total genotype matrix
- B3- Building fitness dataset
- B4- Formatting total genotype matrix and building CT and HD genotype matrixes  
- B5- Building environment dataset
- B6- Building genetic relationship matrixes
- B7- Building emmax linear regression model
- B8- Building adjusted emmax p-values vector and significant p-values dataframe
- B9- Building significant genes table for G22-Environment
- B10- Performing a Fisher´s exact test on the gene category specified of significant genes for G22-Environment
- B11- Performing a permutation test for the category specified on significant genes of G22-Environment

#### <span id="all">[`allFunction.R`]</span>

-   1- emmax() - building a linear regression model accounting for relatedness among individuals 
-   2- gc() - building a genotype matrix for emmax() from a .vcf file 
-   3- grmx() - building a genetic relationship matrix for emmax()
-   4- adjpvgn() -  building a vector with p-values from emmax() adjusted with 'fdr'  method, and a dataframe with the SNP ID and the respective GENE of the significant adjusted p-values (-log10 pval > 1.3) 
-   5- ftest() - performing a Fisher´s exact test on the gene category specified
-   6- permtest() - performing a permutation test on the gene category specified
-   7- freq() - extracting the allele frequencies


### Full(raw) scripts

#### `G1_A1.R`  

- G1 - ANALYSIS 1 - CT~FT + HD~FT
- FUNCTION ------------------------------------
- gc               Genotype Characterization  
- forch            Format and Check
- grmx             Genetic Relationship MatriX
- snpa             Single SNP Analysis
- emmax            EMMAX analysis
- plotting         Plots building

#### `G1_A2_MERGED.R`

- GWAS MERGED DATASETS ( ENVIRONMENT || FITNESS ) WITH UNIQUE SNPS
-   1- Creating Emmax Linear Regression model
-   2- FUNCTION adjpvgn - adjusted significant p-values
-   3- Creating Plots (EmmaxLM_p-value, adjusted_p-value, EmmaxLM_vs_LM, Beta_LM, r2_LM, qqplot_EmmaxLM)
-   4- Beta and Emmax_p-values correlation ENV-FTFT
-   5- FUNCTION  mantx - Manhattan Plot 
-   6- FUNCTION annot - SNPs of Interest
-   7- FUNCTION genet -Gene of Interest

#### `G21_A1.R`

- G21 -- ANALYSIS
- 
- BASEMENT 
-   1 - Building GENOTYPE MATRIX GTG21
-   2 - Building FITNESS DATASET
-   3 - Building corrected GENOTYPE MATRIXES
-   4 - Building ENVIRONMENT DATASET (0 = CT, 1 = HD)
-   5 - Building GRM - Genetic relationship matrix
- ANALYSIS  
-   6 - Building Emmax Linear Model
-   7 - RESULTS - adjpvgn - plotting - mantx - annot -genet





