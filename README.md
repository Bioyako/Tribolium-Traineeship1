
#  - EVOLUTION OF PHENOTYPIC PLASTICITY OF GENE EXPRESSION -

This repository contains scripts and results from my Master's thesis project focused on transcriptome data analysis. The project is divided into two main directories: `GWAS` and `DE_analysis`.

## TABLE OF CONTENTS

- [Project Overview](#project-overview)
- [Directory Structure](#directory-structure)
  - [`GWAS`](./1_GWAS)
    - R Scripts
    - [`data`](./1_GWAS/1_data)
    - [`Results-Output`](./1_GWAS/2_Results-Output)
    - [`Plots`](./1_GWAS/3_Plots)
    - [`Manhattan.html`](./1_GWAS/4_Manhattan.html)
    - [`gProfiler`](./1_GWAS/5_gProfiler)
  - [`DE_analysis`](./2_DE_analysis)
    - R Scripts
    - [`data`](./2_DE_analysis/1a_data)
    - [`Results-Output`](./2_DE_analysis/2a_Results-Output)
    - [`Plots`](./2_DE_analysis/3a_Plots)
    - [`gProfiler`](./2_DE_analysis/4a_gProfiler)


## Project Overview

This project aims to analyze experimental whole-transcriptome RNA-seq data of <i>Tribolium castaneum</i> (red flour beetle) to gain insight into the evolution of phenotypic plasticity in gene expression during adaptation for 20 generations to a stressful environment (Hot-Dry).  
Firstly, multiple GWASs in generations 1 and 21 are established to explore the potential genetic variants associated with fitness, environment, and ultimately with plasticity and evolutionary changes. Gene-Ontology analysis is integrated to identify biological processes, cellular locations and molecular functions of the genes involved in the genetic differences.  
Secondly, plasticity and evolutionary changes are inferred by performing Different Expression Analysis of meaningful comparisons between Control and Hot-Dry conditions. Further Gene-Ontology analysis is performed.  
Then the genetic selection coefficient is calculated for both generations as covariance between expression level and relative fitness, to estimate the response profile of the plasticity (adaptive or maladaptive) and overall both short- and long- term responses to selection.  
Lastly, statistical tests, mostly permutation and Fisher tests, are performed to obtain significant information from the results of the above analysis. Additionally, statistical analysis are used to investigate the relationship between the evolution of plasticity and network position of the genes under selection in functional gene co-expression networks.  
Many plots are integrated with the analysis.

## Directory Structure

This project is organized as follows:

### [`GWAS`](./1_GWAS) - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - -  - - - - - -

Assessment of genetic variants associated with fitness (offspring numbers), and environment in generations 1 and 21.  
Association tests are conducted using *EMMAX* linear model (https://doi.org/10.1038/ng.548) to account for the genetic structure of the groups.  
For each generation 4 GWASs are established:
1. Association of Control environment with fitness
2. Association of Hot-Dry environment with fitness
3. Association of all samples (Control and Hot-Dry together) with fitness
4. Association of all samples with environment (CT = 0; HD = 1)

#### Folder Content 

- [`data`](./1_GWAS/1_data)
- [`Results-Output`](./1_GWAS/2_Results-Output)
- [`Plots`](./1_GWAS/3_Plots)
- [`Manhattan.html`](./1_GWAS/4_Manhattan.html)
- [`gProfiler`](./1_GWAS/5_gProfiler)

### [`DE_analysis`](./2_DE_analysis) + Genetic selection coefficient + Statistical tests  - - - - - - - - - - - - - - - - - - - -

Different expression analysis between categories within generations 1 and 21 using *limma* (https://doi.org/10.1093/nar/gkv007):  

#### DE Categories
1. Generation 1
   - [**PC1**] Plastic genes: CT vs HD
2. Generation 21
   - [**PC21**] Plastic genes: HD.HD vs HD.CT
   - [**EC21**] Evolved genes: CT.HD vs HD.HD

#### Folder Content
- [`data`](./2_DE_analysis/1a_data)
- [`Results-Output`](./2_DE_analysis/2a_Results-Output)
- [`Plots`](./2_DE_analysis/3a_Plots)
- [`gProfiler`](./2_DE_analysis/4a_gProfiler)



