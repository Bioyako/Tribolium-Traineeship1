## OUTPUT CONTENT 

`DE_table.csv` - Results table of DE with all the categories, logFC and selection coefficients  
`PC1_table.csv` - logFC table of genes after p-value correction from DE of generation 1, CT vs HD  
`PC21_table.csv` - logFC table of  genes after p-value correction from DE of generation 21, HD.HD vs HD.CT  
`EC21_table.csv` - logFC table of  genes after p-value correction from DE of generation 21, HD.HD vs CT.HD  




## GENE_LIST CONTENT
List of significant DE genes + GWAS.7 + HUB_genes

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

## DE TESTS SUMMARY
Results of different statistical tests:
- is the category X under stronger selection in generation 1 and 21?
- are the HUB_genes enriched with the category X?
- are the genes with observed changes at the expression level (EC21, NWPC, PCEC) enriched of genes with observed changes at the genetic level (GWAS.7- ENV22)?
- are the plastic and evolutionary responses of the category X in genration 1 and 21 adaptive or maladaptive?