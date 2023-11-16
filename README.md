# OverOptimism_in_GeneSetAnalysis
This repository allows you to 
- reproduce the results from our over-optimism study generated in R 
- inspect the documentation for the web-based applications

@Moritz: Kann man annehmen, dass der Leser weiß, dass das working directory dann genau den Aufbau haben soll wie das Repo auf Github?

 -------------------------------  
**Preparation of the gene expression data for the optimizations**

This is the first script to run when reproducing the results.

- **Random_Phenotype_Permutations.R**: *R* code to obtain
  1. both gene expression data sets and true labels of the respective samples 
  2. the ten random permutations of the true sample labels for each for both gene expression data sets
  
  The generated random phenotype permutations and the gene expression data set "bottomly" are additionally provided in the folder **GeneExpression_data** in this repository
 
 -------------------------------  
**Preparation of the gene expression data for the GSA methods**

 - **PreProcessing_Functions.R**: contains functions provided in the preprocessing of (almost) all of the investigated GSA methods and is therefore sourced in each of the respective R scripts
 - **RNASeq_Transformation.R**: contains two functions to transform the gene expression measurements to match the characteristics for microarray data, as needed for the methods *PADOG* and *GSEA* (web-based application). It is sourced in the corresponding scripts for the methods
 
 -------------------------------
  **Run the optimizations**
The following optimization scripts were generated for all methods implemented in *R*, i.e. *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, and *clusterProfiler*'s GSEA
  
- **n_DEGS_optim_... .R**: *R* code for the optimization of the number of differentially enriched gene sets for the respective computational GSA method (optimization goal 1).
- **rank_p_optim... .R**: *R* code for the optimization of the adjusted p-value and rank of the specific (optimization goals 2 and 3).

For the web-based applications *DAVID*, *GSEA*, and *GSEAPreranked*, we have prepared separate folders (of the same name as the respective method). For their structures see below. 

 -------------------------------
**DAVID**


-------------------------------
**GSEA (web-based application)**

-------------------------------
**GSEAPreranked**



-------------------------------
 **Generate the results illustrations as shown in the paper**
- folder **Results_illustrations**: contains the *R* codes to generate the results images from the paper (and entails the images themselves)
    - **Struktur muss noch beschrieben werden, wenn alles steht**
 
    - 








