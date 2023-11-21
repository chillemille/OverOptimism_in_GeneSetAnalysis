# OverOptimism_in_GeneSetAnalysis
<<<<<<< HEAD
Repository for code and documentation in our analysis on over-optimism on gene set analysis.
=======
This repository allows you to 
- reproduce the results from our over-optimism study generated in R 
- inspect the documentation for the web-based applications

@Moritz: Kann man annehmen, dass der Leser weiß, dass das working directory dann genau den Aufbau haben soll wie das Repo auf Github?

**Note** that all of the scripts are based on the working directory you have to specifiy at the beginning.

## **Preparation of the gene expression data for the optimizations**

This is the first script to run when reproducing the results.

- **Random_Phenotype_Permutations.R**: *R* code to obtain
  1. both gene expression data sets and true labels of the respective samples 
  2. the ten random permutations of the true sample labels for both gene expression data sets
  
  The generated random phenotype permutations and the gene expression data set "bottomly" are additionally provided in the folder **GeneExpression_data** in this repository. These are then sourced in the subsequent functions. 
 
## **Preparation of the gene expression data for the GSA methods**

 - **PreProcessing_Functions.R**: contains functions provided in the preprocessing of (almost) all of the investigated GSA methods and is therefore sourced in each of the respective R scripts.
 - **RNASeq_Transformation.R**: contains two functions to transform the gene expression measurements to match the characteristics for microarray data, as needed for the methods *PADOG* and *GSEA* (web-based application). It is sourced in the corresponding scripts for the methods.
 
## **Define functions required to run the optimisations **
The following scripts were generated for all methods implemented in *R*, i.e. *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, and *clusterProfiler*'s GSEA. These contain all functions required to perform the optimization for the respective computational GSA method. The scripts are then sourced when running the optimisations in the next step. 
  
- **n_DEGS_OptimisationFunctions_... .R**: Functions for the optimization of the number of differentially enriched gene sets for the respective computational GSA method (optimization goal 1).
- **rank_p_OptimisationFunctions_... .R**: Functions for the optimization of the adjusted p-value and rank of the specific (optimization goals 2 and 3).

For the web-based applications *DAVID*, *GSEA*, and *GSEAPreranked*, we have prepared separate folders (of the same name as the respective method). For their structures see below. 

# **Run the optimisations **
Run optimisations for *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, and *clusterProfiler*'s GSEA based on the defined functions: 

- **Run_n_DEGS_optimisations.R**: Run optimisations of the number of differentially enriched gene sets
- **Run_pvalue_rank_optimisations.R**: Run optimisations of the adjusted p-value and rank of the specific gene sets


## **DAVID**


## **GSEA (web-based application)**

The application can be downloaded from (https://www.gsea-msigdb.org/gsea/index.jsp), for which an account must be created.

The documentation for the web-based application is structured by both gene expression data sets (folders **Pickrell** and **Bottomly**). Within each folder, you will find a folder 
- **n_DEGS**: Contains data and documentation for the maximization of the number of differentially enriched gene sets
- **p_adj**: Contains data and documentation for the minimization of the adjusted p-values of two gene sets 

In the folders **n_DEGS** and **p_adj**, there are the following subfolders: 
- **Data**: Contains for the true and the 10 permuted sample labels the ...
  - ... raw, i.e. exported from R, data sets in folder **Raw**
  - ... prepared data sets in the format required by the web application in folder **Prep**
- **Screenshots**: Contains the screenshots of the progression of the optimizatino process for the true sample labels and the 10 permutations 

Additionally, both folders (Raw and Prep) contain a subfolder **Phenotypes** in which the sample labels in their raw (i.e. as exported from *R*) and prepared formats, respectively. 

The preparation of the gene expression data sets was carried out in **Excel** and according to (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats), namely as in Sections 
 - **GCT: Gene Cluster Text file format (*.gct)** for the gene expression measurements
 - **RNK: Ranked list file format (*.rnk)** for the sample labels

Note that the optimization of GSEA was carried out over several months and by hand. In this period, the application was updated as well as the gene set database GO (with subontology biological process). The versions that were current in the respective optimization processes can be found in the screenshots containing the name "Param", namely 
- in the top left corner for the version of the overall web application
- in the specified gene set database (tab **Gene Sets database**) for the version of the gene set database

  In addition to the screenshots, the optimization processes are documented in the respective *R* scripts. 


## **GSEAPreranked**


## **Generate the results illustrations as shown in the paper**
- folder **Results_illustrations**: contains the *R* codes to generate the results images from the paper (and entails the images themselves)
    - **Struktur muss noch beschrieben werden, wenn alles steht**
 
    - 






>>>>>>> cb9285866d9a885d30850ebb0dd258feee0701c6


