# OverOptimism_in_GeneSetAnalysis
<<<<<<< HEAD
Repository for code and documentation in our analysis on over-optimism in gene set analysis.
=======
This repository allows you to 
- reproduce the results from our over-optimism study generated in *R* 
- inspect the documentation for the web-based applications

**Important note 1:** 

Reproducing all of the results takes a long time and can take up more memory than is available. To prevent *R* from crashing, we therefore recommend limiting yourself to reproducing a part of the results.

**Important note 2:**

Some of the considered gene set analysis methods are web-based applications (*DAVID*, *GSEA*, and *GSEAPreranked*). The optimisation processes for these methods were performed by hand (and documented with screenshots). Fully reproducing the results would therefore take months. 

**Important note 3:**

Following 'important note 2', I produced the results for these three web-based application over a period of several months, including the part in which I (partly) prepared the required input objects in *R*. Unfortunately, I was not aware of reproducible environments such as *renv* and therefore naively proceeded over the months without ensuring exact reproducibility by documenting the current versions of all packages needed in the process. However, I compared **many, but of course not all** of the *R* outputs generated using *renv* to those I generated at the time and they were very similar (for instance, for the rankings required as input to *GSEAPreranked* differed only from the third decimal place).  

## Reproduce the figures (in *R*)
To reproduce a figure from the paper, run the corresponding script from folder `R/Figures`.

Note that, while the *R* scripts source the intermediate results from the GSA methods implemented in *R* internally, the results for the web-based applications *GSEA*, *GSEAPreranked*, and *DAVID* were transferred from the corresponding screenshots **by hand** since the optimisation processes could not be run in *R*. 

## Rerun the full experiment (folder *R*)
Note that this takes several days or weeks, depending on the available resources. The following *R* scripts are stored in the folder *R*.

#### Preparation of the gene expression data for the optimisations (folder `R/Prepare_data_and_permutations`)
This is the first script to run when reproducing the results.

- **Random_Phenotype_Permutations.R**: *R* code to obtain
  1. both gene expression data sets and true labels of the respective samples 
  2. the ten random permutations of the true sample labels for both gene expression data sets
  
  The generated random phenotype permutations and the gene expression data set "Bottomly" are additionally provided in the folder **GeneExpression_data** in this repository. These are then sourced in the subsequent functions. 
 
#### Preparation of the gene expression data for the GSA methods (folder `R/Functions`)

 - **PreProcessing_Functions.R**: contains functions provided in the preprocessing of (almost) all of the investigated GSA methods and is therefore sourced in each of the respective R scripts.
 - **RNASeq_Transformation.R**: contains two functions to transform the gene expression measurements to match the characteristics for microarray data, as needed for the methods *PADOG* and *GSEA* (web-based application). It is sourced in the corresponding scripts for the methods.
 
#### Define optimisation functions (folder `R/Functions`)
The following scripts were generated for all methods implemented in *R*, i.e. *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, and *clusterProfiler*'s GSEA. These contain all functions required to perform the optimization for the respective computational GSA method. The scripts are then sourced when running the optimisations in the next step. 
  
- **n_DEGS_OptimisationFunctions_... .R**: Functions for the optimization of the number of differentially enriched gene sets for the respective computational GSA method (optimization goal 1).
- **rank_p_OptimisationFunctions_... .R**: Functions for the optimization of the adjusted p-value and rank of the specific (optimization goals 2 and 3).

For the web-based applications *DAVID*, *GSEA*, and *GSEAPreranked*, the optimisations are structured differently since the optimisations themselves are to be performed in the corresponding web-based application. Therefore, the corresponding *R* scripts to generate the required inputs are defined in folder `R/Functions`, while they do not appear in the *R* scripts in the folder `R/Run_optimisations` 

##### DAVID
For DAVID, a single



#### Run the optimisations for the *R*-based methods (folder `R/Run_optimisations`)
Run optimisations for *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, and *clusterProfiler*'s GSEA in the following scripts:

- **Run_n_DEGS_optimisations.R**: script to run optimisations of the number of differentially enriched gene sets
- **Run_pvalue_rank_optimisations.R**: script to run optimisations of the adjusted p-value and rank of the specific gene sets

Both scripts source all required functions (i.e. preparation and optimisation functions) from the previous sections internally.

#### Generate the results figures (folder `R/Figures`)
To reproduce a results figure from the main document or the supplement from the paper, run the corresponding *R* script.

# Optimisations for the web-based methods *DAVID*, *GSEA*, and *GSEAPreranked*

## DAVID
The web-based application *DAVID* can be accessed via the following link: (https://david.ncifcrf.gov/)
Our analysis was performed with the DAVID Knowledgebase v2023q3.

In your working directory, *generate the folder* *DAVID*

The documentation for DAVID is structured by both gene expression data sets (folders **Pickrell** and **Bottomly**). For both data sets, respectively, the required input objects for DAVID (for the true and the 10 permuted sample labels) are generate in the R scripts
- **GenerateInputs_DAVID_Pickrell.R**
- **GenerateInputs_DAVID_Bottomly.R**

The resulting input objects are then stored in the subfolder **Data** in of the respective folder Pickrell/Bottomly. **Note** that for each data set, the generated input objects are identical for the maximisation of the number of differentially enriched gene sets (**goal 1**) and the minimisation of the adjusted p-value/rank (**goals 2 and 3**). 


**Important:**
The documentation of the number of differentially enriched gene sets (stored in subfolder **n_DEGS**) shows that the number of differentially enriched gene sets could **NOT** be increased for any of the two gene expression data sets and none of the sample labels (neither true nor permuted). Indeed, there was only one case (Bottomly data set, true sample labels) in which there were gene sets with a significant adjusted p-value. 
There was one more optimization step for goal 1 than for goals 2 and 3 (step 3: gene set database KEGG; it never led to an increase in the number of differentially enriched gene set). The remaining optimisation steps were identical between goals 1, 2, and 3. 

We were therefore document the optimization steps for objectives 2 and 3 directly from documentation for objective 1. It (i.e., the Excel files for each optimisation step for goal 1) showed us that both gene expression data sets and for all of the sample labels, the adjusted p-value and rank of the respective gene sets 

- **Demethylation** (GO:0070988) and **t Cell mediated immunity** (GO:0002456) for the Pickrell data set,
- **Metabolic Process** (GO:0008152) and **Cellular Process** (GO:0009987) for the Bottomly data set,
  
could **never** be decreased for their default value of 1. 

## GSEA (web-based application)

The application can be downloaded from (https://www.gsea-msigdb.org/gsea/index.jsp), for which an account must be created.

In your working directory, *generate the folder* *GSEA*

The documentation for GSEA is structured by both gene expression data sets (folders **Pickrell** and **Bottomly**). Within each folder, you will find a folder 
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


## GSEAPreranked

In your working directory, generate the folder* *GSEAPreranked*

## Generate the results illustrations as shown in the paper
- folder **Results_illustrations**: contains the *R* codes to generate the results images from the paper (and entails the images themselves)
    - **Struktur muss noch beschrieben werden, wenn alles steht**
 
    - 






>>>>>>> cb9285866d9a885d30850ebb0dd258feee0701c6


