<!--# OverOptimism_in_GeneSetAnalysis
<<<<<<< HEAD -->
Repository for code and documentation for analysis on over-optimism in GSA.
=======
This repository allows you to reproduce the results from our over-optimism analysis on four *R*-based and three web-based GSA methods. 

### Important notes

**Note 1:** 

Reproducing all results **from scratch** is possible when using 

- identical data sets 
- identical seed(s)
- identical *R* package versions/ versions of the web applications
- identical gene set database versions
  
as used in the original experiments. This is straightforward for the *R*-based methods (using renv), however, it might take up more memory than available and cause your system to crash. For the web-based methods, reproducing the results requires several intermediate steps (see note 4) and is very time-consuming; detailed descriptions are provided below. 

**Note 2:**

The *R* results can be reproduced using *R* package *renv* and command *renv::restore()*. See https://rstudio.github.io/renv/articles/renv.html for information on how to use *renv*. 

**Note 3:**

I produced the results for the three web-based applications (*DAVID*, *GSEA*, and *GSEAPreranked*) over a period of several months, including the part in which I (partly) prepared the required input objects in *R*. Unfortunately, I was not aware of reproducible environments such as *renv* and therefore naively proceeded over the months without ensuring exact reproducibility by documenting the current versions of all packages needed in the process. However, I compared **many, but of course, not all** of the *R* outputs generated using *renv* to those I generated at the time and they were very similar (for instance, for the rankings required as input to *GSEAPreranked* differed only from the third decimal place).  

**Note 4:**

The web-based applications were performed by hand (and documented with screenshots) using the following three steps: 

1. In *R*: Generate input data sets for web application from initial gene expression data 
2. For *GSEA* and *GSEAPreranked*: Further pre-processing steps in Excel (for *GSEA* and *GSEAPreranked*)
3. Data is uploaded to the respective web application and optimisations are performed
4. Results of the optimisation processes are entered manually in *R* (scripts in `R/Code_figures`) for the generation of the results graphs

**Fully reproducing this process would take weeks!**

*** 
## Reproduce figures based on intermediate results
To reproduce a figure from the paper (based on intermediate results), run the corresponding script from folder `R/Code_figures`. The scripts are named after the figure they generate in the manuscript. 

Note that, while the *R* scripts source the intermediate results from the GSA methods implemented in *R* internally, the results for the web-based applications *GSEA*, *GSEAPreranked*, and *DAVID* were transferred from the corresponding screenshots **by hand** since the optimisation processes could not be run in *R*. 


***
## Reproduce figures from scratch
Note that this takes several days or weeks, depending on the available resources. The following *R* scripts are stored in the folder `R/Run_optimisations`.

### 1. Run optimisations for *R*-based GSA methods *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, and *clusterProfiler*'s GSEA in the following scripts:

- **Run_task1_optimisations.R**: required for Figures 3 and S2, 
- **Run_task2_3_optimisations.R**: required for Figures 4, 5, S3, and S4,  

where task 1, 2, and 3 are the optimisation tasks based on which over-optimism is assessed. Both scripts source all required functions and preprocess the initial gene expression data sets internally. 

For the remaining, web-based applications *DAVID*, *GSEA*, and *GSEAPreranked*, the optimisations were performed in the corresponding web-based application (i.e., by hand). To generate the results figures from the manuscript, the optimisation results were therefore also transferred manually to the corresponding *R* scripts (`R/Code_figures`). See below for instructions on how to reproduce the optimisation results for the web-based applications. 

### 2. Generate the results figures (folder `R/Code_figures`)
The *R* scripts to generate the results figures are named directly after the figure.  

*** 

## Detailed: reproduce analysis for web-based methods *DAVID*, *GSEA*, *GSEAPreranked*

Note that for the over-optimism analysis of the web-applications, we strongly recommend reading the manuscript as it requires several steps, most of which must be performed by hand. These steps are structured according to note 4 above and described in detail in the following. Note that for the web-based applications, only tasks 1 and 2 are pursued: 

- task 1: maximise number of differentially enriched gene sets
- task 2: minimise adjusted p-value of a specific gene set

### DAVID
The web-based application *DAVID* can be accessed via the following link: https://david.ncifcrf.gov/
Our analysis was performed with *DAVID* version 6.8 and the DAVID Knowledgebase v2023q3.

#### 1. Generation of input data sets (folder `R/Optimisation_functions`)
The *R* script to generate the input data sets for the Pickrell **and** the Bottomly data set is stored in the file  **generate_Inputs_DAVID.R**. 
  
The generated input data sets are then stored in the folders
- `Results/Intermediate_results/DAVID/Pickrell`
- `Results/Intermediate_results/DAVID/Bottomly`

and in there, in subfolders that corresponds to the given phenotypes (original phenotypes or phenotype permutations 1-10). 

Note that for the Pickrell and the Bottomly data set respectively, the input objects for the web application are identical across both tasks 1 and 2 since the individual optimisation steps are (almost) identical. 

#### 2. Run DAVID optimisation 
Access the link https://david.ncifcrf.gov/. An input list generated in step 1 can be uploaded by clicking on **Start Analysis** and submitting the input gene list under `Step 1: Enter Gene List`. Select the identifier as **Ensembl_Gene_ID** and set the list type as **Gene List**. 

The optimisation steps are the following (see corresponding screenshots for more detailed illustrations): 

1. Upload default list of differentially expressed genes (generated using DESeq2; gene set database GO (BP)): **DEGs_DESeq2_... .txt**
2. Change method to detect differentially expressed genes: upload **DEGs_limma_... .txt**
3. **For task 1 only**: change gene set database to KEGG 
4. Change set of background genes: upload **universe_DESeq2_... .txt** or **universe_limma_... .txt** according to whether step 1 or step 2 leads to better results w.r.t the optimisation task).

In each step, you compare the results w.r.t. the given task to the results from the previous step and choose the option yielding better results. 

#### 3. Inspect documentation of the results (folder `Results/Screenshots_WebApplications`)
For task 1, you find the documentation screenshots (and Excel results files) in the folder 

- **Screenshots_DAVID_Pickrell** for the Pickrell data set 
- **Screenshots_DAVID_Bottomly** for the Bottomly data set 

**Important:**
The documentation shows that the number of differentially enriched gene sets could **NOT** be increased for any of the two gene expression data sets and none of the sample labels (neither true nor permuted). Indeed, there was only one case (Bottomly data set, true sample labels) in which there were gene sets with a significant adjusted p-value. 
There was one more optimisation step for task 1 than for task 2 (step 3: gene set database KEGG; it never led to an increase in the number of differentially enriched gene sets). The remaining optimisation steps were identical between tasks 1 and 2. 

We were therefore able to document the optimisation steps for task 2 directly from the documentation for task 1. It (i.e., the Excel files for each optimisation step for task 1) showed us that both gene expression data sets and for all of the sample labels, the adjusted p-value and rank of the respective gene sets.

- **Demethylation** (GO:0070988) and **t Cell mediated immunity** (GO:0002456) for the Pickrell data set,
- **Metabolic Process** (GO:0008152) and **Cellular Process** (GO:0009987) for the Bottomly data set,
  
could **never** be decreased for their default value of 1. 

***
### GSEA (web-based application)

The application can be downloaded from https://www.gsea-msigdb.org/gsea/index.jsp, for which an account must be created. I performed the experiment on versions 4.2.2/4.2.3 (the exact version can be found on the corresponding screenshot). 

#### 1. Generation of inputs in *R* (folder `R/Optimisation_functions/GSEA_Web`)
For the Pickrell and the Bottomly data set each, you will find the following three *R* scripts:

- **task1_optimisation_GSEAWeb_... .R**: input generation for task 1 
- **task2_optimisation_Demethylation_GSEAWeb_... .R**: input generation task 2; gene set **Demethylation** 
- **task2_optimisation_tCell_GSEAWeb_... .R**: Input generation task 2;  gene set **t Cell mediated immunity**.  

The input data generated from each of the scripts is stored in folder `Results/Intermediate_results/GSEA_Web/...Raw`. These contain the gene expression measurements as well as the phenotype assignments (which are stored in subfolder **Phenotypes**). 

#### 2. Further preprocessing in Excel

These files stored in `Results/Intermediate_results/GSEA_Web/...Raw` require further processing according to https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats. For this, inspect Sections 

- **GCT: Gene Cluster Text file format (*.gct)** for the preprocessing of the gene expression data set in Excel.
- **CLS: Categorical (e.g. tumor vs normal) class file format (*.cls)** for the preprocessing of the corresponding phenotype assignments.

The corresponding preprocessed files are then stored in subfolder `Prep`.

#### 3. Upload data and make optimisations according to task 1 or 2 
The application must be downloaded from https://www.gsea-msigdb.org/gsea/index.jsp. In the application, the required data sets are uploaded in the tab **Load data**, after which you need to proceed to the tab **Run GSEA**. Further information on the necessary fields to be clicked and filled out can be obtained from the screenshots in folder `Results/Intermediate_results/GSEA_Web`. Note that in each step, the random seed is set to **149**. 


Optimisation steps for tasks 1 and 2 : 

1. Upload pre-processed gene expression data set and phenotype assignments: upload **exprdat_default_... .gct**
2. Change RNA-Seq transformation method: upload **exprdat_vst_... .gct**
3. Change pre-filtering method: upload **exprdat_vst_filterByExpr_... .gct**
4. **For task 1 only**: change gene set database to KEGG ("c2.cp.kegg...")
5. Change metric for ranking genes 
    - change to **tTest**
    - change to **Diff_of_classes**
6. Change enrichment statistic 
    - change to **classic**
    - change to **weighted_p1.5**
    - change to **weighted_p2**
  
For each optimisation step, the results of the GSEA run can be found in the provided report tab: 

- For task 1, you add up the number of gene sets reported as "significant at FDR < 0.25" under both tabs "Enrichment in phenotype: ..."
- For task 2, you get the adjusted p-value of the given gene set by transferring both "Detailed enrichment results in html format" to Excel and inspecting column "FDR_qval"

In each optimisation step, you compare the results w.r.t. the given task to the results from the previous step and choose the option yielding better results. 

#### 4. Inspect documentation of the results (folder `Results/Screenshots_WebApplications`)

Note that the optimisation of GSEA was carried out over several months and by hand. In this period, the application was updated as well as the gene set database GO (with subontology biological process). The versions that were current in the respective optimisation processes can be found in the screenshots containing the name "Param", namely 
- in the top left corner for the version of the overall web application
- in the specified gene set database (tab **Gene Sets database**) for the version of the gene set database

In addition to the screenshots, the optimisation processes are documented in the respective *R* scripts as comments. 

***
### GSEAPreranked

*GSEAPreranked* is a variant of the above-described web-based method GSEA and can therefore be accessed via the same application (download from https://www.gsea-msigdb.org/gsea/index.jsp). The experiment was performed on versions 4.2.2 and 4.2.3


#### 1. Generation of inputs (folder `R/Optimisation_functions/GSEAPreranked`)
For the Pickrell and the Bottomly data set each, you will find the following *R* scripts:

- **generateInputs_optimisation_GSEAPreranked_task1_... .R**: Input generation for task 1 for the Pickrell and Bottomly data set + additional documentation of the optimisation procedure as comments (documentation could not be placed in separate .txt files since for some optimisation steps, the options depend on the previous step(s)). 
- **generateInputs_GSEAPreranked_task2_... .R**: Input generation for task 2 (additionally, the corresponding optimisation documentations are stored in text files in this folder).

#### 2. Further preprocessing in Excel 

**Note** that further preprocessing must be performed in Excel according to section **RNK: Ranked list file format (*.rnk)** in https://www.gsea-msigdb.org/gsea/index.jsp. 

#### 3. Upload data and make optimisations according to task 1 or 2 
The required data sets are uploaded in the tab **Load data**, after which you need to proceed to the tab **Run GSEAPreranked**. Further information on the necessary fields to be clicked and filled out can be obtained from the screenshots in folder `Results/Screenshots_WebApplications`. Note that in each step, the random seed is set to **149**. 

Optimisation steps for tasks 1 and 2 : 

1. Upload default gene ranking: upload **DESeq2_ranking_ ... .rnk**
2. Change method to rank the genes: upload **limma_ranking_ ... .rnk**
3. **For task 1 only**: change gene set database to KEGG ("c2.cp.kegg...")
4. Change enrichment statistic 
    - change to **classic**
    - change to **weighted_p1.5**
    - change to **weighted_p2**
  
In each step, you compare the results w.r.t. the given task to the results from the previous step and choose the option yielding better results. 

#### 4. Transfer results to *R* files to generate results figures 

The screenshot documentation is placed in folders 
- **Screenshots_GSEAPreranked_..._task1** for task 1 ("...": Bottomly or Pickrell data set)
- **Screenshots_GSEAPreranked_... task2_GENESET** or task 2 ("...": Bottomly or Pickrell data set)

The results from these screenshots are then transferred to the corresponding *R* scripts in folder **R/Code_figures**. 

Note that the optimisation of GSEAPreranked was carried out over several months and by hand. In this period, the application was updated as well as the gene set database GO (with subontology biological process). The versions that were current in the respective optimisation processes can be found in the screenshots containing the name "Param", namely 
- in the top left corner for the version of the overall web application
- in the specified gene set database (tab **Gene Sets database**) for the version of the gene set database

In addition to the screenshots, the optimisation processes are documented in the respective *R* scripts as comments for task 1 (**task1_optimisation_GSEAPreranked_ ... .R** in folder `R/Optimisation_functions/GSEAPreranked`) and in .txt files **task2_documentation_GSEAPreranked_ ... .txt** for task 2).  


<!--

## Description of additional *R* scripts 

### Folder `R/Functions`

 - **PreProcessing_Functions.R**: contains functions required in the preprocessing of (almost) all of the investigated GSA methods and is therefore sourced in each of the respective *R* scripts.
 - **RNASeq_Transformation.R**: contains two functions to transform the gene expression measurements to match the characteristics for microarray data, as needed for the methods *PADOG* and *GSEA* (web-based application). It is sourced in the corresponding scripts for the optimisations of the methods. 
 - **task1_OptimisationFunctions_... .R**: Functions for the optimisation of the number of differentially enriched gene sets for the respective computational GSA method (task 1; only for *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, *clusterProfiler*'s GSEA).
- **task2_3_OptimisationFunctions_... .R**: Functions for the optimisation of the adjusted p-value and rank of the specific (tasks 2 and 3; only for *GOSeq*, *clusterProfiler*'s ORA, *PADOG*, *clusterProfiler*'s GSEA).

-->




