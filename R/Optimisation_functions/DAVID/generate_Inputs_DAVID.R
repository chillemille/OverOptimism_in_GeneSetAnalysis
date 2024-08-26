################################################################################
### Generate the input objects (i.e. lists of differentially enriched gene sets)
### based on the Pickrell data set for the web-based application DAVID #########
################################################################################


library(DESeq2)
library(limma)
library(edgeR) # for filterByExpr()
library(org.Hs.eg.db)
library(dplyr)



# the link to run DAVID is the following: https://david.ncifcrf.gov/

# If not already done, specify working directory using setwd()

# load gene expression data set with true phenotype randomly permuted phenotype assignments
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")

# load functions required for data preprocessing
source("./R/Help_functions/PreProcessing_Functions.R")


### Generate the folders required to store the inputs

name_data <- c("Pickrell", "Bottomly")
for(j in 1:2){

# create folder to contain everything associated with the Pickrell data set
dir.create(paste0("./R/Optimisation_functions/DAVID/", name_data[j]))

# create folder to contain all data
dir.create("./R/Optimisation_functions/DAVID/", name_data[j], "/Data")

# create folder to contain the data for the true sample labels
dir.create("./R/Optimisation_functions/DAVID/", name_data[j], "/Data/Phen_Original")

for(i in 1:10){

  path_data_permutations <- paste0("./R/Optimisation_functions/DAVID/", name_data[j], "/Data/Permutation", i)
  dir.create(path_data_permutations)

}

}


#########################################
###generate necessary input for DAVID ###
#########################################

DAVID_input_preparation <- function(DE_results){

  #required input for clusterProfiler function: vector of entrez gene ID
  #-> need to pre-process results table DE_results

  #vector of differentially expressed genes
  #DEG_vec serves as input vector for ORA performed by clusterProfiler



  # classify those genes as DE that have an adjusted p-value < 0.05
  DEG_vec<-rownames(DE_results[(DE_results$p_adj < 0.05) & (!is.na(DE_results$p_adj)),])


  # return vector of differentially expressed genes
  return(DEG_vec)

}


################################################################################
### Generate gene lists (input object) for original Phenotype ##################
################################################################################

# j=1: Pickrell data set, j=2: Bottomly data set
for(j in 1:2){

  # get unprocessed gene expression data set
  raw_data <- eval(parse(text = paste0("Biobase::exprs(", tolower(name_data[j]), ".eset)")))

  # building block to access the sample conditions
  condition <- ifelse(j==1, "gender", "strain")

  # access the true sample conditions
  true_phen_labels <- eval(parse(text = paste0(tolower(name_data[j]), ".eset$", condition)))
  # access the permuted sample conditions
  perm_phen_labels <- eval(parse(text = paste0("phen_", tolower(name_data[j]))))


# generate DESeq2 results and rename column for adjusted p-values
DESeq2_results_phenorig <- pre_filt(raw_data, threshold=10) %>%
  deseq_preprocess(phenotype_labels = true_phen_labels) %>%
  DESeq() %>% results() %>%
  as.data.frame() %>% dplyr::rename(p_adj=padj)

# from DESeq2 results, generate list of differentially expressed genes
DEGs_DESeq2_phenorig <- DAVID_input_preparation(DESeq2_results_phenorig)

# count number of differentially expressed genes
n_DEGs_DESeq2_phenorig <- length(DEGs_DESeq2_phenorig)



# create path for storage of DESeq2 results
path_DEGs_deseq2_phenorig <- paste0("./R/Optimisation_functions/DAVID/", name_data[j], "/Data/Original_Phenotype/DEGs_DESeq2_phenorig.txt")

# export DESeq2 results
write.table(DEGs_DESeq2_phenorig,
            file = path_DEGs_deseq2_phenorig,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


### Universe ###



# generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
universe_DESeq2_phenorig <- rownames(DESeq2_results_phenorig[!is.na(DESeq2_results_phenorig$p_adj),])


# create path for storage of DESeq2 results
path_universe_deseq2_phenorig <- paste0("./R/Optimisation_functions/DAVID/",name_data[j],"/Data/Original_Phenotype/universe_DESeq2_phenorig.txt")

# export alternative universe
write.table(universe_DESeq2_phenorig,
            file = path_universe_deseq2_phenorig,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


############
### limma
############

### limma results ###


# generate list of differentially expressed genes using limma

# (i.1) generate pre-filtering indicator using edgeR's filterByExpr
ind_filt_phenorig <- DGEList(raw_data, group= true_phen_labels) %>% filterByExpr()

# (i.2) generate design matrix
mm_phenorig <-  model.matrix(~ true_phen_labels)

# (ii) run limma workflow and rename column containing adjusted p-values
limma_results_phenorig <- DGEList(raw_data[ind_filt_phenorig,], group= true_phen_labels) %>% calcNormFactors() %>%
  voom(design=mm_phenorig) %>% lmFit(design=mm_phenorig) %>% eBayes() %>%
  topTable(coef=ncol(mm_phenorig), number=60000) %>%
  as.data.frame() %>% dplyr::rename(p_adj=adj.P.Val)



# from limma results, generate list of differentially expressed genes
DEGs_limma_phenorig <- DAVID_input_preparation(limma_results_phenorig)

# count number of differentially expressed genes
n_DEGs_limma_phenorig <- length(DEGs_limma_phenorig)



# create path for storage of limma results
path_DEGs_limma_phenorig <- paste0("./R/Optimisation_functions/DAVID/", name_data[j], "/Data/Original_Phenotype/DEGs_limma_phenorig.txt")

# export DESeq2 results
write.table(DEGs_limma_phenorig,
            file = path_DEGs_limma_phenorig,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


# generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
universe_limma_phenorig <- rownames(limma_results_phenorig[!is.na(limma_results_phenorig$p_adj),])


# create path for storage of limma results
path_universe_limma_phenorig <- paste0("./R/Optimisation_functions/DAVID/",name_data[j],"/Data/Original_Phenotype/universe_limma_phenorig.txt")

# export alternative universe
write.table(universe_limma_phenorig,
            file = path_universe_limma_phenorig,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


################################################################################
### Random Phenotype Permutations ##############################################
################################################################################

# for double-checking: store number of differentially enriched gene set for each
# permutation of the true sample labels

n_DEGs_limma_phenperm <- c()
n_DEGs_DESeq2_phenperm <- c()

for(i in 1:ncol(perm_phen_labels)){


  ############
  ### DESeq2
  ############

  ### DESeq2 results ###


  # generate DESeq2 results and rename column for adjusted p-values
  DESeq2_results <- pre_filt(raw_data, threshold=10) %>%
    deseq_preprocess(phenotype_labels = perm_phen_labels[,i]) %>%
    DESeq() %>% results() %>%
    as.data.frame() %>% dplyr::rename(p_adj=padj)

  # from DESeq2 results, generate list of differentially expressed genes
  DEGs_DESeq2 <- DAVID_input_preparation(DESeq2_results)

  # count number of differentially expressed genes
  n_DEGs_DESeq2_phenperm[i] <- length(DEGs_DESeq2)



  # create path for storage of DESeq2 results
  path_DEGs_deseq2 <- paste0("./R/Optimisation_functions/DAVID/",name_data[j],"/Data/Permutation", i, "/DEGs_DESeq2_permutation",i,".txt")

  # # export DESeq2 results
  write.table(DEGs_DESeq2,
              file = path_DEGs_deseq2,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)


  ### Universe ###

  # generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
  universe_DESeq2 <- rownames(DESeq2_results[!is.na(DESeq2_results$p_adj),])


  # create path for storage of DESeq2 results
  path_universe_deseq2 <- paste0("./R/Optimisation_functions/DAVID/",name_data[j],"/Data/Permutation", i, "/universe_DESeq2_permutation",i,".txt")

  # export alternative universe
  write.table(universe_DESeq2,
              file = path_universe_deseq2,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)


  ############
  ### limma
  ############

  ### limma results ###


  # generate list of differentially expressed genes using limma

  # (i.1) generate pre-filtering indicator using edgeR's filterByExpr
  ind_filt <- DGEList(raw_data, group= perm_phen_labels[,i]) %>%
    filterByExpr()

  # (i.2) generate design matrix
  mm <-   mm <- model.matrix(~ perm_phen_labels[,i])

  # (ii) run limma workflow and rename column containing adjusted p-values
  limma_results <- DGEList(counts = raw_data[ind_filt,],group= perm_phen_labels[,i]) %>%
    calcNormFactors() %>% voom(design=mm) %>% lmFit(design=mm) %>% eBayes() %>%
    topTable(coef=ncol(mm), number=60000) %>%
    as.data.frame() %>% dplyr::rename(p_adj=adj.P.Val)



  # from limma results, generate list of differentially expressed genes
  DEGs_limma <- DAVID_input_preparation(limma_results)

  # count number of differentially expressed genes
  n_DEGs_limma_phenperm[i] <- length(DEGs_limma)



  # create path for storage of limma results
  path_DEGs_limma <- paste0("./R/Optimisation_functions/DAVID/",name_data[j],"/Data/Permutation", i, "/DEGs_limma_permutation",i,".txt")

  # export DESeq2 results
  write.table(DEGs_limma,
              file = path_DEGs_limma,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)


  # generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
  universe_limma <- rownames(limma_results[!is.na(limma_results$p_adj),])


  # create path for storage of limma results
  path_universe_limma <- paste0("./R/Optimisation_functions/DAVID/name_data[j]/Data/Permutation", i, "/universe_limma_permutation",i,".txt")

  # export alternative universe
  write.table(universe_limma,
              file = path_universe_limma,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)





}

}











