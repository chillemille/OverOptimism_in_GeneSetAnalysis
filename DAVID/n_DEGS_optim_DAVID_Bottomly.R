#####################################################################################################
### Maximization of number of differentially enriched gene sets obtained by DAVID for Bottomly data #
#####################################################################################################

library(DESeq2)
library(limma)
library(edgeR) # for filterByExpr()
library(dplyr)



# the link to run DAVID is the following: https://david.ncifcrf.gov/

# set working directory
setwd("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/Assessment_OverOptimism")

# load gene expression data set with true phenotype randomly permuted phenotype assignments 
source("./Random_Phenotype_Permutations.R")

# load functions required for data preprocessing 
source("./PreProcessing_Functions.R")


############################################
###generate necessary input for ORA tool ###
############################################

# (i) for all Differential expression results created all DE techniques apart from NOISeq 
DAVID_input_preparation <- function(DE_results){
  
  #required input for clusterProfiler function: vector of entrez gene ID
  #-> need to pre-process results table DE_results
  
  #vector of differentially expressed genes 
  #DEG_vec serves as input vector for ORA performed by clusterProfiler 
  
  
  
  # classify those genes as DE that have an adjusted p-value < 0.05
  DEG_vec<-rownames(DE_results[(DE_results$p_adj<0.05) & (!is.na(DE_results$p_adj)),])
  
  
  # return vector of differentially expressed genes 
  return(DEG_vec)
  
}


################################################################################
### Original Phenotype Assignment ##############################################
################################################################################


# for both DE methods, create vector which contains number of DIFFERENTIALLY EXPRESSED GENES for each random phenotype permutation


  
  ############
  ### DESeq2
  ############
  
  ### DESeq2 results ###
  
  
  # generate DESeq2 results and rename column for adjusted p-values
  DESeq2_results_phenorig <- pre_filt(Biobase::exprs(bottomly.eset), threshold=10) %>% 
                             deseq_preprocess(phenotype_labels = bottomly.eset$strain) %>%
                             DESeq() %>% results() %>% 
                             as.data.frame() %>% dplyr::rename(p_adj=padj)
  
  # from DESeq2 results, generate list of differentially expressed genes 
  DEGs_DESeq2_phenorig <- DAVID_input_preparation(DESeq2_results_phenorig)
  
  # for an overview: count number of differentially expressed genes 
  n_DEGs_DESeq2_phenorig <- length(DEGs_DESeq2_phenorig)

  # create path for storage of DESeq2 results 
  path_DEGs_deseq2_phenorig <- paste0("./DAVID/Bottomly/n_DEGS/Data/Phen_Original/DEGs_DESeq2_phenOriginal.txt")
  
  # # export DESeq2 results 
  write.table(DEGs_DESeq2_phenorig,
              file = path_DEGs_deseq2_phenorig,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  
  ### Universe ###
  
  
  
  # generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
  universe_DESeq2_phenorig <- rownames(DESeq2_results_phenorig[!is.na(DESeq2_results_phenorig$p_adj),])
  
  
  # create path for storage of DESeq2 results 
  path_universe_deseq2_phenorig <- paste0("./DAVID/Bottomly/n_DEGS/Data/Phen_Original/universe_DESeq2_phenOrig.txt")
  
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
  ind_filt_phenorig <- DGEList(Biobase::exprs(bottomly.eset), group= bottomly.eset$strain) %>% filterByExpr()
  
  # (i.2) generate design matrix 
  mm_phenorig <- model.matrix(~ bottomly.eset$strain)
  
  # (ii) run limma workflow and rename column containing adjusted p-values 
  limma_results_phenorig <- DGEList(counts = Biobase::exprs(bottomly.eset)[ind_filt_phenorig,], group= bottomly.eset$strain) %>% 
                            calcNormFactors() %>% voom(design=mm_phenorig) %>% lmFit(design=mm_phenorig) %>% eBayes() %>% 
                            topTable(coef=ncol(mm_phenorig), number=60000) %>% 
                            as.data.frame() %>% dplyr::rename(p_adj=adj.P.Val)
  
  
  
  # from limma results, generate list of differentially expressed genes 
  DEGs_limma_phenorig <- DAVID_input_preparation(limma_results_phenorig)
  
  # for an overview: count number of differentially expressed genes 
  n_DEGS_limma_phenorig <- length(DEGs_limma_phenorig)
  
  
  # create path for storage of limma results 
  path_DEGs_limma_phenorig <- paste0("./DAVID/Bottomly/n_DEGS/Data/Phen_Original/DEGs_limma_phenOrig.txt")
  
  # export DESeq2 results 
  write.table(DEGs_limma_phenorig,
              file = path_DEGs_limma_phenorig,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  
  # generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
  universe_limma_phenorig <- rownames(limma_results_phenorig[!is.na(limma_results_phenorig$p_adj),])
  
  
  # create path for storage of limma results 
  path_universe_limma_phenorig <- paste0("./DAVID/Bottomly/n_DEGS/Data/Phen_Original/universe_limma_phenOrig.txt")
  
  # export alternative universe 
  write.table(universe_limma_phenorig,
              file = path_universe_limma_phenorig,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  
  
##################
# Step 1: Default 
##################
  
# -> adj p-value = 1 
  
  
  
#################
# Step 2: Gene List obtained with limma 
#################
  
# -> adj. p-value = 1 
# ->> return to default list of differntially expressed genes generated using DESeq2
  
  

#####################################
# Step 3: Upload alternative universe
#####################################
  
# -> 0 DEGS 
  
### final results: 2 DEGS -> default configuration corresponds to optimal
# configuration 
  
  
################################################################################
### Random Phenotype Permutations ##############################################
################################################################################


# for both DE methods, create vector which contains number of DIFFERENTIALLY EXPRESSED GENES for each random phenotype permutation
n_DEGs_DESeq2_phenperm  <- c()
n_DEGs_limma_phenperm <- c()

for(i in 1:ncol(phen_bottomly)){
  
  
  ############
  ### DESeq2
  ############
  
  ### DESeq2 results ###
  
  
  # generate DESeq2 results and rename column for adjusted p-values
  DESeq2_results <- pre_filt(Biobase::exprs(bottomly.eset), threshold=10) %>% 
                    deseq_preprocess(phenotype_labels = phen_bottomly[,i]) %>%
                    DESeq() %>% results() %>% 
                    as.data.frame() %>% dplyr::rename(p_adj=padj)
  
  # from DESeq2 results, generate list of differentially expressed genes 
  DEGs_DESeq2 <- DAVID_input_preparation(DESeq2_results)
  
  # count number of differentially expressed genes 
  n_DEGs_DESeq2_phenperm[i] <- length(DEGs_DESeq2)
  
  
  
  # create path for storage of DESeq2 results 
  path_DEGs_deseq2 <- paste0("./DAVID/Bottomly/n_DEGS/Data/Permutation", i, "/DEGs_DESeq2_permutation",i,".txt")
  
  # export DESeq2 results
  write.table(DEGs_DESeq2,
              file = path_DEGs_deseq2,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  
  ### Universe ###
  
  # generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
  universe_DESeq2 <- rownames(DESeq2_results[!is.na(DESeq2_results$p_adj),])
  
  
  # create path for storage of DESeq2 results 
  path_universe_deseq2 <- paste0("./DAVID/Bottomly/n_DEGS/Data/Permutation", i, "/universe_DESeq2_permutation",i,".txt")
  
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
  ind_filt <- DGEList(Biobase::exprs(bottomly.eset), group= phen_bottomly[,i]) %>% filterByExpr()
  
  # (i.2) generate design matrix 
  mm <-   mm <- model.matrix(~ phen_bottomly[,i])
  
  # (ii) run limma workflow and rename column containing adjusted p-values 
  limma_results <- DGEList(counts = Biobase::exprs(bottomly.eset)[ind_filt,], group= phen_bottomly[,i]) %>% 
                   calcNormFactors() %>% voom(design=mm) %>% lmFit(design=mm) %>% eBayes() %>% 
                   topTable(coef=ncol(mm), number=60000) %>% 
                   as.data.frame() %>% dplyr::rename(p_adj=adj.P.Val)
  
  
  
  # from limma results, generate list of differentially expressed genes 
  DEGs_limma <- DAVID_input_preparation(limma_results)
  
  # count number of differentially expressed genes 
  n_DEGs_limma_phenperm[i] <- length(DEGs_limma)
  
  
  
  # create path for storage of limma results 
  path_DEGs_limma <- paste0("./DAVID/Bottomly/n_DEGS/Data/Permutation", i, "/DEGs_limma_permutation",i,".txt")
  
  # export DESeq2 results
  write.table(DEGs_limma,
              file = path_DEGs_limma,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  
  # generate alternative universe (i.e. all those genes from the experiment with a non-NA adjusted p-value)
  universe_limma <- rownames(limma_results[!is.na(limma_results$p_adj),])
  
  
  # create path for storage of limma results 
  path_universe_limma <- paste0("./DAVID/Bottomly/n_DEGS/Data/Permutation", i, "/universe_limma_permutation",i,".txt")
  
 #  export alternative universe
  write.table(universe_limma,
              file = path_universe_limma,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)


}




################################################################################
### Phenotype Permutations 1-9 #################################################
################################################################################

# None of the lists of differentially expressed genes, whether generated with 
# DESeq2 nor limma, contained any genes 

# -> DAVID would therefore always return 0 differentially enriched gene sets 

################################################################################
### Phenotype Permutation 10 ###################################################
################################################################################



##################
# Step 1: Default 
##################

# -> 0 DEGS 



#################
# Step 2: Gene List obtained with limma 
#################

# the list of differentially expressed genes contained 0 differentially expressed 
# gene sets when being generated with limma 

# -> skip step 2 and proceed directly to step 3 


#####################################
# Step 3: Upload alternative universe
#####################################

# -> 0 DEGS 
# ->> final results: 0 DEGS 







