#######################################################################################
### Maximize number of differentially enriched gene sets obtained with GSEAPreranked ##
#######################################################################################


library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(apeglm)

# load pickrell data set with sample labels (true and permuted)
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")

# load required pre-processing functions
source("./R/Help_functions/PreProcessing_Functions.R")

# load help functions for GSEAPreranked
source("./R/Help_functions/helpfunctions_GSEAPreranked.R")


######################################
### generate required folders ########
######################################

dir.create("./Results/Intermediate_results/GSEAPreranked")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Raw")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Prep")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Prep/Original_Phenotype")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Raw/Original_Phenotype")

for(i in 1:10){

  path_raw <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Raw/Phenotype_Permutation",
                     i)

  path_prep <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Prep/Phenotype_Permutation",
                      i)

  dir.create(path_raw)
  dir.create(path_prep)


}


################################################################################
### Export gene rankings #######################################################
################################################################################


# (I) Original Phenotype Assignment

################
# DESeq2 ranking
################

# create DESeq2 results and rank by p-value
DESeq2_ranking_phenorig <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% geneID_conversion_SYMBOL_pickrell(dupl_removal_method = 1) %>%
  deseq_preprocess(phenotype_labels = pickrell.eset$gender ) %>% DESeq() %>%
  lfcShrink(coef = "condition_treated_vs_untreated", type = "apeglm") %>% as.data.frame() %>%
  rankedList_cP(rankby = "p_value", method = "DESeq2")

# create path for storage of DESeq2 ranking
path_DESeq2_phenorig <- "./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Raw/Original_Phenotype/DESeq2_ranking_phenOrig.txt"

# export
write.table(DESeq2_ranking_phenorig,
            file = path_DESeq2_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = FALSE)


##############
#limma ranking
##############

# filtering indicator
keep_phenorig <- DGEList(Biobase::exprs(pickrell.eset), group = pickrell.eset$gender) %>% filterByExpr()

# design matrix
mm_phenorig <- model.matrix( ~ pickrell.eset$gender)

# create limma results and rank by p-value
limma_ranking_phenorig <- geneID_conversion_SYMBOL_pickrell(Biobase::exprs(pickrell.eset)[keep_phenorig, ], dupl_removal_method = 1) %>%
  DGEList(group = pickrell.eset$gender) %>% calcNormFactors() %>%
  voom(design = mm_phenorig) %>% lmFit(design = mm_phenorig) %>% eBayes() %>% topTable(coef = ncol(mm_phenorig), number = 100000) %>%
  rankedList_cP(rankby =  "p_value", method = "limma")

# Create path for storage of limma ranking
path_limma_phenorig <- "./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Raw/Original_Phenotype/limma_ranking_phenOrig.txt"


# export
write.table(limma_ranking_phenorig,
            file = path_limma_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = FALSE)



# (II) Random Phenotype Permutations of Original Phenotypes

for(i in 1:ncol(phen_pickrell)){

  ################
  # DESeq2 ranking
  ################

  # create DESeq2 results and rank by p-value
  DESeq2_ranking_phenperm <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% geneID_conversion_SYMBOL_pickrell(dupl_removal_method = 1) %>%
    deseq_preprocess(phenotype_labels = phen_pickrell[, i] ) %>% DESeq() %>%
    lfcShrink(coef = "condition_treated_vs_untreated", type = "apeglm") %>% as.data.frame() %>%
    rankedList_cP(rankby = "p_value", method = "DESeq2")

  # create path for storage of DESeq2 ranking
  path_DESeq2_phenperm <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Raw/Phenotype_Permutation", i, "/DESeq2_ranking_permutation", i, ".txt")

  # export
  write.table(DESeq2_ranking_phenperm,
              file = path_DESeq2_phenperm,
              quote = FALSE,
              row.names = TRUE,
              col.names = FALSE)


  ##############
  #limma ranking
  ##############

  # filtering indicator
  keep_phenperm <- DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, i]) %>% filterByExpr()

  # design matrix
  mm_phenperm <- model.matrix( ~ phen_pickrell[, i])

  # create limma results and rank by p-value
  limma_ranking_phenperm <- geneID_conversion_SYMBOL_pickrell(Biobase::exprs(pickrell.eset)[keep_phenperm, ], dupl_removal_method = 1) %>%
    DGEList(group = phen_pickrell[, i]) %>% calcNormFactors() %>%
    voom(design = mm_phenperm) %>% lmFit(design = mm_phenperm) %>% eBayes() %>% topTable(coef = ncol(mm_phenperm), number = 100000) %>%
    rankedList_cP(rankby =  "p_value", method = "limma")

  # Create path for storage of limma ranking
  path_limma_phenperm <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task2/Raw/Phenotype_Permutation", i, "/limma_ranking_permutation", i, ".txt")


  # export
  write.table(limma_ranking_phenperm,
              file = path_limma_phenperm,
              quote = FALSE,
              row.names = TRUE,
              col.names = FALSE)



}


