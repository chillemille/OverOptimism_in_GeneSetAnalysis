#######################################################################################
### Maximize number of differentially enriched gene sets obtained with GSEAPreranked ##
#######################################################################################


library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(apeglm)


### run script to obtain gene expression data sets, true and randomly permuted phenotype assignments
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")

# load preprocessing functions
source("./R/Help_functions/PreProcessing_Functions.R")

# load help functions for GSEAPreranked
source("./R/Help_functions/helpfunctions_GSEAPreranked.R")


######################################
### generate required folders ########
######################################

dir.create("./Results/Intermediate_results/GSEAPreranked")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Raw")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Prep")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Prep/Original_Phenotype")
dir.create("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Raw/Original_Phenotype")

for(i in 1:10){

  path_raw <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Raw/Phenotype_Permutation",
                     i)

  path_prep <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Prep/Phenotype_Permutation",
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
DESeq2_ranking_phenorig <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>%
  geneID_conversion_SYMBOL_pickrell(dupl_removal_method = 1) %>%
  deseq_preprocess(phenotype_labels = pickrell.eset$gender ) %>% DESeq() %>%
  lfcShrink(coef="condition_treated_vs_untreated", type="apeglm") %>% as.data.frame() %>%
  rankedList_cP(rankby = "p_value", method = "DESeq2")

# create path for storage of DESeq2 ranking
path_DESeq2_phenorig <- "./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Raw/Original_Phenotype/DESeq2_ranking_phenOrig.txt"

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
keep_phenorig <- DGEList(Biobase::exprs(pickrell.eset), group = pickrell.eset$gender) %>%
  filterByExpr()

# design matrix
mm_phenorig <- model.matrix( ~ pickrell.eset$gender)

# create limma results and rank by p-value
limma_ranking_phenorig <- geneID_conversion_SYMBOL_pickrell(Biobase::exprs(pickrell.eset)[keep_phenorig, ], dupl_removal_method = 1) %>%
  DGEList(group = pickrell.eset$gender) %>% calcNormFactors() %>%
  voom(design=mm_phenorig) %>% lmFit(design=mm_phenorig) %>%
  eBayes() %>% topTable(coef=ncol(mm_phenorig), number=100000) %>%
  rankedList_cP(rankby= "p_value", method="limma")

# Create path for storage of limma ranking
path_limma_phenorig <- "./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Raw/Original_Phenotype/limma_ranking_phenOrig.txt"


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
  DESeq2_ranking_phenperm <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>%
    geneID_conversion_SYMBOL_pickrell(dupl_removal_method = 1) %>%
    deseq_preprocess(phenotype_labels = phen_pickrell[, i] ) %>%
    DESeq() %>%
    lfcShrink(coef="condition_treated_vs_untreated", type="apeglm") %>%
    as.data.frame() %>%
    rankedList_cP(rankby = "p_value", method = "DESeq2")

  # create path for storage of DESeq2 ranking
  path_DESeq2_phenperm <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Raw/Phenotype_Permutation",
                                 i,
                                 "/DESeq2_ranking_permutation",
                                 i,
                                 ".txt")

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
  keep_phenperm <- DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, i]) %>%
    filterByExpr()

  # design matrix
  mm_phenperm <- model.matrix( ~ phen_pickrell[, i])

  # create limma results and rank by p-value
  limma_ranking_phenperm <- geneID_conversion_SYMBOL_pickrell(Biobase::exprs(pickrell.eset)[keep_phenperm, ], dupl_removal_method = 1) %>%
    DGEList(group = phen_pickrell[, i]) %>% calcNormFactors() %>%
    voom(design=mm_phenperm) %>% lmFit(design=mm_phenperm) %>%
    eBayes() %>% topTable(coef=ncol(mm_phenperm), number=100000) %>%
    rankedList_cP(rankby= "p_value", method="limma")

  # Create path for storage of limma ranking
  path_limma_phenperm <- paste0("./Results/Intermediate_results/GSEAPreranked/Pickrell/Data_task1/Raw/Phenotype_Permutation",
                                i,
                                "/limma_ranking_permutation",
                                i,
                                ".txt")


  # export
  write.table(limma_ranking_phenperm,
              file = path_limma_phenperm,
              quote = FALSE,
              row.names = TRUE,
              col.names = FALSE)

}


################################################################################
### (I) Full Optimization Process (Pre-Processing and Internal Parameters) #####
################################################################################


################################################################################
### Original Phenotypes ########################################################
################################################################################



#########
# 1. step: Default
#########

# -> 194 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 541 DEGS
# ->> proceed with limma as DE method to generate ranking of genes


#########
# 3. step: change geneset database to KEGG
#########

# -> 6 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 59 + 316 = 375 DEGS
# alternative 2: exponent 1.5 -> 365 DEGS
# alternative 3: exponent 2 -> 0 DEGS

# ->> return to default exponent 1

### -> final results: 541 DEGS
# achieved with ALTERNATIVE DE method limma



################################################################################
### Random Phenotype Permutations ##############################################
################################################################################



################################################################################
### Phenotype Permutation 1  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 18 + 30 = 48 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 220 + 34 = 254 DEGS
# -> proceed with limma to generate ranked list of genes


#########
# 3. step: change geneset database to KEGG
#########

# -> 11 DEGS
# return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 297 + 54 = 351 DEGS
# alternative 2: exponent 1.5 -> 53 + 25 = 78 DEGS
# alternative 3: exponent 2 -> 117 + 10 = 127 DEGS


### -> final results: 351 DEGS
# achieved with
# ALTERNATIVE DE method limma
# ALTERNATIVE exponent 0



################################################################################
### Phenotype Permutation 2  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 79 + 172 = 251 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 255 + 173 = 428 DEGS
# ->> proceed with ranked gene list generated with limma

#########
# 3. step: change geneset database to KEGG
#########

# -> 14 DEGS
# -> return to default geneset database GO (BP)

#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 448 + 191 = 639 DEGS
# alternative 2: exponent 1.5 -> 16 + 99 = 115 DEGS
# alternative 3: exponent 2 -> 118 + 229 DEGS = 347


# ->> final results: 639 DEGS
# achieved with
# ALTERNATIVE DE method limma
# ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 3  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 1 + 37 = 38 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 127 DEGS
# proceed with ranked gene list generated with limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 12 DEGS
# return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 2 + 131 = 133 DEGS
# alternative 2: exponent 1.5 -> 120 DEGS
# alternative 3: exponent 2 -> 67 DEGS

### final results: 133 DEGS
# achieved with
# ALTERNATIVE DE method limma
# ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 4  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 31 + 17 = 48 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 74 + 272 = 346 DEGS
# ->> proceed with ranked gene list limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 14 + 10 = 24 DEGS
# return to default geneset database GO(BP)

#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 101 + 543 = 644 DEGS
# alternative 2: exponent 1.5 -> 30 + 680 = 710 DEGS
# alternative 3: exponent 2 -> 30 + 58 = 88 DEGS

# final results: 710 DEGS
# achieved with
# ALTERNATIVE DE method limma
# ALTERNATIVE exponent 1.5


################################################################################
### Phenotype Permutation 5  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 12 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 28 DEGS
# ->> proceed with ranked gene list created with limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 6 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 21 + 175 DEGS = 196 DEGS
# alternative 2: exponent 1.5 -> 143 DEGS
# alternative 3: exponent 2 -> 3 DEGS


### final results: 196 DEGS
# achieved with
# ALTERNATIVE DE method limma
# ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 6  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 1 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 1 DEGS
# return to default gene ranking generated with DESeq2


#########
# 3. step: change geneset database to KEGG
#########

# -> 3 DEGS
# -> proceed with geneset database KEGG


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 5 DEGS
# alternative 2: exponent 1.5 -> 10 DEGS
# alternative 3: exponent 2 -> 0 DEGS


################################################################################
### Phenotype Permutation 7  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 33 + 77 = 110 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 70 + 101 = 171 DEGS
# ->> proceed with gene ranking generated with limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 9 DEGS
# ->> return to default geneset database GO (BP)



#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 233 + 119 = 352 DEGS
# alternative 2: exponent 1.5 -> 11 + 33 = 44DEGS
# alternative 3: exponent 2 -> 36 + 38 = 74 DEGS


# final results: 352 DEGS
# achieved with
# ALTERNATIVE DE method limma
# ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 8  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 1


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 1 + 2 = 3 DEGS
## ->> proceed with gene ranking generated with limma



#########
# 3. step: change geneset database to KEGG
#########

# -> 9 DEGS
## -> proceed with geneset database KEGG


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 4 DEGS
# alternative 2: exponent 1.5 -> 3 DEGS
# alternative 3: exponent 2 -> 1 DEGS


### final results: 9 DEGS
# achieved with
# ALTERNATIVE ranking genrated with limma
# ALTERNATIVE geneset database KEGG


################################################################################
### Phenotype Permutation 9  ###################################################
################################################################################


#########
# 1. step: Default
#########

# ->  27 + 2 = 29 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 7 DEGS
# -> return to default DE method DESeq2 to generate ranked list of genes



#########
# 3. step: change geneset database to KEGG
#########

# ->  4 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 140 + 33 = 173 DEGS
# alternative 2: exponent 1.5 -> 8 DEGS
# alternative 3: exponent 2 -> 1 DEGS

### ->> final results: 173 DEGS
# achieved with
# ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 10  ##################################################
################################################################################


#########
# 1. step: Default
#########

# ->  19 + 6 = 25 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> 71 + 1 = 72 DEGS
# ->> proceed with ranked gene list generated with limma

#########
# 3. step: change geneset database to KEGG
#########

# -> 9 + 19 = 28 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 78 + 47 = 125 DEGS
# alternative 2: exponent 1.5 -> 22 DEGS
# alternative 3: exponent 2 -> 29 DEGS


### -> final results: 125 DEGS
# achieved with
# ALTERNATIVE DE method limma
# ALTERNATIVE exponent 0


