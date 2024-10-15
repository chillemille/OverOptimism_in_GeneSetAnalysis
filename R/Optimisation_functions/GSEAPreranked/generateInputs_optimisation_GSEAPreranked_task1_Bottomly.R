#######################################################################################
### Maximize number of differentially enriched gene sets obtained with GSEAPreranked ##
#######################################################################################


library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(apeglm)

##### run script to convert mouse ensembl IDs to human gene Symbols
source("./R/Help_functions/Mapping_MouseID_to_HumanID.R")

### run script to obtain gene expression data sets, true and randomly permuted phenotype assignments
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")

# load preprocessing functions
source("./R/Help_functions/PreProcessing_Functions.R")


######################################
### generate required folders ########
######################################

dir.create("./Results/Intermediate_results/GSEAPreranked")
dir.create("./Results/Intermediate_results/GSEAPreranked/Bottomly")
dir.create("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1")
dir.create("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Raw")
dir.create("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Prep")
dir.create("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Prep/Original_Phenotype")
dir.create("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Raw/Original_Phenotype")

for(i in 1:10){

  path_raw <- paste0("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Raw/Phenotype_Permutation",
                     i)

  path_prep <- paste0("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Prep/Phenotype_Permutation",
                      i)

  dir.create(path_raw)
  dir.create(path_prep)


}





##################################################################################
##create ranked list from DE results##############################################
##################################################################################

#for DESeq2 (method = "DESeq2") and edgeR (method = "edgeR")
#rankby must be in c("p_value", "lfc") to perform ranking based on
#(i) p-value (rank = sign(lfc)*(-1)*log10(unadjusted_pvalue)
#(ii) log fold changes

rankedList_cP <- function(DE_results, rankby, method){

  if(method  ==  "DESeq2"){#create ranking based on DESeq2 results table

    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))

    DE_results$pvalue[ DE_results$pvalue  ==  0] <- min(DE_results$pvalue[ DE_results$pvalue > 0 ]) / 10

    #DE_results <- edgeR_results
    if(rankby  ==  "lfc"){ #ranking by log2 fold change
      #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
      rankvec <- as.vector(DE_results[!is.na(DE_results$pvalue), ]$log2FoldChange)
      names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue), ])
      rankvec <- sort(rankvec, decreasing  =  TRUE)

    }else if (rankby  ==  "p_value"){#ranking by p-value
      #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
      rankvec <- as.vector(sign(DE_results[!is.na(DE_results$pvalue), ]$log2FoldChange)*(-1)*log10(DE_results[!is.na(DE_results$pvalue), ]$pvalue))
      names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue), ])
      rankvec <- sort(rankvec, decreasing  =  TRUE)
    }
  }


  else if(method  ==  "edgeR"){#create ranking based on edgeR results table

    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))

    DE_results$table$PValue[ DE_results$table$PValue == 0 ] <- min(DE_results$table$PValue[ DE_results$table$PValue > 0 ]) /10

    if(rankby  ==  "lfc"){#ranking based on log2 fold change
      rankvec <- as.vector(DE_results$table$logFC)
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing  =  TRUE)
    }

    else if(rankby  ==  "p_value"){#ranking based on p-value
      rankvec <- as.vector(sign(DE_results$table$logFC)*(-1)*log10(DE_results$table$PValue))
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing  =  TRUE)
    }
  }

  else if(method  ==  "limma"){#create ranking based on edgeR results table

    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))

    DE_results$P.Value[ DE_results$P.Value  ==  0 ] <- min(DE_results$P.Value[ DE_results$P.Value > 0 ]) / 10

    if(rankby  ==  "lfc"){#ranking based on log2 fold change
      rankvec <- as.vector(DE_results$logFC)
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing  =  TRUE)
    }

    else if(rankby  ==  "p_value"){#ranking based on p-value
      rankvec <- as.vector(sign(DE_results$logFC)*(-1)*log10(DE_results$P.Value))
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing  =  TRUE)
    }
  }

  # return gene ranking (vector of all genes from differential expression experiment
  # ranked according to the ranking metric)
  return(rankvec)
}


################################################################################
### Export gene rankings #######################################################
################################################################################


# (I) Original Phenotype Assignment

################
# DESeq2 ranking
################

# create DESeq2 results and rank by p-value
DESeq2_ranking_phenorig <- pre_filt(Biobase::exprs(bottomly.eset), threshold  =  10) %>%
  conversion_mouseEnsembl_HumanSymbol(dupl_removal_method  =  1) %>%
  deseq_preprocess(phenotype_labels = bottomly.eset$strain ) %>%
  DESeq() %>%
  lfcShrink(coef = "condition_treated_vs_untreated", type = "apeglm") %>%
  as.data.frame() %>%
  rankedList_cP(rankby = "p_value", method  = "DESeq2")

# create path for storage of DESeq2 ranking
path_DESeq2_phenorig <- "./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Raw/Original_Phenotype/DESeq2_ranking_phenOrig.txt"

# export
write.table(DESeq2_ranking_phenorig,
            file  =  path_DESeq2_phenorig,
            quote  =  FALSE,
            row.names  =  TRUE,
            col.names  =  FALSE)


##############
#limma ranking
##############

# filtering indicator
keep_phenorig <- DGEList(Biobase::exprs(bottomly.eset), group  =  bottomly.eset$strain) %>%
  filterByExpr()

# design matrix
mm_phenorig <- model.matrix( ~ bottomly.eset$strain)

# create limma results and rank by p-value
limma_ranking_phenorig <- conversion_mouseEnsembl_HumanSymbol(Biobase::exprs(bottomly.eset)[keep_phenorig, ], dupl_removal_method  =  1) %>%
  DGEList(group  =  bottomly.eset$strain) %>%
  calcNormFactors() %>%
  voom(design = mm_phenorig) %>% lmFit(design = mm_phenorig) %>%
  eBayes() %>% topTable(coef = ncol(mm_phenorig), number = 100000) %>%
  rankedList_cP(rankby =  "p_value", method = "limma")

# Create path for storage of limma ranking
path_limma_phenorig <- "./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Raw/Original_Phenotype/limma_ranking_phenOrig.txt"


# export
write.table(limma_ranking_phenorig,
            file  =  path_limma_phenorig,
            quote  =  FALSE,
            row.names  =  TRUE,
            col.names  =  FALSE)



# (II) Random Phenotype Permutations of Original Phenotypes

for(i in 1:ncol(phen_bottomly)){

  ################
  # DESeq2 ranking
  ################

  # create DESeq2 results and rank by p-value
  DESeq2_ranking_phenperm <- pre_filt(Biobase::exprs(bottomly.eset), threshold  =  10)  %>%
    conversion_mouseEnsembl_HumanSymbol(dupl_removal_method  =  1) %>%
    deseq_preprocess(phenotype_labels  =  phen_bottomly[, i] ) %>% DESeq() %>%
    lfcShrink(coef = "condition_treated_vs_untreated", type = "apeglm") %>%
    as.data.frame() %>%
    rankedList_cP(rankby  =  "p_value", method  =  "DESeq2")

  # create path for storage of DESeq2 ranking
  path_DESeq2_phenperm <- paste0("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Raw/Phenotype_Permutation",
                                 i,
                                 "/DESeq2_ranking_permutation", i, ".txt")

  # export
  write.table(DESeq2_ranking_phenperm,
              file  =  path_DESeq2_phenperm,
              quote  =  FALSE,
              row.names  =  TRUE,
              col.names  =  FALSE)


  ##############
  #limma ranking
  ##############

  # filtering indicator
  keep_phenperm <- DGEList(Biobase::exprs(bottomly.eset), group  =  phen_bottomly[, i]) %>%
    filterByExpr()

  # design matrix
  mm_phenperm <- model.matrix( ~ phen_bottomly[, i])

  # create limma results and rank by p-value
  limma_ranking_phenperm <-  conversion_mouseEnsembl_HumanSymbol(Biobase::exprs(bottomly.eset)[keep_phenperm, ], dupl_removal_method  =  1) %>%
    DGEList(group  =  phen_bottomly[, i]) %>% calcNormFactors() %>%
    voom(design = mm_phenperm) %>% lmFit(design = mm_phenperm) %>%
    eBayes() %>% topTable(coef = ncol(mm_phenperm), number = 100000) %>%
    rankedList_cP(rankby =  "p_value", method = "limma")

  # Create path for storage of limma ranking
  path_limma_phenperm <- paste0("./Results/Intermediate_results/GSEAPreranked/Bottomly/Data_task1/Raw/Phenotype_Permutation",
                                i,
                                "/limma_ranking_permutation",
                                i,
                                ".txt")


  # export
  write.table(limma_ranking_phenperm,
              file  =  path_limma_phenperm,
              quote  =  FALSE,
              row.names  =  TRUE,
              col.names  =  FALSE)


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

# -> 0 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 0 DEGS
# ->> return to default gene list generated with DESeq2

#########
# 3. step: change geneset database to KEGG
#########

# -> 0 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 ->  16 + 24  =  40 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 ->

# ->> return to default exponent 1

### -> final results: 40 DEGS
# achieved WITH alternative exponent 0



################################################################################
### Random Phenotype Permutations ##############################################
################################################################################



################################################################################
### Phenotype Permutation 1  ###################################################
################################################################################

#########
# 1. step: Default
#########

# -> 65 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 225 DEGS
# ->> proceed with alternative gene ranking generated with limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 0 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 148 + 60  =  208 DEGS
# alternative 2: exponent 1.5 -> 20 + 4  =  24 DEGS
# alternative 3: exponent 2 -> 43 + 1  =  44 DEGS

# final results: 225 DEGS
# achieved with ALTERNATIVE ranking metric limma




################################################################################
### Phenotype Permutation 2  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 0 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 28 DEGS
# ->> proceed with alternative gene ranking generated with limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 8 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 30 DEGS
# alternative 2: exponent 1.5 -> 10 DEGS
# alternative 3: exponent 2 -> 3 DEGS

# final results: 30 DEGS
# achieved with
# ALTERNATIVE ranking generated using limma
# ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 3  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 17 + 1  =  18 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 24 + 9  =  33 DEGS
# ->> proceed with alternative gene ranking generated with limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 1 DEGS
# return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 77 + 91  =  168 DEGS
# alternative 2: exponent 1.5 -> 11 + 36  =  47 DEGS
# alternative 3: exponent 2 -> 5 DEGS

### -> final results: 168 DEGS
# achieved with
# ALTERNATIVE gene ranking using limma
# ALTERNATIVE exponent 0



################################################################################
### Phenotype Permutation 4  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 0 DEGS



#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 8 DEGS
# ->> proceed with alternative ranking generated using limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 0 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 3 + 32  =  35 DEGS
# alternative 2: exponent 1.5 -> 1 DEGS
# alternative 3: exponent 2 -> 0 DEGS

################################################################################
### Phenotype Permutation 5  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 239 + 41  =  280 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 52 + 56  =  108 DEGS
# return to default ranking generated with DESeq2


#########
# 3. step: change geneset database to KEGG
#########

# -> 6 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 239 + 142  =  381 DEGS
# alternative 2: exponent 1.5 -> 14 + 3  =  17 DEGS
# alternative 3: exponent 2 -> 104+ 4  =  108 DEGS

### -> final results: 381 DEGS
# achieved with ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 6  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 124 DEGS

#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 6 DEGS
# ->> return to default gene ranking generated with DESeq2


#########
# 3. step: change geneset database to KEGG
#########

# -> 0 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 40 DEGS
# alternative 2: exponent 1.5 -> 29 DEGS
# alternative 3: exponent 2 -> 0 DEGS

# ->> final results: 124 DEGS
# optimal parameter configuration coincides with default configuration
# ->> number of differentially enriched gene sets could not be increased


################################################################################
### Phenotype Permutation 7  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 114 + 17  =  131 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 146 + 56  =  202 DEGS
# ->> proceed with alternative ranking generated using limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 10 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 382 + 117  =  499 DEGS
# alternative 2: exponent 1.5 -> 24 DEGS
# alternative 3: exponent 2 -> 8 DEGS

# -> final results: 499 DEGS
# achieved with
# ALTERNATIVE gene ranking generated using limma
# ALTERNATIVE exponent 0




################################################################################
### Phenotype Permutation 8  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 7 + 1  =  8 DEGS

#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 273 + 1  =  274 DEGS
# ->> proceed with alternative gene ranking generated using limma


#########
# 3. step: change geneset database to KEGG
#########

# -> 1 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 246 + 70  =  316 DEGS
# alternative 2: exponent 1.5 -> 8 + 5  =  13 DEGS
# alternative 3: exponent 2 -> 119 DEGS





################################################################################
### Phenotype Permutation 9  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> 262 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 238 + 2  =  240 DEGS
# return to default gene ranking generated with DESeq2


#########
# 3. step: change geneset database to KEGG
#########

# -> 1 DEGS
# ->> return to default geneset database GO (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 535 + 19  =  554 DEGS
# alternative 2: exponent 1.5 -> 1 DEGS
# alternative 3: exponent 2 -> 17 DEGS

# ->> final results: 554 DEGS
# achieved with ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 10  ##################################################
################################################################################


#########
# 1. step: Default
#########

# -> 63 + 661  =  724 DEGS


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> 52 + 651  =  703 DEGS
# ->> return to default gene ranking generated with DESeq2


#########
# 3. step: change geneset database to KEGG
#########

# -> 4 + 23  =  27 DEGS
# ->> return to default geneset database GP (BP)


#########
# 4. step: change exponent
#########

# alternative 1: exponent 0 -> 198 + 764  =  962 DEGS
# alternative 2: exponent 1.5 -> 757 + 7  =  764 DEGS
# alternative 3: exponent 2 -> 12 + 474  =  486 DEGS

# final results: 962 DEGS
# -> achieved with ALTERNATIVE exponent 0



