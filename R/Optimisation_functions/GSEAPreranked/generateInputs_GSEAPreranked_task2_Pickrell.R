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




#######################################################################################
###pre-filtering function #############################################################
#######################################################################################

pre_filt <- function(expression_data, threshold){

  expression_data_filt <- expression_data[rowSums(expression_data) >= threshold, ]

  # return filtered gene expression data set
  return(expression_data_filt)


}

##################################################################################
##Pre-Process expression data for DESeq2##########################################
##################################################################################


deseq_preprocess <- function(expression_data, phenotype_labels){

  # check whether the number of samples in expression_data equals the length of the phenotype vector
  if(ncol(expression_data)!= length(phenotype_labels)) stop("Error: number of samples in expression expression_data and p
                                                phenotype labels do not match")

  # generate expression_data frame containing sample labels of respective samples and column name "condition"
  coldata <- data.frame(phenotype_labels,
                      row.names = colnames(expression_data))
  colnames(coldata) <- "condition"
  coldata$condition <- factor(coldata$condition,
                              labels = c("untreated", "treated"))

  # generate DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = expression_data, #submit expression data set of counts
    colData = coldata, #as generated above
    design = ~ condition)

  return(dds)

}


#######################################################################################
### Gene ID conversion and duplicate gene ID removal ##################################
#######################################################################################

# particularity of GSEAPreranked: Recommended/ required gene ID format is Hugo symbols
# and not Entrez gene IDs

# background: GSEA typically run using genesets from MSigDB which consists of human gene symbols
# If input data contain other identifiers, the IDs need to be converted to gene symbols
# option "Collapse/ Remap to gene symbols" performs conversion which handles the case of
# several feature identifiers mapping to same gene identifier.
# This method was developed and tuned for gene expression data, however, the ranked list
# of genes in GSEAPreranked was created using unspecified ranking procedure outside of GSEA

# -> recommended to provide ranked list with genes already converted to gene SYMBOLS and
# select parameter "NO_Collapse"

# function arguments
# - expression_data: gene expression data set
# - dupl_removal_method: manner of removing duplicated gene IDs that result from
#                         gene ID conversion

geneID_conversion_SYMBOL <- function(expression_data, dupl_removal_method){


  # identify organism for which gene expression is measured from the format of the gene IDs (all are Ensembl)
  # ENSEMBL gene ID ENSMUSGXXXXXXXXXXX corresponds to mouse; identifiable by substring "ENSMUSG"
  # ENSEMBL gene ID ENSGXXXXXXXXXXX corresponds to homo sapiens; identifiable by string "ENSG"

  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl,
                         X = c("ENSG", "ENSMUSG"),
                         x = rownames(expression_data)[1])
  # choose suitable organism (required for function bitr)
  organism <- unlist(c(org.Hs.eg.db , org.Mm.eg.db)[ind_organism])[[1]]




  # Gene ID conversion via clusterProfiler::bitr()
  bitr_enstoentr <- bitr(rownames(expression_data) ,
                         fromType = "ENSEMBL",
                         toType = "SYMBOL",
                         OrgDb = organism)
  # note: not all ENSEMBL IDs could be converted to a corresponding ENTREZ Gene ID
  dim(bitr_enstoentr)

  # results:
  # - not all Ensembl gene IDs can be mapped to a corresponding Entrez gene ID
  # - some single Ensembl gene IDs were mapped to multiple distinct Entrez gene IDs
  # - some distincts Ensembl gene IDs were mapped to an identical Entrez gene ID

  ### merge

  # this step is independent of sample IDs
  # merge by row names of expression data set and ENSEMBL ID of conversion data set
  expression_data <- merge(expression_data,
                           bitr_enstoentr,
                           by.x = 0,
                           by.y = "ENSEMBL",
                           all.y = TRUE,
                           sort = TRUE)
  dim(expression_data)

  ### take closer look at duplicates

  # CASE 1: single ENSEMBL IDs are mapped to multiple ENTREZ IDs
  # View(bitr_toKEGG[(duplicated(bitr_toKEGG$ENSEMBL)), ])
  sum(duplicated(bitr_enstoentr$ENSEMBL)) #number of times an ENSEMBL gene ID was converted to several ENTREZ IDs
  # determine all duplicated ENSEMBL gene IDS
  dupl_ensembl <- unique(bitr_enstoentr$ENSEMBL[duplicated(bitr_enstoentr$ENSEMBL)])
  # number of ENSEMBL IDs that have at least one duplicate
  length(dupl_ensembl)
  # display of conversion scheme of duplicated ENSEMBL IDs
  duplicated_conversion_ens <- bitr_enstoentr[bitr_enstoentr$ENSEMBL %in% dupl_ensembl, ]
  dim(duplicated_conversion_ens)


  # CASE 2: multiple ENSEMBL IDs are mapped to single entrez ID
  sum(duplicated(bitr_enstoentr$ENTREZ)) #number of times several ENSEMBL gene IDs were converted to a single ENTREZ ID
  # determine all duplicated ENTREZ IDs
  dupl_entrez <- unique(bitr_enstoentr$SYMBOL[duplicated(bitr_enstoentr$SYMBOL)])
  # number of ENTREZ IDs that have at least one duplicate
  length(dupl_entrez)
  # display of conversion scheme of duplicated ENTREZ IDs
  duplicated_conversion_entrez <- bitr_enstoentr[bitr_enstoentr$SYMBOL %in% dupl_entrez, ]
  dim(duplicated_conversion_entrez)




  #######
  #2. step: Removal of Duplicated gene IDs
  #######


  if(dupl_removal_method ==  1){

    ### 1. option: keep first subscript among duplicates #########################

    #1. remove duplicated ENTREZ gene IDs
    exprdat_dupl <- expression_data[!duplicated(expression_data$SYMBOL), ]
    dim(expression_data)

    #2. remove duplicated ENSEMBL gene IDs
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names), ]
    dim(exprdat_dupl)

    #3. ENTREZ IDs as row names and
    rownames(exprdat_dupl) <- exprdat_dupl$SYMBOL
    #Remove columns containing ENSEMBL and ENTREZ IDs
    exprdat_dupl <- subset(exprdat_dupl, select = -c(Row.names, SYMBOL))
    dim(exprdat_dupl)

  } else if(dupl_removal_method == 2){

    ### 2. option: keep mean expression value of all duplicated gene IDs  #########################



    #generate matrix to contain (rounded) mean expression values of all rows that
    #have same ENTREZ gene ID
    #ncol = ncol(expression_data)-2 since data set contains 2 columns with IDs at this point
    mean_entrez <- matrix(, nrow = 0, ncol = ncol(expression_data)-2)


    # -> There are cases where no ENSEMBL gene IDs are mapped to an identical Entrez gene ID
    # -> nrow(dupl_entrez) = 0
    # -> add if-condition since the following part of code would otherwise produce NaN values which would in turn lead to errors
    # when pre-processing the gene expression data set

    # -> if no distinct ENSEMBL IDs are mapped to an identical Entrez gene ID then mean_entrez remains an empty data frame


    if(length(dupl_entrez) != 0){


      #1.remove duplicated ENTREZ gene IDs (case 2)
      #i.e. multiple different ENSEMBL IDs that are mapped to the same single ENTREZ ID

      for(i in 1:length(dupl_entrez)){#go through each ENTREZ IDs which occurs multiple times
        #determine all rows whose ENTREZ IDs correspond to current ENTREZ ID
        counts_dupl <- expression_data[expression_data$SYMBOL %in% unique(dupl_entrez)[i], ]
        #for rows duplicated ENTREZ ID compute (rounded) mean expression value
        dupl_id <- round(colMeans(counts_dupl[, c(2:(ncol(expression_data)-1))]))
        #store rounded mean expression value in matrix
        mean_entrez <- rbind(mean_entrez, dupl_id)
      }
    }


    #test whether the number of rows in mean_entrez corresponds to the number ENTREZ IDs
    #that occur more than once
    nrow(mean_entrez) == length(dupl_entrez)

    #remove all rows from the expression data whose ENTREZ ID has at least one duplicate
    exprdat_dupl <- expression_data[!expression_data$SYMBOL %in% dupl_entrez, ]
    #test whether number of rows in resulting data set equals nrow of inital data set
    #minus number of genes with at least one duplicate
    nrow(exprdat_dupl) == nrow(expression_data)-nrow(duplicated_conversion_entrez)
    dim(exprdat_dupl)

    #set corresponding ENTREZ gene IDs as rownames
    rownames(mean_entrez) <- unique(dupl_entrez)




    #2. remove duplicated ENSEMBL IDs
    #caution: single ENSEMBL IDs that are mapped to multiple ENTREZ ID naturally generate
    #identical count data for all corresponding ENTREZ IDs
    #->pointless to compute mean expression values
    #verifiable by looking at data set only containing those ENSEMBL IDs that are
    #mapped by multiple ENTREZ IDs:
    #test_dupl_ensembl <- expression_data[expression_data$Row.names %in% dupl_ensembl, ]
    #View(test_dupl_ensembl)

    #therefore: proceed as in option 1 and use ENTREZ ID that occurs first, remove the rest
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names), ]
    dim(exprdat_dupl)
    #set ENTREZ ID as rownames
    rownames(exprdat_dupl) <- exprdat_dupl$SYMBOL
    #remove any columns containing IDs
    exprdat_dupl <- subset(exprdat_dupl, select =  -c(Row.names, SYMBOL))
    #add rows to data set that contain mean expression values of duplicate ENTREZ IDs
    exprdat_dupl <- rbind(exprdat_dupl, mean_entrez)
    #dimension of remaining expression data set:
    #dim(exprdat_dupl)

  }else if(dupl_removal_method == 3){

    ###option 3: among duplicates, keep row with highest overall expression values (i.e highest counts across all samples)

    #intuition: row with highest counts values has  highest power of detecting
    #differential expression later on
    #as in option 2, this applies only to duplicates that result from multiple ENSEMBL IDs
    #that are mapped to the same ENTREZ ID


    #case 2: (case 1 below) multiple ENSEMBL IDs that are converted to the same single ENTREZ ID

    #generate matrix to later contain row with highest count values among ID duplicates
    highest_count_entrez <- matrix(, nrow = 0, ncol = ncol(expression_data))
    #go through each ENTREZ ID that occurs multiple times
    for(i in 1:length(dupl_entrez)){
      #determine all rows with specific ENTREZ ID which occurs multiple times
      counts_dupl <- expression_data[expression_data$SYMBOL %in% unique(dupl_entrez)[i], ]
      #detect row with highest count values and order in decreasing manner
      order_rowsums <- order(rowSums(counts_dupl[, 2:(ncol(counts_dupl)-1)]), decreasing = TRUE)
      dupl_id <- counts_dupl[order_rowsums == 1, ]
      #store rounded mean expression value in matrix
      highest_count_entrez <- rbind(highest_count_entrez, dupl_id)
      #View(highest_count_entrez)
      #remove rows in counts_dupl from count data set successively
    }

    #Remove all initial values with ENTREZ duplicates from the dataset
    exprdat_dupl <- expression_data[! expression_data$SYMBOL %in% unique(dupl_entrez), ]


    #case 1: single ENSEMBL ID that is mapped to multiple ENTREZ gene IDs
    #as in option 2, pointless to detect row with highest count values as all rows
    #corresponding to the same ENSEMBL ID naturally contain identical count data
    #therefore: remove duplicate ENSEMBL ID that occurs first
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names), ]

    #Add all rows contain initially duplicate ENTREZ IDs but contain highest
    #count values among those
    exprdat_dupl <- rbind(exprdat_dupl, highest_count_entrez )

    #Set ENTREZ IDs as rownames remove all columns containing any ID info and
    rownames(exprdat_dupl) <- exprdat_dupl$SYMBOL
    #Remove any column that contains gene IDs
    exprdat_dupl <- subset(exprdat_dupl, select = -c(Row.names, SYMBOL))
    #dim(exprdat_dupl)
  }

  # store resulting gene expression data sets (all with duplicated gene IDs removed
  # in different manners)
  return(exprdat_dupl)



}

##################################################################################
##create ranked list from DE results##############################################
##################################################################################

# function arguments:
# - DE_results: results table of differential expression analysis
# - rankby: must be in c("p_value", "lfc") to perform ranking based on
#           (i) p-value (rank = sign(lfc)*(-1)*log10(unadjusted_pvalue)
#           (ii) log fold changes
# - method: method by which differential expression analysis was performed
#           -> "edgeR" or "DESeq2"

rankedList_cP <- function(DE_results, rankby, method){


  if(method == "DESeq2"){#create ranking based on DESeq2 results table

    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))

    DE_results$pvalue[ DE_results$pvalue ==  0 ] <- min(DE_results$pvalue[ DE_results$pvalue > 0 ]) / 10

    #DE_results <- edgeR_results
    if(rankby == "lfc"){ #ranking by log2 fold change
      #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
      rankvec <- as.vector(DE_results[!is.na(DE_results$pvalue), ]$log2FoldChange)
      names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue), ])
      rankvec <- sort(rankvec, decreasing = TRUE)

    }else if (rankby == "p_value"){#ranking by p-value
      #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
      rankvec <- as.vector(sign(DE_results[!is.na(DE_results$pvalue), ]$log2FoldChange)*(-1)*log10(DE_results[!is.na(DE_results$pvalue), ]$pvalue))
      names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue), ])
      rankvec <- sort(rankvec, decreasing = TRUE)
    }
  }


  else if(method == "edgeR"){#create ranking based on edgeR results table

    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))

    DE_results$table$PValue[ DE_results$table$PValue ==  0 ] <- min(DE_results$table$PValue[ DE_results$table$PValue > 0 ]) /10

    if(rankby == "lfc"){#ranking based on log2 fold change
      rankvec <- as.vector(DE_results$table$logFC)
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing = TRUE)
    }

    else if(rankby == "p_value"){#ranking based on p-value
      rankvec <- as.vector(sign(DE_results$table$logFC)*(-1)*log10(DE_results$table$PValue))
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing = TRUE)
    }
  }

  else if(method == "limma"){#create ranking based on edgeR results table

    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))

    DE_results$P.Value[ DE_results$P.Value ==  0 ] <- min(DE_results$P.Value[ DE_results$P.Value > 0 ]) / 10

    if(rankby == "lfc"){#ranking based on log2 fold change
      rankvec <- as.vector(DE_results$logFC)
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing = TRUE)
    }

    else if(rankby == "p_value"){#ranking based on p-value
      rankvec <- as.vector(sign(DE_results$logFC)*(-1)*log10(DE_results$P.Value))
      names(rankvec) <- rownames(DE_results)
      rankvec <- sort(rankvec, decreasing = TRUE)
    }
  }
  # return gene ranking (ranking based on results table of differential expression
  # analysis and ranking metric)
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
DESeq2_ranking_phenorig <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% geneID_conversion_SYMBOL(dupl_removal_method = 1) %>%
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
limma_ranking_phenorig <- geneID_conversion_SYMBOL(Biobase::exprs(pickrell.eset)[keep_phenorig, ], dupl_removal_method = 1) %>%
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
  DESeq2_ranking_phenperm <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% geneID_conversion_SYMBOL(dupl_removal_method = 1) %>%
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
  limma_ranking_phenperm <- geneID_conversion_SYMBOL(Biobase::exprs(pickrell.eset)[keep_phenperm, ], dupl_removal_method = 1) %>%
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


