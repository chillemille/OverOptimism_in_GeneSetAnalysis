################################################################################
### help functions for GSEAPreranked ###########################################
################################################################################

#######################################################################################
### Gene ID conversion and duplicate gene ID removal ##################################
#######################################################################################

# particularity of GSEAPreranked: Recommended/ required gene ID format is Hugo symbols
# and not human HGNC symbols

# background: GSEA typically run using genesets from MSigDB which consists of human gene symbols
# If input data contain other identifiers, the IDs need to be converted to gene symbols
# option "Collapse/ Remap to gene symbols" performs conversion which handles the case of
# several feature identifiers mapping to same gene identifier.
# This method was developed and tuned for gene expression data, however, the ranked list
# of genes in GSEAPreranked was created using unspecified ranking procedure outside of GSEA

# -> recommended to provide ranked list with genes already converted to gene SYMBOLS and
# select parameter "NO_Collapse"

# function arguments:
# - expression_data: gene expression data set
# - dupl_removal_method: manner of removing duplicated gene IDs resulting
#                        from gene ID conversion

geneID_conversion_SYMBOL_bottomly <- function(expression_data, dupl_removal_method){

  # Conversion of the gene IDs from Mouse mouse ENSEMBL IDs to human HGNC symbols using
  # biomaRt


  # IMPORTANT NOTE: using the default host (i.e. current version of the ensembl website) returned an error
  # we therefore had to work with archived versions: we went down the list in listEnsemblArchives() and
  # selected the most recent archive that works


  # connect to Biomart database for mouse
  mouse <- useMart(biomart = "ensembl",
                   dataset = "mmusculus_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org")
  # connect to Biomart database for human
  human <- useMart(biomart = "ensembl",
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org")

  mouseEnsembl_to_humanSymbol <- getLDS(
    mart = mouse,
    attributes = c('ensembl_gene_id'),
    martL = human,
    attributesL = c('hgnc_symbol'),
    # filters = 'ensembl_gene_id',
    values = rownames(expression_data)
  )


  # results:
  # - not all mouse mouse ENSEMBL IDs can be mapped to a corresponding human HGNC symbol
  # - some single mouse mouse ENSEMBL IDs were mapped to multiple distinct human HGNC symbols
  # - some distincts mouse mouse ENSEMBL IDs were mapped to an identical human HGNC symbol

  ### merge

  # this step is independent of sample IDs
  # merge by row names of expression data set and mouse ENSEMBL ID of conversion data set
  expression_data <- merge(expression_data,
                           mouseEnsembl_to_humanSymbol,
                           by.x = 0,
                           by.y = "Gene.stable.ID",
                           all = FALSE)
  dim(expression_data)

  ###take closer look at duplicates

  # CASE 1: single mouse ENSEMBL IDs are mapped to multiple human HGNC symbols
  # View(bitr_toKEGG[(duplicated(bitr_toKEGG$ENSEMBL)), ])
  sum(duplicated(mouseEnsembl_to_humanSymbol$Gene.stable.ID)) #number of times an mouse ENSEMBL ID was converted to several human HGNC symbols
  # determine all duplicated mouse ENSEMBL IDS
  dupl_ensembl <- unique(mouseEnsembl_to_humanSymbol$Gene.stable.ID[duplicated(mouseEnsembl_to_humanSymbol$Gene.stable.ID)])
  # number of mouse ENSEMBL IDs that have at least one duplicate
  length(dupl_ensembl)
  # display of conversion scheme of duplicated mouse ENSEMBL IDs
  duplicated_conversion_ens <- mouseEnsembl_to_humanSymbol[mouseEnsembl_to_humanSymbol$Gene.stable.ID %in% dupl_ensembl, ]
  dim(duplicated_conversion_ens)


  # CASE 2: multiple mouse ENSEMBL IDs are mapped to single human HGNC symbol
  sum(duplicated(mouseEnsembl_to_humanSymbol$HGNC.symbol)) #number of times several mouse ENSEMBL IDs were converted to a single human HGNC symbol
  # determine all duplicated human HGNC symbols
  dupl_entrez <- unique(mouseEnsembl_to_humanSymbol$HGNC.symbol[duplicated(mouseEnsembl_to_humanSymbol$HGNC.symbol)])
  # number of human HGNC symbols that have at least one duplicate
  length(dupl_entrez)
  # display of conversion scheme of duplicated human HGNC symbolmouseEnsembl_to_humanSymbol$HGNC.symbols
  duplicated_conversion_entrez <- mouseEnsembl_to_humanSymbol[mouseEnsembl_to_humanSymbol$HGNC.symbol %in% dupl_entrez, ]
  dim(duplicated_conversion_entrez)




  #######
  #2. step: Removal of Duplicated gene IDs
  #######


  if(dupl_removal_method ==  1){

    ### 1. option: keep first subscript among duplicates #########################

    #1. remove duplicated human HGNC symbols
    exprdat_dupl <- expression_data[!duplicated(expression_data$HGNC.symbol), ]
    dim(expression_data)

    #2. remove duplicated mouse ENSEMBL IDs
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names), ]
    dim(exprdat_dupl)

    #3. human HGNC symbols as row names and
    rownames(exprdat_dupl) <- exprdat_dupl$HGNC.symbol
    # Remove columns containing ENSEMBL and human HGNC symbols
    exprdat_dupl <- subset(exprdat_dupl, select = -c(Row.names, HGNC.symbol))
    dim(exprdat_dupl)

  } else if(dupl_removal_method ==  2){

    ### 2. option: keep mean expression value of all duplicated gene IDs  #########################



    # generate matrix to contain (rounded) mean expression values of all rows that
    # have same human HGNC symbol
    # ncol = ncol(expression_data)-2 since data set contains 2 columns with IDs at this point
    mean_entrez <- matrix(, nrow = 0, ncol = ncol(expression_data)-2)


    # -> There are cases where no mouse ENSEMBL IDs are mapped to an identical human HGNC symbol
    # -> nrow(dupl_entrez) = 0
    # -> add if-condition since the following part of code would otherwise produce NaN values which would in turn lead to errors
    # when pre-processing the gene expression data set

    # -> if no distinct mouse ENSEMBL IDs are mapped to an identical human HGNC symbol then mean_entrez remains an empty data frame


    if(length(dupl_entrez) != 0){


      #1. remove duplicated human HGNC symbols (case 2)
      #i.e. multiple different mouse ENSEMBL IDs that are mapped to the same single human HGNC symbol

      for(i in 1:length(dupl_entrez)){#go through each human HGNC symbols which occurs multiple times
        # determine all rows whose human HGNC symbols correspond to current human HGNC symbol
        counts_dupl <- expression_data[expression_data$HGNC.symbol %in% unique(dupl_entrez)[i], ]
        # for rows duplicated human HGNC symbol compute (rounded) mean expression value
        dupl_id <- round(colMeans(counts_dupl[, c(2:(ncol(expression_data)-1))]))
        # store rounded mean expression value in matrix
        mean_entrez <- rbind(mean_entrez, dupl_id)
      }
    }


    # test whether the number of rows in mean_entrez corresponds to the number human HGNC symbols
    # that occur more than once
    nrow(mean_entrez) ==  length(dupl_entrez)

    # remove all rows from the expression data whose human HGNC symbol has at least one duplicate
    exprdat_dupl <- expression_data[!expression_data$HGNC.symbol %in% dupl_entrez, ]
    # test whether number of rows in resulting data set equals nrow of inital data set
    # minus number of genes with at least one duplicate
    nrow(exprdat_dupl) ==  nrow(expression_data) - nrow(duplicated_conversion_entrez)
    dim(exprdat_dupl)

    # set corresponding human HGNC symbols as rownames
    rownames(mean_entrez) <- unique(dupl_entrez)




    #2. remove duplicated mouse ENSEMBL IDs
    # caution: single mouse ENSEMBL IDs that are mapped to multiple human HGNC symbol naturally generate
    # identical count data for all corresponding human HGNC symbols
    # ->pointless to compute mean expression values
    # verifiable by looking at data set only containing those mouse ENSEMBL IDs that are
    # mapped by multiple human HGNC symbols:
    # test_dupl_ensembl <- expression_data[expression_data$Row.names %in% dupl_ensembl, ]
    # View(test_dupl_ensembl)

    # therefore: proceed as in option 1 and use human HGNC symbol that occurs first, remove the rest
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names), ]
    dim(exprdat_dupl)
    # set human HGNC symbol as rownames
    rownames(exprdat_dupl) <- exprdat_dupl$HGNC.symbol
    # remove any columns containing IDs
    exprdat_dupl <- subset(exprdat_dupl, select =  -c(Row.names, HGNC.symbol))
    # add rows to data set that contain mean expression values of duplicate human HGNC symbols
    exprdat_dupl <- rbind(exprdat_dupl, mean_entrez)
    # dimension of remaining expression data set:
    # dim(exprdat_dupl)

  }else if(dupl_removal_method == 3){

    ### option 3: among duplicates, keep row with highest overall expression values (i.e highest counts across all samples)

    # intuition: row with highest counts values has  highest power of detecting
    # differential expression later on
    # as in option 2, this applies only to duplicates that result from multiple mouse ENSEMBL IDs
    # that are mapped to the same human HGNC symbol


    # case 2: (case 1 below) multiple mouse ENSEMBL IDs that are converted to the same single human HGNC symbol

    # generate matrix to later contain row with highest count values among ID duplicates
    highest_count_entrez <- matrix(, nrow = 0, ncol = ncol(expression_data))
    # go through each human HGNC symbol that occurs multiple times
    for(i in 1:length(dupl_entrez)){
      # determine all rows with specific human HGNC symbol which occurs multiple times
      counts_dupl <- expression_data[expression_data$HGNC.symbol %in% unique(dupl_entrez)[i], ]
      # detect row with highest count values and order in decreasing manner
      order_rowsums <- order(rowSums(counts_dupl[, 2:(ncol(counts_dupl)-1)]), decreasing = TRUE)
      dupl_id <- counts_dupl[order_rowsums == 1, ]
      # store rounded mean expression value in matrix
      highest_count_entrez <- rbind(highest_count_entrez, dupl_id)
      # View(highest_count_entrez)
      # remove rows in counts_dupl from count data set successively
    }

    # Remove all initial values with ENTREZ duplicates from the dataset
    exprdat_dupl <- expression_data[! expression_data$HGNC.symbol %in% unique(dupl_entrez), ]


    # case 1: single mouse ENSEMBL ID that is mapped to multiple human HGNC symbols
    # as in option 2, pointless to detect row with highest count values as all rows
    # corresponding to the same mouse ENSEMBL ID naturally contain identical count data
    # therefore: remove duplicate mouse ENSEMBL ID that occurs first
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names), ]

    # Add all rows contain initially duplicate human HGNC symbols but contain highest
    # count values among those
    exprdat_dupl <- rbind(exprdat_dupl, highest_count_entrez )

    # Set human HGNC symbols as rownames remove all columns containing any ID info and
    rownames(exprdat_dupl) <- exprdat_dupl$HGNC.symbol
    # Remove any column that contains gene IDs
    exprdat_dupl <- subset(exprdat_dupl, select = -c(Row.names, HGNC.symbol))
    # dim(exprdat_dupl)
  }

  # store resulting gene expression data sets in list (each list element
  # with duplicated gene IDs removed in different manner)
  return(exprdat_dupl)



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

geneID_conversion_SYMBOL_pickrell <- function(expression_data, dupl_removal_method){


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

# for DESeq2 (method = "DESeq2") and edgeR (method = "edgeR")
# rankby must be in c("p_value", "lfc") to perform ranking based on
# (i) p-value (rank = sign(lfc)*(-1)*log10(unadjusted_pvalue)
# (ii) log fold

rankedList_cP <- function(DE_results, rankby, method){

  if(method == "DESeq2"){#create ranking based on DESeq2 results table

    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))

    DE_results$pvalue[ DE_results$pvalue ==  0 ] <- min(DE_results$pvalue[ DE_results$pvalue > 0 ]) / 10

    # DE_results <- edgeR_results
    if(rankby == "lfc"){ #ranking by log2 fold change
      # remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
      rankvec <- as.vector(DE_results[!is.na(DE_results$pvalue), ]$log2FoldChange)
      names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue), ])
      rankvec <- sort(rankvec, decreasing = TRUE)

    }else if (rankby == "p_value"){#ranking by p-value
      # remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
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
  return(rankvec)
}



