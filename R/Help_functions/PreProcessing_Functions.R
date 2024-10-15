## Definition of pre-processing functions required for (almost) all of the computational
# GSA methods

# load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(dplyr)

#######################################################################################
###input preparation for DESeq2 #######################################################
#######################################################################################

# generate input required for DESeq2

# function inputs:
# - data: gene expression data to be processed
# - phenotype_labels: vector of binary labels indicating the status of each sample
#       in the gene expression data set

deseq_preprocess <- function(data, phenotype_labels){

  #check whether the number of samples in data equals the length of the phenotype vector
  if(ncol(data) != length(phenotype_labels)) stop("Error: number of samples in expression data and p
  phenotype labels do not match")

  #generate data frame containing sample labels of respective samples and column name "condition"
  coldata <- data.frame(phenotype_labels,
        row.names = colnames(data))

  colnames(coldata) <- "condition"
  coldata$condition <- factor(coldata$condition, labels = c("untreated", "treated"))

  #generate DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = data, #submit expression data set of counts
            colData = coldata, #as generated above
            design = ~ condition)

  # return gene expression data in correct format for DESeq2
  return(dds)

}


#######################################################################################
###pre-filtering function #############################################################
#######################################################################################

# manual pre-filtering: Remove all genes with a number of counts across all samples
# lower than the specified threshold

# function arguments:
# expression_data: gene expression data set
# threshold: minimum number of counts across all samples below which a gene is
#            removed from the data set
pre_filt <- function(expression_data, threshold){

  expression_data_filt <- expression_data[rowSums(expression_data) >= threshold, ]

  # return filtered gene expression data set
  return(expression_data_filt)


}





#######################################################################################
### Gene ID conversion and duplicate gene ID removal ##################################
#######################################################################################

# function inputs:
# - expression_data: gene expression data set with genes in Ensembl ID format
# - dupl_removal_method: manner of removing duplicates that result from the conversion
#           "1": keep row with lowest subscript
#           "2": keep mean expression value of all duplicated gene IDs
geneID_conversion <- function(expression_data, dupl_removal_method){

  # identify organism for which gene expression is measured from the format of the gene IDs (all are Ensembl)
  # ENSEMBL gene ID ENSMUSGXXXXXXXXXXX corresponds to mouse; identifiable by substring "ENSMUSG"
  # ENSEMBL gene ID ENSGXXXXXXXXXXX corresponds to homo sapiens; identifiable by string "ENSG"

  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl,
         X = c("ENSG", "ENSMUSG"),
         x = rownames(expression_data)[1])

  # choose corresponding organism (required for function bitr)
  organism <- unlist(c(org.Hs.eg.db , org.Mm.eg.db)[ind_organism])[[1]]


  #Gene ID conversion via clusterProfiler's function bitr()
  # Conversion from ENSEMBL to Entrez gene ID
  bitr_enstoentr <- bitr(rownames(expression_data),
         fromType = "ENSEMBL",
         toType = "ENTREZID",
         OrgDb = organism)

  #note: not all ENSEMBL IDs could be converted to a corresponding ENTREZ Gene ID:
  dim(bitr_enstoentr)

  #results:
  # - not all Ensembl gene IDs can be mapped to a corresponding Entrez gene ID
  # - some single Ensembl gene IDs were mapped to multiple distinct Entrez gene IDs
  # - some distincts Ensembl gene IDs were mapped to an identical Entrez gene ID

  ### merge

  #this step is independent of sample IDs
  #merge by row names of expression data set and ENSEMBL ID of conversion data set
  expression_data <- merge(expression_data, bitr_enstoentr,
         by.x = 0,
         by.y = "ENSEMBL",
         all.y = TRUE,
         sort = TRUE)
  # inspect dimension
  dim(expression_data)

  ###take closer look at duplicates

  #CASE 1: single ENSEMBL IDs are mapped to multiple ENTREZ IDs
  #View(bitr_toKEGG[(duplicated(bitr_toKEGG$ENSEMBL)), ])
  sum(duplicated(bitr_enstoentr$ENSEMBL)) #number of times an ENSEMBL gene ID was converted to several ENTREZ IDs
  #determine all duplicated ENSEMBL gene IDS
  dupl_ensembl <- unique(bitr_enstoentr$ENSEMBL[duplicated(bitr_enstoentr$ENSEMBL)])
  #number of ENSEMBL IDs that have at least one duplicate
  length(dupl_ensembl)
  #display of conversion scheme of duplicated ENSEMBL IDs
  duplicated_conversion_ens <- bitr_enstoentr[bitr_enstoentr$ENSEMBL %in% dupl_ensembl, ]
  dim(duplicated_conversion_ens)


  #CASE 2: multiple ENSEMBL IDs are mapped to single entrez ID
  sum(duplicated(bitr_enstoentr$ENTREZ)) #number of times several ENSEMBL gene IDs were converted to a single ENTREZ ID
  #determine all duplicated ENTREZ IDs
  dupl_entrez <- unique(bitr_enstoentr$ENTREZID[duplicated(bitr_enstoentr$ENTREZID)])
  #number of ENTREZ IDs that have at least one duplicate
  length(dupl_entrez)
  #display of conversion scheme of duplicated ENTREZ IDs
  duplicated_conversion_entrez <- bitr_enstoentr[bitr_enstoentr$ENTREZID %in% dupl_entrez, ]
  dim(duplicated_conversion_entrez)




  #######
  #2. step: Removal of Duplicated gene IDs
  #######


  if(dupl_removal_method == 1){

  ### 1. option: keep first subscript among duplicates #########################

  #1. remove duplicated ENTREZ gene IDs
  exprdat_dupl <- expression_data[!duplicated(expression_data$ENTREZID), ]
  dim(expression_data)

  #2. remove duplicated ENSEMBL gene IDs
  exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names), ]
  dim(exprdat_dupl)

  #3. ENTREZ IDs as row names and
  rownames(exprdat_dupl) <- exprdat_dupl$ENTREZID
  #Remove columns containing ENSEMBL and ENTREZ IDs
  exprdat_dupl <- subset(exprdat_dupl, select = -c(Row.names, ENTREZID))
  dim(exprdat_dupl)

  } else if(dupl_removal_method == 2){

  ### 2. option: keep mean expression value of all duplicated gene IDs #########################

  #1.remove duplicated ENTREZ gene IDs (case 2) IF THERE ARE ANY

  if(length(dupl_entrez) != 0){
  #i.e. multiple different ENSEMBL IDs that are mapped to the same single ENTREZ ID

  #generate matrix to contain (rounded) mean expression values of all rows that
  #have same ENTREZ gene ID
  #ncol = ncol(expression_data)-2 since data set contains 2 columns with IDs at this point
  mean_entrez <- matrix(, nrow = 0, ncol = ncol(expression_data)-2)

  for(i in 1:length(dupl_entrez)){#go through each ENTREZ IDs which occurs multiple times
    #determine all rows whose ENTREZ IDs correspond to current ENTREZ ID
    counts_dupl <- expression_data[expression_data$ENTREZID %in% unique(dupl_entrez)[i], ]
    #for rows duplicated ENTREZ ID compute (rounded) mean expression value
    dupl_id <- round(colMeans(counts_dupl[, c(2:(ncol(expression_data)-1))]))
    #store rounded mean expression value in matrix
    mean_entrez <- rbind(mean_entrez, dupl_id)
  }


  #test whether the number of rows in mean_entrez corresponds to the number ENTREZ IDs
  #that occur more than once
  nrow(mean_entrez) == length(dupl_entrez)

  #remove all rows from the expression data whose ENTREZ ID has at least one duplicate
  exprdat_dupl <- expression_data[!expression_data$ENTREZID %in% dupl_entrez, ]

  #test whether number of rows in resulting data set equals nrow of inital data set
  #minus number of genes with at least one duplicate
  nrow(exprdat_dupl) == nrow(expression_data) - nrow(duplicated_conversion_entrez)
  dim(exprdat_dupl)

  #set corresponding ENTREZ gene IDs as rownames
  rownames(mean_entrez) <- unique(dupl_entrez)

  }else{

  # if there are no duplicated Entrez IDs then we start with the "initial" entered as
  # a function argument at this point

  exprdat_dupl <- expression_data
  }


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
  rownames(exprdat_dupl) <- exprdat_dupl$ENTREZID

  #remove any columns containing IDs
  exprdat_dupl <- subset(exprdat_dupl, select = -c(Row.names, ENTREZID))

  if(length(dupl_entrez) != 0){
  # If there were any duplicated Entrez IDs removed in step 1, then we now add
  # the corresponding mean expression values here

  #add rows to data set that contain mean expression values of duplicate ENTREZ IDs
  exprdat_dupl <- rbind(exprdat_dupl, mean_entrez)

  }
  #dimension of remaining expression data set:
  dim(exprdat_dupl)

  }

  #store resulting gene expression data sets in list (with converted gene IDs and
  # duplicates removed according to function argument dupl_removal_method)
  return(exprdat_dupl)



}

