#######################################################################################
### Maximize number of differentially enriched gene sets obtained with GSEAPreranked ##
#######################################################################################

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(apeglm)
library(ashr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(biomaRt) #note that here, we have to work with biomaRt to convert the gene 
# IDs instead of clusterProfiler's bitr(). The reason for this is that we have 
# to convert between organisms, namely from mouse mouse ENSEMBL IDs to Human HGNC gene symbols, 
# which does not work with bitr().


# set working directory 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/Assessment_OverOptimism")

# load gene expression data set with true phenotype randomly permuted phenotype assignments 
source("./Random_Phenotype_Permutations.R")

# load required pre-processing functions
source("./PreProcessing_Functions.R")

# load functions for RNA-Seq transformation
source("./RNASeq_Transformation.R")


#######################################################################################
###pre-filtering function #############################################################
#######################################################################################

pre_filt<-function(expression_data, threshold){
  
  expression_data_filt<-expression_data[rowSums(expression_data)>=threshold,]
  
  return(expression_data_filt)
  
  
}

#######################################################################################
### Gene ID conversion and duplicate gene ID removal ##################################
#######################################################################################

# particularity of GSEAPreranked: Recommended/ required gene ID format is Hugo symbols 
# and not human HGNC symbols 

# background: GSEA typically run using genesets from MSigDB which consists of human gene symbols 
#If input data contain other identifiers, the IDs need to be converted to gene symbols
# option "Collapse/ Remap to gene symbols" performs conversion which handles the case of 
# several feature identifiers mapping to same gene identifier.
# This method was developed and tuned for gene expression data, however, the ranked list
# of genes in GSEAPreranked was created using unspecified ranking procedure outside of GSEA 

# -> recommended to provide ranked list with genes already converted to gene SYMBOLS and
# select parameter "NO_Collapse"

geneID_conversion_SYMBOL <- function(expression_data, dupl_removal_method){
  
  # Conversion of the gene IDs from Mouse mouse ENSEMBL IDs to human HGNC symbols using 
  # biomaRt 
  

  # IMPORTANT NOTE: using the default host (i.e. current version of the ensembl website) returned an error
  # we therefore had to work with archived versions: we went down the list in listEnsemblArchives() and 
  # selected the most recent archive that works 
  
  
  # connect to Biomart database for mouse 
  mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  # connect to Biomart database for human 
  human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  
  mouseEnsembl_to_humanSymbol <- getLDS(
    mart = mouse,
    attributes = c('ensembl_gene_id'),
    martL = human,
    attributesL = c('hgnc_symbol'),
    #filters = 'ensembl_gene_id',
    values = rownames(expression_data)
  )
  
  
  #results: 
  #-not all mouse mouse ENSEMBL IDs can be mapped to a corresponding human HGNC symbol
  #-some single mouse mouse ENSEMBL IDs were mapped to multiple distinct human HGNC symbols
  #-some distincts mouse mouse ENSEMBL IDs were mapped to an identical human HGNC symbol
  
  ### merge
  
  #this step is independent of sample IDs
  #merge by row names of expression data set and mouse ENSEMBL ID of conversion data set
  expression_data <- merge(expression_data, mouseEnsembl_to_humanSymbol, by.x=0, by.y="Gene.stable.ID", all = FALSE)
  dim(expression_data)
  
  ###take closer look at duplicates 
  
  #CASE 1: single mouse ENSEMBL IDs are mapped to multiple human HGNC symbols
  #View(bitr_toKEGG[(duplicated(bitr_toKEGG$ENSEMBL)),])
  sum(duplicated(mouseEnsembl_to_humanSymbol$Gene.stable.ID)) #number of times an mouse ENSEMBL ID was converted to several human HGNC symbols
  #determine all duplicated mouse ENSEMBL IDS
  dupl_ensembl<-unique(mouseEnsembl_to_humanSymbol$Gene.stable.ID[duplicated(mouseEnsembl_to_humanSymbol$Gene.stable.ID)])
  #number of mouse ENSEMBL IDs that have at least one duplicate
  length(dupl_ensembl)
  #display of conversion scheme of duplicated mouse ENSEMBL IDs
  duplicated_conversion_ens<-mouseEnsembl_to_humanSymbol[mouseEnsembl_to_humanSymbol$Gene.stable.ID %in% dupl_ensembl,]
  dim(duplicated_conversion_ens)
  
  
  #CASE 2: multiple mouse ENSEMBL IDs are mapped to single human HGNC symbol
  sum(duplicated(mouseEnsembl_to_humanSymbol$HGNC.symbol)) #number of times several mouse ENSEMBL IDs were converted to a single human HGNC symbol
  #determine all duplicated human HGNC symbols
  dupl_entrez<-unique(mouseEnsembl_to_humanSymbol$HGNC.symbol[duplicated(mouseEnsembl_to_humanSymbol$HGNC.symbol)])
  #number of human HGNC symbols that have at least one duplicate
  length(dupl_entrez)
  #display of conversion scheme of duplicated human HGNC symbolmouseEnsembl_to_humanSymbol$HGNC.symbols
  duplicated_conversion_entrez<-mouseEnsembl_to_humanSymbol[mouseEnsembl_to_humanSymbol$HGNC.symbol %in% dupl_entrez,]
  dim(duplicated_conversion_entrez)
  
  
  
  
  #######
  #2. step: Removal of Duplicated gene IDs 
  #######
  
  
  if(dupl_removal_method == 1){
    
    ### 1. option: keep first subscript among duplicates #########################
    
    #1. remove duplicated human HGNC symbols
    exprdat_dupl<-expression_data[!duplicated(expression_data$HGNC.symbol),]
    dim(expression_data)
    
    #2. remove duplicated mouse ENSEMBL IDs
    exprdat_dupl<-exprdat_dupl[!duplicated(exprdat_dupl$Row.names),]
    dim(exprdat_dupl)
    
    #3. human HGNC symbols as row names and 
    rownames(exprdat_dupl)<-exprdat_dupl$HGNC.symbol
    #Remove columns containing ENSEMBL and human HGNC symbols
    exprdat_dupl<-subset(exprdat_dupl, select=-c(Row.names,HGNC.symbol))
    dim(exprdat_dupl)
    
  } else if(dupl_removal_method == 2){
    
    ### 2. option: keep mean expression value of all duplicated gene IDs  #########################
    
    
    
    #generate matrix to contain (rounded) mean expression values of all rows that 
    #have same human HGNC symbol
    #ncol=ncol(expression_data)-2 since data set contains 2 columns with IDs at this point
    mean_entrez<-matrix(, nrow=0, ncol=ncol(expression_data)-2)
    
    
    # -> There are cases where no mouse ENSEMBL IDs are mapped to an identical human HGNC symbol
    # -> nrow(dupl_entrez) = 0 
    # -> add if-condition since the following part of code would otherwise produce NaN values which would in turn lead to errors
    # when pre-processing the gene expression data set 
    
    # -> if no distinct mouse ENSEMBL IDs are mapped to an identical human HGNC symbol then mean_entrez remains an empty data frame
    
    
    if(length(dupl_entrez) !=0){
      
      
      #1.remove duplicated human HGNC symbols (case 2)
      #i.e. multiple different mouse ENSEMBL IDs that are mapped to the same single human HGNC symbol
      
      for(i in 1:length(dupl_entrez)){#go through each human HGNC symbols which occurs multiple times
        #determine all rows whose human HGNC symbols correspond to current human HGNC symbol 
        counts_dupl<-expression_data[expression_data$HGNC.symbol %in% unique(dupl_entrez)[i],]
        #for rows duplicated human HGNC symbol compute (rounded) mean expression value 
        dupl_id<-round(colMeans(counts_dupl[,c(2:(ncol(expression_data)-1))]))
        #store rounded mean expression value in matrix 
        mean_entrez<-rbind(mean_entrez,dupl_id)
      }
    }
    
    
    #test whether the number of rows in mean_entrez corresponds to the number human HGNC symbols
    #that occur more than once
    nrow(mean_entrez)==length(dupl_entrez)
    
    #remove all rows from the expression data whose human HGNC symbol has at least one duplicate
    exprdat_dupl<-expression_data[!expression_data$HGNC.symbol %in% dupl_entrez,]
    #test whether number of rows in resulting data set equals nrow of inital data set 
    #minus number of genes with at least one duplicate
    nrow(exprdat_dupl)==nrow(expression_data)-nrow(duplicated_conversion_entrez)
    dim(exprdat_dupl)
    
    #set corresponding human HGNC symbols as rownames
    rownames(mean_entrez)<-unique(dupl_entrez)
    
    
    
    
    #2. remove duplicated mouse ENSEMBL IDs
    #caution: single mouse ENSEMBL IDs that are mapped to multiple human HGNC symbol naturally generate
    #identical count data for all corresponding human HGNC symbols
    #->pointless to compute mean expression values
    #verifiable by looking at data set only containing those mouse ENSEMBL IDs that are
    #mapped by multiple human HGNC symbols:
    #test_dupl_ensembl<-expression_data[expression_data$Row.names %in% dupl_ensembl,]
    #View(test_dupl_ensembl)
    
    #therefore: proceed as in option 1 and use human HGNC symbol that occurs first, remove the rest
    exprdat_dupl<-exprdat_dupl[!duplicated(exprdat_dupl$Row.names),]
    dim(exprdat_dupl)
    #set human HGNC symbol as rownames
    rownames(exprdat_dupl)<-exprdat_dupl$HGNC.symbol
    #remove any columns containing IDs
    exprdat_dupl<-subset(exprdat_dupl,select= -c(Row.names,HGNC.symbol))
    #add rows to data set that contain mean expression values of duplicate human HGNC symbols
    exprdat_dupl<-rbind(exprdat_dupl,mean_entrez)
    #dimension of remaining expression data set:
    #dim(exprdat_dupl)
    
  }else if(dupl_removal_method ==3){
    
    ###option 3: among duplicates, keep row with highest overall expression values (i.e highest counts across all samples)
    
    #intuition: row with highest counts values has  highest power of detecting
    #differential expression later on
    #as in option 2, this applies only to duplicates that result from multiple mouse ENSEMBL IDs
    #that are mapped to the same human HGNC symbol
    
    
    #case 2: (case 1 below) multiple mouse ENSEMBL IDs that are converted to the same single human HGNC symbol
    
    #generate matrix to later contain row with highest count values among ID duplicates
    highest_count_entrez<-matrix(, nrow=0, ncol=ncol(expression_data))
    #go through each human HGNC symbol that occurs multiple times
    for(i in 1:length(dupl_entrez)){
      #determine all rows with specific human HGNC symbol which occurs multiple times
      counts_dupl<-expression_data[expression_data$HGNC.symbol %in% unique(dupl_entrez)[i],]
      #detect row with highest count values and order in decreasing manner
      order_rowsums<-order(rowSums(counts_dupl[,2:(ncol(counts_dupl)-1)]),decreasing=TRUE)
      dupl_id<-counts_dupl[order_rowsums==1,]
      #store rounded mean expression value in matrix 
      highest_count_entrez<-rbind(highest_count_entrez,dupl_id)
      #View(highest_count_entrez)
      #remove rows in counts_dupl from count data set successively
    }
    
    #Remove all initial values with ENTREZ duplicates from the dataset
    exprdat_dupl<-expression_data[! expression_data$HGNC.symbol %in% unique(dupl_entrez),]
    
    
    #case 1: single mouse ENSEMBL ID that is mapped to multiple human HGNC symbols 
    #as in option 2, pointless to detect row with highest count values as all rows
    #corresponding to the same mouse ENSEMBL ID naturally contain identical count data
    #therefore: remove duplicate mouse ENSEMBL ID that occurs first 
    exprdat_dupl<-exprdat_dupl[!duplicated(exprdat_dupl$Row.names),]
    
    #Add all rows contain initially duplicate human HGNC symbols but contain highest
    #count values among those 
    exprdat_dupl<-rbind(exprdat_dupl,highest_count_entrez )
    
    #Set human HGNC symbols as rownames remove all columns containing any ID info and
    rownames(exprdat_dupl)<-exprdat_dupl$HGNC.symbol
    #Remove any column that contains gene IDs
    exprdat_dupl<-subset(exprdat_dupl, select=-c(Row.names,HGNC.symbol))
    #dim(exprdat_dupl)
  }
  
  #store resulting gene expression data sets in list 
  return(exprdat_dupl)
  
  
  
}

##################################################################################
##create ranked list from DE results##############################################
##################################################################################

#for DESeq2 (method="DESeq2") and edgeR (method="edgeR")
#rankby must be in c("p_value", "lfc") to perform ranking based on 
#(i) p-value (rank=sign(lfc)*(-1)*log10(unadjusted_pvalue)
#(ii) log fold changes 
rankedList_cP<-function(DE_results, rankby, method){
  
  if(method=="DESeq2"){#create ranking based on DESeq2 results table
    
    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))
    
    DE_results$pvalue[ DE_results$pvalue == 0 ] <- min(DE_results$pvalue[ DE_results$pvalue > 0 ]) / 10 
    
    #DE_results<-edgeR_results
    if(rankby=="lfc"){ #ranking by log2 fold change
      #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
      rankvec<-as.vector(DE_results[!is.na(DE_results$pvalue),]$log2FoldChange)
      names(rankvec)<-rownames(DE_results[!is.na(DE_results$pvalue),])
      rankvec<-sort(rankvec, decreasing=TRUE)
      
    }else if (rankby=="p_value"){#ranking by p-value
      #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
      rankvec<-as.vector(sign(DE_results[!is.na(DE_results$pvalue),]$log2FoldChange)*(-1)*log10(DE_results[!is.na(DE_results$pvalue),]$pvalue))
      names(rankvec)<-rownames(DE_results[!is.na(DE_results$pvalue),])
      rankvec<-sort(rankvec, decreasing=TRUE)
    }
  }
  
  
  else if(method=="edgeR"){#create ranking based on edgeR results table
    
    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))
    
    DE_results$table$PValue[ DE_results$table$PValue == 0 ] <- min(DE_results$table$PValue[ DE_results$table$PValue > 0 ]) /10
    
    if(rankby=="lfc"){#ranking based on log2 fold change
      rankvec<-as.vector(DE_results$table$logFC)
      names(rankvec)<-rownames(DE_results)
      rankvec<-sort(rankvec, decreasing=TRUE)
    }
    
    else if(rankby=="p_value"){#ranking based on p-value
      rankvec<-as.vector(sign(DE_results$table$logFC)*(-1)*log10(DE_results$table$PValue))
      names(rankvec)<-rownames(DE_results)
      rankvec<-sort(rankvec, decreasing=TRUE)
    }
  }
  
  else if(method=="limma"){#create ranking based on edgeR results table
    
    # first step: replace p-values of 0 with the smallest representable positive number
    # in R (necessary since one term of ranking metric is equal to log10(p-value))
    
    DE_results$P.Value[ DE_results$P.Value == 0 ] <- min(DE_results$P.Value[ DE_results$P.Value > 0 ]) / 10 
    
    if(rankby=="lfc"){#ranking based on log2 fold change
      rankvec<-as.vector(DE_results$logFC)
      names(rankvec)<-rownames(DE_results)
      rankvec<-sort(rankvec, decreasing=TRUE)
    }
    
    else if(rankby=="p_value"){#ranking based on p-value
      rankvec<-as.vector(sign(DE_results$logFC)*(-1)*log10(DE_results$P.Value))
      names(rankvec)<-rownames(DE_results)
      rankvec<-sort(rankvec, decreasing=TRUE)
    }
  }
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
DESeq2_ranking_phenorig <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% geneID_conversion_SYMBOL(dupl_removal_method = 1) %>% 
  deseq_preprocess(phenotype_labels = bottomly.eset$strain ) %>% DESeq() %>%
  lfcShrink(coef="condition_treated_vs_untreated", type="apeglm") %>% as.data.frame() %>%
  rankedList_cP(rankby = "p_value", method = "DESeq2")

# create path for storage of DESeq2 ranking
path_DESeq2_phenorig <- "/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/GSEAPreranked/Bottomly/p_value/Data/Raw/Original_Phenotype/DESeq2_ranking_phenOrig.txt"

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
keep_phenorig <- DGEList(Biobase::exprs(bottomly.eset), group = bottomly.eset$strain) %>% filterByExpr()

# design matrix 
mm_phenorig <- model.matrix( ~ bottomly.eset$strain)

# generate limma results
limma_results <- geneID_conversion_SYMBOL(Biobase::exprs(bottomly.eset)[keep_phenorig,], dupl_removal_method = 1) %>% 
  DGEList(group = bottomly.eset$strain) %>% calcNormFactors() %>%
  voom(design=mm_phenorig) %>% lmFit(design=mm_phenorig) %>% eBayes() %>% topTable(coef=ncol(mm_phenorig), number=100000)

# create limma results and rank by p-value 
limma_ranking_phenorig <- rankedList_cP(limma_results, rankby= "p_value", method="limma")

# Create path for storage of limma ranking            
path_limma_phenorig <- "/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/GSEAPreranked/Bottomly/p_value/Data/Raw/Original_Phenotype/limma_ranking_phenOrig.txt"


# export 
write.table(limma_ranking_phenorig, 
            file = path_limma_phenorig, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = FALSE)   



# (II) Random Phenotype Permutations of Original Phenotypes 

for(i in 1:ncol(phen_bottomly)){
  
  ################
  # DESeq2 ranking
  ################

  # create DESeq2 results and rank by p-value
  DESeq2_ranking_phenperm <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% geneID_conversion_SYMBOL(dupl_removal_method = 1) %>%
    deseq_preprocess(phenotype_labels = phen_bottomly[,i] ) %>% DESeq() %>%
    lfcShrink(coef="condition_treated_vs_untreated", type="apeglm") %>% as.data.frame() %>%
    rankedList_cP(rankby = "p_value", method = "DESeq2")

  # create path for storage of DESeq2 ranking
  path_DESeq2_phenperm <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/GSEAPreranked/Bottomly/p_value/Demethylation/Data/Raw/Phenotype_Permutation",i,"/DESeq2_ranking_permutation",i,".txt")

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
  keep_phenperm <- DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,i]) %>% filterByExpr()
  
  # design matrix 
  mm_phenperm <- model.matrix( ~ phen_bottomly[,i])
  
  # create limma results 
  limma_results_phenperm <- geneID_conversion_SYMBOL(Biobase::exprs(bottomly.eset)[keep_phenperm,], dupl_removal_method = 1) %>% 
    DGEList(group = phen_bottomly[,i]) %>% calcNormFactors() %>%
    voom(design=mm_phenperm) %>% lmFit(design=mm_phenperm) %>% eBayes() %>% topTable(coef=ncol(mm_phenperm), number=100000) 
    
  
  # create ranking by p-value
  limma_ranking_phenperm <- rankedList_cP(limma_results_phenperm, rankby= "p_value", method="limma")
  
  # Create path for storage of limma ranking            
  path_limma_phenperm <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/GSEAPreranked/Bottomly/p_value/Data/Raw/Phenotype_Permutation",i,"/limma_ranking_permutation",i,".txt")
  
  
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

# -> adj. p-value = 1


#########
# 2. step: change method to generate ranking
#########

# alternative 1: limma -> adj. p-value = 1

# -> proceed with the default ranking generated using DESeq2 


#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 1

### -> final results: adj. p-value = 1
# optimal setting coincides with default setting 



################################################################################
### Random Phenotype Permutations ##############################################
################################################################################


################################################################################
### Phenotype Permutation 1  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> ad. p-value = 0.6604594



#########
# 2. step: change method to generate ranking
#########

# alternative: limma, ranking by p-value -> adj. p-value = 0.759344
# -> return to alternative ranking generated using DESeq2 

#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 0.93330926
# alternative 2: exponent 1.5 -> adj. p-value = 0.8523406
# alternative 3: exponent 2 -> adj. p-value = 0.7235578


# ->> final results: adj. p-value = 0.6604594
# default configuration coincides with optimal configuration


################################################################################
### Phenotype Permutation 2  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj. p-value = 1


#########
# 2. step: change method to generate ranking
#########

# alternative 1: limma, ranking by p-value -> adj. p-value = 0.7542026

# -> proceed with alternative ranking generated using limma

#########
# 4. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 0.97179914
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.9234213



### -> final results: adj. p-value = 0.7542026
# achieved with ALTERNATIVE ranking generated using limma




################################################################################
### Phenotype Permutation 3  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj. p.value = 0.9798141


#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> p.adj = 0.9576295

# -> proceed with alternative ranking generated using limma 



#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 0.7589818
# alternative 2: exponent 1.5 -> adj. p-value = 0.9908181
# alternative 3: exponent 2 -> adj. p-value = 0.9933294

### final results: adj. p-value = 0.7589818
# achieved with 
# ALTERNATIVE ranking generated using limma 
# ALTERNATIVE exponent 0 



################################################################################
### Phenotype Permutation 4  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj p-value = 0.57627034

#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> adj. p-value = 0.9729733

# -> proceed with default ranking generated using DESeq2


#########
# 4. step: change exponent 
######### 

# alternative 1: exponent 0 -> adj. p-value = 0.8624485
# alternative 2: exponent 1.5 -> adj. p-value = 0.6021893
# alternative 3: exponent 2 -> adj. p-value = 0.7074654

# final results: adj. p-value = 0.57627034
# optimal configuration coincides with default configuration of parameters 





################################################################################
### Phenotype Permutation 5  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj. p-value = 0.822709


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> adj. p-value = 0.929875

# -> return to default ranking generated using DESeq2 



#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 0.9414946
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.9689939


### final results: adj. p-value = 0.822709
# optimal configuration of coincides with default configuration


################################################################################
### Phenotype Permutation 6  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj. p-value = 0.937159


#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> adj. p-value = 1

# -> return to default ranking generated using DESeq2


#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.83406526

# final results: adj. p-value =  0.83406526
# achieved with ALTERNATIVE exponent 2 



################################################################################
### Phenotype Permutation 7  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj. p-value = 0.88564014

#########
# 2. step: change method to generate ranking
#########

# alternative : limma -> adj. p-value = 0.7940431

# proceed with alternative ranking generated using limma 

#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 0.82804126
# alternative 2: exponent 1.5 -> adj. p-value = 0.7937083
# alternative 3: exponent 2 -> adj. p-value = 0.9930696



# final results: adj. p-value = 0.7937083
# achieved with 
# ALTERNATIVE ranking using limma 
# ALTERNATIVE exponent 1.5 




################################################################################
### Phenotype Permutation 8  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj. p-value = 0.62759393


#########
# 2. step: change method to generate ranking
#########

# alternative 1: limma -> adj. p-value = 0.6549921

# proceed with default ranking generated using DESeq2

#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 0.8552458
# alternative 2: exponent 1.5 -> adj. p-value = 0.75988406
# alternative 3: exponent 2 -> adj. p-value = 0.43375164


### final results: adj. p-value = 0.43375164
# achieved with ALTERNATIVE exponent 2 


################################################################################
### Phenotype Permutation 9  ###################################################
################################################################################


#########
# 1. step: Default
#########

# -> adj. p-value = 0.99037695

#########
# 2. step: change method to generate ranking
#########

# alternative: limma -> adj. p-value = 0.45262173

# -> proceed with alternative ranking generated using limma 

#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 -> adj. p-value = 0.21605846
# alternative 2: exponent 1.5 -> adj. p-value = 0.058096815
# alternative 3: exponent 2 -> adj. p-value = 0.468571

### ->> final results: final results: adj. p-value = 0.058096815
# achieved with
# ALTERNATIVE ranking generated using limma 
# ALTERNATIVE exponent 1.5



################################################################################
### Phenotype Permutation 10  ##################################################
################################################################################


#########
# 1. step: Default
#########

# ->  adj. p-value = 0.8831213



#########
# 2. step: change method to generate ranking
#########

# alternative : limma, ranking by p-value  -> adj. p-value = 0.8108312

# -> proceed with alternative ranking generated using limma 


#########
# 3. step: change exponent 
#########

# alternative 1: exponent 0 ->  adj. p-value = 0.98778003
# alternative 2: exponent 1.5 -> adj. p-value = 0.93723977
# alternative 3: exponent 2 -> adj. p-value = 0.64864224


### -> final results: adj. p-value = 0.64864224
# achieved with 
# ALTERNATIVE ranking generated using limma 
# ALTERNATIVE exponent 2







