#Pre-Processing of the expression data set for FCS methods **PADOG and GSEA (web application)**

#In this file, two functions are defined to perform voom and the variance
#stabilizing transformation on the expression data set(s)

#-> required for PADOG and GSEA


library(edgeR)
library(limma)
library(DESeq2)
library(dplyr)



################################################################################
#option 1: voom transformation ################################################# 
################################################################################

#argument normmethod is in c("TMM","upperquartile")
voom_trans <- function(expression_data, phenotype_labels, normmethod = "TMM"){

  
#step 1: generate DGEList object from the expression data 
counts_voom <- DGEList(expression_data, 
                         group = phenotype_labels)

#step 2: perform normalization (calculate normalization/scaling factor)
counts_voom <- calcNormFactors(counts_voom, 
                                 method = normmethod)

#step 3: perform voom function (do NOT make use of precision weights)
counts_voom <- voom(counts_voom)

#return transformed expression data set 
return(counts_voom$E)
}



################################################################################
#option 2: varianceStabilizingTransformation via DESeq2#########################
################################################################################
# perform a variance stabilizing transformation so that the resulting 
# transformed data are homoscedastic (and can be modeled as Microarray Data,
# according to EnrichmentBrowser Package)

variancetransform <- function(expression_data, phenotype_labels){
  
  #generate data frame that contains information on samples (required format for DESeq2)
  coldata <- data.frame(phenotype_labels,
                        row.names = colnames(expression_data))
  colnames(coldata) <- "condition"
  coldata$condition <- factor(coldata$condition,labels = c("untreated","treated"))
  
  #create DESeqDataSet (required format to later output normalized count data)
  dds <- DESeqDataSetFromMatrix(countData = expression_data, 
                                colData = coldata,
                                design = ~ condition)
  
  #perform variance stabilizing transformation
  #function assay() retrieves transformed counts 
  dds_transform <- vst(dds, 
                       blind = TRUE, 
                       fitType = "parametric") %>% 
                    assay()
  
  return(dds_transform)
  
}

