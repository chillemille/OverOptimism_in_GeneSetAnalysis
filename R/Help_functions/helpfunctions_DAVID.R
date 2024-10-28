################################################################################
### help function DAVID ########################################################
################################################################################

#########################################
###generate necessary input for DAVID ###
#########################################

DAVID_input_preparation <- function(DE_results){

  #required input for clusterProfiler function: vector of entrez gene ID
  #-> need to pre-process results table DE_results

  #vector of differentially expressed genes
  #DEG_vec serves as input vector for ORA performed by clusterProfiler



  # classify those genes as DE that have an adjusted p-value < 0.05
  DEG_vec <- rownames(DE_results[(DE_results$p_adj < 0.05) & (!is.na(DE_results$p_adj)), ])


  # return vector of differentially expressed genes
  return(DEG_vec)

}
