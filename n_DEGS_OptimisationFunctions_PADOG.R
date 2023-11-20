#PADOG optimization


library(PADOG)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(edgeR) # for pre-filtering function filterByExpr()


# load gene expression data set with true phenotype randomly permuted phenotype assignments 
source("./Random_Phenotype_Permutations.R")

# load require pre-processing functions 
source("./PreProcessing_Functions.R")

# load functions to perform RNA-Seq transformation for (approximate) alignment with the normal distribution
source("./RNASeq_Transformation.R")


####################################################################################
###phenotype label preparation######################################################
####################################################################################
#padog requires character vector with class labels of the samples
#-> can only contain "c" for control samples or "d" for disease samples

phenotype_prep <- function(phenotype_labels){
  
  if(class(phenotype_labels) == "factor"){
    
    # change the levels to required format "c" and "d"
    levels(phenotype_labels)  <-  c("c", "d")
    phenotype_labels <- as.character(phenotype_labels)
    
    return(phenotype_labels)
    
  } else{
  
  # from initial phenotype labels, create vector with levels "c" and "d"
  return(factor(phenotype_labels, 
                levels = c(0,1), 
                labels = c("c","d")))
  
        }
  }


####################################################################################
###lossfunction#####################################################################
####################################################################################

lossfunction_padog <- function(gsa_padog_result){
  
  #count number of differentially enriched gene sets by adjusted p-value
  n_DEGS <- sum(gsa_padog_result$p_adj < 0.05)
  
  return(n_DEGS)
}


####################################################################################
###PADOG Optimization Function#####################################################
####################################################################################

PADOG_optim <- function(expression_data, phenotype_labels){
  
  # identify organism based on gene IDs 
  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism  <- sapply(FUN = grepl, X = c("ENSG", "ENSMUSG"), 
                          x = rownames(expression_data)[1])
  
  # (ii) specification of organism for KEGG
  # hsa: homo sapiens (human)
  # mmu: mus musculus (mouse)
  organism  <- (c("hsa", "mmu"))[ind_organism][[1]]
  
  
  # prepare phenotype permutations to have format as required by function padog()
  phen_prep <- phenotype_prep(phenotype_labels)
  

  ##########
  # 1. step: Choose optimal pre-filtering threshold
  #########
  
  #default threshold: 10 read counts across all samples 
  #alternatives: 
  #(i) thresholds of 0, 20, 50 read counts (0 corresponds to NO pre-filtering)
  #(ii) pre-filtering performed using edgeR's builtin function filterByExpr 
  
  
  #generate documentation frame 
  
  doc <- data.frame(step=c("Default", "Pre-Filtering Threshold", "Duplicate Gene ID removal", "RNA-Seq Data Transformation Method"),
                    optimal_parameter = NA, 
                    n_DEGS = NA)
  
  
  # (i) Default: pre-filtering approach: manual pre-filtering with threshold 10 
  
  exprdat_prefilt_alt <- pre_filt(expression_data, 
                                  threshold = 10)
                                
  
  # (ii) pre-filtering using edgeR's builtin function filterByExpr()
  ind_genes1 <- filterByExpr(expression_data, 
                             group=phenotype_labels)
  #check whether the order of the genes is aligned
  all(names(ind_genes1) == rownames(expression_data))
 
  #add pre-filtered expression data set to the list 
  exprdat_list_prefilt <- list(exprdat_prefilt_alt, 
                               expression_data[ind_genes1,]) 
            
  # Run Padog with both approaches to pre-filtering
  PADOG_prefilt_list  <-  list()
  for(i in 1:length(exprdat_list_prefilt)){
    
    
    # (i) transform ENSEMBL gene IDs to Entrez gene IDs (choose default removal manner) AND
    # transform RNA-Seq data using voom
    
    expression_data_transf  <-  geneID_conversion(exprdat_list_prefilt[[i]], 
                                                  dupl_removal_method = 1) %>% 
                                voom_trans(phenotype_labels = phenotype_labels, 
                                           normmethod = "TMM")

    # (ii) run PADOG
    seed <- 1007
    PADOG_prefilt_list[[i]] <-  padog(expression_data_transf, 
                                      group=phen_prep, 
                                      dseed=seed, 
                                      organism = organism)
    
    # (iii) perform multiple test adjustment using Benjamini and Hochberg 
    PADOG_prefilt_list[[i]]$p_adj  <-  p.adjust(PADOG_prefilt_list[[i]]$Ppadog, 
                                                method="BH")
    
    
  }
  
  # get number of differentially enriched gene sets for each pre-filtering approach
  n_DEGS_prefilt  <-  unlist(lapply(X=PADOG_prefilt_list, 
                                    FUN=lossfunction_padog))
  
  #get index of optimal pre-filtering threshold, i.e. threshold that maximizes the number
  #of differentially enriched gene sets 
  #-> in case of a tie: choose lower index 
  ind_opt_prefilt  <-  min(which(n_DEGS_prefilt == max(n_DEGS_prefilt), 
                                 arr.ind=TRUE))
  
  #update documentation frame 
  #update documentation frame 
  doc[1, "n_DEGS"]  <-  n_DEGS_prefilt[1] #default result
  #optimal result:
  ind_opt_prefilt  <-  min(which(n_DEGS_prefilt == max(n_DEGS_prefilt), 
                                 arr.ind = TRUE)) #index of optimal pre-filtering manner
  
  doc[2, "optimal_parameter"]  <-  c("manual pre-filtering with threshold 10", "filterByExpr")[ind_opt_prefilt] #optimal pre-filtering manner 
  doc[2, "n_DEGS"]  <-  max(n_DEGS_prefilt)
  
  #set current optimal number of differentially enriched gene sets 
  n_DEGS  <-  max(n_DEGS_prefilt)
  #set optimal PRE-FILTERED gene expression data set 
  exprdat_prefilt_opt  <- exprdat_list_prefilt[[ind_opt_prefilt]]
  
  
  
  ##########
  #2. step: Choose optimal manner of removing duplicated gene IDs caused by gene ID conversion 
           #from ENSEMBL to ENTREZ gene ID 
  ##########

  #perform 2 manners of gene ID conversion for the optimally pre-filtered gene expression data set 
  exprdat_convID  <-  lapply(FUN = geneID_conversion, 
                             expression_data = exprdat_prefilt_opt, 
                             X = 1:2)
  
  #generate list containing GSEA results for each "converted" expression data set
  PADOG_results_IDconv_list <- list()
  #perform default DESeq2 workflow for each of the expression data sets
  for(i in 1:length(exprdat_convID)){
    
    # (i) prepare input as required for clusterProfiler's ORA tool
    expression_data_transf  <-  voom_trans(exprdat_convID[[i]], 
                                           phenotype_labels, "TMM")
    
    # (ii) run default ORA (geneset database: KEGG)
    seed <- 1007
    PADOG_results_IDconv_list[[i]] <-  padog(expression_data_transf, 
                                             group= phen_prep, 
                                             dseed=seed,
                                             organism = organism)
    
    # (iii) perform multiple test adjustment using Benjamini and Hochberg 
    PADOG_results_IDconv_list[[i]]$p_adj  <-  p.adjust(PADOG_results_IDconv_list[[i]]$Ppadog, 
                                                       method="BH")
    
  }
  
    #calculate lossfunction for each result and store in vector
    n_DEGS_vec_convID <- unlist(lapply(X= PADOG_results_IDconv_list, 
                                       FUN = lossfunction_padog))
  
  #update documentation
  
  #index of optimal gene ID removal technique 
  ind_opt_conv_ID  <-  min(which(n_DEGS_vec_convID == max(n_DEGS_vec_convID), 
                                 arr.ind = TRUE))
  doc[3, "optimal_parameter"]  <-  ind_opt_conv_ID
  doc[3, "n_DEGS"]  <-  max(n_DEGS_vec_convID)
  
  n_DEGS  <-  max(n_DEGS_vec_convID)
  
  #get optimal PRE-FILTERED gene expression data with converted gene IDs 
  expression_data_opt  <-  exprdat_convID[[ind_opt_conv_ID]]
  
  
  
  ##########
  # 3. step: Choose optimal RNA-Seq transformation method
  
  #alternative to voom transformation (default): DESeq2's varianceStabilizingTransformation
  seed  <-  1007
  PADOG_vst  <-  variancetransform(expression_data_opt, phenotype_labels) %>% 
                                   padog(group = phenotype_prep(phenotype_labels), 
                                         dseed = seed, 
                                         organism = organism) %>%
                                   mutate(p_adj = p.adjust(Ppadog, method="BH"))
   
   #update documentation frame
  doc[4, "optimal_parameter"] <- ifelse(lossfunction_padog(PADOG_vst) > n_DEGS , 
                                        "varianceStabilizingTransformation", 
                                        "voom")  
  
  doc[4, "n_DEGS"]  <-  max(lossfunction_padog(PADOG_vst), n_DEGS)
  
  #set optimal PADOG results which will then be returned by this optimization function 
  if(lossfunction_padog(PADOG_vst) > n_DEGS){
    
    # if alternative transformation method vst leads to an increased number of differentially 
    # enriched gene sets, then update the optimal PADOG results 
    PADOG_opt  <-  PADOG_vst
    
  } else{
    # if vst does not lead to an increased number of differentially enriched gene sets, 
    # then we return to voom (i.e. the optimal PADOG results from the previous step)
    PADOG_opt  <-  PADOG_results_IDconv_list[[ind_opt_conv_ID]]}
  
  
  #return optimal documentation frame 
  return(list(default = PADOG_prefilt_list[[1]], #default results
              optim = PADOG_opt, #optimal results
              documentation = doc)) # documentation frame 
  
  
  
}

################################################################################
### Run Optimization Functions #################################################
################################################################################



# # create list of vectors from phen_pickrell 
# # -> for use of function lapply() since apply() returns an error 
# phen_pickrell_list <- list()
# 
# for(i in 1:ncol(phen_pickrell)){
#   
#   phen_pickrell_list[[i]] <- phen_pickrell[,i]
#   
# }
# 
# 
# # create list of vectors from phen_bottomly 
# # -> for use of function lapply() since apply() returns an error 
# phen_bottomly_list <- list()
# 
# for(i in 1:ncol(phen_bottomly)){
#   
#   phen_bottomly_list[[i]] <- phen_bottomly[,i]
#   
# }
# 
# #############
# ### Pickrell 
# #############
# 
# # original phenotype assignment 
# optim_PADOG_results_Pickrell_originalphenotype <- PADOG_optim(Biobase::exprs(pickrell.eset),
#                                                                                pickrell.eset$gender)
# 
# # save results
# save(optim_PADOG_results_Pickrell_originalphenotype,
# file = "./Results/PADOG_Results_Pickrell_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optim_PADOG_results_Pickrell_phenotypepermutation <- lapply(FUN = PADOG_optim, expression_data = Biobase::exprs(pickrell.eset), 
#                                                              X = phen_pickrell_list)
# 
# # save results
# save(optim_PADOG_results_Pickrell_phenotypepermutation,
# file = "./Results/PADOG_Results_Pickrell_PhenotypePermutations.RData")
# 
# #############
# ### Bottomly 
# #############
# 
# # original phenotype assignment 
# optim_PADOG_results_Bottomly_originalphenotype <- PADOG_optim(Biobase::exprs(bottomly.eset),
#                                                                                bottomly.eset$strain)
# 
# # save results
# save(optim_PADOG_results_Bottomly_originalphenotype,
# file = "./Results/PADOG_Results_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 random permutations of the sample labels 
# optim_PADOG_results_Bottomly_phenotypepermutation <- lapply(FUN = PADOG_optim, expression_data = Biobase::exprs(bottomly.eset), 
#                                                              X = phen_bottomly_list)
# 
# 
# # save results
# save(optim_PADOG_results_Bottomly_phenotypepermutation,
# file = "./Results/PADOG_Results_Bottomly_PhenotypePermutations.RData")
# 
# 
# 
# 
# 
