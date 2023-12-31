#PADOG: Minimize adjusted p-value or rank of a specific gene set

library(PADOG)
library(dplyr)
library(clusterProfiler)
library(edgeR)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


# load gene expression data set with true phenotype randomly permuted phenotype assignments 
source("./Random_Phenotype_Permutations.R")

# load required pre-processing functions
source("./PreProcessing_Functions.R")

# load functions for RNA-Seq transformation
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


#######################################################################################
###pre-filtering function #############################################################
#######################################################################################

pre_filt <- function(expression_data, threshold){
 
 expression_data_filt <- expression_data[rowSums(expression_data) >= threshold, ]
 
 return(expression_data_filt)
 
 
}


################################################################################
### get gene's adjusted p-value or rank (based on what is needed) ##############
################################################################################

############### NOTE: 
# encountered problem in the optimization process of the rank or (adjusted) pvalue of a given gene set:
#clusterProfiler only reports those gene sets which contain at least a single gene from the input list
#of differentially expressed genes. 

# this specifically means that, if a gene set does not contain any genes from the input list, function 
#pvalue_rank returns integer(0) which cannot be compared to numerical number. 
#To bypass this issue we set the rank resp. adjusted p-value to 0 if a gene set is not reported in the results.
#p-value of a given gene set to Inf 


###detect whether gene set does not contain any genes from the input list of differentially expressed genes 
is.integer0 <- function(x){
 is.integer(x) && length(x) == 0L
}
##############


#in the case that a gene's
#(i) (adjusted) p-value is needed: enter "p_adj" in argument metric 
#(ii) relative rank among the remaining gene sets is needed: enter "rank" in argument metric

#note: function works for GO as well as KEGG as results table "padog_results" 
#contains relevant information automatically
#argument term must be in the form of a GO ID resp. KEGG ID 
pvalue_rank_padog <- function(term, padog_results, metric){

 
 #check whether there are any gene sets reported in the results 
 if (nrow(padog_results) == 0){ #if there are no gene sets then the value 1 shall be returned
  #-> no gene sets also means that the gene set of interest does not contain any 
  #of the DE genes 
  
  return(1)
  
 } else {


 if(metric == "rank"){
   
   # prepare the relative ranks such that all gene sets with the same adjusted p-value
   # are given the same (i.e. average) rank
   # We divide by the maximum rank such that those gene sets with the highest 
   # adjusted p-values are given the (relative) rank 1, which is the worst rank
   
   padog_results$ranks <- rank(padog_results$p_adj) /max(rank(padog_results$p_adj)) 
   
   # get the relative rank of the gene set of interest 
   rank <- padog_results$ranks[grep(term, padog_results$ID)]
   
   
  # if the gene set is not contained in the results, return the rank 1.2
  # this shall serve as an indicator in the results that the gene set 
  # was not contained in the results 
  return(ifelse(!is.integer0(grep(term, padog_results$ID)), 
         rank, 
         1))
   
  #note: in the case that a gene set is not reported in the results table of padog_results,
  #ifelse() in combination with !is.integer0() then ensures that a rank of Inf is returned,
  #meaning that each adaption leading to any infinite rank is considered an improvement
  
 }else if(metric == "p_adj"){
  
  #identify row number of respective gene set 
  ind_row <- grep(term, padog_results$ID)
  
  #return respective adjusted p-value
  return(ifelse(!is.integer0(ind_row), 
         padog_results$p_adj[ind_row], 
         1))
  
  #note: in the case that a gene set is not reported in the results table of padog_results, 
  #ifelse() in combination with !is.integer0() then ensures that an adjusted p-value of 1.2 is returned,
  #meaning that each adaption leading to a an adjusted p-value in (0,1] is considered an improvement
 }
 }
}





####################################################################################
###PADOG Optimization Function#####################################################
####################################################################################

#expression data sets must be provided in list
#default transformation: voom 
PADOG_rankp_optim <- function(geneset, expression_data, phenotype_labels, metric){
 
 
 #####################################
 # identify organism based on gene IDs 
 #####################################
 
 # indicate which of the two strings can be found in the gene IDs of the expression data at hand
 ind_organism  <-  sapply(FUN = grepl, 
                         X = c("ENSG", "ENSMUSG"), 
                         x = rownames(expression_data)[1])
           
 # (ii) specification of organism for KEGG
 # hsa: homo sapiens (human)
 # mmu: mus musculus (mouse)
 organism  <-  (c("hsa", "mmu"))[ind_organism][[1]]
 
 
 
 #############################
 #generate documentation frame 
 #############################
 
 if(metric == "rank"){
 
 doc <- data.frame(step = c("Default", "Pre-Filtering Threshold", "Duplicate Gene ID removal", "RNA-Seq Data Transformation Method"),
                          optimal_parameter = NA, 
                          rank = NA)
 
 }else if(metric == "p_adj"){
  
  doc <- data.frame(step = c("Default", "Pre-Filtering Threshold", "Duplicate Gene ID removal", "RNA-Seq Data Transformation Method"),
          optimal_parameter = NA, 
          p_adj = NA)
  
 }
 
 
 ###########################
 # transform phenotype labels to format required by function padog()
 ###########################
 phen_prep <- phenotype_prep(phenotype_labels)
 
 
 ##########
 # 1. step: Choose optimal pre-filtering threshold
 #########
 
 #default threshold: 10 read counts across all samples 
 #alternatives: pre-filtering performed using edgeR's builtin function filterByExpr 
 
 filt_threshold_alt <- 10
 
 exprdat_list_prefilt <- lapply(X = filt_threshold_alt, 
                                FUN = pre_filt, 
                                expression_data = expression_data)
 
 #further option: use built-in function from edgeR to determine which genes have sufficiently
 #large counts to be retained
 
 #pre-filtering using edgeR's builtin function filterByExprs()
 ind_genes1 <- filterByExpr(expression_data, 
                            group = phenotype_labels)
 #check whether the order of the genes is alignes 
 all(names(ind_genes1) == rownames(expression_data))
 
 
 #add pre-filtered expression data set to the list 
 exprdat_list_prefilt <- append(exprdat_list_prefilt, 
                                list(expression_data[ind_genes1, ]))
 
 
 PADOG_prefilt_list <- list()
 
 for(i in 1:length(exprdat_list_prefilt)){
  
  
  # (i) transform ENSEMBL gene IDs to Entrez gene IDs (choose default removal manner) AND
  # transform RNA-Seq data using voom 
  
  expression_data_transf <- geneID_conversion(exprdat_list_prefilt[[i]], dupl_removal_method = 1) %>% 
                                              voom_trans(phenotype_labels = phenotype_labels)
  
  # (ii) run PADOG
  seed <- 1007
  PADOG_prefilt_list[[i]] <- padog(expression_data_transf, 
                                   group = phen_prep, 
                                   dseed = seed, 
                                   organism = organism)
  
  # 4. step: perform multiple test adjustment using Benjamini and Hochberg 
  PADOG_prefilt_list[[i]]$p_adj <- p.adjust(PADOG_prefilt_list[[i]]$Ppadog, 
                                            method = "BH")
  
  
 }
 
 metric_vec_prefilt <- unlist(lapply(X = PADOG_prefilt_list, 
                                     FUN = pvalue_rank_padog, 
                                     term = geneset, 
                                     metric = metric))
 
 #get index of optimal pre-filtering threshold, i.e. threshold that maximizes the number
 #of differentially enriched gene sets 
 #-> in case of a tie: choose lower index 
 ind_opt_prefilt <- min(which(metric_vec_prefilt == min(metric_vec_prefilt), arr.ind = TRUE))
 
 #update documentation frame 
 #update documentation frame 
 doc[1, metric] <- metric_vec_prefilt[1] #default result
 #optimal result:
 ind_opt_prefilt <- min(which(metric_vec_prefilt == min(metric_vec_prefilt), arr.ind = TRUE)) #index of optimal pre-filtering manner
 # update documentation frame
 doc[2, "optimal_parameter"] <- c(filt_threshold_alt, 
                                  "filterByExpr", 
                                  "by cpm>1 in at least 2 samples")[ind_opt_prefilt] #optimal pre-filtering manner 
 doc[2, metric] <- min(metric_vec_prefilt)
 
 #set current optimal number of differentially enriched gene sets 
 metric_value <- min(metric_vec_prefilt)
 #set optimal PRE-FILTERED gene expression data set 
 exprdat_prefilt_opt <-  exprdat_list_prefilt[[ind_opt_prefilt]]
 
 
 
 ##########
 #2. step: Choose optimal manner of removing duplicated gene IDs caused by gene ID conversion 
      #from ENSEMBL to ENTREZ gene ID 
 ##########

 #perform 3 manners of gene ID conversion for the optimally pre-filtered gene expression data set 
 exprdat_convID <- lapply(FUN = geneID_conversion, 
                          expression_data = exprdat_prefilt_opt, 
                          X = 1:2)
 
 #generate list containing GSEA results for each expression data set
 PADOG_results_IDconv_list <- list()
 #perform default DESeq2 workflow for each of the expression data sets
 for(i in 1:length(exprdat_convID)){
  
  # (i) prepare input as required for clusterProfiler's ORA tool
  expression_data_transf <- voom_trans(exprdat_convID[[i]], 
                                       phenotype_labels, 
                                       "TMM")
  
  # (ii) run PADOG
  seed <- 1007
  PADOG_results_IDconv_list[[i]] <- padog(expression_data_transf, 
                                          group = phen_prep, 
                                          dseed = seed, 
                                          organism = organism)
  
  # 4. step: perform multiple test adjustment using Benjamini and Hochberg 
  PADOG_results_IDconv_list[[i]]$p_adj <- p.adjust(PADOG_results_IDconv_list[[i]]$Ppadog, 
                                                   method = "BH")
  
 }
 
 
 #calculate lossfunction for each result and store in vector
 metric_vec_convID <- unlist(lapply(X = PADOG_results_IDconv_list, 
                                    FUN = pvalue_rank_padog, 
                                    term = geneset, 
                                    metric = metric))
 
 #update documentation
 
 #index of optimal gene ID removal technique 
 ind_opt_conv_ID <- min(which(metric_vec_convID == min(metric_vec_convID), arr.ind = TRUE))
 doc[3, "optimal_parameter"] <- ind_opt_conv_ID
 doc[3, metric] <- min(metric_vec_convID)
 
 metric_value <- min(metric_vec_convID)
 
 #get optimal PRE-FILTERED gene expression data with converted gene IDs 
 expression_data_opt <- exprdat_convID[[ind_opt_conv_ID]]
 
 
 
 ##########
 # 3. step: Choose optimal RNA-Seq transformation method
 
 #alternative to voom transformation (default): DESeq2's varianceStabilizingTransformation
  seed <- 1007
  PADOG_vst <- variancetransform(expression_data_opt, phenotype_labels) %>% 
                                 padog(group = phenotype_prep(phenotype_labels), 
                                       dseed = seed, 
                                       organism = organism) %>%
                                 mutate(p_adj = p.adjust(Ppadog, method = "BH"))
  
  #update documentation frame
 doc[4, "optimal_parameter"] <- ifelse(pvalue_rank_padog(geneset,PADOG_vst, metric) < metric_value , 
                                       "varianceStabilizingTransformation", 
                                       "voom") 
 
 doc[4, metric] <- min(pvalue_rank_padog(geneset,PADOG_vst, metric), metric_value )
 
 #set optimal PADOG results which will then be returned by this optimization function 
 if(pvalue_rank_padog(geneset,PADOG_vst, metric) < metric_value){
  
  PADOG_opt <- PADOG_vst
  
 } else{PADOG_opt <- PADOG_results_IDconv_list[[ind_opt_conv_ID]]}
 
 
 #return optimal documentation frame 
 return(list(default = PADOG_prefilt_list[[1]], #default results
       optim = PADOG_opt, #optimal results
       documentation = doc)) # documentation frame 
 
 
 
}



################################################################################
### Run Optimization Functions #################################################
################################################################################

# phen_pickrell_list <- list()
# 
# for(i in 1:ncol(phen_pickrell)){
#  
#  phen_pickrell_list[[i]] <- phen_pickrell[,i]
#  
# }
# 
# phen_bottomly_list <- list()
# 
# for(i in 1:ncol(phen_bottomly)){
#  
#  phen_bottomly_list[[i]] <- phen_bottomly[,i]
#  
# }
# 
# # optimize the adjusted p-values 
# 
# 
# #############
# ### Pickrell 
# #############
# 
# ### (I) Gene Set Primary Immunodeficiency -> "05340"
# 
# # original phenotype assignment 
# optimP_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype <- 
#   PADOG_rankp_optim(geneset = "05340",
#                     Biobase::exprs(pickrell.eset), 
#                     pickrell.eset$gender, 
#                     metric = "p_adj")
# 
# # save results
# save(optimP_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype, 
#    file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Pickrell_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations <- 
#                       lapply(FUN = PADOG_rankp_optim, 
#                             geneset =  "05340", 
#                             expression_data = Biobase::exprs(pickrell.eset), 
#                             metric = "p_adj",
#                             X = phen_pickrell_list)
# 
# # save results
# save(optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations, 
#    file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Pickrell_PhenotypePermutations.RData")
# 
# 
# ### (II) Gene Set "Graft vs host disease" -> "05332"
# 
# # original phenotype assignment 
# optimP_PADOG_GraftvsHost_Pickrell_originalphenotype <- PADOG_rankp_optim(geneset = "05332",
#                                                                          Biobase::exprs(pickrell.eset), 
#                                                                          pickrell.eset$gender, metric = "p_adj")
# 
# # save results
# save(optimP_PADOG_GraftvsHost_Pickrell_originalphenotype, 
#    file = "./Results/optimP_PADOG_GraftvsHost_Pickrell_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations <- lapply(FUN = PADOG_rankp_optim, 
#                                                                  geneset =  "05332", 
#                                                                  expression_data = Biobase::exprs(pickrell.eset), 
#                                                                  metric = "p_adj",
#                                                                  X = phen_pickrell_list)
# 
# save(optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations, 
#    file = "./Results/optimP_PADOG_GraftvsHost_Pickrell_PhenotypePermutations.RData")
# 
# 
# 
# 
# #############
# ### Bottomly 
# #############
# 
# # We use the same gene sets as for the pickrell data set as the gene sets provided 
# # for human and mouse are almost identical 
# 
# ### (I) Gene Set Primary Immunodeficiency -> " 	05340"
# 
# # original phenotype assignment 
# optimP_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype <- 
#                         PADOG_rankp_optim(geneset = "05340",
#                                         Biobase::exprs(bottomly.eset), 
#                                         bottomly.eset$strain, 
#                                         metric = "p_adj")
# 
# # save results
# save(optimP_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype, 
#    file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimP_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations <- 
#                               lapply(FUN = PADOG_rankp_optim, 
#                                      geneset =  "05340", 
#                                      expression_data = Biobase::exprs(bottomly.eset), 
#                                      metric = "p_adj",
#                                      X = phen_bottomly_list)
# 
# # save results
# save(optimP_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations, 
#    file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_PhenotypePermutations.RData")
# 
# 
# ### (II) Gene Set "Graft vs host disease" -> "05332"
# 
# # original phenotype assignment 
# optimP_PADOG_GraftvsHost_Bottomly_originalphenotype <- 
#   PADOG_rankp_optim(geneset = "05332",
#                     Biobase::exprs(bottomly.eset), 
#                     bottomly.eset$strain, 
#                     metric = "p_adj")
# 
# # save results
# save(optimP_PADOG_GraftvsHost_Bottomly_originalphenotype, 
#    file = "./Results/optimP_PADOG_GraftvsHost_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimP_PADOG_GraftvsHost_Bottomly_phenotypepermutations <- 
#                         lapply(FUN = PADOG_rankp_optim, 
#                                geneset =  "05332", 
#                                expression_data = Biobase::exprs(bottomly.eset), 
#                                metric = "p_adj",
#                                X = phen_bottomly_list)
# 
# save(optimP_PADOG_GraftvsHost_Bottomly_phenotypepermutations, 
#    file = "./Results/optimP_PADOG_GraftvsHost_Bottomly_PhenotypePermutations.RData")
# 
# 
# 
# # optimize the ranks 
# 
# #############
# ### Pickrell 
# #############
# 
# ### (I) Gene Set Primary Immunodeficiency -> "05340"
# 
# # original phenotype assignment 
# optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype <- 
#   PADOG_rankp_optim(geneset = "05340",
#                     Biobase::exprs(pickrell.eset), 
#                     pickrell.eset$gender, 
#                     metric = "rank")
# 
# # save results
# save(optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype, 
#      file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations <- 
#   lapply(FUN = PADOG_rankp_optim, 
#          geneset =  "05340", 
#          expression_data = Biobase::exprs(pickrell.eset), 
#          metric = "rank",
#          X = phen_pickrell_list)
# 
# # save results
# save(optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations, 
#      file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_PhenotypePermutations.RData")
# 
# 
# ### (II) Gene Set "Graft vs host disease" -> "05332"
# 
# # original phenotype assignment 
# optimRank_PADOG_GraftvsHost_Pickrell_originalphenotype <- PADOG_rankp_optim(geneset = "05332",
#                                                                          Biobase::exprs(pickrell.eset), 
#                                                                          pickrell.eset$gender, metric = "rank")
# 
# # save results
# save(optimRank_PADOG_GraftvsHost_Pickrell_originalphenotype, 
#      file = "./Results/optimRank_PADOG_GraftvsHost_Pickrell_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimRank_PADOG_GraftvsHost_Pickrell_phenotypepermutations <- lapply(FUN = PADOG_rankp_optim, 
#                                                                   geneset =  "05332", 
#                                                                   expression_data = Biobase::exprs(pickrell.eset), 
#                                                                   metric = "rank",
#                                                                   X = phen_pickrell_list)
# 
# save(optimRank_PADOG_GraftvsHost_Pickrell_phenotypepermutations, 
#      file = "./Results/optimRank_PADOG_GraftvsHost_Pickrell_PhenotypePermutations.RData")
# 
# 
# 
# 
# #############
# ### Bottomly 
# #############
# 
# # We use the same gene sets as for the pickrell data set as the gene sets provided 
# # for human and mouse are almost identical 
# 
# ### (I) Gene Set Primary Immunodeficiency -> " 	05340"
# 
# # original phenotype assignment 
# optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype <- 
#   PADOG_rankp_optim(geneset = "05340",
#                     Biobase::exprs(bottomly.eset), 
#                     bottomly.eset$strain, 
#                     metric = "rank")
# 
# # save results
# save(optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype, 
#      file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations <- 
#   lapply(FUN = PADOG_rankp_optim, 
#          geneset =  "05340", 
#          expression_data = Biobase::exprs(bottomly.eset), 
#          metric = "rank",
#          X = phen_bottomly_list)
# 
# # save results
# save(optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations, 
#      file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_PhenotypePermutations.RData")
# 
# 
# ### (II) Gene Set "Graft vs host disease" -> "05332"
# 
# # original phenotype assignment 
# optimRank_PADOG_GraftvsHost_Bottomly_originalphenotype <- 
#   PADOG_rankp_optim(geneset = "05332",
#                     Biobase::exprs(bottomly.eset), 
#                     bottomly.eset$strain, 
#                     metric = "rank")
# 
# # save results
# save(optimRank_PADOG_GraftvsHost_Bottomly_originalphenotype, 
#      file = "./Results/optimRank_PADOG_GraftvsHost_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimRank_PADOG_GraftvsHost_Bottomly_phenotypepermutations <- 
#   lapply(FUN = PADOG_rankp_optim, 
#          geneset =  "05332", 
#          expression_data = Biobase::exprs(bottomly.eset), 
#          metric = "rank",
#          X = phen_bottomly_list)
# 
# save(optimRank_PADOG_GraftvsHost_Bottomly_phenotypepermutations, 
#      file = "./Results/optimRank_PADOG_GraftvsHost_Bottomly_PhenotypePermutations.RData")
# 
# 
# 









