#### clusterProfiler's ORA tool: minimization of a given gene set's 
#(i) (adjusted) p-value / 
#(ii) rank among the remaining gene sets 

library(clusterProfiler)
library(dplyr)
library(DESeq2)
library(edgeR)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


#load phenotype permutations and gene expression data sets 
source("./Random_Phenotype_Permutations.R")

# load the preprocessing functions
source("./PreProcessing_Functions.R")



################################################################################
### get gene's adjusted p-value or rank (based on what is needed) ##############
################################################################################

############### NOTE: 
# encountered problem in the optimization process of the rank or (adjusted) pvalue of a given gene set:
#clusterProfiler only reports those gene sets which contain at least a single gene from the input list
#of differentially expressed genes. 

# this specifically means that, if a gene set does not contain any genes from the input list, function 
#pvalue_rank_ORA returns integer(0) which cannot be compared to numerical number. 
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
# -> the rank is then relative to the total number of gene sets in the GSA results such that a 
# rank of 1 indicates that the gene set has the highest adjusted p-value among the remaining gene sets


#note: function works for GO as well as KEGG as results table "ora_results" 
#contains relevant information automatically
#argument term must be in the form of a GO ID resp. KEGG ID 

# note that clusterProfiler's ORA only includes those gene sets from the results that 
# have at least one gene from the input as a member. 
# This means that when a gene set does not appear in the GSA results, we do not know
# whether it is not differentially enriched or it does not appear in the gene set database at all

pvalue_rank_ORA <- function(term, ora_results, metric){
 
if(metric == "rank"){
 # if the ora_result is NULL, this means that no gene set is reported as differentially enriched
 # note that ORA only reports those gene sets from the geneset database that have 
 # at least one differentially enriched gene sets from the gene set database 
 # We therefore cannot know whether the gene set was contained by the gene set database or not
 # we then set the rank of the given gene set to the worst possible value, i.e. 1.2
  
  # there are also cases when none of the differentially expressed genes (i.e. input)
  # can be matched to a gene set (in this case, differential enrichment cannot be 
  # precluded either) -> set rank to 1.2 
 
 if(is.null(ora_results)){
  
 return(1.2)
  
 } else if(nrow(ora_results) == 0){


  return(1.2)

 } else if(nrow(ora_results) != 0){
  # the ORA results are non-empty, i.e. there is at least one differentially 
   # enriched gene set
   
   
   # prepare the relative ranks such that all gene sets with the same adjusted p-value
   # are given the same (i.e. average) rank
   # We divide by the maximum rank such that those gene sets with the highest 
   # adjusted p-values are given the (relative) rank 1, which is the worst rank
   
   ora_results$ranks <- rank(ora_results$p.adjust) / max(rank(ora_results$p.adjust)) 
   
   # get the relative rank of the gene set of interest 
   rank <- ora_results$ranks[grep(term, ora_results$ID)]
   
   
  # if the gene set is not contained in the results, return the rank 1.2
  # this shall serve as an indicator in the results that the gene set 
  # was not contained in the results 

 #return row number of respective gene set
 return(ifelse(!is.integer0(grep(term, ora_results$ID)), 
               rank, 
               1.2))
 # note: in the case that a gene set is not reported in the results table of ora_results,
 # ifelse() in combination with !is.integer0() then ensures that a rank of 1.2 is returned,
 # meaning that each adaption resulting in the gene set appearing in the results leads
 # to an improvement of the results
  
 } 
 
}else if(metric == "p_adj"){
 
 if(is.null(ora_results)){
  
  return(1.2)
  
 } else if(nrow(ora_results) == 0){
  
  
  return(1.2)
 }
 
 if(nrow(ora_results) != 0){
  
 #for doublecheck reasons: order results by adjusted p-value in ascencing manner 
 ora_results <- ora_results[order(ora_results$p.adjust, decreasing = FALSE),]
 
 #identify row number of respective gene set 
 ind_row <- grep(term, ora_results$ID)
 
 #return respective adjusted p-value
 return(ifelse(!is.integer0(ind_row), 
               ora_results$p.adjust[ind_row], 
               1.2))
 #note: in the case that a gene set is not reported in the results table of ora_results, 
 #ifelse() in combination with !is.integer0() then ensures that an adjusted p-value of 1.2 is returned,
 #meaning that each adaption leading to a an adjusted p-value in (0,1] is considered an improvement
 
} # if no gene sets are reported in the ORA results (which is the case when the input
 # list of differentially expressed genes is empty) return the rank Inf for the given gene set 
}
 
}



#######################################################################################
###generate necessary input for ORA tool ##############################################
#######################################################################################

cP_ORA_input_preparation <- function(DE_results){
 
 #required input for clusterProfiler function: vector of entrez gene ID
 #-> need to pre-process results table DE_results
 
 #vector of differentially expressed genes 
 #DEG_vec serves as input vector for ORA performed by clusterProfiler 
 DEG_vec <- rownames(DE_results[(DE_results$p_adj<0.05) & (!is.na(DE_results$p_adj)),])
 
 return(DEG_vec)
 
}

################################################################################
### default ORA pipelines ######################################################
################################################################################
 
#(i) for geneset database KEGG 

ORA_pipeline_default <- function(DE_results, geneset_database, organism){
 
 
 # note: set pvalueCutoff and qvalueCutoff to 1 to ensure that all gene sets are displayed in the results 
 
 ############
 ### option 1: geneset database KEGG 
 ###########
 
 if(geneset_database == "KEGG"){
  
 
 ORA_results <- cP_ORA_input_preparation(DE_results) %>%
                enrichKEGG(organism = organism, 
                          keyType = "kegg", 
                          pvalueCutoff = 1, 
                          qvalueCutoff = 1) %>% 
                as.data.frame()
 
  ############
  ### option 2: geneset database GO
  ##########

 } else if(geneset_database == "GO"){# default ORA pipeline for geneset database GO
  

   
   ORA_results <- cP_ORA_input_preparation(DE_results) %>%
                  enrichGO(OrgDb = organism, 
                          ont = "BP", 
                          pvalueCutoff = 1, 
                          qvalueCutoff = 1) %>% 
                  as.data.frame()
   

 }
 
 return(ORA_results)
}






################################################################################
### optimization for KEGG or GO ################################################
################################################################################

#minimize either rank or adjusted p-value of the gene set of choice 
#for all runs of enrichGO, we set p-value cutoff and q-value cutoff values to 1
#to ensure that ALL gene sets (whether significantly differentially enriched or not)
#are reported

# optimization function can be run for KEGG as well as GO (gene set database of choice
# can be adjusted using argument geneset_database)
# -> format of gene set geneset must be adjusted accordingly since GO and KEGG offer
# different gene sets 

ORA_rank_pvalue_optim <- function(geneset, metric, expression_data, phenotype_labels, geneset_database){

 
 # identify organism for which gene expression is measured from the format of the gene IDs (all are Ensembl)
 # ENSEMBL gene ID ENSMUSGXXXXXXXXXXX corresponds to mouse; identifiable by substring "ENSMUSG"
 # ENSEMBL gene ID ENSGXXXXXXXXXXX corresponds to homo sapiens; identifiable by string "ENSG"
 
 
 # for GO and KEGG, the correct organism must be identified, however, they must be 
 # specified differently for KEGG and GO 
 
 
 if(geneset_database == "GO"){
  
  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl, 
                         X = c("ENSG", "ENSMUSG"), 
                         x = rownames(expression_data)[1])
  
  # (i) specification of organism for GO
  organism <- c(org.Hs.eg.db , org.Mm.eg.db)[ind_organism][[1]]
  
 } else if(geneset_database == "KEGG"){
  
  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl, 
                         X = c("ENSG", "ENSMUSG"), 
                         x = rownames(expression_data)[1])
  
  # (ii) specification of organism for KEGG
  # hsa: homo sapiens (human)
  # mmu: mus musculus (mouse)
  organism <- (c("hsa", "mmu"))[ind_organism][[1]]
  
  
 }
 
 
 ###### generate documentation data frame to document development 
 #we want the column names of the documentation frame to depend on the 
 #metric value to optimize (i.e. rank vs. adjusted p-value)
 if(metric == "rank"){
 doc <- data.frame(step = c("Default", "Differential Expression Technique", "Pre-filtering Threshold", "Duplicate Gene ID Removal"),
         optimal_parameter = NA, 
         rank = NA)

 
 }else if(metric == "p_adj"){
  doc <- data.frame(step = c("Default", "Differential Expression Technique", "Pre-filtering Threshold", "Duplicate Gene ID Removal"),
          optimal_parameter = NA, 
          p_adj = NA)
  

  
 }
 

 
 
 ##############
 #### 1.step: choose optimal DE technique
 #############
 
 de_techniques <- c("DESeq2", "limma")
 metricvalue_detechniques <- c() #to store adjusted p-value / rank of each DE technique 
 
 
 ### default DE technique: DESeq2 
 
 # (i) run default DESeq2 workflow (including default pre-filtering using threshold of 10 and default removal of duplicated gene IDs)
 deseq2_results <- pre_filt(expression_data, threshold = 10) %>% 
                          geneID_conversion(dupl_removal_method = 1) %>%
                          deseq_preprocess(phenotype_labels) %>% 
                          DESeq() %>% results() %>% 
                          as.data.frame() %>% dplyr::rename( p_adj = padj)
 
 # (ii) run ORA 
  ORA_results_deseq2 <- ORA_pipeline_default(deseq2_results, 
                                            geneset_database = geneset_database, 
                                            organism = organism)
 
 # (iii) store metric value of geneset 
 metricvalue_detechniques[1] <- pvalue_rank_ORA(geneset, 
                                            ORA_results_deseq2, 
                                            metric)
 
 # update documentation frame
 doc[1, metric] <- metricvalue_detechniques[1]
 
 
  
 ### alternative 2: limma 
   
   
  # (i.1) generate pre-filtering indicator using filterByExpr()
   keep <- DGEList(expression_data, 
                   group = phenotype_labels) %>% 
          filterByExpr() 
   
   
  # (i.2) generate design matrix 
   mm <- model.matrix(~ phenotype_labels)
   
  # (i.3) run default limma pipeline (reuse filtering indicator keep from edgeR)
   limma_results <- geneID_conversion(expression_data[keep,], dupl_removal_method = 1) %>%
                                      DGEList(group = phenotype_labels) %>% 
                                      calcNormFactors() %>% voom(design = mm) %>% 
                                      lmFit(design = mm) %>% eBayes() %>% 
                                      topTable(coef = ncol(mm), number = 100000) %>%
                                      as.data.frame() %>% dplyr::rename(p_adj = adj.P.Val)
   
  # (ii) run ORA
   ORA_results_limma <- ORA_pipeline_default(limma_results, 
                        geneset_database = geneset_database, 
                        organism)
   
  # (iii) store metricvalue of geneset
   metricvalue_detechniques[2] <- pvalue_rank_ORA(geneset, 
                         ORA_results_limma, 
                         metric)
   
   
  
 ####################  
   
 # get index of optimal DE technique
 ind_opt_DE_technique <- min(which( metricvalue_detechniques == min( metricvalue_detechniques)))
 
 # get optimal DE technique
 opt_DE_technique <- de_techniques[ind_opt_DE_technique]
 
 
 #update documentation frame
 doc[2, "optimal_parameter"] <- opt_DE_technique
 doc[2, metric] <- min(metricvalue_detechniques)
 
 #for the specific gene sets tested here, edgeR never lead to a decreased RANK of the given gene set among 
 #the remaining. This can be substantiated by the fact that edgeR does not deal with outlier values, resulting
 #in genes containing outlier values having small (adjusted) pvalues and superseding other genes in the lowest
 #RANKS
 


 
 ##############
 ##### 2. step: Choose optimal pre-filtering threshold 
 #############
 
 # depending on the optimal DE technique established in the previous optimization step, different alternatives of pre-filtering are 
 # applied to the gene expression data set 
 
 if(opt_DE_technique == "DESeq2"){# for DESeq2, different manual pre-filtering thresholds are applied the the gene expression data set 
  
  # different pre-filtering threshold 
  filt_methods <- c(10, 50)
  
  # differently pre-filtered gene expression data sets 
  exprdat_list_prefilt <- lapply(X = filt_methods, 
                                FUN = pre_filt, 
                                expression_data = expression_data)
  
  # run DESeq2 default pipeline for each of the pre-filtered gene expression data sets 
  DESeq2_results_prefilt <- lapply(X = exprdat_list_prefilt, FUN = geneID_conversion, dupl_removal_method = 1) %>% 
                            lapply(FUN = deseq_preprocess, phenotype_labels = phenotype_labels) %>% 
                            lapply(FUN = DESeq) %>% lapply(FUN = results) %>% 
                            lapply(FUN = as.data.frame) %>% lapply(FUN = dplyr::rename, p_adj = padj)
  
  # run default ORA pipeline for each of the DE results created with DESeq2
  ORA_results_prefilt <- lapply(X = DESeq2_results_prefilt, 
                                FUN = ORA_pipeline_default, 
                                geneset_database = geneset_database, 
                                organism = organism)
  
  
  
 } else if(opt_DE_technique  == "limma"){
  # for edgeR and limma, pre-filtering is performed differently, namely using filterByExpr (default) OR using cpm()
  
  # DEFAULT pre-filtering using filterByExpr
  keep1 <- DGEList(counts = expression_data, group = phenotype_labels) %>% 
       filterByExpr()
  # 1. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 1 count per million in at least 2 samples
  keep2 <- rowSums(cpm(expression_data) > 1) >= 2 

  
  filt_methods <- c("by FilterbyExpr()", "by cpm>1 in at least 2 samples")
  
  
  # for all methods of pre-filtering: store resulting gene expression data sets in list 
  exprdat_list_prefilt <- list()
  # pre-filtered gene expression data set using keep1 
  exprdat_list_prefilt[[1]] <- expression_data[keep1,]
  # pre-filtered gene expression data set using keep2
  exprdat_list_prefilt[[2]] <- expression_data[keep2,]
 
  

   
   # run limma pipeline for each of the pre-filtered gene expression data sets 
   limma_results_prefilt <- lapply(FUN = geneID_conversion, X = exprdat_list_prefilt, dupl_removal_method = 1) %>%
                            lapply(FUN = DGEList,group = phenotype_labels) %>% 
                            lapply(FUN = calcNormFactors) %>% lapply(FUN = voom, design = mm) %>% 
                            lapply(FUN = lmFit, design = mm) %>% lapply(FUN = eBayes) %>% 
                            lapply(FUN = topTable,coef = ncol(mm), number = 100000) %>%
                            lapply(FUN = as.data.frame) %>% lapply(FUN = dplyr::rename, p_adj = adj.P.Val)
   
   # run ORA for each of the DE results created with limma
   ORA_results_prefilt <- lapply(FUN = ORA_pipeline_default, 
                  X = limma_results_prefilt, 
                  geneset_database = geneset_database, 
                  organism = organism)
   
 
 } 
 
 # get metric value for each of the pre-filtering thresholds
 metricvalue_prefilt <- unlist(lapply(FUN = pvalue_rank_ORA, 
                                      X = ORA_results_prefilt, 
                                      term = geneset, 
                                      metric = metric))
 
 # get index of optimal pre-filtering threshold (choose lower index in case of ties)
 ind_opt_prefilt <- min(which(metricvalue_prefilt == min(metricvalue_prefilt)))
 
 # get optimal pre-filtering threshold
 opt_prefilt <- filt_methods[ind_opt_prefilt]
 
 # update documentation frame
 doc[3, "optimal_parameter"] <- opt_prefilt
 doc[3, metric] <- min(metricvalue_prefilt)
 

  
  # get optimally pre-filtered gene expression data set 
  exprdat_prefilt_opt <- exprdat_list_prefilt[[ind_opt_prefilt]] 
  

 


 
 
 ###########
 # 3. step: Choose optimal manner of duplicate gene ID removal
 ###########
 
 if(opt_DE_technique == "DESeq2"){
  
  # run DESeq2 pipeline for each of the duplicate gene ID removal schemes 
  DE_results_convID <- lapply(FUN = geneID_conversion, X = 1:2, expression_data = exprdat_prefilt_opt) %>% 
                       lapply(FUN = deseq_preprocess, phenotype_labels = phenotype_labels) %>% 
                       lapply(FUN = DESeq) %>% lapply(FUN = results) %>% 
                       lapply(FUN = as.data.frame) %>% lapply(FUN = dplyr::rename, p_adj = padj)
  
  # run ORA pipeline for each of the DE results created with DESeq2
  ORA_results_convID <- lapply(FUN = ORA_pipeline_default, 
                              X = DE_results_convID, 
                              geneset_database = geneset_database, 
                              organism = organism) 
             
  
  
 }else if(opt_DE_technique == "limma"){
  
  # run limma pipeline for each of the duplicate gene id removal schemes 
  DE_results_convID <- lapply(FUN = geneID_conversion, expression_data = exprdat_prefilt_opt, X = 1:2) %>%
                       lapply(FUN = DGEList,group = phenotype_labels) %>% 
                       lapply(FUN = calcNormFactors) %>% lapply(FUN = voom,design = mm) %>% 
                       lapply(FUN = lmFit,design = mm) %>% lapply(FUN = eBayes) %>% 
                       lapply(FUN = topTable,coef = ncol(mm), number = 100000) %>%
                       lapply(FUN = as.data.frame) %>% lapply(FUN = dplyr::rename,p_adj = adj.P.Val)
  
  # run ORA for each of the DE results created with limma 
  ORA_results_convID <- lapply(FUN = ORA_pipeline_default, 
                              X = DE_results_convID, 
                              geneset_database = geneset_database, 
                              organism = organism) 

 } 
 
 # get metric value for each of the ORA results 
 metricvalue_convID <- unlist(lapply(FUN = pvalue_rank_ORA, 
                                    X = ORA_results_convID, 
                                    term = geneset, 
                                    metric = metric))
 
 # index of optimal duplicate gene ID removal scheme (in case of ties choose lower index)
 ind_opt_convID <- min(which(metricvalue_convID == min(metricvalue_convID)))
 
 # update documentation frame
 doc[4, "optimal_parameter"] <- ind_opt_convID
 doc[4, metric] <- min(metricvalue_convID) 
 
 
 # current optimal metric value
 metricvalue <- min(metricvalue_convID)
 

 

 
 # get optimally pre-filtered gene expression data set with optimally converted gene IDs 
 #exprdat_prefilt_convID_opt <- geneID_conversion(expression_data = exprdat_prefilt_opt, dupl_removal_method = ind_opt_convID)
 
 # get optimal results from differential expression analysis
 # -> corresponds to DE results resulting from the optimal manner of duplicate gene ID removal 
 opt_DE_results <- DE_results_convID[[ind_opt_convID]]

 
 
 #############
 ### 4. step: optimize internal parameter set of KEGG 
 ############
 #-> only parameter to be optimized within KEGG: universe (i.e. set of background genes)
 #-> all genes annotated to KEGG (default) vs. all genes measured in experiment (i.e. all 
 #genes with a non-NA adjusted p-value present in the optimal DE result from cP_ORA_optim)
 
 #alternative unviverse: extract all genes from the optimal DE results that have an adjusted p-value (i.e. whose gene expression
 # was actually tested )
 # extract universe from optimal limma results 
 univ_alt <- rownames(opt_DE_results[!is.na(opt_DE_results$p_adj), ])
 
 # now run ORA based on optimal results stored in cP_ORA_optim and the alternative universe
 # -> different functions for gene set databases KEGG and GO 
 
 if(geneset_database == "KEGG"){
  
  ORA_univ <- enrichKEGG(gene = cP_ORA_input_preparation(DE_results_convID[[ind_opt_convID]]), 
                        organism = organism, keyType = "kegg", 
                        universe = univ_alt , pvalueCutoff = 1, qvalueCutoff = 1) %>%
              as.data.frame()
  
  
 }else if(geneset_database == "GO"){
  
  ORA_univ <- enrichGO(gene = cP_ORA_input_preparation(DE_results_convID[[ind_opt_convID]]), 
                       OrgDb = organism,ont = "BP", universe = univ_alt, 
                       pvalueCutoff = 1, qvalueCutoff = 1) %>%
              as.data.frame()
  
  
 }
 
 
 
 

 #run final optimal ORA result depending on if the alternative universe leads to
 #an im provement of the metric value 
 
 if(pvalue_rank_ORA(geneset, ORA_univ, metric) < metricvalue){
  
  #set final optimal ORA run as the one with adapted universe 
  ORA_final<- ORA_univ
  #update documentation frame 
  doc[nrow(doc) + 1, "step"] <- "Universe"
  doc[nrow(doc), "optimal_parameter"] <- "Adapted"
  doc[nrow(doc), metric] <- pvalue_rank_ORA(geneset, ORA_univ, metric)
  
  
 }else if(pvalue_rank_ORA(geneset, ORA_univ, metric) >= metricvalue){
  
  #set final optimal ORA run as the one with the original universe 
  ORA_final <- ORA_results_convID[[ind_opt_convID]]
  #update documentation frame 
  doc[nrow(doc) + 1, "step"] <- "Universe"
  doc[nrow(doc), "optimal_parameter"] <- "Default"
  doc[nrow(doc), metric] <- metricvalue
  
}


 
 return(list(default_ORA = ORA_results_deseq2, #default ORA result
     optim_ORA = ORA_final, #optimal DE result
     documentation = doc)) #documentation frame 
 
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
# 
# ###### optimization (i.e. minimization) of adjusted p-values 
# 
# 
# #############
# ### Pickrell 
# #############
# 
# 
# 
# ############################################################
# ### (I) Gene Set "t-cell mediated immunity" -> GO:0002456
# ############################################################
# 
# # original phenotype assignment 
# optimP_cP_ORA_tcell_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0002456",
#                                  metric = "p_adj", 
#                                  Biobase::exprs(pickrell.eset), 
#                                  pickrell.eset$gender, 
#                                  "GO")
# 
# # save results
# save(optimP_cP_ORA_tcell_Pickrell_originalphenotype, 
#    file = "./Results/optimP_ORA_tCell_Pickrell_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimP_cP_ORA_tcell_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                             geneset =  "GO:0002456", 
#                             metric = "p_adj",
#                             expression_data = Biobase::exprs(pickrell.eset), 
#                             geneset_database = "GO",
#                             X = phen_pickrell_list)
# 
# # save results
# save(optimP_cP_ORA_tcell_Pickrell_phenotypepermutations, 
#    file = "./Results/optimP_ORA_tCell_Pickrell_PhenotypePermutations.RData")
# 
# 
# ############################################################
# ### (II) Gene Set "Demethylation" -> GO:0070988
# ############################################################
# 
# # Note: in other GSA tools, demethylation might also be identified by GO:0080111
# 
# # original phenotype assignment 
# optimP_cP_ORA_Demethylation_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0070988",
#                                      metric = "p_adj", 
#                                      Biobase::exprs(pickrell.eset), 
#                                      pickrell.eset$gender, 
#                                      "GO")
# 
# # save results
# save(optimP_cP_ORA_Demethylation_Pickrell_originalphenotype, 
#    file = "./Results/optimP_ORA_Demethylation_Pickrell_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                                geneset =  "GO:0070988", 
#                                metric = "p_adj",
#                                expression_data = Biobase::exprs(pickrell.eset), 
#                                geneset_database = "GO",
#                                X = phen_pickrell_list)
# 
# # save results
# save(optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations, 
#    file = "./Results/optimP_ORA_Demethylation_Pickrell_PhenotypePermutations.RData")
# 
# 
# 
# #############
# ### Bottomly 
# #############
# 
# ### (I) Gene set "Metabolic Process: GO:0008152"
# 
# # original phenotype assignment 
# optimP_cP_ORA_MetabolicProcess_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0008152","p_adj", Biobase::exprs(bottomly.eset), bottomly.eset$strain, "GO")
# 
# # save results
# save(optimP_cP_ORA_MetabolicProcess_Bottomly_originalphenotype, 
#    file = "./Results/optimP_ORA_MetabolicProcess_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimP_cP_ORA_MetabolicProcess_Bottomly_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                                  geneset = "GO:0008152", 
#                                  metric = "p_adj",
#                                  expression_data = Biobase::exprs(bottomly.eset), 
#                                  geneset_database = "GO",
#                                  X = phen_bottomly_list)
# 
# # save results 
# save(optimP_cP_ORA_MetabolicProcess_Bottomly_phenotypepermutations, 
#    file = "./Results/optimP_ORA_MetabolicProcess_Bottomly_PhenotypePermutations.RData")
# 
# 
# 
# 
# ### (II) Gene set "Cellular Process: GO:0009987"
# 
# # original phenotype assignment 
# optimP_cP_ORA_CellularProcess_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0009987",metric = "p_adj", Biobase::exprs(bottomly.eset), bottomly.eset$strain, "GO")
# 
# # save results
# save(optimP_cP_ORA_CellularProcess_Bottomly_originalphenotype, 
#    file = "./Results/optimP_ORA_CellularProcess_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimP_cP_ORA_CellularProcess_Bottomly_Phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                                                                       geneset = "GO:0009987", 
#                                                                       metric = "p_adj",
#                                                                       expression_data = Biobase::exprs(bottomly.eset), 
#                                                                       geneset_database = "GO",
#                                                                       X = phen_bottomly_list)
# # save results 
# save(optimP_cP_ORA_CellularProcess_Bottomly_Phenotypepermutations, 
#    file = "./Results/optimP_ORA_CellularProcess_Bottomly_PhenotypePermutations.RData")
# 
# 
# # optimization (i.e. minimization of the ranks)
# 
# 
# 
# #############
# ### Pickrell 
# #############
# 
# 
# 
# ############################################################
# ### (I) Gene Set "t-cell mediated immunity" -> GO:0002456
# ############################################################
# 
# # original phenotype assignment 
# optimRank_cP_ORA_tcell_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0002456",
#                                                                         metric = "rank", 
#                                                                         Biobase::exprs(pickrell.eset), 
#                                                                         pickrell.eset$gender, 
#                                                                         "GO")
# 
# # save results
# save(optimRank_cP_ORA_tcell_Pickrell_originalphenotype, 
#      file = "./Results/optimRank_ORA_tCell_Pickrell_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimRank_cP_ORA_tcell_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                                                              geneset =  "GO:0002456", 
#                                                              metric = "rank",
#                                                              expression_data = Biobase::exprs(pickrell.eset), 
#                                                              geneset_database = "GO",
#                                                              X = phen_pickrell_list)
# 
# # save results
# save(optimRank_cP_ORA_tcell_Pickrell_phenotypepermutations, 
#      file = "./Results/optimRank_ORA_tCell_Pickrell_PhenotypePermutations.RData")
# 
# 
# ############################################################
# ### (II) Gene Set "Demethylation" -> GO:0070988
# ############################################################
# 
# # Note: in other GSA tools, demethylation might also be identified by GO:0080111
# 
# # original phenotype assignment 
# optimRank_cP_ORA_Demethylation_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0070988",
#                                                                                 metric = "rank", 
#                                                                                 Biobase::exprs(pickrell.eset), 
#                                                                                 pickrell.eset$gender, 
#                                                                                 "GO")
# 
# # save results
# save(optimRank_cP_ORA_Demethylation_Pickrell_originalphenotype, 
#      file = "./Results/optimRank_ORA_Demethylation_Pickrell_OriginalPhenotype.RData")
# 
# # 10 random phenotype permutations
# optimRank_cP_ORA_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                                                                      geneset =  "GO:0070988", 
#                                                                      metric = "rank",
#                                                                      expression_data = Biobase::exprs(pickrell.eset), 
#                                                                      geneset_database = "GO",
#                                                                      X = phen_pickrell_list)
# 
# # save results
# save(optimRank_cP_ORA_Demethylation_Pickrell_phenotypepermutations, 
#      file = "./Results/optimRank_ORA_Demethylation_Pickrell_PhenotypePermutations.RData")
# 
# 
# 
# #############
# ### Bottomly 
# #############
# 
# ### (I) Gene set "Metabolic Process: GO:0008152"
# 
# # original phenotype assignment 
# optimRank_cP_ORA_MetabolicProcess_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0008152","rank", Biobase::exprs(bottomly.eset), bottomly.eset$strain, "GO")
# 
# # save results
# save(optimRank_cP_ORA_MetabolicProcess_Bottomly_originalphenotype, 
#      file = "./Results/optimRank_ORA_MetabolicProcess_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimRank_cP_ORA_MetabolicProcess_Bottomly_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                                                                         geneset = "GO:0008152", 
#                                                                         metric = "rank",
#                                                                         expression_data = Biobase::exprs(bottomly.eset), 
#                                                                         geneset_database = "GO",
#                                                                         X = phen_bottomly_list)
# 
# # save results 
# save(optimRank_cP_ORA_MetabolicProcess_Bottomly_phenotypepermutations, 
#      file = "./Results/optimRank_ORA_MetabolicProcess_Bottomly_PhenotypePermutations.RData")
# 
# 
# 
# 
# ### (II) Gene set "Cellular Process: GO:0009987"
# 
# # original phenotype assignment 
# optimRank_cP_ORA_CellularProcess_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0009987",metric = "rank", Biobase::exprs(bottomly.eset), bottomly.eset$strain, "GO")
# 
# # save results
# save(optimRank_cP_ORA_CellularProcess_Bottomly_originalphenotype, 
#      file = "./Results/optimRank_ORA_CellularProcess_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 random phenotype permutations
# optimRank_cP_ORA_CellularProcess_Bottomly_Phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
#                                                                        geneset = "GO:0009987", 
#                                                                        metric = "rank",
#                                                                        expression_data = Biobase::exprs(bottomly.eset), 
#                                                                        geneset_database = "GO",
#                                                                        X = phen_bottomly_list)
# # save results 
# save(optimRank_cP_ORA_CellularProcess_Bottomly_Phenotypepermutations, 
#      file = "./Results/optimRank_ORA_CellularProcess_Bottomly_PhenotypePermutations.RData")
# 
# 
# 
# 
# 
