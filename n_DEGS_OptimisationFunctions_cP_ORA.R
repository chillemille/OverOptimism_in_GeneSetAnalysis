# Investigation of Potential for over-optimistic results for clusterProfiler's ORA tool

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(dplyr)
library(edgeR)

# get gene expression datasets, phenotype assignments and random phenotype permutations
source("./Random_Phenotype_Permutations.R")

# load required pre-processing functions 
source("./PreProcessing_Functions.R")




#introduction of  necessary functions to perform optimization of GSA results 

######################################################################################
###lossfunction for ORA tool #########################################################
######################################################################################

#generate lossfunction to count the number of differentially enriched gene sets 
lossfunction_cPORA <- function(gsa_results){
  
  #number of Differentially Enriched Gene Sets 
  n_DEGS <- sum(as.data.frame(gsa_results)$p.adjust < 0.05)
  
  return(n_DEGS)
  
}


#############################################
###generate necessary input for ORA tool  ###
#############################################

cP_ORA_input_preparation <- function(DE_results){
  
  # required input for clusterProfiler function: vector of entrez gene ID
  # that contains all differentially expressed genes 
  #-> need to pre-process results table DE_results 
  
 
  # classify those genes as DE that have an EXISTING adjusted p-value < 0.05 in the results
  # table of differential expression analysis 
  DEG_vec <- rownames(DE_results[(DE_results$p_adj < 0.05) & (!is.na(DE_results$p_adj)), ])
  
  # return vector of differentially expressed genes 
  return(DEG_vec)
  
}


#######################################################################################
### (I) optimization of preparatory parameters ########################################
#######################################################################################
#optimization of input list of differentially expressed genes to maximize number of 
#differentially enriched gene sets 



#######################################################################################
### optimization of pre-processing ####################################################
#######################################################################################

cP_ORA_preprocess_optim <- function(expression_data, phenotype_labels){
  
  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl, X = c("ENSG", "ENSMUSG"), 
                         x = rownames(expression_data)[1])
  # choose suitable organism (required for function bitr)
  organism <- c(org.Hs.eg.db , org.Mm.eg.db)[ind_organism][[1]]
  
  
  
  ###### generate documentation data frame to document development of number of differentially enriched gene sets
  doc <- data.frame(step=c("Default", "Differential Expression Technique", "Pre-filtering Threshold", "Duplicate Gene ID Removal Technique"),
                    optimal_parameter=NA, 
                    n_DEGS=NA)
  
  
  
  ##############
  ##### 1. step: choose optimal Differential Expression analysis technique (based on default parameters
  #############
  
  # note that for all techniques individually, a suiting default pre-filtering technique is chosen and based on which 
  # DE technique then leads to the highest n_DEGS, the respective pre-filtering methods as optimized in the next step 
  
  DE_techniques <- c("DESeq2", "limma")
  n_DEGS_detechniques <- c()
  
  
  # (I) DESeq2 (using manual pre-filtering and a default pre-filtering threshold of 10 read counts across all samples)
  ############
  
  # run DESeq2 with default settings and a pre-filtering threshold of 10 
  deseq2_results_default <- pre_filt(expression_data, threshold = 10) %>%  
                            geneID_conversion(dupl_removal_method = 1) %>% 
                            deseq_preprocess(phenotype_labels = phenotype_labels) %>%
                            DESeq() %>% results() %>% 
                            as.data.frame() %>% 
                            dplyr::rename(p_adj=padj)
  
  # run ORA tool with DESeq2 result as input and internal parameters in their default configuration 
  cP_ORA_results_DESeq2 <- cP_ORA_input_preparation(deseq2_results_default) %>%
                           enrichGO(OrgDb=organism, ont = "BP") %>% 
                           as.data.frame()
  
  # store number of differentially enriched gene sets 
  n_DEGS_detechniques[1] <- lossfunction_cPORA(cP_ORA_results_DESeq2)
  
  # insert default results into documentation frame 
  doc[1, "n_DEGS"] <- n_DEGS_detechniques[1]
  
  
  # (II) limma (pre-filtering using filterByExpr)
  ############
  
  # get indicator of pre-filtering using filterByExpr
  # -> apply to  expression data set where gene IDs have not been converted yet as 
  # filtering is performed prior to gene ID conversion 
  keep <- DGEList(expression_data, group = phenotype_labels) %>% filterByExpr()
  
  # generate design matrix 
  mm <- model.matrix(~ phenotype_labels)
  
  # as before, pre-filtering is performed prior to the conversion of gene IDs 
  limma_results <- geneID_conversion(expression_data[keep, ], dupl_removal_method = 1) %>% 
                   DGEList(group = phenotype_labels) %>% calcNormFactors() %>%
                   voom(design = mm) %>% lmFit(design = mm) %>% eBayes() %>% topTable(coef = ncol(mm), number = 60000) %>%
                   as.data.frame() %>% dplyr::rename(p_adj = adj.P.Val)
  
  
  #run clusterProfiler's ORA tool (with defaut parameters) for each of the resulting lists of differentially
  #expressed gene sets -> default gene set database: GO
  
  cP_ORA_results_limma<- cP_ORA_input_preparation(limma_results) %>%
                         enrichGO(OrgDb = organism, ont = "BP") %>% 
                         as.data.frame()
  
  # store number of differentially enriched gene sets 
  n_DEGS_detechniques[2] <- lossfunction_cPORA(cP_ORA_results_limma)
  
  
  #################################
  
  # indicator of optimal DE technique
  ind_opt_DE <- min(which(n_DEGS_detechniques == max(n_DEGS_detechniques)))
  # optimal DE technique 
  opt_DE_technique <- DE_techniques[ind_opt_DE]
  
  # get current optimal number of DEGS 
  n_DEGS <- max(n_DEGS_detechniques)
  
  # update documentation frame 
  doc[2,"optimal_parameter"] <- opt_DE_technique
  doc[2,"n_DEGS"] <- n_DEGS
  
  
  
  ##############
  ##### 2. step: choose optimal pre-filtering threshold 
  #############
  
  # note: in this step, different options for pre-filtering are used for the 
  #different DE techniques (in accordance with the corresponding user manuals)
  # for DESeq2, manual pre-filtering is performed using different filtering thresholds 
  # for limma, pre-filtering is performed using 
  # (i) filterByExpr (default) and 
  # (ii) cpm-transformed values 
  
  
  # 1. case: DESeq2 is optimal DE technique
  
  if(opt_DE_technique == "DESeq2"){
    
    # alternative pre-filtering thresholds of X read counts across all samples 
    filt_threshold_alt <- c(10, 50)
    
    # perform pre-filtering according to both (default and alternative) thresholds
    exprdat_list_prefilt <- lapply(X = filt_threshold_alt, 
                                   FUN = pre_filt, 
                                   expression_data = expression_data)
    
    
    # (i) perform differential expression analysis with default parameters
    
    # run DESeq2 workflow
    deseq2_results_prefilt <- lapply(X = exprdat_list_prefilt, FUN = geneID_conversion, dupl_removal_method = 1) %>%  
                              lapply(FUN=deseq_preprocess, phenotype_labels = phenotype_labels) %>%
                              lapply(FUN = DESeq) %>% lapply(FUN = results) %>% 
                              lapply(as.data.frame) %>% lapply(dplyr::rename, p_adj = padj)
    
    #(ii) run clusterProfiler's ORA tool with default parameters 
    cP_ORA_results_prefilt <- lapply(X = deseq2_results_prefilt,FUN = cP_ORA_input_preparation) %>%
                              lapply(FUN = enrichGO, OrgDb=organism, ont = "BP") %>% 
                              lapply(as.data.frame)
    
    #count number of differentially enriched gene sets resulting from the alternative prefiltering thresholds
    n_DEGS_prefilt <- unlist(lapply(X = cP_ORA_results_prefilt, 
                                    FUN = lossfunction_cPORA))
    
    # get index of optimal pre-filtering threshold 
    ind_opt_prefilt <- min(which( n_DEGS_prefilt == max( n_DEGS_prefilt)))
    # get optimal pre-filtering threshold 
    opt_prefilt <-  filt_threshold_alt[ind_opt_prefilt]
    
    # get current optimal number of differentially enriched gene sets 
    n_DEGS <- max(n_DEGS_prefilt)
    
    # update documentation frame 
    doc[3, "n_DEGS"] <- n_DEGS
    doc[3, "optimal_parameter"] <- opt_prefilt
    
    # optimal pre-filtered gene expression data set (gene IDs have not been converted yet)
    exprdat_prefilt_opt <- exprdat_list_prefilt[[ind_opt_prefilt]]
    
    
    
    
    
  } else if(opt_DE_technique  == "limma"){ 
    # for limma, pre-filtering is performed differently, namely using 
    # filterByExpr (default) OR using cpm()
    
    # DEFAULT pre-filtering using filterByExpr
    keep1 <- DGEList(counts=expression_data, 
                     group=phenotype_labels) %>% 
             filterByExpr()
    
    # 1. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 1 count per million in at least 2 samples
    keep2 <- rowSums(cpm(expression_data) > 1) >=2 
   
    
    filt_methods <- c("by FilterbyExpr()", "by cpm>1 in at least 2 samples")
    
    
    # for all methods of pre-filtering: store resulting gene expression data sets in list 
    exprdat_list_prefilt <- list()
    # pre-filtered gene expression data set using keep1 
    exprdat_list_prefilt[[1]] <-  expression_data[keep1, ]
    # pre-filtered gene expression data set using keep2
    exprdat_list_prefilt[[2]] <-  expression_data[keep2, ]
    
      
    limma_resuls <- list()
      
    limma_results <- lapply(FUN = geneID_conversion, X = exprdat_list_prefilt, dupl_removal_method = 1) %>% 
                     lapply(FUN = DGEList, 
                             group = phenotype_labels) %>% 
                     lapply(FUN = calcNormFactors) %>%
                     lapply(FUN = voom,design = mm) %>% lapply(FUN = lmFit,design = mm) %>% lapply(FUN = eBayes) %>% 
                     lapply(FUN = topTable, 
                             coef = ncol(mm), 
                             number = 60000) %>%
                     lapply(FUN = as.data.frame) %>% 
                     lapply(FUN = dplyr::rename, p_adj = adj.P.Val)
      
      
    # perform ORA
    cP_ORA_results_prefilt <- lapply(FUN=cP_ORA_input_preparation, X=limma_results) %>%
                              lapply(FUN=enrichGO,OrgDb=organism, ont = "BP") %>% 
                              lapply(FUN=as.data.frame)
      
    # get number of differentially enriched gene sets for each pre-filtering approach
    n_DEGS_prefilt <- unlist(lapply(X = cP_ORA_results_prefilt, 
                                    FUN = lossfunction_cPORA))
      
    # get index optimal pre-filtering method (in case of tie, choose default filterByExpr)
    ind_opt_prefilt <- min(which(n_DEGS_prefilt == max(n_DEGS_prefilt)))
    # get optimal pre-filtering method 
    opt_prefilt <-  filt_methods[ind_opt_prefilt ]
      
    # current optimal number of differentially enriched gene sets 
    n_DEGS <- max(n_DEGS_prefilt)
      
    
    
    
    # update documentation frame 
    doc[3, "n_DEGS"] <- n_DEGS
    doc[3, "optimal_parameter"] <- opt_prefilt
    
    # optimal pre-filtered gene expression data set (gene IDs have not been converted yet)
    exprdat_prefilt_opt <- exprdat_list_prefilt[[ind_opt_prefilt]]
    
    
  }
  
  
  
  ##############
  ##### 3. step: choose optimal gene expression data resulting from duplicate gene ID removal 
  #############
  # -> use DE technique and manner of pre-filtering as established in first optimization step 
  
  
  
  if(opt_DE_technique == "DESeq2"){
    
    # run DESeq2 workflow for each of the three gene ID conversion schemes 
    DE_results_convID_list <- lapply(FUN = geneID_conversion, X=1:2, expression_data = exprdat_prefilt_opt) %>% 
                              lapply(FUN = deseq_preprocess,phenotype_labels = phenotype_labels) %>%
                              lapply(FUN = DESeq) %>% lapply(FUN = results) %>% 
                              lapply(as.data.frame) %>% lapply(FUN = dplyr::rename,p_adj = padj)
    
    ### perform ORA for both DE results resulting from duplicate gene ID removal 
    cP_ORA_results_convID <- lapply(FUN = cP_ORA_input_preparation, X = DE_results_convID_list) %>%
                             lapply(FUN = enrichGO,OrgDb=organism, ont = "BP") %>% 
                             lapply(FUN = as.data.frame)
    
  } else if(opt_DE_technique == "limma"){
    
    # perform edgeR workflow for each of the three gene ID conversion schemes 
    DE_results_convID_list <- lapply(geneID_conversion, expression_data=exprdat_prefilt_opt, X = 1:2) %>% 
                              lapply(FUN = DGEList,  group = phenotype_labels) %>% 
                              lapply(calcNormFactors) %>% lapply(FUN=voom,design= mm) %>% 
                              lapply(FUN=lmFit,design= mm ) %>% lapply(FUN = eBayes) %>% 
                              lapply(FUN = topTable, coef=ncol(mm), number = 60000) %>%
                              lapply(as.data.frame) %>% lapply(FUN = dplyr::rename,p_adj = adj.P.Val)
    
    ### perform ORA for all of the 3 DE results resulting from duplicate gene ID removal 
    cP_ORA_results_convID <- lapply(FUN = cP_ORA_input_preparation, X = DE_results_convID_list) %>%
                             lapply(FUN = enrichGO,OrgDb=organism, ont = "BP") %>% 
                             lapply(FUN = as.data.frame)
    
  }
  
  
  # for each of the 3 ORA results, get n_DEGS 
  n_DEGS_convID <- unlist(lapply(FUN = lossfunction_cPORA, X = cP_ORA_results_convID))
  
  # index of optimal gene ID removal technique 
  ind_convID <- min(which(n_DEGS_convID == max(n_DEGS_convID)))
  
  
  # update documentation frame 
  doc[4, "optimal_parameter"] <- ind_convID
  doc[4, "n_DEGS"] <- max(n_DEGS_convID)
  
  # generate CONVERTED and PRE-FILTERED gene expression data set as usual 
  exprdat_conv_prefilt_opt <- geneID_conversion(expression_data = exprdat_prefilt_opt, 
                                                dupl_removal_method =  ind_convID) 
  
  
  # get optimal list of differentially expressed genes 
  
  optim_DEGs <- cP_ORA_input_preparation(DE_results_convID_list[[ind_convID]])
  
  # generate optimal ORA results
  ORA_opt <- enrichGO(optim_DEGs, 
                      OrgDb = organism, 
                      ont = "BP") %>% as.data.frame()
  
  #return final results
  return(list(default_ORA = cP_ORA_results_DESeq2, #default results table
              optim_ORA = ORA_opt, #optimal results table
              optim_DE_results = DE_results_convID_list[[ind_convID]], # optimal differential expression analysis results
              ind_organism = ind_organism, # indicate organism 
              documentation = doc)) #documentation frame 
  
}





#######################################################################################
### (II) Optimization of internal parameters ##########################################
#######################################################################################



#######################################################################################
### Create input for ORA in default manner ############################################
#######################################################################################

get_DEGs_default <- function(expression_data, phenotype_labels){
  
  
  # identify organism for which gene expression is measured from the format of the gene IDs (all are Ensembl)
  # ENSEMBL gene ID ENSMUSGXXXXXXXXXXX corresponds to mouse; identifiable by substring "ENSMUSG"
  # ENSEMBL gene ID ENSGXXXXXXXXXXX corresponds to homo sapiens; identifiable by string "ENSG"
  
  
  # for GO and KEGG, the correct organism must be identified, however, they must be 
  # specified differently for KEGG and GO 
  
  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl, 
                         X = c("ENSG", "ENSMUSG"), 
                        x = rownames(expression_data)[1])
  

  # (i) generate input list of differentially expressed genes (default parameters)
  DE_results <- pre_filt(expression_data, threshold = 10) %>% 
                geneID_conversion(dupl_removal_method = 1) %>%
                deseq_preprocess( phenotype_labels = phenotype_labels) %>% 
                DESeq() %>% results() %>% 
                as.data.frame() %>% dplyr::rename(p_adj = padj)
  
  # Return: 
  # (i) differential expression analysis results 
  # (ii) indication of organism (1 stands for human, 2 stands for mouse)
  return(list(DE_results = DE_results, 
              ind_organism = ind_organism))
  
  
  
}




#######################################################################################
###regular KEGG optimization###########################################################
#######################################################################################

#optimize number of differentially enriched DEG with regular (r) KEGG 
#input as provided by function global_DE_optim: results table from DE analysis
cP_ORA_rkegg_optim <- function(DE_results, oragnism_kegg){
  
  # create vector of differentially expressed genes from DE results 
  DEG_vec <- cP_ORA_input_preparation(DE_results)
  
  
  # argument "universe" is left out here since "if missing, the all genes listed 
  #in the database (eg TERM2GENE table) will be used as background.
  
  cp_ora_kegg_default <- enrichKEGG(DEG_vec,
                                    organism = organism_kegg, #homo sapiens
                                    keyType = "kegg")
  
  #determine number of significantly enriched gene sets with current configuration of parameters
  n_DEGS <- lossfunction_cPORA(cp_ora_kegg_default)
  n_DEGS_doc <- n_DEGS
  opt_arg_doc <- "KEGG"
  
  ##########
  #change universe to all genes measured in the experiment (i.e. all genes with a non-missing adjusted
  # p-value )
  univ_adapt <- rownames(DE_results[!is.na(DE_results$p_adj), ])
  #rerun function
  cp_ora_kegg_univ <- enrichKEGG(DEG_vec,
                               organism= organism_kegg, 
                               keyType="kegg",
                               universe=univ_adapt)
  
  #evaluate whether number of DEGS has increased
  if(lossfunction_cPORA(cp_ora_kegg_univ) > n_DEGS){
    #case 1: specification of the alternative universe leads to an increased 
    # number of differentially enriched gene sets 
    
    #update current number of differentially enriched gene sets
    n_DEGS<-lossfunction_cPORA(cp_ora_kegg_univ)
    #for documentation
    n_DEGS_doc<-c(n_DEGS_doc, n_DEGS)
    opt_arg_doc<-c(opt_arg_doc, "Adapted")
    
    
    
    return(list(optim=as.data.frame(cp_ora_kegg_univ), # optimized ORA results
                parameters=opt_arg_doc, # documentation frame 
                n_DEGS=n_DEGS_doc)) # optimal number of differentially enriched gene sets 
    
    
    
  }
  
  else{#case 2: changing universe does not increase number of differentially enriched gene sets
    #means: n_DEGS does not change from default configuration
    
    # update documentation: 
    # add current optimal (i.e. default) number of differentially enriched gene sets 
    n_DEGS_doc<-c(n_DEGS_doc, n_DEGS)
    # add default universe as optimal choice of the universe 
    opt_arg_doc<-c(opt_arg_doc, "Default")
    
    
    
    return(list(optim=as.data.frame(cp_ora_kegg_default), # initial ORA results
                parameters=opt_arg_doc, # documentation frame 
                n_DEGS=n_DEGS_doc)) # number of differentially enriched gene sets 
    
  }
  
}


#######################################################################################
### GO optimization ###################################################################
#######################################################################################

cPORA_GO_optim <- function(DE_results, organism_go){

  #obtain input vector for KEGG analysis with input preparation function
  DEG_vec <- cP_ORA_input_preparation(DE_results)

  
  # extract adapted universe from DE results
  # -> all genes whose differential expression was tested 
  
  
  # adapted universe consists of all genes with a nonNA adjusted p-value  
  univ_adapt<-rownames(DE_results[!is.na(DE_results$p_adj), ])
  

  
  #run function with default configuration 
  cp_ora_go_default<-enrichGO(DEG_vec,
                              OrgDb= organism_go,
                              ont="BP")
  
  
  
  #current number of differentially enriched gene sets
  n_DEGS<-lossfunction_cPORA(cp_ora_go_default)
  #store for documentation of optimization process 
  n_DEGS_doc<-lossfunction_cPORA(cp_ora_go_default)
  opt_arg_doc<-"GO"
  
  ##############
  #change universe to all genes measured in the experiment 
  cp_ora_go_univ<-enrichGO(DEG_vec,
                           OrgDb= organism_go,
                           ont="BP",
                           universe = univ_adapt)
  
  
  #if the number of differentially enriched gene sets increases with updated 
  #universe, proceed with new universe 
  if(lossfunction_cPORA(cp_ora_go_univ) > n_DEGS){
    
    # count number of differentially enriched gene sets 
    n_DEGS<-lossfunction_cPORA(cp_ora_go_univ)
    
    # update documentation frame 
    n_DEGS_doc<-c(n_DEGS_doc, n_DEGS)
    opt_arg_doc<-c(opt_arg_doc, "Adapted")
    
    
    return(list(optim=as.data.frame(cp_ora_go_univ),
                parameters=opt_arg_doc, 
                n_DEGS=n_DEGS_doc))
    
  } else{ #case 1: updating universe does not lead to an increase in the number 
    #of differentially enriched gene sets 
    
    # update documentation 
    # add (unchanged but current optimal) number of differentially enriched gene sets 
    n_DEGS_doc <- c(n_DEGS_doc, n_DEGS)
    # add information that default universe is the optimal choice in terms of the number 
    # of differentially enriched gene sets 
    opt_arg_doc <- c(opt_arg_doc, "Default")
    
    
    
    return(list(optim=as.data.frame(cp_ora_go_default), # (default but current optimal ORA results)
                parameters = opt_arg_doc, # documentation frame of parameter choice 
                n_DEGS=n_DEGS_doc)) # documentation frame of number of differentially enriched gene sets 
    
  }
  
  
}

##########################################################################################
###global ORA optimization################################################################
##########################################################################################
#stepwise analysis: from GO and KEGG choose gene set database which yields highest
#number of DEGS with DEFAULT parameters
#-> choose method yielding highest number of DEGS and perform "Stellschraubenanalyse" with
#respective function 

#note: as initial gene expression data set is in Ensembl gene ID format,

cP_ORA_internparam_optim <- function(DE_results, ind_organism){
  
  
  # (i) specification of organism for GO
  organism_go <- c(org.Hs.eg.db , org.Mm.eg.db)[ind_organism][[1]]
  
  # (ii) specification of organism for KEGG
  # hsa: homo sapiens (human)
  # mmu: mus musculus (mouse)
  organism_kegg <- (c("hsa", "mmu"))[ind_organism][[1]]
  
  
  #obtain input vector for KEGG analysis with input preparation function
  DEG_vec <- cP_ORA_input_preparation(DE_results)
  
  
  

  #prepare dataframe to document stepwise optimization
  doc <- data.frame(step=c("Default", "Gene Set Database","Universe"),
                    optimal_parameter=NA, 
                    n_DEGS=NA)
  
  ##########
  ### step 1: choose optimal geneset database
  ##########

  
  #(i) default: GO 
  cP_default_go <- enrichGO(DEG_vec,
                          OrgDb = organism_go,
                          ont="BP")
  
  #default result (incl. default unranked list)
  doc[1,"n_DEGS"] <- lossfunction_cPORA(cP_default_go)
  
  
  
  # (ii) alternative 1: KEGG 
  cP_default_kegg <- enrichKEGG(DEG_vec,
                                organism= organism_kegg, #homo sapiens
                                keyType="kegg",
                                pAdjustMethod = "BH")
  
  
  #compute lossfunction for all 3 ORA results and store store them in a matrix
  n_DEGS_vec <- c()
  n_DEGS_vec[1] <- lossfunction_cPORA(cP_default_go)
  n_DEGS_vec[2] <- lossfunction_cPORA(cP_default_kegg)
  #indicate which gene set database yields highest number of DEGS
  ind_db  <- c("GO","KEGG")
  
  #depending on which gene set database yields highest number of DEGS with default
  #configuration, perform optimization with function corresponding function 
  #in case of a tie, procees with gene set database with lower index 
  best_db <- ind_db[min(which(n_DEGS_vec == max(n_DEGS_vec)))]
  
  #perform optimization for gene set database which yields highest 
  #number of DEGS with default configuration
  
  if(best_db == "GO"){
    
    #perform GO "Stellschraubenanalyse"
    cP_ora_optim <- cPORA_GO_optim(DE_results, 
                                   organism_go)
    
  }else if(best_db == "KEGG"){
    #perform KEGG "Stellschraubenanalyse"
    cP_ora_optim <- cP_ORA_rkegg_optim(DE_results, 
                                   organism_kegg)
    
  }
  
  doc[2:3, "n_DEGS"] <- cP_ora_optim$n_DEGS
  doc[2:3, "optimal_parameter"] <- cP_ora_optim$parameters
  
  #return final result
  return(list(default_list = as.data.frame(cP_default_go), #results with default input and ORA's parameters in default configuration
              optim = cP_ora_optim$optim, # optimal ORA result
              documentation = doc)) #documentation frame 
  
}



################################################################################
### (III) Joint optimization of Pre-Processing and Internal Parameters #########
################################################################################

cP_ORA_joint_optimization <- function(expression_data, phenotype_labels){
  
  # run optimization of pre-processing steps
  optim_preprocess <- cP_ORA_preprocess_optim(expression_data, 
                                              phenotype_labels)
  
  # run optimization of internal parameters of ORA (with optimal results from
  # the optimization of the pre-processing)
  optim_internalparam <- cP_ORA_internparam_optim(optim_preprocess$optim_DE_results, 
                                                  optim_preprocess$ind_organism)
  
  
  # merge documentation frames 
  doc <- rbind(optim_preprocess$documentation, 
               optim_internalparam$documentation[optim_internalparam$documentation$step != "Default", ])
  
  
  
  return(list(ORA_default = optim_preprocess$default, 
              ORA_optim = optim_internalparam$optim, 
              documentation = doc))
  
  
}




################################################################################
### Run Optimization Functions #################################################
################################################################################


# #############
# ### Pickrell 
# #############
# 
# # true sample labels
# optim_ORA_results_Pickrell_originalphenotype <- cP_ORA_joint_optimization(Biobase::exprs(pickrell.eset), 
#                                                                           pickrell.eset$gender)
# # save results
# save(optim_ORA_results_Pickrell_originalphenotype, 
#      file = "./Results/ORA_Results_Pickrell_OriginalPhenotype.RData")
# 
# 
# 
# # 10 permuted sample labels
# optim_ORA_results_Pickrell_phenotypepermutation <- apply(FUN = cP_ORA_joint_optimization, 
#                                                          expression_data = Biobase::exprs(pickrell.eset), 
#                                                          X = phen_pickrell, MARGIN = 2)
# 
# # save results
# save(optim_ORA_results_Pickrell_phenotypepermutation, 
#      file = "./Results/ORA_Results_Pickrell_PhenotypePermutations.RData")
# 
# #############
# ### Bottomly 
# #############
# 
# # true sample labels
# optim_ORA_results_Bottomly_originalphenotype <- cP_ORA_joint_optimization(Biobase::exprs(bottomly.eset), 
#                                                                           bottomly.eset$strain)
# 
# # save results
# save(optim_ORA_results_Bottomly_originalphenotype, 
#      file = "./Results/ORA_Results_Bottomly_OriginalPhenotype.RData")
# 
# 
# # 10 permuted sample labels
# optim_ORA_results_Bottomly_phenotypepermutation <- apply(FUN = cP_ORA_joint_optimization, 
#                                                          expression_data = Biobase::exprs(bottomly.eset), 
#                                                          X = phen_bottomly, 
#                                                          MARGIN = 2)
# 
# 
# # save results
# save(optim_ORA_results_Bottomly_phenotypepermutation, 
#      file = "./Results/ORA_Results_Bottomly_PhenotypePermutations.RData")
# 
