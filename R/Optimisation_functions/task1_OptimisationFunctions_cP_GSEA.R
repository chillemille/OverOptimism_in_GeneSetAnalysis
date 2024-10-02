#clusterProfiler GSEA optimization process

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(apeglm)
library(ggplot2)
library(dplyr)


##################################################################################
##create ranked list from DE results##############################################
##################################################################################

#for DESeq2 (method = "DESeq2") and limma (method = "limma")
#rankby must be in c("p_value", "lfc") to perform ranking based on
#(i) p-value (rank = sign(lfc)*(-1)*log10(unadjusted_pvalue)
#(ii) log fold changes

rankedList_cP_GSEA <- function(DE_results, rankby, method){

 if(method == "DESeq2"){#create ranking based on DESeq2 results table

   if(rankby == "lfc") {
     #ranking by log2 fold change
     #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
     rankvec <-
       as.vector(DE_results[!is.na(DE_results$pvalue),]$log2FoldChange)
     names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue),])
     # sort in decreasing manner
     rankvec <- sort(rankvec, decreasing = TRUE)

   }else if (rankby == "p_value") {
     #ranking by p-value
     #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
     rankvec <-
       as.vector(sign(DE_results[!is.na(DE_results$pvalue),]$log2FoldChange) *
                   (-1) * log10(DE_results[!is.na(DE_results$pvalue),]$pvalue))
     names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue),])
     # sort in decreasing manner
     rankvec <- sort(rankvec, decreasing = TRUE)
   }
      }


 else if(method == "limma"){#create ranking based on limma's results table

 # first step: replace p-values of 0 with the smallest representable positive number
 # in R (necessary since one term of ranking metric is equal to log10(p-value))

   if(rankby == "lfc") {
     #ranking based on log2 fold change
     rankvec <- as.vector(DE_results$logFC)
     names(rankvec) <- rownames(DE_results)
     rankvec <- sort(rankvec, decreasing = TRUE)
   }

   else if (rankby == "p_value") {
     #ranking based on p-value
     rankvec <-
       as.vector(sign(DE_results$logFC) * (-1) * log10(DE_results$P.Value))
     names(rankvec) <- rownames(DE_results)
     rankvec <- sort(rankvec, decreasing = TRUE)
   }
 }
  return(rankvec)
}


##################################################################################
###lossfunction###################################################################
##################################################################################
lossfunction_cP_GSEA <- function(gsa_result){

 #count number of gene sets with adjusted p-value smaller than 0.05
 n_DEGS <- sum(as.data.frame(gsa_result)$p.adjust < 0.05)

 return(n_DEGS)
}


#######################################################################################
### Default GSEA (for KEGG or GO) #####################################################
#######################################################################################

GSEA_pipeline_default <- function(gene_ranking, geneset_database, exp = 1, organism ){


 if(geneset_database == "KEGG"){

 # set seed
 set.seed(713)

 # run default GSEA
 GSEA <- gseKEGG(gene_ranking,
                 organism = organism,
                 eps = 0, # to ensure that all p-values are computed exactly
                 seed = TRUE, # ensure that results are reproducible
                 exponent = exp)

 } else if(geneset_database == "GO"){

 # set seed
 set.seed(713)

 GSEA <- gseGO(gene_ranking,
               OrgDb = organism,
               ont = "BP", # subontology Molecular Function
               eps = 0, # to ensure that all p-values are computed exactly
               seed = TRUE, # ensure that results are reproducible
               exponent = exp)
 }


 # return GSEA results
 return(GSEA %>% as.data.frame())



}


#######################################################################################
### (I) optimization of preparatory parameters ########################################
#######################################################################################
#optimization of required input for clusterProfiler's GSEA tool: list of all genes in
#experiment ranked by their magnitude of differential expression


cP_GSEA_preprocess_optim <- function(expression_data, geneset_database, phenotype_labels){

 # identify organism for which gene expression is measured from the format of the gene IDs (all are Ensembl)
 # ENSEMBL gene ID ENSMUSGXXXXXXXXXXX corresponds to mouse; identifiable by substring "ENSMUSG"
 # ENSEMBL gene ID ENSGXXXXXXXXXXX corresponds to homo sapiens; identifiable by string "ENSG"


 # for GO and KEGG, the correct organism must be identified, however, they must be
 # specified differently for KEGG and GO


 if(geneset_database == "GO"){

 # indicate which of the two strings can be found in the gene IDs of the expression data at hand
 ind_organism <- sapply(FUN = grepl, X = c("ENSG", "ENSMUSG"), x = rownames(expression_data)[1])

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



 #create frame to document progression of the optimization workflow
 doc <- data.frame(step = c("Default", "Ranking Metric", "Pre-Filtering Threshold", "Duplicate Gene ID Removal"),
                   optimal_parameter = NA,
                   n_DEGS = NA)



 ###########
 # 1. step: Choose optimal manner of creating ranking metric (default method is DESeq2)
 ###########

 ranking_metrics <- c("DESeq2", "limma")


 gsea_results_rankmetrics <- list()

 n_DEGS_rankmetrics <- c()


 ### DEFAULT option: ranking using DESeq2 ###

 # create default DESeq2 results (incl. default pre-filtering and duplicate gene ID removal technique)
 DESeq2_results <- pre_filt(expression_data, threshold = 10) %>% #
                   geneID_conversion(dupl_removal_method = 1) %>%
                   deseq_preprocess(phenotype_labels = phenotype_labels) %>%
                   DESeq() %>%
                   lfcShrink(coef = "condition_treated_vs_untreated", type = "apeglm") %>%
                   as.data.frame()

 # create ranked list from DESEq2 results
 ranked_list_DESeq2 <- rankedList_cP_GSEA(DESeq2_results,
                                     rankby = "p_value",
                                     method = "DESeq2")

 # run GSEA (using GO)
 gsea_results_rankmetrics[[1]] <- GSEA_pipeline_default(ranked_list_DESeq2,
                                                        geneset_database,
                                                        organism = organism)

 # store n_DEGS to compare to remaining rankings
 n_DEGS_rankmetrics[1] <- lossfunction_cP_GSEA(gsea_results_rankmetrics[[1]])


 # update documentation frame
 doc[1, "n_DEGS"] <- lossfunction_cP_GSEA(gsea_results_rankmetrics[[1]])



 ### alternative option 1 : ranking using limma ###

 # pre-filtering:
 keep <- DGEList(expression_data, group = phenotype_labels) %>%
         filterByExpr()

 # design matrix
 mm <- model.matrix(~phenotype_labels)

 # as before, pre-filtering is performed prior to the conversion of gene IDs
 # filtering indicator keep recycled from alternative 1 edgeR
 limma_results <- geneID_conversion(expression_data[keep, ], dupl_removal_method = 1) %>%
                  DGEList(group = phenotype_labels) %>%
                  calcNormFactors() %>%
                  voom(design = mm) %>% lmFit(design = mm) %>%
                  eBayes() %>% topTable(coef = ncol(mm), number = 100000)


 #perform ranking by p-value
 ranked_list_limma <- rankedList_cP_GSEA(limma_results,
                                    rankby = "p_value",
                                    method = "limma")
 #perform default GSEA
 gsea_results_rankmetrics[[2]] <- GSEA_pipeline_default(ranked_list_limma,
                                                        geneset_database,
                                                        organism = organism)


 n_DEGS_rankmetrics[2] <- lossfunction_cP_GSEA( gsea_results_rankmetrics[[2]])




 ###choose optimal ranking metric and update documentation frame
 ind_opt_rankmetric <- min(which(n_DEGS_rankmetrics == max(n_DEGS_rankmetrics), arr.ind = TRUE))
 opt_rankmetric <- ranking_metrics[ind_opt_rankmetric]
 doc[2, "optimal_parameter"] <- opt_rankmetric
 doc[2, "n_DEGS"] <- max(n_DEGS_rankmetrics)



 ###########
 # 2. step: Choose optimal pre-filtering threshold
 ###########

 # NOTE: different methods of pre-filtering are used for the different methods to create the
 # ranking of the genes
 # DESeq2 and t-Test: manual pre-filtering that exclude all genes with less than X counts across all samples
 # edgeR and limma: pre-filtering is performed using functions filterByExpr() and cpm()


 if(opt_rankmetric == "DESeq2"){

 # create pre-filtered gene expression data sets using manual pre-filtering with different
 # filtering thresholds

 # pre-filtering thresholds
 prefilt_thresh_list <- c(10,50)
 # vector to store resulting number of n_DEGS for the different pre-filtering thresholds

 exprdat_prefilt_list <- lapply(X = prefilt_thresh_list,
                                FUN = pre_filt,
                                expression_data = expression_data)

 GSEA_results_prefilt_list <- list()


  # run DESEq2 workflow for each of the pre-filtered gene expression data sets
  DESeq2_results_prefilt <- lapply(FUN = geneID_conversion, X = exprdat_prefilt_list, dupl_removal_method = 1) %>%
                            lapply(FUN = deseq_preprocess, phenotype_labels = phenotype_labels) %>%
                            lapply(FUN = DESeq) %>%
                            lapply(FUN = lfcShrink,coef = "condition_treated_vs_untreated", type = "apeglm") %>%
                            lapply(FUN = as.data.frame)

  # generate ranking for each of the DESeq2 results
  ranked_lists_DESeq2 <- lapply(FUN = rankedList_cP_GSEA,
                                X = DESeq2_results_prefilt,
                                rankby = "p_value",
                                method = "DESeq2")

  # run GSEA (using GO)
  GSEA_results_prefilt_list <- lapply(FUN = GSEA_pipeline_default,
                                      X = ranked_lists_DESeq2,
                                      geneset_database = geneset_database,
                                      organism = organism)





 # count number of differentially enriched gene sets resulting from each pre-filtering threshold
 n_DEGS_prefilt <- unlist(lapply(FUN = lossfunction_cP_GSEA,
                                 X = GSEA_results_prefilt_list))

 # get index of optimal pre-filtering threshold (in case of ties, choose lower index)
 ind_opt_prefilt <- min(which( n_DEGS_prefilt == max( n_DEGS_prefilt)))

 # update documentation frame
 doc[3, "optimal_parameter"] <- prefilt_thresh_list[ind_opt_prefilt]
 doc[3, "n_DEGS"] <- max( n_DEGS_prefilt)



 } else if(opt_rankmetric == "limma"){

 # create pre-filtered gene expression data sets using functions filterByExpr() and cpm()
 # filtering thresholds

 keep1 <- DGEList(counts = expression_data, group = phenotype_labels) %>% filterByExpr()
 # 1. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 1 count per million in at least 2 samples
 keep2 <- rowSums(cpm(expression_data) > 1) >= 2
 # 2. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 2 count per million in at least 3 samples


 # specify both (default and alternative) pre-filtering approaches
 filt_methods <- c("by FilterbyExpr()", "by cpm>1 in at least 2 samples")


 # for all methods of pre-filtering: store resulting gene expression data sets in list
 exprdat_prefilt_list <- list()
 # pre-filtered gene expression data set using keep1
 exprdat_prefilt_list[[1]] <- expression_data[keep1, ]
 # pre-filtered gene expression data set using keep2
 exprdat_prefilt_list[[2]] <- expression_data[keep2, ]


  # run limma workflow for both pre-filtered gene expression data sets
  limma_results_prefilt <- lapply(FUN = geneID_conversion,
                                  X = exprdat_prefilt_list,
                                  dupl_removal_method = 1) %>%
                           lapply(FUN = DGEList,group = phenotype_labels) %>%
                           lapply(FUN = calcNormFactors) %>%
                           lapply(FUN = voom,design = mm) %>%
                           lapply(FUN = lmFit,design = mm) %>%
                           lapply(FUN = eBayes) %>%
                           lapply(FUN = topTable,coef = ncol(mm), number = 100000)


  #perform ranking by p-value for each of the limma results
  ranked_lists_limma <- lapply(FUN = rankedList_cP_GSEA,
                               X = limma_results_prefilt,
                               rankby = "p_value",
                               method = "limma")

  # run GSEA (GO) for each of the rankings generated with limma
  GSEA_results_prefilt_list <- lapply(FUN = GSEA_pipeline_default,
                                      X = ranked_lists_limma,
                                      geneset_database = geneset_database,
                                      organism = organism)


 # count number of differentially enriched gene sets resulting from each pre-filtering threshold
 n_DEGS_prefilt <- unlist(lapply(FUN = lossfunction_cP_GSEA,
                                 X = GSEA_results_prefilt_list))

 # get index of optimal pre-filtering threshold (in case of ties, choose lower index)
 ind_opt_prefilt <- min(which(n_DEGS_prefilt == max( n_DEGS_prefilt)))

 # update documentation frame
 doc[3, "optimal_parameter"] <- filt_methods[ind_opt_prefilt ]
 doc[3, "n_DEGS"] <- max(n_DEGS_prefilt)


 }

 # get OPTIMAL PRE-FILTERED gene expression data set
 # gene IDs are in Ensembl gene ID format as gene id conversion is performed in next step

 exprdat_prefilt_opt <- exprdat_prefilt_list[[ind_opt_prefilt]]




 ###########
 # 3. step: Choose optimal manner of duplicate gene ID removal
 ###########

 if(opt_rankmetric == "DESeq2"){

 # run DESEq2 workflow for each of the manners of duplicate gene ID removal
 DESeq2_results_convID <- lapply(FUN = geneID_conversion,
                                 X = 1:2,
                                 expression_data = exprdat_prefilt_opt) %>%
                          lapply(FUN = deseq_preprocess,
                                 phenotype_labels = phenotype_labels) %>%
                          lapply(FUN = DESeq) %>%
                          lapply(FUN = lfcShrink,
                                 coef = "condition_treated_vs_untreated",
                                 type = "apeglm") %>%
                          lapply(FUN = as.data.frame)

 # generate ranking for each of the DESeq2 results
 ranked_lists_convID <- lapply(FUN = rankedList_cP_GSEA,
                               X = DESeq2_results_convID,
                               rankby = "p_value" ,
                               method = "DESeq2")

 # run GSEA (using GO)
 GSEA_results_convID_list <- lapply(FUN = GSEA_pipeline_default,
                                    X = ranked_lists_convID,
                              geneset_database = geneset_database,
                              organism = organism)



 } else if(opt_rankmetric == "limma"){

 # run limma workflow for each of the manners of duplicate gene ID removal
 limma_results_convID <- lapply(FUN = geneID_conversion,
                                X = 1:2,
                                expression_data = exprdat_prefilt_opt) %>%
                         lapply(FUN = DGEList,group = phenotype_labels) %>%
                         lapply(FUN = calcNormFactors) %>%
                         lapply(FUN = voom,design = mm) %>%
                         lapply(FUN = lmFit,design = mm) %>%
                         lapply(FUN = eBayes) %>%
                         lapply(FUN = topTable,coef = ncol(mm), number = 100000)

 #perform ranking by p-value for each of the limma results
 ranked_lists_convID <- lapply(FUN = rankedList_cP_GSEA,
                               X = limma_results_convID,
                               rankby = "p_value",
                               method = "limma")

 GSEA_results_convID_list <- lapply(FUN = GSEA_pipeline_default,
                                    X = ranked_lists_convID,
                                    geneset_database = geneset_database,
                                    organism = organism)


 }

 # count number of differentially enriched gene sets for each of the duplicate gene ID removal methods
 n_DEGS_convID <- unlist(lapply(FUN = lossfunction_cP_GSEA,
                                X = GSEA_results_convID_list))
 # get index of optimal duplicate gene ID removal method (choose lower index in case of tie)
 ind_convID_opt <- min(which( n_DEGS_convID == max( n_DEGS_convID)))

 # update documentation frame
 doc[4, "optimal_parameter"] <- ind_convID_opt
 doc[4, "n_DEGS"] <- max( n_DEGS_convID)

 # optimal PRE-FILTERED and CONVERTED gene expression data set
 exprdat_prefilt_convID_opt <- geneID_conversion(expression_data = exprdat_prefilt_opt,
                                                 dupl_removal_method = ind_convID_opt)






 return(list(default_GSEA = gsea_results_rankmetrics[[1]], #default GSEA results
    optim_GSEA = GSEA_results_convID_list[[ind_convID_opt]], #optimal GSEA results
    optim_ranking = ranked_lists_convID[[ind_convID_opt]], # optimal ranking
    ind_organism = ind_organism, # indicate organism
    documentation = doc))


}




##################################################################################
#####GSEA - optimization of internal parameters ##################################
##################################################################################

cP_GSEA_go_optim <- function(ranking, organism_GO){

 #for documentation of optimization process
 doc <- data.frame(step = c("Gene Set Database", "Exponent"),
                    optimal_parameter = NA,
                    n_DEGS = NA)

 #default GSEA
 set.seed(713)
 go_default <- gseGO(ranking,
                     exponent = 1, #default exponent
                     ont = "BP",
                     OrgDb = organism_GO,
                     seed = TRUE) #default multiple test adjustment

 #count number of DEGS
 n_DEGS <- lossfunction_cP_GSEA(go_default)
 #documentation
 doc[1, "n_DEGS"] <- n_DEGS
 doc[1, "optimal_parameter"] <- "GO"


 ##########
 #change 1: exponent
 exp_options <- c(1,0,1.5,2) # possible options for exponent

 n_DEGS_exp <- rep(NA, times = length(exp_options))
 n_DEGS_exp[1] <- n_DEGS #number of DEGS resulting from exponent = 1

 #for each of the possible exponent (except for exp = 1 for which GSEA has already been run )
 # run GSEA and store the resulting number of differentially enriched gene sets
 for(i in 2:length(exp_options)){

 #perform gsea with respective exponent
 set.seed(713)
 go_exp <- gseGO(ranking,
                 exponent = exp_options[i], #provide respective component
                 OrgDb = organism_GO,
                 seed = TRUE,
                 ont = "BP")

 #store n_DEGS in vector to compare to all exponents
 n_DEGS_exp[i] <- lossfunction_cP_GSEA(go_exp)

 }

 #determine which exponent leads to highest number of DEGS
 #in case of tie choose lower index
 ind_exp_max <- min(which(n_DEGS_exp == max(n_DEGS_exp), arr.ind = TRUE))

 #update exponent
 exp <- exp_options[ind_exp_max]
 #update number of DEGS
 n_DEGS <- max(n_DEGS_exp)
 #for documentation
 doc[2, "n_DEGS"] <- n_DEGS
 doc[2, "optimal_parameter"] <- exp



 #return optimal GSEA results and optimal parameters
 set.seed(713)
 go_final <- gseGO(ranking,
                   exponent = 1,
                   ont = "BP",
                   OrgDb = organism_GO,
                   seed = TRUE)

 #return optimal result GSEA and optimal parameters
 return(list(optim_gsea = as.data.frame(go_final),
             documentation = doc))
}









##################################################################################
#####GSEA KEGG optimization#######################################################
##################################################################################

cP_GSEA_kegg_optim <- function(ranking, organism_KEGG){

 #for documentation of optimization process
 doc <- data.frame(step = c("Gene Set Database", "Exponent"),
                   optimal_parameter = NA,
                   n_DEGS = NA)

 #default GSEA
 set.seed(713)
 kegg_default <- gseKEGG(ranking,
                         exponent = 1, #default exponent
                         organism = organism_KEGG,
                         keyType = "kegg",
                         seed = TRUE) #default multiple test adjustment

 #count number of DEGS
 n_DEGS <- lossfunction_cP_GSEA(kegg_default)
 #documentation
 doc[1, "n_DEGS"] <- n_DEGS
 doc[1, "optimal_parameter"] <- "KEGG"


 ##########
 #change 1: exponent
 exp_options <- c(1,0,1.5,2) #possible options for exponent

 n_DEGS_exp <- rep(NA, times = length(exp_options))
 n_DEGS_exp[1] <- n_DEGS #number of DEGS resulting from exponent = 1

 for(i in 2:length(exp_options)){

 #perform gsea with respective exponent
 set.seed(713)
 kegg_exp <- gseKEGG(ranking,
                     exponent = exp_options[i],
                     organism = organism_KEGG,
                     keyType = "kegg",
                     seed = TRUE)

 #store n_DEGS in vector to compare to all exponents
 n_DEGS_exp[i] <- lossfunction_cP_GSEA(kegg_exp)

 }

 #determine which exponent leads to highest number of DEGS
 #in case of tie choose lower index
 ind_exp_max <- min(which(n_DEGS_exp == max(n_DEGS_exp), arr.ind = TRUE))

 #update exponent
 exp <- exp_options[ind_exp_max]
 #update number of DEGS
 n_DEGS <- max(n_DEGS_exp)
 doc[2, "n_DEGS"] <- n_DEGS
 doc[2, "optimal_parameter"] <- exp





 #return optimal GSEA results and optimal parameters
 set.seed(713)
 kegg_final <- gseKEGG(ranking,
                      exponent = exp,
                      organism = organism_KEGG,
                      keyType = "kegg",
                      seed = TRUE)

 #return optimal result GSEA and optimal parameters
 return(list(optim_gsea = as.data.frame(kegg_final),
    documentation = doc))
}



##################################################################################
#global optimization##############################################################
##################################################################################

#optimization of the internal parameters of clusterProfiler's GSEA
# ind_organism: indicates from which organism the gene expression data originates
# 1: human
# 2: mouse
cP_GSEA_internparam_optim <- function(ranking, ind_organism){

 ### generate organism indicator in required format for GO and KEGG

 # (i) specification of organism for GO
 organism_GO <- c(org.Hs.eg.db , org.Mm.eg.db)[ind_organism][[1]]


 # (ii) specification of organism for KEGG
 # hsa: homo sapiens (human)
 # mmu: mus musculus (mouse)
 organism_KEGG <- (c("hsa", "mmu"))[ind_organism][[1]]




 ##########
 # 1. step: Optimization of the choice of gene set database
 ##########



 #documentation of optimization process
 doc <- data.frame(step = "Default",
                   optimal_parameter = NA,
                   n_DEGS = NA)

 #storage for number of DEGS of all three options
 geneset_db <- c("GO","KEGG")
 n_DEGS_geneset_db <- c()

 #change 1: Gene Set Database
 #option 1: GO (default)
 #set seed
 set.seed(713)
 gsea_go <- gseGO(ranking,
                 ont = "BP",
                 OrgDb = organism_GO,
                 seed = TRUE)

 #count current number of DEGS
 n_DEGS_geneset_db[1] <- lossfunction_cP_GSEA(gsea_go)
 doc[1, "n_DEGS"] <- n_DEGS_geneset_db[1]


 #option 2: KEGG
 set.seed(713)
 gsea_kegg <- gseKEGG(ranking,
                      organism = organism_KEGG,
                      keyType = "kegg",
                      seed = TRUE)

 n_DEGS_geneset_db[2] <- lossfunction_cP_GSEA(gsea_kegg)


 #determine which gene set database leads to highest number of DEGS with
 #remaining parameters in default configuration
 ind_max <- min(which(n_DEGS_geneset_db == max(n_DEGS_geneset_db), arr.ind = TRUE))

 #optimal gene set database
 geneset_db_opt <- geneset_db[ind_max]
 n_DEGS <- n_DEGS_geneset_db[ind_max]



 #proceed with optimal gene set database
 if(geneset_db_opt == "GO"){

 #perform GO GSEA optimization
 optim_GSEA_results <- cP_GSEA_go_optim(ranking, organism_GO)

 }else if(geneset_db_opt == "KEGG"){

 #perform KEGG GSEA optimization
 optim_GSEA_results <- cP_GSEA_kegg_optim(ranking,
                                         organism_KEGG)

 }

 #complete documentation frame
 doc_complete <- rbind(doc,
                       optim_GSEA_results$documentation)



 #return optimal GSEA result and optimal parameters
 return(list(default_GSEA = gsea_go, #default GSEA
             optimal_GSEA = optim_GSEA_results$optim_gsea, #optimal GSEA
             optim_ranking = optim_GSEA_results$optim_ranking,
             documentation = doc_complete)) #optimal GSEA parameters

}

##################################################################################
### Joint optimization of pre-processing and internal parameters #################
##################################################################################

cP_GSEA_joint_optimization <- function(expression_data, phenotype_labels){

 # optimize ranked input list of genes to maximize number of differentially enriched gene sets
 optim_preprocess <- cP_GSEA_preprocess_optim(expression_data, "GO", phenotype_labels)

 # optimize internal parameters of GSEA tools to maximize number of differentially enriched gene sets
 optim_internalparam <- cP_GSEA_internparam_optim(optim_preprocess$optim_ranking,
                        optim_preprocess$ind_organism)

 # merge documentation frames
 doc <- rbind(optim_preprocess$documentation,
    optim_internalparam$documentation[optim_internalparam$documentation$step != "Default", ])


 return(list(default_GSEA = optim_preprocess$default_GSEA,
    optim_GSEA = optim_internalparam$optimal_GSEA,
    documentation = doc))


}


