#clusterProfiler's GSEA: Minimization of the adjusted p-value or rank of a gene set of choice


library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(apeglm)



##################################################################################
##create ranked list from DE results##############################################
##################################################################################

#for DESeq2 (method = "DESeq2") and limma (method = "limma")
#rankby must be in c("p_value", "lfc") to perform ranking based on
#(i) p-value (rank = sign(lfc)*(-1)*log10(unadjusted_pvalue)
#(ii) log fold changes
rankedList_cP <- function(DE_results, rankby, method){

 if(method == "DESeq2"){#create ranking based on DESeq2 results table

 # first step: replace p-values of 0 with the smallest representable positive number
 # in R (necessary since one term of ranking metric is equal to log10(p-value))

 DE_results$pvalue[ DE_results$pvalue == 0 ] <- min(DE_results$pvalue[ DE_results$pvalue > 0 ]) / 10

 if(rankby == "lfc") {
   #ranking by log2 fold change
   #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
   rankvec <-
     as.vector(DE_results[!is.na(DE_results$pvalue), ]$log2FoldChange)
   names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue), ])
   rankvec <- sort(rankvec, decreasing = TRUE)

 }else if (rankby == "p_value") {
   #ranking by p-value
   #remove rows containing NA p-values (relevant if Cook's outlier detection turned on)
   rankvec <-
     as.vector(sign(DE_results[!is.na(DE_results$pvalue), ]$log2FoldChange) *
                 (-1) * log10(DE_results[!is.na(DE_results$pvalue), ]$pvalue))
   names(rankvec) <- rownames(DE_results[!is.na(DE_results$pvalue), ])
   rankvec <- sort(rankvec, decreasing = TRUE)
 }

 }else if(method == "limma"){#create ranking based on limma results table

 # first step: replace p-values of 0 with the smallest representable positive number
 # in R (necessary since one term of ranking metric is equal to log10(p-value))

 DE_results$P.Value[ DE_results$P.Value == 0 ] <- min(DE_results$P.Value[ DE_results$P.Value > 0 ]) / 10

 if(rankby == "lfc") {
   #ranking based on log2 fold change
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


################################################################################
### get gene's adjusted p-value or rank (based on what is needed) ##############
################################################################################

############### NOTE:
# encountered problem in the optimization process of the rank or (adjusted) pvalue of a given gene set:
#clusterProfiler only reports those gene sets which contain at least a single gene from the input list
#of differentially expressed genes.

# this specifically means that, if a gene set does not contain any genes from the input list, function
#pvalue_rank_GSEA returns integer(0) which cannot be compared to numerical number.
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
# -> note that the relative ranks are relative to the total number of gene sets in the results
# such that a rank of 1 means that the gene set has the highest adjusted p-value (i.e. of least
# relevance) in the GSA results

#note: function works for GO as well as KEGG as results table "GSEA_results"
#contains relevant information automatically
#argument term must be in the form of a GO ID resp. KEGG ID
pvalue_rank_GSEA <- function(term, GSEA_results, metric){

 if(metric == "rank"){

   # prepare the relative ranks such that all gene sets with the same adjusted p-value
   # are given the same (i.e. average) rank
   # We divide by the maximum rank such that those gene sets with the highest
   # adjusted p-values are given the (relative) rank 1, which is the worst rank

   GSEA_results$ranks <- rank(GSEA_results$p.adjust) /max(rank(GSEA_results$p.adjust))

   # get the relative rank of the gene set of interest
   rank <- GSEA_results$ranks[grep(term, GSEA_results$ID)]


   # if the gene set is not contained in the results, return the rank 1.2
   # this shall serve as an indicator in the results that the gene set
   # was not contained in the results
    return(ifelse(!is.integer0(grep(term, GSEA_results$ID)),
               rank,
               1))

 }else if(metric == "p_adj"){

 #identify row number of respective gene set
 ind_row <- grep(term, GSEA_results$ID)

 #return respective adjusted p-value
 return(ifelse(!is.integer0(ind_row),
               GSEA_results$p.adjust[ind_row],
               1))
 #note: in the case that a gene set is not reported in the results table of GSEA_results,
 #ifelse() in combination with !is.integer0() then ensures that an adjusted p-value of 1.2 is returned,
 #meaning that each adaption leading to a an adjusted p-value in (0,1] is considered an improvement
 }
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
                 pvalueCutoff = 1, # ensure that ALL gene sets are included in results, not just differentially enriched gene sets
                 exponent = exp)

 } else if(geneset_database == "GO"){

 # set seed
 set.seed(713)

 GSEA <- gseGO( gene_ranking,
                OrgDb = organism,
                ont = "BP", # subontology Molecular Function
                pvalueCutoff = 1, # ensure that All gene sets are included in results, not just differentially enriched gene sets
                eps = 0, # to ensure that all p-values are computed exactly
                seed = TRUE, # ensure that results are reproducible
                exponent = exp)
 }


 # return GSEA results
 return(GSEA %>% as.data.frame())



}




########################################################################################
########################################################################################
### (I) Minimization of a Gene Set's rank or adjusted p-value ##########################
########################################################################################
########################################################################################

#required inputs:
#(i) Gene set from KEGG gene set database
#(ii) expression data set (in this case in the Ensembl gene ID format)
#-> gene ID conversion to Entrez is performed within the function
#(iii) phenotype labels to assign the individual samples to the opposing conditions
#(iv) metric: indicate whether a genes adjusted p-value or its rank among the
#remaining gene sets shall be minimized

cP_GSEA_rankp_optim <- function(geneset,geneset_database, expression_data, phenotype_labels, metric){


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



 #create frame to document progression of the optimization workflow
 if(metric == "rank"){

 # raw documentation frame to be filled successively if the rank of a given gene set is to be optimized
  doc <- data.frame(step = c("Default", "Ranking Metric", "Pre-Filtering Threshold", "Duplicate gene ID Removal"),
                    optimal_parameter = NA,
                    rank = NA)

 }else if(metric == "p_adj"){

 # raw documentation frame to be filled successively if adjusted p-value of given gene set is to be optimized
  doc <- data.frame(step = c("Default", "Ranking Metric", "Pre-Filtering Threshold", "Duplicate gene ID Removal"),
                    optimal_parameter = NA,
                    p_adj = NA)


 }



 ###########
 # 1. step: Choose optimal manner of creating ranking metric
 ###########

 # set of options to create ranking of genes
 ranking_metrics <- c("DESeq2" ,"limma")

 # vector to store ranking metrics
 metricvalue_rankingmetrics <- c()


 ### default: create ranked list of genes using DESeq2

 # (i) run default DESeq2 workflow (incl. default pre-filtering, gene ID conversion and shrinkage)
 DESeq2_results <- pre_filt(expression_data, threshold = 10) %>%
                            geneID_conversion(dupl_removal_method = 1) %>%
                            deseq_preprocess(phenotype_labels = phenotype_labels) %>%
                            DESeq() %>% DESeq2::lfcShrink(coef = "condition_treated_vs_untreated") %>%
                            as.data.frame()

 # (ii) create ranked list of the genes (default ranking)
 ranked_list_DESeq2 <- rankedList_cP(DESeq2_results,
                                     "p_value",
                                     "DESeq2")

 # (iii) run default GSEA
 GSEA_results_DESeq2ranking <- GSEA_pipeline_default(ranked_list_DESeq2,
                                                     geneset_database,
                                                     organism = organism)

 # (iv) get metric value of given gene set
 metricvalue_rankingmetrics[1] <- pvalue_rank_GSEA(geneset,
                                              GSEA_results_DESeq2ranking,
                                              metric)

 # (v) update documentation frame
 doc[1, metric] <- metricvalue_rankingmetrics[1]


 ### alternative 2: limma (with internal parameters in default configuration)

 # (i.1) generate design matrix
 mm <- model.matrix(~phenotype_labels)

 # (i.2) generate pre-filtering indicator
 keep <- DGEList(expression_data, group = phenotype_labels) %>%
                 filterByExpr()

 # (i.3) run limma pipeline (filtering indicator keep recycled from previous alternative edgeR)
 limma_results <- geneID_conversion(expression_data[keep,], dupl_removal_method = 1) %>%
                                    DGEList(group = phenotype_labels) %>% calcNormFactors() %>%
                                    voom(design = mm) %>% lmFit(design = mm) %>%
                                    eBayes() %>% topTable(coef = ncol(mm), number = 100000)


 # (ii) generate ranking based on p-value
 ranked_list_limma <- rankedList_cP(limma_results,
                                    rankby = "p_value",
                                    method = "limma")

 # (iii) run GSEA
 GSEA_results_limma <- GSEA_pipeline_default(ranked_list_limma,
                                             geneset_database,
                                             organism = organism)

 # (iv) store metric value of given geneset
 metricvalue_rankingmetrics[2] <- pvalue_rank_GSEA(geneset,
                                              GSEA_results_limma,
                                              metric)


 ###########################################################################################################

 ###choose optimal ranking metric and update documentation frame
 ind_opt_rankmetric <- min(which(metricvalue_rankingmetrics == min(metricvalue_rankingmetrics)))
 opt_rankmetric <- ranking_metrics[ind_opt_rankmetric] #optimal manner of generating ranked gene list
 doc[2, "optimal_parameter"] <- opt_rankmetric
 doc[2, metric] <- min(metricvalue_rankingmetrics)


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


 exprdat_prefilt_list <- lapply(X = prefilt_thresh_list,
                                FUN = pre_filt,
                                expression_data = expression_data)

 GSEA_results_prefilt_list <- list()

 # run DESEq2 workflow for each of the pre-filtered gene expression data sets
 DESeq2_results_prefilt <- lapply(FUN = geneID_conversion,X = exprdat_prefilt_list, dupl_removal_method = 1) %>%
                           lapply(FUN = deseq_preprocess, phenotype_labels = phenotype_labels) %>%
                           lapply(FUN = DESeq) %>% lapply(FUN = lfcShrink,coef = "condition_treated_vs_untreated", type = "apeglm") %>%
                           lapply(FUN = as.data.frame)

 # generate ranking for each of the DESeq2 results
 ranked_lists_DESeq2 <- lapply(FUN = rankedList_cP,
                               X = DESeq2_results_prefilt,
                               rankby = "p_value",
                               method = "DESeq2")

 # run GSEA
 GSEA_results_prefilt_list <- lapply(FUN = GSEA_pipeline_default,
                                     X = ranked_lists_DESeq2,
                                     geneset_database = geneset_database,
                                     organism = organism)


 # count number of differentially enriched gene sets resulting from each pre-filtering threshold
 metricvalue_prefilt <- unlist(lapply(FUN = pvalue_rank_GSEA,
                                      X = GSEA_results_prefilt_list,
                                      term = geneset, metric = metric))

 # get index of optimal pre-filtering threshold (in case of ties, choose lower index)
 ind_opt_prefilt <- min(which(metricvalue_prefilt == min(metricvalue_prefilt)))

 # update documentation frame
 doc[3, "optimal_parameter"] <- prefilt_thresh_list[ind_opt_prefilt ]
 doc[3, metric] <- min(metricvalue_prefilt)



 } else if(opt_rankmetric == "limma"){

 # create pre-filtered gene expression data sets using functions filterByExpr() and cpm()
 # filtering thresholds

 keep1 <- DGEList(counts = expression_data, group = phenotype_labels) %>%
                  filterByExpr()
 # 1. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 1 count per million in at least 2 samples
 keep2 <- rowSums(cpm(expression_data) > 1) >= 2
 # 2. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 2 count per million in at least 3 samples


 filt_methods <- c("by FilterbyExpr()", "by cpm>1 in at least 2 samples")


 # for all methods of pre-filtering: store resulting gene expression data sets in list
 exprdat_prefilt_list <- list()
 # pre-filtered gene expression data set using keep1
 exprdat_prefilt_list[[1]] <- expression_data[keep1,]
 # pre-filtered gene expression data set using keep2
 exprdat_prefilt_list[[2]] <- expression_data[keep2,]


 limma_results_prefilt <- lapply(FUN = geneID_conversion,X = exprdat_prefilt_list, dupl_removal_method = 1) %>%
                          lapply(FUN = DGEList,group = phenotype_labels) %>%
                          lapply(FUN = calcNormFactors) %>%
                          lapply(FUN = voom,design = mm) %>%
                          lapply(FUN = lmFit,design = mm) %>%
                          lapply(FUN = eBayes) %>%
                          lapply(FUN = topTable,coef = ncol(mm), number = 100000)


 #perform ranking by p-value for each of the limma results
 ranked_lists_limma <- lapply(FUN = rankedList_cP,
                              X = limma_results_prefilt,
                              rankby = "p_value",
                              method = "limma")

 # run GSEA (GO) for each of the rankings generated with limma
 GSEA_results_prefilt_list <- lapply(FUN = GSEA_pipeline_default,
                                     X = ranked_lists_limma,
                                     geneset_database = geneset_database,
                                     organism = organism)

 # get metric value resulting from each pre-filtering threshold
 metricvalue_prefilt <- unlist(lapply(FUN = pvalue_rank_GSEA,
                                      X = GSEA_results_prefilt_list,
                                      term = geneset,
                                      metric = metric))

 # get index of optimal pre-filtering threshold (in case of ties, choose lower index)
 ind_opt_prefilt <- min(which(metricvalue_prefilt == min(metricvalue_prefilt)))


 # update documentation frame
 doc[3, "optimal_parameter"] <- filt_methods[ind_opt_prefilt ]
 doc[3, metric] <- min(metricvalue_prefilt)


 }

 # get OPTIMAL PRE-FILTERED gene expression data set
 # gene IDs are in Ensembl gene ID format as gene id conversion is performed in next step

 exprdat_prefilt_opt <- exprdat_prefilt_list[[ind_opt_prefilt]]




 ###########
 # 3. step: Choose optimal manner of duplicate gene ID removal
 ###########

 if(opt_rankmetric == "DESeq2"){

 # (i) run DESeq2 pipeline for each of the duplicate gene ID removal schemes
 deseq2_results_convID <- lapply(FUN = geneID_conversion, X = 1:2, expression_data = exprdat_prefilt_opt) %>%
                          lapply(FUN = deseq_preprocess, phenotype_labels = phenotype_labels) %>%
                          lapply(FUN = DESeq) %>%
                          lapply(FUN = lfcShrink,coef = "condition_treated_vs_untreated", type = "apeglm") %>%
                          lapply(FUN = as.data.frame)


 # (ii) generate ranked list for each of the gene expression data sets resulting from the different manners of duplicate
 # gene ID removal

 rankedLists_duplremoval <- lapply(FUN = rankedList_cP,
                                   X = deseq2_results_convID,
                                   rankby = "p_value",
                                   method = "DESeq2")

 # (iii) run default GSEA for each of the rankings resulting from the different manners of duplicate gene ID removal
 GSEA_results_convID_list <- lapply(FUN = GSEA_pipeline_default,
                                    X = rankedLists_duplremoval,
                                    geneset_database = geneset_database,
                                    organism = organism)




 }else if(opt_rankmetric == "limma"){

 # (i) run limma pipeline for each of the duplicate gene id removal schemes
 limma_results_convID <- lapply(FUN = geneID_conversion, expression_data = exprdat_prefilt_opt, X = 1:2) %>%
                         lapply(FUN = DGEList,group = phenotype_labels) %>%
                         lapply(FUN = calcNormFactors) %>% lapply(FUN = voom,design = mm) %>%
                         lapply(FUN = lmFit,design = mm) %>% lapply(FUN = eBayes) %>%
                         lapply(FUN = topTable,coef = ncol(mm), number = 100000)

 # (ii) generate ranked list for each of the gene expression data sets resulting from the different manners of duplicate gene
 # ID removal
 rankedLists_duplremoval <- lapply(FUN = rankedList_cP,
                                   X = limma_results_convID,
                                   rankby = "p_value",
                                   method = "limma")

 # (iii) run default GSEA for each of the rankings resulting from the different manners of duplicate gene ID removal
 GSEA_results_convID_list <- lapply(FUN = GSEA_pipeline_default,
                                    X = rankedLists_duplremoval,
                                    geneset_database = geneset_database,
                                    organism = organism)
 }

 # get metric value for each of the ORA results
 metricvalue_convID <- unlist(lapply(FUN = pvalue_rank_GSEA,
                                     X = GSEA_results_convID_list,
                                     term = geneset,
                                     metric = metric))

 # index of optimal duplicate gene ID removal scheme (in case of ties choose lower index)
 ind_opt_convID <- min(which(metricvalue_convID == min(metricvalue_convID)))

 # update documentation frame
 doc[4, "optimal_parameter"] <- ind_opt_convID
 doc[4, metric] <- min(metricvalue_convID)

 # get optimal gene ranking
 optim_ranking <- rankedLists_duplremoval[[ind_opt_convID]]



 ##########
 # last step: Change exponent values
 ##########


 exp_options <- c(1, 0, 1.5, 2) #possible options for exponent


 #perform GSEA with respective exponent
 GSEA_results_exp_list <- lapply(FUN = GSEA_pipeline_default,
                                 X = exp_options,
                                 gene_ranking = optim_ranking,
                                 geneset_database = geneset_database,
                                 organism = organism)


 #store metric_value in vector to compare to all exponents
 metricvalue_exp <- unlist(lapply(X = GSEA_results_exp_list,
                                  FUN = pvalue_rank_GSEA,
                                  term = geneset,
                                  metric = metric))


 #determine which exponent leads to highest number of DEGS
 #in case of tie choose lower index
 ind_exp_opt <- min(which(metricvalue_exp == min(metricvalue_exp)))

 #update exponent
 exp_opt <- exp_options[ind_exp_opt]
 # extend documentation frame by row corresponding to exponent choice
 doc[nrow(doc) + 1, "step"] <- "Exponent"
 doc[nrow(doc), "optimal_parameter"] <- exp_opt
 doc[nrow(doc), metric] <- min(metricvalue_exp)


 return(list(default = GSEA_results_DESeq2ranking, #default results
 optim = as.data.frame(GSEA_results_exp_list[[ind_exp_opt]]), #optimal results
 documentation = doc))


}


