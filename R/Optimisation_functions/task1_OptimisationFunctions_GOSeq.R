#Investigation of Potential for over-optimistic results w.r.t. GOSeq

#packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)
library(edgeR)
library(goseq)

#introduction of  necessary functions to perform optimization of parameters

######################################################################################
###lossfunction for GOSeq ############################################################
######################################################################################

#generate lossfunction to count the number of differentially enriched gene sets
lossfunction_GOSeq <- function(gsa_results){

  #count number of differentially expressed gene sets which are over-represented
  #according to the adjusted p-value of over-representation
  n_DEGS <- sum(gsa_results$p_adj_overrep < 0.05)


  return(n_DEGS)

}



#######################################################################################
###generate necessary input for GOSeq #################################################
#######################################################################################
GOSeq_input_preparation <- function(DE_results){

  # prepare required input vector for DE results generated with DESeq2 or limma


  #remove all genes with adjusted p-value set to NA in DE analysis -> especially concerns DE results generated with DESeq2
  DE_results_nona <- DE_results[!is.na(DE_results$p_adj),]
  #create named binary vector for all genes with adjusted p-value NOT set to NA
  DEG_vec_bin <- ifelse((DE_results_nona$p_adj < 0.05) & (!is.na(DE_results_nona$p_adj)), 1,0)
  names(DEG_vec_bin)<-rownames(DE_results_nona)

  # return binary names vector
  return(DEG_vec_bin)

}


########################################################################################
### GOSeq pipeline #####################################################################
########################################################################################

GOSeq_pipeline <- function(DE_results, geneset_database = "GO:BP", bias = NULL,
                           extended_universe = FALSE, calc_method = "Wallenius", geneID, organism){


  # create input vector of required format from the DESeq2/ limma results
  input_vector <- GOSeq_input_preparation(DE_results)


  GOSeq_results <- nullp(DEgenes = input_vector,
                         genome = organism, id = geneID, plot.fit = FALSE,
                         bias.data = bias) %>%
    goseq(genome = organism, id = geneID, test.cats = geneset_database,
          use_genes_without_cat = FALSE, method = calc_method)  %>%
    as.data.frame() %>%
    mutate(p_adj_overrep = p.adjust(over_represented_pvalue)) # external conduct of multiple test adjustment


  return(GOSeq_results)


}




#######################################################################################
### (I) optimization of preparatory parameters ########################################
#######################################################################################
# optimization of input list of differentially expressed genes to maximize number of
# differentially enriched gene sets




#######################################################################################
###optimization of pre-processing #####################################################
#######################################################################################

GOSeq_preprocess_optim <- function(expression_data, phenotype_labels){


  ###### generate documentation data frame to document development of number of differentially enriched gene sets
  doc <- data.frame(step = c("Default", "Differential Expression Technique", "Pre-Filtering Threshold"),
                    optimal_parameter = NA,
                    n_DEGS=NA)

  # required for GOSeq workflow: get identifier of gene IDs
  gene_ID <- ifelse(grepl("ENS",
                          rownames(expression_data)[1]),
                    "ensGene",
                    "knownGene")


  ### identify organism whose gene expression was measured in experiment
  # human (homo sapiens): ENSEMBL gene ID starts with "ENSG"  -> identified by GOSeq as "hg19"
  # mouse (mus musculus): "ENSEMBL gene ID starts with "ENSMUSG" -> identifued by GOSeq as "mm9"

  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl,
                         X = c("ENSG", "ENSMUSG"),
                         x = rownames(expression_data)[1])
  # choose suitable organism
  organism <- c("hg19", "mm9")[ind_organism][[1]]


  #############
  ##### 1. step: choose optimal differential expression technique
  #############

  # choose between DESeq2 amd limma and to classify genes as DE/ not DE
  de_techniques <- c("DESeq2", "limma")
  n_DEGS_detechniques <- c()

  #######################
  ### default: DESeq2 ###
  #######################

  # (i) perform differential expression analysis with default parameters
  deseq2_results_default <- pre_filt(expression_data, threshold=10) %>%
    deseq_preprocess(phenotype_labels = phenotype_labels) %>%
    DESeq() %>% results() %>%
    as.data.frame() %>% rename( p_adj = padj)

  # (ii) run GOSeq
  GOSeq_results_deseq2 <- GOSeq_pipeline(deseq2_results_default,
                                         geneID = gene_ID,
                                         organism = organism)

  # (iii) update documentation frame
  doc[1,"n_DEGS"] <- lossfunction_GOSeq(GOSeq_results_deseq2)

  # count number of differentially enriched gene sets
  n_DEGS_detechniques[1] <- lossfunction_GOSeq(GOSeq_results_deseq2)



  ############################
  ### alternative 2: limma ###
  ############################

  # (i.2) perform pre-filtering using builtin function filterByExpr()
  # -> results in indication on which genes are included and which are omitted from analysis
  keep <- DGEList(expression_data, group = phenotype_labels) %>%
    filterByExpr()

  # (i.2) generate design matrix
  mm <- model.matrix( ~ phenotype_labels)

  # (i.3) run standard limma workflow
  # for pre-filtering: reuse indicator keep from preceding alternative edgeR
  limma_results_default <- DGEList(expression_data[keep,], group=phenotype_labels) %>%
    calcNormFactors() %>% voom(design=mm) %>% lmFit(design=mm) %>%
    eBayes() %>% topTable(coef=ncol(mm), number=100000) %>%
    as.data.frame() %>% rename(p_adj=adj.P.Val)

  # (ii) run GOSeq
  GOSeq_results_limma <- GOSeq_pipeline(limma_results_default,
                                        geneID = gene_ID,
                                        organism = organism)

  # (iii) count number of differentially enriched gene sets
  n_DEGS_detechniques[2] <- lossfunction_GOSeq(GOSeq_results_limma)



  #####################################################################################

  # get index of optimal DE technique (in case of ties choose lower index)
  ind_opt_DEtechnique <- min(which(n_DEGS_detechniques == max(n_DEGS_detechniques)))
  #choose optimal gene expression technique
  opt_DE_technique <- de_techniques[ind_opt_DEtechnique]


  # update documentation frame
  doc[2, "optimal_parameter"] <- opt_DE_technique
  doc[2,"n_DEGS"] <- max( n_DEGS_detechniques )





  ##############
  ##### 2. step: choose optimal pre-filtering threshold
  #############
  #provide default result to documentation frame

  # note: in this step, different options of pre-filtering are used for the different DE techniques (in accordance with the
  # corresponding user manuals)
  # for DESeq2, manual pre-filtering is performed using different filtering thresholds
  # for limma, pre-filtering is performed using (i) filterByExpr (default) and (ii) cpm-transformed values


  # 1. case: DESeq2 is optimal DE technique

  if(opt_DE_technique == "DESeq2"){

    # alternative pre-filtering thresholds of X read counts across all samples
    filt_threshold_alt <- c(10, 50)
    # perform pre-filtering according to the alternative thresholds
    exprdat_list_prefilt <- lapply(X=filt_threshold_alt,
                                   FUN=pre_filt,
                                   expression_data=expression_data)


    # (i) perform differential expression analysis with default parameters

    DE_results_prefilt <- lapply(FUN=deseq_preprocess, X= exprdat_list_prefilt, phenotype_labels=phenotype_labels) %>%
      lapply(FUN=DESeq) %>% lapply(FUN=results) %>%
      lapply(as.data.frame) %>% lapply(rename, p_adj=padj)

    #(ii) run GOSeq with default parameters for each of the pre-filtered gene expression data sets
    GOSeq_results_prefilt <- lapply(FUN = GOSeq_pipeline,
                                    X = DE_results_prefilt,
                                    geneID= gene_ID,
                                    organism = organism)


    #count number of differentially enriched gene sets resulting from the alternative prefiltering thresholds
    n_DEGS_prefilt <- unlist(lapply(X= GOSeq_results_prefilt,
                                    FUN=lossfunction_GOSeq))

    # get index of optimal pre-filtering threshold (in case of ties choose lower index)
    ind_opt_prefilt <- min(which( n_DEGS_prefilt == max( n_DEGS_prefilt)))
    # get optimal pre-filtering threshold
    opt_prefilt <-  filt_threshold_alt[ind_opt_prefilt]

    # get current optimal number of differentially enriched gene sets
    n_DEGS <- max(n_DEGS_prefilt)

    # update documentation frame
    doc[3, "optimal_parameter"] <- opt_prefilt
    doc[3, "n_DEGS"] <-  max(n_DEGS_prefilt)





  } else if(opt_DE_technique == "limma"){ # for edgeR and limma, pre-filtering is performed differently, namely using filterByExpr (default) OR using cpm()

    # DEFAULT pre-filtering using filterByExpr
    keep1 <- DGEList(counts=expression_data, group=phenotype_labels) %>% filterByExpr()
    # 1. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 1 count per million in at least 2 samples
    keep2 <- rowSums(cpm(expression_data) > 1) >=2



    filt_methods <- c("by FilterbyExpr()", "by cpm>1 in at least 2 samples")


    # for all methods of pre-filtering: store resulting gene expression data sets in list
    exprdat_list_prefilt <- list()
    # pre-filtered gene expression data set using keep1
    exprdat_list_prefilt[[1]] <- expression_data[keep1,]
    # pre-filtered gene expression data set using keep2
    exprdat_list_prefilt[[2]] <- expression_data[keep2,]




    DE_results_prefilt <- lapply(FUN=DGEList,X= exprdat_list_prefilt,  group=phenotype_labels) %>%
      lapply(FUN=calcNormFactors) %>% lapply(FUN=voom,design=mm) %>%
      lapply(FUN=lmFit,design=mm) %>% lapply(FUN=eBayes) %>%
      lapply(FUN=topTable, coef=ncol(mm), number=100000) %>%
      lapply(FUN=as.data.frame) %>% lapply(FUN=rename,p_adj=adj.P.Val)


    # perform GOSeq for each of the pre-filtered gene expression data sets
    GOSeq_results_prefilt<- lapply(FUN = GOSeq_pipeline,
                                   X = DE_results_prefilt,
                                   geneID = gene_ID,
                                   organism = organism)


    # store number of differentially enriched gene sets for each GOSeq result
    n_DEGS_prefilt <- unlist(lapply(X= GOSeq_results_prefilt, FUN=lossfunction_GOSeq))

    # get index optimal pre-filtering method (in case of tie, choose default filterByExpr)
    ind_opt_prefilt <- min(which(n_DEGS_prefilt == max(n_DEGS_prefilt)))
    # get optimal pre-filtering method
    opt_prefilt <-  filt_methods[ind_opt_prefilt]

    # current optimal number of differentially enriched gene sets
    n_DEGS <- max( n_DEGS_prefilt )


    # update documentation frame
    doc[3, "optimal_parameter"] <- opt_prefilt
    doc[3, "n_DEGS"] <- max(n_DEGS_prefilt)

  }


  ### optimally pre-filtered gene expression data set
  exprdat_prefilt_opt <- exprdat_list_prefilt[[ind_opt_prefilt]]


  # optimal differential expression analysis results
  optim_DE_results <- DE_results_prefilt[[ind_opt_prefilt]]

  # optimal GOSeq results
  optim_GOSeq_results <- GOSeq_pipeline(optim_DE_results,
                                        geneID = gene_ID,
                                        organism = organism)



  #return final results
  return(list(default_GOSeq = GOSeq_results_deseq2[GOSeq_results_deseq2$p_adj_overrep < 0.05, ], #default results table
              optim_GOSeq = optim_GOSeq_results , #optimal results table
              optim_DE_results = optim_DE_results, # optimal differential expression analysis results
              documentation = doc)) #documentation frame

}



#######################################################################################
### (II) Optimization of internal parameters ##########################################
#######################################################################################



##########################################################################################
###global Goseq optimization##############################################################
##########################################################################################
# stepwise analysis: from GO, KEGG and MKEGG, choose gene set database which yields highest
# number of DEGS with DEFAULT parameters
# -> choose method yielding highest number of DEGS and perform optimization with
# respective function

GOSeq_optim <- function(DE_results, expression_data){

  # obtain default input vector for GOSeq in required format
  # obtain optimized input vector for GOSeq in required format


  ### identify organism whose gene expression was measured in experiment
  # human (homo sapiens): ENSEMBL gene ID starts with "ENSG"  -> identified by GOSeq as "hg19"
  # mouse (mus musculus): ENSEMBL gene ID starts with "ENSMUSG" -> identified by GOSeq as "mm9"

  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl,
                         X = c("ENSG", "ENSMUSG"),
                         x = rownames(DE_results)[1])

  # choose suitable organism
  organism <- c("hg19", "mm9")[ind_organism][[1]]


  # from the differential expression analysis results provided as an argument,
  # generate the required input object for GOSeq
  DEG_vec_bin <- GOSeq_input_preparation(DE_results)


  # determine gene Identifier
  # If gene IDs begin with "ENS" -> set gene identifier to ensGene (ENSEMBL IDs)
  # else: set gene identifier to ENTREZ gene ID
  gene_ID <- ifelse(all(grepl("ENS", rownames(DE_results))), "ensGene", "knownGene")
  #-> this parameter remains untouched throughout the following analysis


  # prepare dataframe to document stepwise optimization
  doc<-data.frame(step=c("Default", "Bias", "Gene Set Database","Universe: Use Genes w/o Gene Set", "p-Value Calculation"),
                  optimal_parameter=NA,
                  n_DEGS=NA)

  #################################
  ### step 1: Default GOSeq results
  #################################

  GOSeq_default<- GOSeq_pipeline(DE_results,
                                 geneset_database = "GO:BP",
                                 geneID = gene_ID,
                                 organism = organism)


  #count number of differentially enriched gene sets
  n_DEGS <- lossfunction_GOSeq(GOSeq_default)
  doc[1,"n_DEGS"] <- lossfunction_GOSeq(GOSeq_default)


  #######################
  ### step 2: change bias
  #######################

  # account for all possible biases present in the data (i.e. read count bias),
  # not just length bias
  # calculate read count bias for all genes present in input vector DEG_vec_bin
  countbias <- rowSums(expression_data[rownames(expression_data) %in% names(DEG_vec_bin),])
  # note: non of the row sums are 0 since lowly expressed gened
  # were filtered out in "Step1_Preprocessing"

  #supply read count bias to the probability weighting function
  GOSeq_allbias <- GOSeq_pipeline(DE_results,
                                  geneID = gene_ID,
                                  bias = countbias, # supply read count bias
                                  organism = organism)

  #set current pwf_results
  #if considering all possible biases (countbias) leads to higher number of differentially
  #enriched gene sets then update the pwf function from now on

  if(lossfunction_GOSeq(GOSeq_allbias) > n_DEGS){#case1: accounting for all possible
    #biases leads to higher number of DEGs
    bias <- countbias
    #count current number of DEGS
    n_DEGS <- lossfunction_GOSeq(GOSeq_allbias)
    #define parameter which will indicate in the results whether bias has been updated
    doc[2,"optimal_parameter"] <- "All Biases"
    doc[2, "n_DEGS"] <- lossfunction_GOSeq(GOSeq_allbias)



  } else { #case2: accounting for all possible biases does not lead to higher
    #number of DEGs -> keep default version

    bias <- NULL
    #define parameter which will indicate in the results whether bias has been updated
    doc[2,"optimal_parameter"] <-"Length Bias Only"
    doc[2, "n_DEGS"] <- n_DEGS

  }

  ############################################
  ### step 3: change gene set database to KEGG
  ############################################

  GOSeq_kegg <- GOSeq_pipeline(DE_results,
                               geneset_database = "KEGG", # change gene set database to KEGG
                               geneID = gene_ID,
                               bias = countbias,
                               organism = organism)


  #update gene set database in case of an increase in the number of DEGS
  gene_set_db <- ifelse(lossfunction_GOSeq(GOSeq_kegg) > n_DEGS, "KEGG", "GO:BP")
  #update current number of differentiallys enriched gene sets
  n_DEGS <- max(lossfunction_GOSeq(GOSeq_kegg), n_DEGS)
  #for documentation
  doc[3,"optimal_parameter"] <- gene_set_db
  doc[3, "n_DEGS"] <- n_DEGS




  #############################
  ### step 4: Change Background
  #############################

  # change background: genes without category shall now be counted towards the total number of
  # genes outside a category (are otherwise excluded from analysis by default)

  GOSeq_genesnocat<- GOSeq_pipeline(DE_results,
                                    geneset_database = gene_set_db,
                                    bias = bias,
                                    geneID = gene_ID,
                                    extended_universe = TRUE, #change background
                                    organism = organism)



  #update gene set database in case of an increase in the number of DEGS
  genes_wio_cat <- ifelse(lossfunction_GOSeq(GOSeq_genesnocat) > n_DEGS, TRUE, FALSE)
  #update current number of differentially enriched gene sets
  n_DEGS <- max(lossfunction_GOSeq(GOSeq_genesnocat), n_DEGS)
  #for documentation
  doc[4,"optimal_parameter"] <- genes_wio_cat
  doc[4, "n_DEGS"] <- n_DEGS


  ##########################################
  ### step 5: Change calculation method ####
  ##########################################

  #change calculation method to Random Sampling (which is the exact method)
  # argument repcnt (number of random samples) is kept at its default value
  GOSeq_calcmethod <- GOSeq_pipeline(DE_results,
                                     geneset_database = gene_set_db,
                                     bias = bias,
                                     geneID = gene_ID,
                                     extended_universe = genes_wio_cat,
                                     calc_method = "Sampling", # change calculation method to random sampling
                                     organism = organism)



  #update calculation method in case of an increased number of differentially
  #enriched gene sets
  method <- ifelse(lossfunction_GOSeq(GOSeq_calcmethod) > n_DEGS, "Sampling", "Wallenius")
  n_DEGS <- max(lossfunction_GOSeq(GOSeq_calcmethod), n_DEGS)
  #for documentation
  doc[5,"optimal_parameter"] <- method
  doc[5, "n_DEGS"] <- n_DEGS


  # run final GOSeq analysis with optimal parameter configuration ( multiple test adjustment
  # method is optimized in the very last step)

  GOSeq_final <- GOSeq_pipeline(DE_results,
                                geneset_database = gene_set_db,
                                bias = bias,
                                geneID = gene_ID,
                                extended_universe = genes_wio_cat,
                                calc_method = method,
                                organism = organism)





  #return optimal result with only those gene sets detected as differetially enriched
  return(list(default_GOSeq =GOSeq_default[GOSeq_default$p_adj_overrep< 0.05,], # default GOSeq results
              optim_GOSeq =GOSeq_final[GOSeq_final$p_adj_overrep< 0.05,], # optimal GOSeq results
              documentation=doc))

}


################################################################################
### (III) Joint optimization of Pre-Processing and Internal Parameters #########
################################################################################

GOSeq_joint_optimization <- function(expression_data, phenotype_labels){


  # (i) optimize binary input list of DE vs. non-DE to maximize number of differentially enriched gene sets
  optim_preprocess <- GOSeq_preprocess_optim(expression_data, phenotype_labels)

  # (ii) based on optimal binary input list from (i), optimize internal parameters of GOSeq to further maximize
  # number of differentially enriched gene sets
  optim_internalparam <- GOSeq_optim(optim_preprocess$optim_DE_results,
                                     expression_data)

  # for control purposes
  doc <- rbind(optim_preprocess$documentation, optim_internalparam$documentation[optim_internalparam$documentation$step != "Default",])



  return(list(default_GOSeq = optim_preprocess$default_GOSeq, # default GOSeq results
              optim_GOSeq = optim_internalparam$optim_GOSeq, # optimal GOSeq results
              documentation = doc)) # documentation frame

}

