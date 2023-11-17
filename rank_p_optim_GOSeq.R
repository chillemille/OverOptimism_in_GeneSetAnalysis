#### GOSeql: minimization of a given gene set's 
#(i) (adjusted) p-value / 
#(ii) rank among the remaining gene sets 

library(clusterProfiler)
library(dplyr)
library(DESeq2)
library(edgeR)
library(org.Hs.eg.db)
library(rlang)
library(goseq)

# load gene expression data set with true phenotype randomly permuted phenotype assignments 
source("./Random_Phenotype_Permutations.R")

# load required pre-processing functions
source("./PreProcessing_Functions.R")




################################################################################
### get gene's adjusted p-value or rank (based on what is needed) ##############
################################################################################

###function used to detect whether gene set is not containted in the GOSeq 
# results 
is.integer0  <-  function(x){
  is.integer(x) && length(x) == 0L
}
##############


#in the case that a gene's
#(i) (adjusted) p-value is needed: enter "p_adj" in argument metric 
#(ii) relative rank among the remaining gene sets is needed: enter "rank" in argument metric

#note: function works for GO as well as KEGG as results table 
#contains relevant information automatically
#argument term must be in the form of a GO ID resp. KEGG ID 

# note that the GOSeq results are automatically ordered by over-representation
# so that we do not have to perform any further ordering
pvalue_rank_goseq <- function(term, goseq_results, metric){
  
  #remove the letters hsa from the gene set ID
  term  <-  gsub("hsa", "", term)
  
  if(metric == "rank"){
    
 
    if(nrow(goseq_results) == 0 ){
      # case 1: GOSeq results table is empty
      # -> return rank 1 (i.e. worst rank as no rank can be assigned to the given gene set 
      
      return(1)
      
    } else{ # case: GOSeq results table is NOT empty 
    
    #return row number of respective gene set
    return(ifelse(!is.integer0(grep(term, goseq_results$category)), 
                  grep(term, goseq_results$category)/nrow(goseq_results), 
                  1))
    #note: in the case that a gene set is not reported in the results table of goseq_results,
    #ifelse() in combination with !is.integer0() then ensures that a rank of 1 is returned,
    #meaning that each adaption leading to any infinite rank is considered an improvement
    
    }
    
  }else if(metric == "p_adj"){
    
    if(nrow(goseq_results) == 0) {
      # case 1: GOSeq results table is empty 
      # -> return adjusted p-value 1 as no adjusted p-value can be assigned to the given gene set 
      
      return(1)
      
    } else{ # case 2: GOSeq results table is NOT empty 
    
    #identify row number of respective gene set 
    ind_row <- grep(term, goseq_results$category)
    
    #return respective adjusted p-value
    return(ifelse(!is.integer0(ind_row),goseq_results$p_adj_overrep[ind_row], 1))
    #note: in the case that a gene set is not reported in the results table of goseq_results, 
    #ifelse() in combination with !is.integer0() then ensures that an adjusted p-value/rank of 1 is returned,
    #meaning that each adaption leading to a an adjusted p-value/rank in (0,1) is considered an improvement
    }
  }
}


#######################################################################################
###generate necessary input for GOSeq #################################################
#######################################################################################
goseq_input_preparation  <-  function(DE_results){
  
  # prepare required input vector for DE results generated with DESeq2, or limma (default)

    
    #remove all genes with adjusted p-value set to NA in DE analysis -> especially concerns DE results generated with DESeq2 
    DE_results_nona  <-  DE_results[!is.na(DE_results$p_adj),]
    #create named binary vector for all genes with adjusted p-value NOT set to NA
    DEG_vec_bin  <-  ifelse((DE_results_nona$p_adj < 0.05) & (!is.na(DE_results_nona$p_adj)), 1,0)
    names(DEG_vec_bin) <- rownames(DE_results_nona)
    

  
  return(DEG_vec_bin)
  
}

################################################################################
### default GOSeq pipeline #####################################################
################################################################################

GOSeq_pipeline  <-  function(DE_results, geneset_database, calc_method = "Wallenius", 
                             bias = NULL, genes_wocat = FALSE, geneID, organism){
  
  # create input vector of required format of DESeq2/ limma results 
  input_vector  <-  goseq_input_preparation(DE_results)
    
  # run GOSeq pipeline with input vector 
  GOSeq_results  <-  nullp(DEgenes = input_vector, 
                           genome = organism, 
                           id = geneID, 
                           plot.fit = FALSE, 
                           bias.data = bias) %>% 
                      goseq(genome = organism, 
                            id = geneID, 
                            test.cats = geneset_database, 
                            method = calc_method, 
                            use_genes_without_cat = genes_wocat)  %>% 
                      as.data.frame() %>%
                      mutate(p_adj_overrep = p.adjust(over_represented_pvalue)) # external conduct of multiple test adjustment

  
  return(GOSeq_results)
  
  
}


################################################################################
### optimization of adjusted p-value or rank (for KEGG or GO) ##################
################################################################################

#minimize either rank or adjusted p-value of the gene set of choice 
#for all runs of enrichGO, we set p-value cutoff and q-value cutoff values to 1
#to ensure that ALL gene sets (whether significantly differentially enriched or not)
#are reported

goseq_rank_pvalue_optim  <-  function(geneset, geneset_database, metric,expression_data, 
                                      phenotype_labels){
  

  #default parameters
  #pre-filtering threshold: 10 genes across all samples
  #DE technique: DESeq2 (with all parameters in default configuration)
  
  #NOTE: As GOSeq works with Ensembl as gene ID format, a conversion of the gene IDs and resulting
  #necessary removal of duplicated gene IDs is not necessary
  
  
  #identification of gene ID format 
  gene_ID  <-  ifelse(grepl("ENS", 
                            rownames(expression_data)[1]), 
                            "ensGene", 
                            "knownGene")
  
  
  
  ### identify organism whose gene expression was measured in experiment
  # human (homo sapiens): ENSEMBL gene ID starts with "ENSG"  -> identified by GOSeq as "hg19"
  # mouse (mus musculus): "ENSEMBL gene ID starts with "ENSMUSG" -> identifued by GOSeq as "mm9"
  
  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism  <-  sapply(FUN = grepl, 
                           X = c("ENSG", "ENSMUSG"), 
                           x = rownames(expression_data)[1])
  # choose suitable organism (required for function bitr)
  organism  <-  c("hg19", "mm9")[ind_organism][[1]]
  
  
  ###### generate documentation data frame to document development 
  #we want the column names of the documentation frame to depend on the 
  #metric value to optimize (i.e. rank vs. adjusted p-value)
  if(metric == "rank"){
    doc <- data.frame(step = c("Default", "Differential Expression Technique", "Pre-filtering Threshold"),
                      optimal_parameter = NA, 
                      rank = NA)

    
    
  }else if(metric == "p_adj"){
    doc <- data.frame(step = c("Default", "Differential Expression Technique", "Pre-filtering Threshold"),
                    optimal_parameter = NA, 
                    p_adj = NA)
    
    
  }
  
  
  ##############
  #### 1.step: choose optimal DE technique
  #############
  
  # differential expression techniques 
  DE_techniques  <-  c("DESeq2", "limma")
  metricvalue_detechniques  <-  c()
  
  
  #############################################
  ### default DE technique: DESeq2 ############
  #############################################
  
  # (i) default DESeq2 pipeline (using default pre-filtering) 
  deseq2_results  <-  pre_filt(expression_data, threshold = 10) %>% 
                      deseq_preprocess(phenotype_labels = phenotype_labels) %>% 
                      DESeq() %>% results() %>% 
                      as.data.frame() %>% dplyr::rename(p_adj = padj)
  
  # (ii) default GOSeq 
  goseq_results_deseq2  <-  GOSeq_pipeline(deseq2_results, 
                                           geneset_database = geneset_database, 
                                           geneID = gene_ID, 
                                           organism = organism)
  
  # store metric value for comparison
  metricvalue_detechniques[1]  <-  pvalue_rank_goseq(term = geneset, 
                                                     goseq_results_deseq2, 
                                                     metric)
  
  
  # update documentation frame for default results
  doc[1, metric]  <- metricvalue_detechniques[1]
  
 
  
  ############################
  ### alternative 2: limma ###
  ############################
  
  
  #(i.1) generate pre-filtering indicator using edgeR's builtin function filterByExpr()
  keep  <-  DGEList(expression_data, group = phenotype_labels) %>% 
            filterByExpr()
  
  # (i.2) generate design matrix 
  mm  <-  model.matrix(~phenotype_labels)
  
  
  # (i.3) run default limma pipeline (reuse filtering indicator keep from previous alternative edgeR)
  limma_results  <-  DGEList(counts = expression_data[keep,], group = phenotype_labels) %>% 
                     calcNormFactors() %>%
                     voom(design = mm) %>% lmFit(design = mm) %>% 
                     eBayes() %>% topTable(coef = ncol(mm), number = 100000) %>%
                     as.data.frame() %>% dplyr::rename(p_adj = adj.P.Val)
  
  # (ii) run GOseq
  goseq_results_limma  <-  GOSeq_pipeline(limma_results, 
                                          geneset_database, 
                                          geneID = gene_ID, 
                                          organism = organism)
  
  # store metric value for comparison
  metricvalue_detechniques[2]  <-  pvalue_rank_goseq(term = geneset, 
                                                     goseq_results_limma, 
                                                     metric)
  
  ###########################################################################################
  
  
  # get index of optimal DE technique (choose lower index in case of ties)
  ind_opt_detechnique  <-  min(which( metricvalue_detechniques == min(metricvalue_detechniques)))
  # get optimal DE technique 
  opt_DE_technique  <-  DE_techniques[ind_opt_detechnique]
  
  
  #update documentation frame
  doc[2, "optimal_parameter"]  <-  opt_DE_technique
  doc[2, metric]  <-  min(metricvalue_detechniques)
  
  
  
  ##############
  ##### 2. step: choose optimal pre-filtering threshold 
  #############
  #provide default result to documentation frame
  
  # note: in this step, different options of pre-filtering are used for the different DE techniques (in accordance with the 
  # corresponding user manuals)
  # for DESeq2, manual pre-filtering is performed using different filtering thresholds 
  # for edgeR and limma, pre-filtering is performed using (i) filterByExpr (default) and (ii) cpm-transformed values 
  
  
  # 1. case: DESeq2 is optimal DE technique
  
  if(opt_DE_technique == "DESeq2"){
    
    # alternative pre-filtering thresholds of X read counts across all samples 
    filt_threshold_alt <- c(10,50)
    # perform pre-filtering according to the alternative thresholds 
    exprdat_list_prefilt <- lapply(X = filt_threshold_alt, 
                                   FUN = pre_filt, 
                                   expression_data = expression_data)
    
    
    # (i) perform differential expression analysis with default parameters
    
    DE_results_prefilt <- lapply(FUN = deseq_preprocess, 
                                 X = exprdat_list_prefilt, 
                                 phenotype_labels = phenotype_labels) %>%
                          lapply(FUN = DESeq) %>% lapply(FUN = results) %>% 
                          lapply(as.data.frame) %>% lapply(dplyr::rename, p_adj = padj)
    
    #(ii) run GOSeq with default parameters for each of the pre-filtered gene expression data sets 
    GOSeq_results_prefilt <-  lapply(FUN = GOSeq_pipeline,
                                     X = DE_results_prefilt, 
                                     geneset_database = geneset_database, 
                                     geneID = gene_ID, 
                                     organism = organism)
    
    #count number of differentially enriched gene sets resulting from the alternative prefiltering thresholds
   metricvalue_prefilt <- unlist(lapply(X = GOSeq_results_prefilt, 
                                        FUN = pvalue_rank_goseq, 
                                        term = geneset, 
                                        metric = metric))
    
    # get index of optimal pre-filtering threshold  (in case of ties choose lower index)
    ind_opt_prefilt  <-  min(which(metricvalue_prefilt == min(metricvalue_prefilt)))
    # get optimal pre-filtering threshold 
    opt_prefilt  <-   filt_threshold_alt[ind_opt_prefilt]
    
    # get current optimal number of differentially enriched gene sets 
   metricvalue  <-  min(metricvalue_prefilt)
    
    # update documentation frame 
    doc[3, "optimal_parameter"]  <- opt_prefilt
    doc[3, metric]  <-  metricvalue
    
    ### get optimal pre-filtered gene expression data set 
    
    exprdat_prefilt_opt  <- exprdat_list_prefilt[[ind_opt_prefilt]]
    
    
    
    
  } else if(opt_DE_technique  ==  "limma") { # for edgeR and limma, pre-filtering is performed differently, namely using filterByExpr (default) OR using cpm()
    
    # DEFAULT pre-filtering using filterByExpr
    keep1  <-  DGEList(counts = expression_data, 
                       group = phenotype_labels) %>% 
               filterByExpr()
    
    # 1. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 1 count per million in at least 2 samples
    keep2  <-  rowSums(cpm(expression_data) > 1) >= 2 
    # 2. ALTERNATIVE pre-filtering using cpm: keep those genes that have at least 2 count per million in at least 3 samples 
    
    
    filt_methods  <-  c("by FilterbyExpr()", "by cpm>1 in at least 2 samples")
    
    
    # for all methods of pre-filtering: store resulting gene expression data sets in list 
    exprdat_list_prefilt  <-  list()
    # pre-filtered gene expression data set using keep1 
    exprdat_list_prefilt[[1]]  <-   expression_data[keep1,]
    # pre-filtered gene expression data set using keep2
    exprdat_list_prefilt[[2]]  <-   expression_data[keep2,]


      # run limma for each pre-filtered gene expression data set 
     DE_results_prefilt  <-  lapply(FUN = DGEList,
                                    X = exprdat_list_prefilt,  
                                    group = phenotype_labels) %>% 
                             lapply(FUN = calcNormFactors) %>% lapply(FUN = voom, design = mm) %>% 
                             lapply(FUN = lmFit, design = mm) %>% lapply(FUN = eBayes) %>% 
                             lapply(FUN = topTable, coef = ncol(mm), number=100000) %>%
                             lapply(FUN = as.data.frame) %>% lapply(FUN = dplyr::rename, p_adj = adj.P.Val)
      
      
      # perform GOSeq for each of the pre-filtered gene expression data sets
      GOSeq_results_prefilt <-  lapply(FUN = GOSeq_pipeline, 
                                       X = DE_results_prefilt, 
                                       geneset_database = geneset_database, 
                                       geneID = gene_ID, 
                                       organism = organism)
      
      metricvalue_prefilt  <-  unlist(lapply(X = GOSeq_results_prefilt, 
                                             FUN = pvalue_rank_goseq, 
                                             term = geneset, 
                                             metric = metric))
      
      # get index optimal pre-filtering method (in case of tie, choose default filterByExpr)
      ind_opt_prefilt  <-  min(which(metricvalue_prefilt == min( metricvalue_prefilt)))
      # get optimal pre-filtering method 
      opt_prefilt  <-   filt_methods[ind_opt_prefilt ]
      
      # current optimal number of differentially enriched gene sets 
      metricvalue  <-  min(metricvalue_prefilt)
      

    
    metricvalue_prefilt  <-  unlist(lapply(X = GOSeq_results_prefilt, 
                                           FUN = pvalue_rank_goseq, 
                                           metric = metric, 
                                           term = geneset))
    
    # get index optimal pre-filtering method (in case of tie, choose default filterByExpr)
    ind_opt_prefilt  <-  min(which(metricvalue_prefilt == min(metricvalue_prefilt)))
    # get optimal pre-filtering method 
    opt_prefilt  <-   filt_methods[ind_opt_prefilt ]
    
    # current optimal number of differentially enriched gene sets 
    metricvalue  <-  min(metricvalue_prefilt)
    

    # update documentation frame 
    doc[3, "optimal_parameter"]  <-  opt_prefilt
    doc[3, metric]  <-   min(metricvalue_prefilt)
    
    ### get optimal pre-filtered gene expression data set 
    
    exprdat_prefilt_opt  <-  exprdat_list_prefilt[[ind_opt_prefilt]]
    
  } 

  

 
  
  # set optimal DE results as those resulting from the optimal manner of duplicate gene ID removal 
  optim_DE  <-  DE_results_prefilt[[ind_opt_prefilt]]
  

  
  
  ############################################
  ### optimize internal parameters of GOSeq
  ############################################
  
  # get binary vector of ALL genes whose expression was measured
  DEG_vec_bin  <-  goseq_input_preparation(optim_DE)
  
  
  # generate documentation frame for internal parameters of GOSeq
  # this will be later on merged with the overall documentation frame doc 
  doc_intern  <-  data.frame(step = c("Bias", "Universe", "p-Value Calculation"), 
                             optimal_parameter = NA)
  
  
  ############
  ### 4. step: optimize bias to be accounted for 
  ############
  
  #account for all possible biases present in the data (i.e. read count bias),
  #not just length bias
  #calculate read count bias for all genes present in input vector DEG_vec_bin
  countbias <- rowSums(expression_data[rownames(expression_data) %in% names(DEG_vec_bin),])
  #note: non of the row sums are 0 since lowly expressed gened 
  #were filtered out in "Step1_Preprocessing"
  
  # run GOSeq and considering all biases 
  GOSeq_allbias  <-  GOSeq_pipeline(optim_DE, 
                                    bias = countbias, 
                                    geneID = gene_ID, 
                                    geneset_database = geneset_database, 
                                    organism = organism)
 
  
  # depending on whether accounting for ALL biases leads to a decrease in the metric value, update the bias to account for
  if(pvalue_rank_goseq(geneset, GOSeq_allbias, metric) < metricvalue){ #case 1: accounting for all biases leads to decreases metric
    
    # set bias to be accounted for as all biases in the gene expression data set 
    bias  <-  countbias
    
    # update documentation frame 
    doc_intern[1, "optimal_parameter"]  <-  "All Biases"
    doc_intern[1, metric]  <-  pvalue_rank_goseq(geneset, 
                                                 GOSeq_allbias, metric)
    
    metricvalue  <-  pvalue_rank_goseq(geneset, 
                                       GOSeq_allbias, 
                                       metric)
    

  } else{ # metric value remains unchanged
    
    # set bias to NULL, that way, only length bias is accounted for in the data set
    bias  <-  NULL
    
    # update documentation frame 
    doc_intern[1, "optimal_parameter"]  <-  "Length Bias Only"
    doc_intern[1, metric]  <-  metricvalue
    
    

  }
  
  ############
  ### 5.step: change universe as ALL genes from experiment, not only those that can be annotated to gene set database 
  ############
  
  # run GOSeq: Include ALL genes in the calculation of the p-value (genes_wocat = TRUE)
  GOSeq_univ  <-  GOSeq_pipeline(optim_DE, 
                                bias = bias, 
                                geneID = gene_ID, 
                                geneset_database = geneset_database, 
                                genes_wocat = TRUE, 
                                organism = organism)
  
  
  
  # update parameter depending on whether the adapted universe leads to a decreased metric value 
  genes_wocat  <-  ifelse(pvalue_rank_goseq(geneset, GOSeq_univ, metric) < metricvalue, TRUE, FALSE)
  
  metricvalue  <-  min(pvalue_rank_goseq(geneset, GOSeq_univ, metric), metricvalue)
  
  # update documentation frame
  doc_intern[2, "optimal_parameter"]  <-  ifelse(pvalue_rank_goseq(geneset, GOSeq_univ, metric) < metricvalue, "Adapted", "Default")
  doc_intern[2, metric]  <-  metricvalue
  
 
   ############
   ### 6. step: change method for calculation of the p-value 
   ############
  
    calc_methods  <-  c("Wallenius", "Sampling")
  
    GOSeq_calcmethods  <-  lapply(FUN = GOSeq_pipeline, 
                                  DE_results = optim_DE, 
                                  geneset_database = geneset_database, 
                                  X = calc_methods, 
                                  bias = bias, 
                                  genes_wocat = genes_wocat, 
                                  geneID = gene_ID, 
                                  organism = organism)
    
    # get metricvalue resulting from both calculation methods
    metricvalue_calcmethod  <-  sapply(FUN = pvalue_rank_goseq, 
                                       term = geneset, 
                                       X = GOSeq_calcmethods, 
                                       metric = metric)
    
    # index of optimal calculation method (in case of ties choose smaller index)
    ind_opt_calcmethod  <-  min(which( metricvalue_calcmethod == min( metricvalue_calcmethod)))
    
    # optimal calculation method 
    opt_calcmethod  <-  calc_methods[ind_opt_calcmethod]
    
    # update documentation frame
    doc_intern[3, "optimal_parameter"]  <-  opt_calcmethod
    doc_intern[3, metric]  <-  min(metricvalue_calcmethod)
    
    # merge documentation frames 
    doc  <-  rbind(doc, doc_intern)


  return(list(default = goseq_results_deseq2, #default DAVID result
              optim = GOSeq_calcmethods[[ind_opt_calcmethod]], # GOSeq results 
              documentation = doc)) #documentation frame 
  
}

################################################################################
### Run optimizations ##########################################################
################################################################################

phen_pickrell_list <- list()

for(i in 1:ncol(phen_pickrell)){
  
  phen_pickrell_list[[i]] <- phen_pickrell[,i]
  
}

phen_bottomly_list <- list()

for(i in 1:ncol(phen_bottomly)){
  
  phen_bottomly_list[[i]] <- phen_bottomly[,i]
  
}

### optimization of adjusted p-values 

#############
### Pickrell 
#############

### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimP_GOSeq_tcell_Pickrell_originalphenotype <-  goseq_rank_pvalue_optim(geneset = "GO:0002456", 
                                                                    "GO:BP",
                                                                    metric = "p_adj", 
                                                                    Biobase::exprs(pickrell.eset), 
                                                                    pickrell.eset$gender)

# save results
save(optimP_GOSeq_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimP_GOSeq_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_GOSeq_tcell_Pickrell_phenotypepermutations <-  lapply(FUN = goseq_rank_pvalue_optim, 
                                                       geneset =  "GO:0002456", 
                                                       geneset_database = "GO:BP",
                                                       metric = "p_adj",
                                                       expression_data = Biobase::exprs(pickrell.eset), 
                                                       X = phen_pickrell_list)

# save results
save(optimP_GOSeq_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_GOSeq_tCell_Pickrell_PhenotypePermutations.RData")


### (II) Gene Set "Demethylation" -> "GO:0070988"

# original phenotype assignment 
optimP_GOSeq_Demethylation_Pickrell_originalphenotype <-  goseq_rank_pvalue_optim(geneset= "GO:0070988", "GO:BP",metric = "p_adj", Biobase::exprs(pickrell.eset), pickrell.eset$gender)

# save results
save(optimP_GOSeq_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimP_GOSeq_Demethylation_Pickrell_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations <-  lapply(FUN = goseq_rank_pvalue_optim, 
                                                                      geneset =  "GO:0070988", 
                                                                      metric = "p_adj",
                                                                      expression_data = Biobase::exprs(pickrell.eset), 
                                                                      geneset_database = "GO:BP",
                                                                      X = phen_pickrell_list)

save(optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_GOSeq_Demethylation_Pickrell_PhenotypePermutations.RData")


#############
### Bottomly 
#############

### (I) Gene set "Metabolic Process: GO:0008152"

# original phenotype assignment 
optimP_GOSeq_MetabolicProcess_Bottomly_originalphenotype <- 
  goseq_rank_pvalue_optim(geneset = "GO:0008152", 
                          "GO:BP","p_adj", 
                          Biobase::exprs(bottomly.eset), 
                          bottomly.eset$strain)

# save results
save(optimP_GOSeq_MetabolicProcess_Bottomly_originalphenotype, 
     file = "./Results/optimP_GOSeq_MetabolicProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                 geneset = "GO:0008152", 
                                                                 geneset_database = "GO:BP",
                                                                 metric = "p_adj",
                                                                 expression_data = Biobase::exprs(bottomly.eset), 
                                                                 X = phen_bottomly_list)
# save results 
save(optimP_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations, 
     file = "./Results/optimP_GOSeq_MetabolicProcess_Bottomly_PhenotypePermutations.RData")




### (II) Gene set "Cellular Process: GO:0009987"

# original phenotype assignment 
optimP_GOSeq_CellularProcess_Bottomly_originalphenotype <- 
  goseq_rank_pvalue_optim(geneset = 
                            "GO:0009987", 
                            "GO:BP",
                            metric = "p_adj", 
                            Biobase::exprs(bottomly.eset), 
                            bottomly.eset$strain)

  # save results
save(optimP_GOSeq_CellularProcess_Bottomly_originalphenotype, 
     file = "./Results/optimP_GOSeq_CellularProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_GOSeq_CellularProcess_Bottomly_Phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                geneset = "GO:0009987", 
                                                                geneset_database = "GO:BP",
                                                                metric = "p_adj",
                                                                expression_data = Biobase::exprs(bottomly.eset), 
                                                                X = phen_bottomly_list)
# save results 
save(optimP_GOSeq_CellularProcess_Bottomly_Phenotypepermutations, 
     file = "./Results/optimP_GOSeq_CellularProcess_Bottomly_PhenotypePermutations.RData")



# optimization of the ranks in the GSA results


#############
### Pickrell 
#############

### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimRank_GOSeq_tcell_Pickrell_originalphenotype <-  goseq_rank_pvalue_optim(geneset= "GO:0002456", 
                                                                       "GO:BP",
                                                                       metric = "rank", 
                                                                       Biobase::exprs(pickrell.eset), 
                                                                       pickrell.eset$gender)

# save results
save(optimRank_GOSeq_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimRank_GOSeq_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_GOSeq_tcell_Pickrell_phenotypepermutations <-  lapply(FUN = goseq_rank_pvalue_optim, 
                                                               geneset =  "GO:0002456", 
                                                               geneset_database = "GO:BP",
                                                               metric = "rank",
                                                               expression_data = Biobase::exprs(pickrell.eset), 
                                                               X = phen_pickrell_list)

  # save results
save(optimRank_GOSeq_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_tCell_Pickrell_PhenotypePermutations.RData")


### (II) Gene Set "Demethylation" -> "GO:0070988"

# original phenotype assignment 
optimRank_GOSeq_Demethylation_Pickrell_originalphenotype <-  
                    goseq_rank_pvalue_optim(geneset = "GO:0070988", 
                                            "GO:BP",
                                            metric = "rank", 
                                            Biobase::exprs(pickrell.eset), 
                                            pickrell.eset$gender)

# save results
save(optimRank_GOSeq_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimRank_GOSeq_Demethylation_Pickrell_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_GOSeq_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                       geneset =  "GO:0070988", 
                                                                       metric = "rank",
                                                                       expression_data = Biobase::exprs(pickrell.eset), 
                                                                       geneset_database = "GO:BP",
                                                                       X = phen_pickrell_list)

save(optimRank_GOSeq_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_Demethylation_Pickrell_PhenotypePermutations.RData")


#############
### Bottomly 
#############

### (I) Gene set "Metabolic Process: GO:0008152"

# original phenotype assignment 
optimRank_GOSeq_MetabolicProcess_Bottomly_originalphenotype <- goseq_rank_pvalue_optim(geneset = "GO:0008152", 
                                                                                       "GO:BP",
                                                                                       "rank", 
                                                                                       Biobase::exprs(bottomly.eset), 
                                                                                       bottomly.eset$strain)

# save results
save(optimRank_GOSeq_MetabolicProcess_Bottomly_originalphenotype, 
     file = "./Results/optimRank_GOSeq_MetabolicProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                         geneset = "GO:0008152", 
                                                                         geneset_database = "GO:BP",
                                                                         metric = "rank",
                                                                         expression_data = Biobase::exprs(bottomly.eset), 
                                                                         X = phen_bottomly_list)
# save results 
save(optimRank_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_MetabolicProcess_Bottomly_PhenotypePermutations.RData")




### (II) Gene set "Cellular Process: GO:0009987"

# original phenotype assignment 
optimRank_GOSeq_CellularProcess_Bottomly_originalphenotype <- goseq_rank_pvalue_optim(geneset = "GO:0009987", 
                                                                                      "GO:BP",
                                                                                      metric = "rank", 
                                                                                      Biobase::exprs(bottomly.eset), 
                                                                                      bottomly.eset$strain)

# save results
save(optimRank_GOSeq_CellularProcess_Bottomly_originalphenotype, 
     file = "./Results/optimRank_GOSeq_CellularProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_GOSeq_CellularProcess_Bottomly_Phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                        geneset = "GO:0009987", 
                                                                        geneset_database = "GO:BP",
                                                                        metric = "rank",
                                                                        expression_data = Biobase::exprs(bottomly.eset), 
                                                                        X = phen_bottomly_list)
# save results 
save(optimRank_GOSeq_CellularProcess_Bottomly_Phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_CellularProcess_Bottomly_PhenotypePermutations.RData")




    

