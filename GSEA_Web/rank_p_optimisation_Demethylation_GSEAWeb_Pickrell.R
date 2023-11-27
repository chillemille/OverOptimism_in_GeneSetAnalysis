#############################################################################################################################
### Optimization of GSEA Web Application for Pickrell data set: Mimimize adjusted p-value of gene set "Demethylation" (GO:BP)
#############################################################################################################################

library(tweeDEseqCountData) # contains Pickrell data set 
library(edgeR) # for builtin pre-filtering 
library(dplyr)


# load Pickrell data set along with the true and permuted sample conditions
source("./Random_Phenotype_Permutations.R")

# load the required pre-processing functions
source("./Random_Phenotype_Permutations.R")

# load functions that perform a transformation of the RNA-Seq data 
source("./RNASeq_Transformation.R")


################################################################################
### Generate required ordner structure #########################################
################################################################################

dir.create("./Pickrell")
dir.create("./Pickrell/p_adj")
dir.create("./Pickrell/p_adj/Data_Demethylation")
dir.create("./Pickrell/p_adj/Data_Demethylation/Prep")
dir.create("./Pickrell/p_adj/Data_Demethylation/Prep/Phen_Original")
dir.create("./Pickrell/p_adj/Data_Demethylation/Prep/Phenotypes")
dir.create("./Pickrell/p_adj/Data_Demethylation/Raw")
dir.create("./Pickrell/p_adj/Data_Demethylation/Raw/Phen_Original")
dir.create("./Pickrell/p_adj/Data_Demethylation/Raw/Phenotypes")



for(i in 1:10){
  
  dir.create(paste0("./Pickrell/p_adj/Data_Demethylation/Prep/Phen_Permutation", i))
  
  dir.create(paste0("./Pickrell/p_adj/Data_Demethylation/Raw/Phen_Permutation", i))
  
  
  
}

################################################################################
### (I) Complete Optimization Process (pre-processing and internal parameters) #
################################################################################



################################################################################
### Original Phenotype Assignment ##############################################
################################################################################


# export phenotype assignments and gene expression data set pre-processed in default manner 

  
  ##############
  ### phenotypes 
  ##############
  
  # convert levels ("female", "male") to 0 and 1 
  phen_orig <- c()
  phen_orig[ pickrell.eset$gender == levels(pickrell.eset$gender)[[1]]] <- 0
  phen_orig[ pickrell.eset$gender == levels(pickrell.eset$gender)[[2]]] <- 1
  
  path_phen <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phenotypes/Phenotype_original.txt")
  
  write.table(phen_orig,
              file = path_phen,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  
  ####################################
  ### default gene expression data set 
  ####################################
  
  # default pre-process pickrell data 
  dat_default_phenorig <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
                          voom_trans(phenotype_labels = pickrell.eset$gender)
  
  # generate path 
  path_dat_phenorig <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Original/exprdat_default_phen_original.txt")
  
  ### export 
  write.table(dat_default_phenorig,
              file = path_dat_phenorig,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)
  
  ############################################
  ### vst transformed gene expression data set
  ############################################
  
  
  # default pre-process pickrell data 
  dat_vst_phenorig <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
                      variancetransform(phenotype_labels = pickrell.eset$gender )
  
  # generate path 
  path_dat_vst_phenorig <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Original/exprdat_vst_phen_original.txt")
  
  
  # export
  write.table(dat_vst_phenorig,
              file = path_dat_vst_phenorig,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)
  







#########
# step 1: default 
#########

# -> adjusted p-value = 0.13153844
# note: the gene set is already detected as differentially enriched in default configurations (FDR < 0.25)

#########
# step 2: vst-transformed gene expression data set 
#########

# -> adjusted p-value = 1
# ->> return to default data set 

#########
# step 3: pre-filtering using edgeR's builtin function filterByExpr (based on default)
#########
  
  
  # generate pre-filtering indicator using filterByExpr()
  prefilt_ind_phenorig <-  DGEList(Biobase::exprs(pickrell.eset), group = pickrell.eset$gender) %>% filterByExpr()
  # perform voom-transformation on accordingly pre-filtered pickrell data set
  exprdat_prefilt_phenorig <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phenorig,], phenotype_labels = pickrell.eset$gender)
  
  ### export 
  path_filterByExpr_phenorig <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Original/exprdat_filterByExpr_phenorig.txt")
  
  write.table(exprdat_prefilt_phenorig, 
              file = path_filterByExpr_phenorig, 
              quote = FALSE, 
              row.names = TRUE, 
              col.names = TRUE)
  

# -> adjusted p-value = 1
# ->> return to default pre-filtering (transformation method voom)



########
# step 4: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.226
# alternative 2: Difference of Classes -> adj. p-value = 1 

# ->> return to gene-level ranking metric Signal2NoiseRatio

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.00176
# alternative 3: exponent 2 -> adj. p-value = 0.001406


### final result: adj. p-value = 0.001406
# achieved with alternative exponent 2

################################################################################
### Random Permutations of Phenotype Assignment ################################
################################################################################


# export phenotype assignments and gene expression data set pre-processed in default manner 
for(i in 1:ncol(phen_pickrell)){
  
  ##############
  ### phenotypes 
  ##############
  
  # convert levels ("female", "male") to 0 and 1 
  phen <- c()
  phen[phen_pickrell[,i] == "female"] <- 0
  phen[phen_pickrell[,i] == "male"] <- 1
  
  path_phen <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phenotypes/Phenotype_Permutation",
                      i, 
                      ".txt")
  
  write.table(phen,
              file = path_phen,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  
  ####################################
  ### default gene expression data set 
  ####################################
  
  # default pre-process pickrell data 
  dat_default <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
                 voom_trans(phenotype_labels = phen_pickrell[,i])
  
  # generate path 
  path_dat <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                     i,
                     "/exprdat_default_phen_permutation", 
                     i, 
                     ".txt")
  
  ### export 
  write.table(dat_default,
              file = path_dat,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)
  
  ############################################
  ### vst transformed gene expression data set
  ############################################
  
  
  # default pre-process pickrell data 
  dat_vst <- pre_filt(Biobase::exprs(pickrell.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
             variancetransform(phenotype_labels = phen_pickrell[,i])
  
  # generate path 
  path_dat_vst <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                         i,
                         "/exprdat_vst_phen_permutation", 
                         i, 
                         ".txt")
  
  
  #export
  write.table(dat_vst,
              file = path_dat_vst,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)

}




################################################################################
### Phenotype Permutation 1 ####################################################
################################################################################

i <- 1

#########
# step 1: default 
#########

# -> adj. p-value = 0.9893133

  

#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 0.6722215
# -> proceed with vst-transformed RNA-Seq data set 


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen1 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,1]) %>% 
                      filterByExpr()
# perform vst-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen1 <- variancetransform(Biobase::exprs(pickrell.eset)[prefilt_ind_phen1,], 
                                           phenotype_labels = phen_pickrell[,1])

### export 
path_filterByExpr_phen1 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, 
                                  ".txt")

write.table(exprdat_prefilt_phen1, 
            file = path_filterByExpr_phen1, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 0.93461543
# -> return to default pre-filtering (vst-transformed RNA-Seq data set)




########
# step 4: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.675105
# alternative 2: Difference of Classes  -> adj. p-value = 0.85543203

# ->> return to default ranking metric Signal2Noise ratio 


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.87602484
# alternative 2: exponent 1.5 -> adj. p-value = 0.39360744
# alternative 3: exponent 2 -> adj. p-value = 0.47316462


### final results: adj. p-value = 0.39360744
# achieved with 
# alternative RNA-Seq transformation vst 
# alternative exponent 1.5 




################################################################################
### Phenotype Permutation 2 ####################################################
################################################################################

i <- 2

#########
# step 1: default 
#########

# -> adj. p-value = 0.9358796


#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 1 
# -> return to voom transformation 


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen2 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,2]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen2 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen2,], 
                                    phenotype_labels = phen_pickrell[,2])

### export 
path_filterByExpr_phen2 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, 
                                  ".txt")

write.table(exprdat_prefilt_phen2, 
            file = path_filterByExpr_phen2, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 0.99565524
# -> return to default pre-fitering (RNA-Seq transformation voom)

########
# step 6: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.94376874
# alternative 2: Difference of Classes -> adj. p-value = 1 

# -> proceed with Signal2Noise ratio as gene-level ranking metric 


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1 
# alternative 2: exponent 1.5 -> adj- p-value = 0.9849069
# alternative 3: exponent 2 -> adj. p-value = 0.8779408


# final results: adj p-value = 0.8779408
# achieved with alternative exponent 2 




################################################################################
### Phenotype Permutation 3 ####################################################
################################################################################

i <- 3


#########
# step 1: default 
#########

# -> adj p-value = 0.9609585



#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj p-value = 1
# ->> return to default RNA-Seq transformation mehod voom 



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr (on vst-transformed gene expression data set)
########


# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen3 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,3]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen3 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen3,], 
                                    phenotype_labels = phen_pickrell[,3])

### export 
path_filterByExpr_phen3 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, 
                                  ".txt")

write.table(exprdat_prefilt_phen3, 
            file = path_filterByExpr_phen3, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 0.996943
# --> return to default pre-filtering (based on transformation method voom)


########
# step 5: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.9687702
# alternative 2: Difference of Classes -> adj. p-value = 1

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.6908496
# alternative 2: exponent 1.5 -> adj. p-value = 1 
# alternative 3: exponent 2 -> adj. p-value = 0.7190371


### final results: adj p-value = 0.6908496
# achieved with alternative exponent 0 




################################################################################
### Phenotype Permutation 4 ####################################################
################################################################################

i <- 4


#########
# step 1: default 
#########

# -> adj. p-value = 1 


#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 0.805615
# ->> proceed with vst as RNA-Seq transformation method 


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen4 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,4]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen4 <- variancetransform(Biobase::exprs(pickrell.eset)[prefilt_ind_phen4,], 
                                           phenotype_labels = phen_pickrell[,4])

### export 
path_filterByExpr_phen4 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, 
                                  ".txt")

write.table(exprdat_prefilt_phen4, 
            file = path_filterByExpr_phen4, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 0.72907937
# ->> proceed with pre-filtering using filterByExpr (based on vst-transformed RNA-Seq data set)




########
# step 4: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.7264533
# alternative 2: Difference of Classes -> adj. p.value = 0.907747


# ->> proceed with t-Test (t-statistic) as gene-level ranking metric



######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.75388193
# alternative 2: exponent 1.5 -> adj. p-value = 0.7826668
# alternative 3: exponent 2 -> adj. p-value = 0.6283782

# final results: adj. p-value = 0.6283782
# achieved with 
# alternative RNA-Seq transformation vst
# alternative pre-filtering method filterByExpr
# alternative gene-level ranking metric tTest
# alternative exponent 2 



################################################################################
### Phenotype Permutation 5 ####################################################
################################################################################

i <- 5


#########
# step 1: default 
#########

# -> adj. p-value = 1

#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 1
# ->> return to default transformation method voom 


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen5 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,5]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen5 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen5,], 
                                    phenotype_labels = phen_pickrell[,5])

### export 
path_filterByExpr_phen5 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, 
                                  ".txt")

write.table(exprdat_prefilt_phen5, 
            file = path_filterByExpr_phen5, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 1 
# -> return to default pre-filtering using treshold 10 (voom - transformation)


########
# step 6: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 1 
# alternative 2: Difference of Classes -> adj. p-value = 1

# -> return to default gene-level ranking metric Signal2Noise ratio 


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.96786606
# alternative 2: exponent 1.5 -> adj. p-value = 0.9941131
# alternative 3: exponent 2 -> adj. p-value = 1 

# final results: adj. p-value = 0.96786606
# -> achieved with alternative exponent 0 


################################################################################
### Phenotype Permutation 6 ####################################################
################################################################################

i <- 6


#########
# step 1: default 
#########

# -> adj. p-value = 0.8138262



#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 0.99553865

# ->> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen6 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,6]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen6 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen6,], 
                                    phenotype_labels = phen_pickrell[,6])

### export 
path_filterByExpr_phen6 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, 
                                  ".txt")

write.table(exprdat_prefilt_phen6, 
            file = path_filterByExpr_phen6, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)


# ->> adj. p-value = 1 
# -> return to default pre-filtering (based on voom-transformed RNA-Seq data set)


########
# step 5: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.8101068
# alternative 2: Difference of Classes -> adj. p-value = 0.9396362


# -> proceed with gene-level ranking metric t-test

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.93282014
# alternative 2: exponent 1.5 -> adj. p-value = 0.9470678
# alternative 3: exponent 2 -> adj. p-value = 0.5283564

### final results: adj. p-value = 0.5283564
# achieved with 
# ALTERNATIVE gene-level ranking metric t-Test
# ALTERNATIVE exponent 2 






################################################################################
### Phenotype Permutation 7 ####################################################
################################################################################

i <- 7

#########
# step 1: default 
#########

# -> adj. p-value = 0.5119471



#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 0.87930465
# -> return to default transformation method voom 


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen7 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,7]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen7 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen7,], 
                                    phenotype_labels = phen_pickrell[,7])

### export 
path_filterByExpr_phen7 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, 
                                  ".txt")

write.table(exprdat_prefilt_phen7, 
            file = path_filterByExpr_phen7, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# ->  adj. p-value = 1
# ->> return to default pre-filtering (based on voom-transformation)


########
# step 5: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.4938039
# alternative 2: Difference of Classes -> adj. p-value = 0.6989737


# -> proceed with gene-level ranking metric t-Test 

######
# step 6: change exponent
######

# alternative 1: exponent 0 ->  adj. p-value = 0.8617825
# alternative 2: exponent 1.5 -> adj. p-value = 0.8713911
# alternative 3: exponent 2 -> adj. p-value = 0.54182476

### end results: adj. p-value = 0.4938039
# Achieved with alternative gene-level ranking metric tTest

################################################################################
### Phenotype Permutation 8 ####################################################
################################################################################

i <- 8

#########
# step 1: default 
#########

# -> adj. p-value = 1



#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 1 
# -> return to default RNA-Seq transformation method voom



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen8 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,8]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen8 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen8,], 
                                    phenotype_labels = phen_pickrell[,8])

### export 
path_filterByExpr_phen8 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", 
                                  i, "
                                  .txt")

write.table(exprdat_prefilt_phen8, 
            file = path_filterByExpr_phen8, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 1
# ->> return to default pre-filtering (based on RNA-Seq transformation method voom)



########
# step 4: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 1
# alternative 2: Difference of Classes -> adj. p-value = 1

# ->> return to default gene-level ranking metric Signal2Noise ratio 

######
# step 5: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1 
# alternative 2: exponent 1.5 -> adj. p-value = 0.92030305
# alternative 3: exponent 2 -> adj. p-value = 1

# final results: adj p-value = 0.92030305
# achieved with alternative exponent 1.5




################################################################################
### Phenotype Permutation 9 ####################################################
################################################################################

i <- 9

#########
# step 1: default 
#########

# -> adj. p-value = 0.7891276

#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 1
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen9 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,9]) %>% 
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen9 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen9,], 
                                    phenotype_labels = phen_pickrell[,9])

### export 
path_filterByExpr_phen9 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen9, 
            file = path_filterByExpr_phen9, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 1
# ->> return to default pre-filtering (based on voom-transformation)



########
# step 4: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 0.7790329
# alternative 2: Difference of Classes -> adj. p-value = 0.85717106

# ->> return to default gene-level ranking metric Signal2Noise ratio




######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.66642624
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.7071824

# final results: adj p-value = 0.66642624
# achieved with alternative exponent 0 





################################################################################
### Phenotype Permutation 10 ####################################################
################################################################################

i <- 10


#########
# step 1: default 
#########

# -> adj. p-value = 1


#########
# step 2: vst-transformed gene expression data set 
#########

# -> adj. p-value = 1
# -> return to default transformation method voom 


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen10 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[,10]) %>% 
                       filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen10 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen10,], 
                                     phenotype_labels = phen_pickrell[,10])

### export 
path_filterByExpr_phen10 <- paste0("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/GSEA_Web/Pickrell/p_value/Data_Demethylation/Raw/Phen_Permutation",
                                   i,
                                   "/exprdat_filterByExpr_phen_permutation", 
                                   i, 
                                   ".txt")

write.table(exprdat_prefilt_phen10, 
            file = path_filterByExpr_phen10, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = TRUE)

# -> adj. p-value = 1
# -> return to default pre-filtering using filterByExpr



########
# step 4: change gene-level ranking metric 
########

# alternative 1: t-Test -> adj. p-value = 1
# alternative 2: Difference of Classes -> adj. p-value = 1

# ->> return to default gene-level ranking metric signal2noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.9977417
# alternative 3: exponent 2 -> adj. p-value = 1


# final results: adj. p-value = 0.9977417
# achieved with alternative exponent 1.5




################################################################################
### (II) Optimization of Internal Parameters ###################################
################################################################################



################################################################################
### Original Phenotype Assignment ##############################################
################################################################################







################################################################################
### Random Permutations of Phenotype Assignment ################################
################################################################################



#########
# step 1: default 
#########

# -> 6 DEGS 

########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS 
# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 5 DEGS 
# alternative 2: Difference of Classes -> 0 DEGS 

# -> return to default gene-level ranking metric Signal2Noise ratio 


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 20 DEGS 
# alternative 2: exponent 1.5 -> 21 DEGS 
# alternative 3: exponent 2 -> 32 DEGS


### final results: 32 DEGS 
# achieved with ALTERNATIVE exponent 2 






################################################################################
### Phenotype Permutation 1 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 1 DEGS 
# -> return to default geneset database GO (BP)


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS 



########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 1 DEGS 
# alternative 2: Difference of Classes -> 0 DEGS 

# ->> return to default gene-level ranking metric Signal2Noise Ratio


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 1 DEGS 
# alternative 3: exponent 2 -> 2 DEGS 

### -> final results: 2 DEGS 
# achieved with alternative exponent 2


################################################################################
### Phenotype Permutation 2 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 79 DEGS 


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS 
# -> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 80 DEGS 
# alternative 2: Difference of Classes -> 61 DEGS 

# ->> proceed with alternative gene-level ranking metric t-Test


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 32 DEGS 
# alternative 2: exponent 1.5 -> 8 DEGS 
# alternative 3: exponent 2 -> 6 DEGS 

### final results: 80 DEGS
# achieved with alternative gene-level ranking metric t-Test 


################################################################################
### Phenotype Permutation 3 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 0 DEGS 

########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS 
# -> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 0 DEGS 
# alternative 2: Difference of Classes -> 0 DEGS 

# -> return to default gene-level ranking metric Signal2Noise Ratio 


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 0 DEGS 
# alternative 3: exponent 2 -> 0 DEGS 

# -> end results: 0 DEGS 
# optimal parameter configuration coincides with default configuration 
# number of differentially enriched gene sets could not be increased 


################################################################################
### Phenotype Permutation 4 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 2 DEGS 



########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS 
# ->> return to default geneset database GO (BP)



########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 3 DEGS 
# alternative 2: Difference of Classes 0 DEGS 

# ->> proceed with alternative gene-level ranking metric t-Test


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 1 DEGS 
# alternative 2: exponent 1.5 -> 0 DEGS 
# alternative 3: exponent 2 -> 0 DEGS 


# final results: 3 DEGS 
# achieved with alternative gene-level ranking metric t-Test 




################################################################################
### Phenotype Permutation 5 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 3 DEGS 


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS 
# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 3 DEGS 
# alternative 2: Difference of Classes -> 0 DEGS 

# -> return to default gene-level ranking metric Signal2Noise Ratio 


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 0 DEGS 
# alternative 3: exponent 2 -> 0 DEGS 

# -> return to default exponent 1

### final results: 3 DEGS 
# optimal parameter configuration coincides with default parameter configuration
# number of differentially enriched gene sets could not be increased 


################################################################################
### Phenotype Permutation 6 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 0 DEGS 


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS 
# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 0 DEGS
# alternative 2: Difference of Classes -> 0 DEGS 

# ->> return to default gene-level ranking metric Signal2Noise Ratio 



######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 0 DEGS 
# alternative 3: exponent 2 -> 1 DEGS 

### final results: 1 DEGS 
# achieved with ALTERNATIVE exponent 2 


################################################################################
### Phenotype Permutation 7 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 155 DEGS 


########
# step 2: change gene set database to KEGG
########

# -> 4 DEGS 
# ->> return to default geneset database GO (BP)



########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 152 DEGS 
# alternative 2: Difference of Classes -> 118 DEGS 

# ->> return to default gene-level ranking metric Signal2Noise Ratio 


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 0 DEGS 
# alternative 3: exponent 2 -> 70 DEGS 

#->> final results: 155 DEGS
# optimal parameter configuration coincides with default parameter configuration
# number of differentially enriched gene sets could not be increased 


################################################################################
### Phenotype Permutation 8 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 0 DEGS 


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS 
# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 0 DEGS 
# alternative 2: Difference of Classes -> 0 DEGS 

# ->> return to default gene-level ranking metric Signal2Noise ratio 


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 40 DEGS 
# alternative 3: exponent 2 -> 1 DEGS 


### final results: 40 DEGS 
# achieved with ALERNATIVE exponent 1.5 


################################################################################
### Phenotype Permutation 9 ####################################################
################################################################################


#########
# step 1: default 
#########

# -> 20 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS 
# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 17 DEGS 
# alternative 2: Difference of Classes -> 6 DEGS 

# ->> return to default gene-level ranking metric Difference of Classes 


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 0 DEGS 
# alternative 3: exponent 2 -> 0 DEGS 

# -> return to default exponent 1 

### final results: 20 DEGS 
# optimal parameter configuration coincides with default parameter configuration
# number of differentially enriched gene sets could not be increased 


################################################################################
### Phenotype Permutation 10 ###################################################
################################################################################


#########
# step 1: default 
#########

# -> 0 DEGS 


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS 
# -> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric 
########

# alternative 1: t-Test -> 0 DEGS 
# alternative 2: Difference of Classes -> 0 DEGS 

# ->> return to default gene-level ranking metric Signal2Noise Ratio 


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS 
# alternative 2: exponent 1.5 -> 0 DEGS 
# alternative 3: exponent 2 -> 0 DEGS

### return to default exponent 1 

### -> final results: 0 DEGS 
# optimal parameter configuration coincides with default paramter configuration 
# number of differentially enriched gene sets could not be increased 























