#############################################################################################################################
### Optimization of GSEA Web Application for Pickrell data set: Mimimize adjusted p-value of gene set "Demethylation" (GO:BP)
#############################################################################################################################

library(edgeR) # for built-in pre-filtering
library(dplyr)



# load Bottomly data set along with the true and permuted sample conditions
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")

# load the required pre-processing functions
source("./R/Functions/PreProcessing_Functions.R")

# load functions that perform a transformation of the RNA-Seq data
source("./R/Functions/RNASeq_Transformation.R")


################################################################################
### Generate required ordner structure #########################################
################################################################################


dir.create("./Results/Intermediate_results/GSEA_Web")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Prep")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Prep/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Prep/Phenotypes")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phenotypes")

for(i in 1:10){

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Prep/Phen_Permutation", i))

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation", i))



}


################################################################################
### Original Phenotype Assignment ##############################################
################################################################################


# export phenotype assignments and gene expression data set pre-processed in default manner


  ##############
  ### phenotypes
  ##############

  # convert levels ("female", "male") to 0 and 1
  phen_orig <- c()
  phen_orig[ bottomly.eset$strain == levels(bottomly.eset$strain)[[1]]] <- 0
  phen_orig[ bottomly.eset$strain == levels(bottomly.eset$strain)[[2]]] <- 1

  path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phenotypes/Phenotype_original.txt")

  write.table(phen_orig,
              file = path_phen,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)


  ####################################
  ### default gene expression data set
  ####################################

  # default pre-process pickrell data
  dat_default_phenorig <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
    voom_trans(phenotype_labels = bottomly.eset$strain)

  # generate path
  path_dat_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Original/exprdat_default_phen_original.txt")

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
  dat_vst_phenorig <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
    variancetransform(phenotype_labels = bottomly.eset$strain )

  # generate path
  path_dat_vst_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Original/exprdat_vst_phen_original.txt")


  # export
  write.table(dat_vst_phenorig,
              file = path_dat_vst_phenorig,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)








#########
# step 1: default
#########

# -> adjusted p-value = 0.96452993


#########
# step 2: vst-transformed gene expression data set
#########

# -> adjusted p-value = 0.9636415
# ->> proceed with alternative transformation vst

#########
# step 3: pre-filtering using edgeR's builtin function filterByExpr (based on default)
#########


  # generate pre-filtering indicator using filterByExpr()
  prefilt_ind_phenorig <-  DGEList(Biobase::exprs(bottomly.eset), group = bottomly.eset$strain) %>% filterByExpr()
  # perform voom-transformation on accordingly pre-filtered pickrell data set
  exprdat_prefilt_phenorig <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phenorig,], phenotype_labels = bottomly.eset$strain)

  ### export
  path_filterByExpr_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Original/exprdat_filterByExpr_phenorig.txt")

  write.table(exprdat_prefilt_phenorig,
              file = path_filterByExpr_phenorig,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)


# -> adjusted p-value = 0.98552865
# ->> return to default pre-filtering approach



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.9746558
# alternative 2: Difference of Classes -> adj. p-value =1

# ->> return to default gene-level statistic Signal2Noise ratio

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.94919527
# alternative 3: exponent 2 -> adj. p-value = 0.86383224


### final result: adj. p-value = 0.86383224
# achieved with
# ALTERNATIVE transformation method vst
# ALTERNATIVE exponent 2
#

################################################################################
### Random Permutations of Phenotype Assignment ################################
################################################################################


# export phenotype assignments and gene expression data set pre-processed in default manner
for(i in 1:ncol(phen_bottomly)){

  ##############
  ### phenotypes
  ##############

  # convert levels ("female", "male") to 0 and 1
  phen <- c()
  phen[phen_bottomly[,i] == levels(bottomly.eset$strain)[[1]]] <- 0
  phen[phen_bottomly[,i] == levels(bottomly.eset$strain)[[2]]] <- 1

  path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phenotypes/Phenotype_Permutation",
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
  dat_default <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
    voom_trans(phenotype_labels = phen_bottomly[,i])

  # generate path
  path_dat <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
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
  dat_vst <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
    variancetransform(phenotype_labels = phen_bottomly[,i])

  # generate path
  path_dat_vst <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
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

# -> adj. p-value = 0.7845396



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.8127381
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen1 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,1]) %>% filterByExpr()
# perform vst-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen1 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen1,], phenotype_labels = phen_bottomly[,1])

### export
path_filterByExpr_phen1 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen1,
            file = path_filterByExpr_phen1,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.6134048
# -> proceed with alternative pre-filtering using filterByExpr()



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.6119803
# alternative 2: Difference of Classes  -> adj. p-value = 0.76970565

# ->> proceed with alternative gene-level statistic t-Test


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.514389
# alternative 3: exponent 2 -> adj. p-value = 0.7757481


### final results: adj. p-value = 0.514389
# achieved with
# ALTERNATIVE pre-filtering using filterByExpr()
# ALTERNATIVE gene-level statistic t-Test
# ALTERNATIVE exponent 1.5




################################################################################
### Phenotype Permutation 2 ####################################################
################################################################################

i <- 2

#########
# step 1: default
#########

# -> adj. p-value = 0.85781413


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.7875908
# -> proceed with alternative transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen2 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,2]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen2 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen2,], phenotype_labels = phen_bottomly[,2])

### export
path_filterByExpr_phen2 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen2,
            file = path_filterByExpr_phen2,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.85839516
# -> return to default pre-filtering approach

########
# step 6: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.78607225
# alternative 2: Difference of Classes -> adj. p-value =  0.9051226

# -> proceed with alternative gene-level statistic t-Test


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj- p-value = 0.9282115
# alternative 3: exponent 2 -> adj. p-value = 0.60153824


# final results: adj p-value = 0.60153824
# achieved with
# ALTERNATIVE transformation vst
# ALTERNATIVE gene-level statistic t-Test
# ALTERNATIVE exponent 2




################################################################################
### Phenotype Permutation 3 ####################################################
################################################################################

i <- 3


#########
# step 1: default
#########

# -> adj p-value = 1



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj p-value = 1
# ->> return to default transformation method voom



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########


# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen3 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,3]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen3 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen3,], phenotype_labels = phen_bottomly[,3])

### export
path_filterByExpr_phen3 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen3,
            file = path_filterByExpr_phen3,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value =1
# --> return to default pre-filtering


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 1
# alternative 2: Difference of Classes -> adj. p-value = 1

# ->> return to default gene-level statistic Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.97550386
# alternative 2: exponent 1.5 -> adj. p-value = 0.6831988
# alternative 3: exponent 2 -> adj. p-value = 0.97805965


### final results: adj p-value = 0.6831988
# achieved with ALTERNATIVE exponent 1.5




################################################################################
### Phenotype Permutation 4 ####################################################
################################################################################

i <- 4


#########
# step 1: default
#########

# -> adj. p-value = 0.803125


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9080778
# ->> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen4 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,4]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen4 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen4,],
                                    phenotype_labels = phen_bottomly[,4])

### export
path_filterByExpr_phen4 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen4,
            file = path_filterByExpr_phen4,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.65209365
# ->> proceed with alternative pre-filtering approach using filterByExpr()




########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.6623734
# alternative 2: Difference of Classes -> adj. p.value = 0.8624298

# -> return to default gene-level statistic Signal2Noise ratio



######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.60986155
# alternative 3: exponent 2 -> adj. p-value = 0.6089982

# final results: adj. p-value = 0.6089982
# achieved with
# ALTERNATIVE pre-filtering using filterByExpr
# ALTERNATIVE exponent 2




################################################################################
### Phenotype Permutation 5 ####################################################
################################################################################

i <- 5


#########
# step 1: default
#########

# -> adj. p-value = 0.8983813

#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.89139265
# ->> proceed with alternative transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen5 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,5]) %>%
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen5 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen5,],
                                           phenotype_labels = phen_bottomly[,5])

### export
path_filterByExpr_phen5 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen5,
            file = path_filterByExpr_phen5,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.8857586
# -> proceed with alternative filtering using filterByExpr()


########
# step 6: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.88453394
# alternative 2: Difference of Classes -> adj. p-value = 0.95386595

# -> proceed with alternative gene-level statistic t-Test


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.9454567
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.8309977

# final results: adj. p-value = 0.8309977
# achieved with
# ALTERNATIVE transformation method vst
# ALTERNATIVE pre-filtering using filterByExpr
# ALTERNATIVE gene-level statistic t-Test
# ALTERNATIVE exponent 2


################################################################################
### Phenotype Permutation 6 ####################################################
################################################################################

i <- 6


#########
# step 1: default
#########

# -> adj. p-value = 1


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9998107
# ->> proceed with alternative transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen6 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,6]) %>%
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen6 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen6,],
                                           phenotype_labels = phen_bottomly[,6])

### export
path_filterByExpr_phen6 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen6,
            file = path_filterByExpr_phen6,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# ->> adj. p-value = 1
# -> return to default pre-filtering


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 1
# alternative 2: Difference of Classes -> adj. p-value = 0.9124334

# -> proceed with alternative gene-level statistic Difference of Classes

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.8656485

### final results: adj. p-value = 0.8656485
# achieved with
# ALTERNATIVE transformation method vst
# ALTERNATIVE gene-level statistic Difference of Classes
# ALTERNATIVE exponent 2






################################################################################
### Phenotype Permutation 7 ####################################################
################################################################################

i <- 7

#########
# step 1: default
#########

# -> adj. p-value = 1



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9995487
# -> proceed with alternative transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen7 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,7]) %>%
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen7 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen7,],
                                           phenotype_labels = phen_bottomly[,7])

### export
path_filterByExpr_phen7 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen7,
            file = path_filterByExpr_phen7,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# ->  adj. p-value = 0.9836634
# ->> proceed with alternative pre-filtering method filterByExpr()


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.9823875
# alternative 2: Difference of Classes -> adj. p-value = 0.9488665


# -> proceed with alternative gene-level statistic Difference of Classes

######
# step 6: change exponent
######

# alternative 1: exponent 0 ->  adj. p-value = 0.9373173
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 1

### end results: adj. p-value = 0.9373173
# achieved with
# ALTERNATIVE transformatio method vst
# ALTERNATIVE pre-filtering using filterByExpr()
# ALTERNATIVE gene-level statistic Difference of Classes
# ALTERNATIVE exponent 0

################################################################################
### Phenotype Permutation 8 ####################################################
################################################################################

i <- 8

#########
# step 1: default
#########

# -> adj. p-value = 0.91667956



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.7369642
# -> proceed with alternative transformation method vst



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen8 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,8]) %>%
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen8 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen8,],
                                           phenotype_labels = phen_bottomly[,8])

### export
path_filterByExpr_phen8 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen8,
            file = path_filterByExpr_phen8,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.6223146
# ->> proceed with alternative pre-filtering using filterByExpr()


########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.62395775
# alternative 2: Difference of Classes -> adj. p-value =0.73889905

# ->> return to default gene-level statistic Signal2Noise ratio

######
# step 5: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.7731031
# alternative 2: exponent 1.5 -> adj. p-value = 0.4990408
# alternative 3: exponent 2 -> adj. p-value = 0.5296165


# final results: adj p-value = 0.4990408
# Achieved with
# ALTERNATIVE transformation method vst
# ALTERNATIVE pre-filtering using filterByExpr
# ALTERNATIVE exponent 1.5




################################################################################
### Phenotype Permutation 9 ####################################################
################################################################################

i <- 9

#########
# step 1: default
#########

# -> adj. p-value = 1

#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9934711
# -> proceed with alternative transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen9 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,9]) %>%
                      filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen9 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen9,],
                                           phenotype_labels = phen_bottomly[,9])

### export
path_filterByExpr_phen9 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen9,
            file = path_filterByExpr_phen9,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.8224765
# ->> proceed with alternative pre-filtering using filterByExpr()



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.8261976
# alternative 2: Difference of Classes -> adj. p-value = 0.7467399

# ->> proceed with alternativ gene-level statistic Difference of Classes




######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.62139297
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.6309108

# final results: adjusted p-value = 0.62139297
# achieved with
# ALTERNATIVE transformation method vst
# ALTERNATIVE pre-filtering using filterByExpr()
# ALTERNATIVE gene-level statistic Difference of Classes
# ALTERNATIVE exponent 1.5





################################################################################
### Phenotype Permutation 10 ####################################################
################################################################################

i <- 10


#########
# step 1: default
#########

# -> adj. p-value = 0.77643055


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9948916
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen10 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[,10]) %>%
                       filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen10 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen10,],
                                     phenotype_labels = phen_bottomly[,10])

### export
path_filterByExpr_phen10 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/p_value/Data_tCell/Raw/Phen_Permutation",
                                   i,
                                   "/exprdat_filterByExpr_phen_permutation",
                                   i,
                                   ".txt")

write.table(exprdat_prefilt_phen10,
            file = path_filterByExpr_phen10,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.85818964
# -> return to default pre-filtering approach



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.81606025
# alternative 2: Difference of Classes -> adj. p-value = 0.9881301

# ->> return to default gene-level statistic Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.84359455
# alternative 2: exponent 1.5 -> adj. p-value = 0.40972874
# alternative 3: exponent 2 -> adj. p-value = 0.89952475


# final results: adj. p-value = adj. p-value = 0.40972874
# achieved with
# ALTERNATIVE exponent 1.5










