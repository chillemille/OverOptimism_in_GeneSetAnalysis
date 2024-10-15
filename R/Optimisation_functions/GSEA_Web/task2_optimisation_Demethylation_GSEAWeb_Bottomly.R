#############################################################################################################################
### Optimization of GSEA Web Application for Pickrell data set: Mimimize adjusted p-value of gene set "Demethylation" (GO:BP)
#############################################################################################################################

library(edgeR) # for built-in pre-filtering
library(dplyr)



# load Bottomly data set along with the true and permuted sample conditions
source("./R/Help_functions/PreProcessing_Functions.R")

# load the required pre-processing functions
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")

# load functions that perform a transformation of the RNA-Seq data
source("./R/Help_functions/RNASeq_Transformation.R")


################################################################################
### Generate required ordner structure #########################################
################################################################################

dir.create("./Results/Intermediate_results/GSEA_Web")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Prep")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Prep/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Prep/Phenotypes")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phenotypes")



for(i in 1:10){

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Prep/Phen_Permutation", i))

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation", i))



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

path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phenotypes/Phenotype_original.txt")

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
path_dat_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Original/exprdat_default_phen_original.txt")

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
path_dat_vst_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Original/exprdat_vst_phen_original.txt")


# export
write.table(dat_vst_phenorig,
            file = path_dat_vst_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)



#########
# step 1: default
#########

# -> adjusted p-value = 0.9703367


#########
# step 2: vst-transformed gene expression data set
#########

# -> adjusted p-value = 0.9636415
# ->> proceed with alternative variancestabilizingtransformation

#########
# step 3: pre-filtering using edgeR's builtin function filterByExpr (based on default)
#########


# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phenorig <-  DGEList(Biobase::exprs(bottomly.eset), group = bottomly.eset$strain) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phenorig <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phenorig, ],
                                              phenotype_labels = bottomly.eset$strain)

### export
path_filterByExpr_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Original/exprdat_filterByExpr_phenorig.txt")

write.table(exprdat_prefilt_phenorig,
            file = path_filterByExpr_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# -> adjusted p-value = 0.98552865
# ->> return to default pre-filtering (vst-transformed data)



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.960218
# alternative 2: Difference of Classes -> adj. p-value = 0.9993506

# ->> proceed using alternative ranking metric t-Test

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.9368053
# alternative 2: exponent 1.5 -> adj. p-value = 0.94885767
# alternative 3: exponent 2 -> adj. p-value = 0.8922197


### final result: adj. p-value = 0.8922197
# achieved with
# alternative transformation method vst
# alternative ranking metric t-statistic
# exponent 2

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
  phen[phen_bottomly[, i] == levels(bottomly.eset$strain)[[1]]] <- 0
  phen[phen_bottomly[, i] == levels(bottomly.eset$strain)[[2]]] <- 1

  path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phenotypes/Phenotype_Permutation",
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
    voom_trans(phenotype_labels = phen_bottomly[, i])

  # generate path
  path_dat <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                     i,
                     "/exprdat_default_phen_permutation", i, ".txt")

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
    variancetransform(phenotype_labels = phen_bottomly[, i])

  # generate path
  path_dat_vst <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
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

# -> adj. p-value = 0.8477939



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 1
# -> return to default data set transformed using voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen1 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 1]) %>%
  filterByExpr()
# perform vst-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen1 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen1, ],
                                    phenotype_labels = phen_bottomly[, 1])

### export
path_filterByExpr_phen1 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen1,
            file = path_filterByExpr_phen1,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.6898383
# -> proceed with alternative filtering using filterByExpr(voom-transformed data set)




########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.68992144
# alternative 2: Difference of Classes  -> adj. p-value = 0.78695667

# ->> proceed with default ranking metric signal2noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.5555019
# alternative 3: exponent 2 -> adj. p-value = 0.7324062


### final results: adj. p-value = 0.5555019
# achieved with
# ALTERNATIVE filtering using filterByExpr
# ALTERNATIVE exponent 1.5




################################################################################
### Phenotype Permutation 2 ####################################################
################################################################################

i <- 2

#########
# step 1: default
#########

# -> adj. p-value = 0.95900315


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.96648276
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen2 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 2]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen2 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen2, ],
                                    phenotype_labels = phen_bottomly[, 2])

### export
path_filterByExpr_phen2 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen2,
            file = path_filterByExpr_phen2,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.99573845
# -> return to default pre-filtering

########
# step 6: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.9598944
# alternative 2: Difference of Classes -> adj. p-value =  0.9930457

# -> return to default ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj- p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.9975129


# final results: adj p-value = 0.95900315
# default results coincide with optimal results




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
# ->> return to default transformation voom



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr (on vst-transformed gene expression data set)
########


# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen3 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 3]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen3 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen3, ],
                                    phenotype_labels = phen_bottomly[, 3])

### export
path_filterByExpr_phen3 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen3,
            file = path_filterByExpr_phen3,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 1
# --> return to default filtering


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 1
# alternative 2: Difference of Classes -> adj. p-value = 1

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.9784683
# alternative 2: exponent 1.5 -> adj. p-value = 0.8153635
# alternative 3: exponent 2 -> adj. p-value = 0.97203565


### final results: adj p-value = 0.8153635
# achieved with alternative exponent 1.5




################################################################################
### Phenotype Permutation 4 ####################################################
################################################################################

i <- 4


#########
# step 1: default
#########

# -> adj. p-value = 0.562401


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.6947677
# ->> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen4 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 4]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen4 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen4, ],
                                    phenotype_labels = phen_bottomly[, 4])

### export
path_filterByExpr_phen4 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                  i, "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen4,
            file = path_filterByExpr_phen4,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.39904216
# ->> proceed with alternative pre-filtering using filterByExpr()




########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.41587377
# alternative 2: Difference of Classes -> adj. p.value = 0.5877484
# -> return to default gene-level ranking metric Signal2Noise ratio


# ->>



######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.9460303
# alternative 2: exponent 1.5 -> adj. p-value = 0.30026412
# alternative 3: exponent 2 -> adj. p-value = 0.63662845

# final results: adj. p-value = 0.30026412
# achieved with
# Alternative pre-filtering
# Alternative exponent 1.5
#
#



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
prefilt_ind_phen5 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 5]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen5 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen5, ],
                                    phenotype_labels = phen_bottomly[, 5])

### export
path_filterByExpr_phen5 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen5,
            file = path_filterByExpr_phen5,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.96100223
# -> proceed with alternative filtering using filterByExpr()


########
# step 6: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value =0.96057487
# alternative 2: Difference of Classes -> adj. p-value = 0.9977999

# -> proceed with alternative gene-level statistic t-Test


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.87200683
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value =  1

# final results: adj. p-value = 0.87200683
# -> achieved with
# ALTERNATIVE pre-filtering using filterByExpr
# ALTERNATIVE exponent 0


################################################################################
### Phenotype Permutation 6 ####################################################
################################################################################

i <- 6


#########
# step 1: default
#########

# -> adj. p-value = 0.95976955


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9025592
# ->> proceed with alternative transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen6 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 6]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen6 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen6, ],
                                           phenotype_labels = phen_bottomly[, 6])

### export
path_filterByExpr_phen6 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen6,
            file = path_filterByExpr_phen6,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# ->> adj. p-value = 0.985849
# -> return to default pre-filtering


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value =  0.9030102
# alternative 2: Difference of Classes -> adj. p-value = 0.9401092


# -> returun to default gene-level statistic Signal2NoiseRatio

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.9824037
# alternative 3: exponent 2 -> adj. p-value = 0.6045143

### final results: adj. p-value = 0.6045143
# achieved with
# ALTERNATIVE transformation method vst
# ALTERNATIVE exponent 2






################################################################################
### Phenotype Permutation 7 ####################################################
################################################################################

i <- 7

#########
# step 1: default
#########

# -> adj. p-value = 0.91576385



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9618282
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen7 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 7]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen7 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen7, ],
                                    phenotype_labels = phen_bottomly[, 7])

### export
path_filterByExpr_phen7 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen7,
            file = path_filterByExpr_phen7,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# ->  adj. p-value = 0.6725676
# ->> proceed with alternative pre-filtering using filterByExpr()


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.6693831
# alternative 2: Difference of Classes -> adj. p-value = 0.736255


# -> proceed with alternative gene-level ranking metric t-Test

######
# step 6: change exponent
######

# alternative 1: exponent 0 ->  adj. p-value = 0.8311909
# alternative 2: exponent 1.5 -> adj. p-value = 0.84496665
# alternative 3: exponent 2 -> adj. p-value = 0.7271227

### end results: adj. p-value = 0.6693831
# Achieved with
# ALTERNATIVE pre-filtering using filterByExpr()
# ALTERNATIVE gene-level statistic t-Test

################################################################################
### Phenotype Permutation 8 ####################################################
################################################################################

i <- 8

#########
# step 1: default
#########

# -> adj. p-value = 0.7584467



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.7210969
# -> proceed with alternative transformation method vst



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen8 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 8]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen8 <- variancetransform(Biobase::exprs(bottomly.eset)[prefilt_ind_phen8, ],
                                           phenotype_labels = phen_bottomly[, 8])

### export
path_filterByExpr_phen8 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen8,
            file = path_filterByExpr_phen8,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.91942537
# ->> return to default pre-filtering


########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.7208875
# alternative 2: Difference of Classes -> adj. p-value =0.85748076

# ->> proceed with alternative gene-level statistic t-Test

######
# step 5: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.91090083
# alternative 2: exponent 1.5 -> adj. p-value = 0.5801879
# alternative 3: exponent 2 -> adj. p-value = 0.6889628

# final results: adj p-value = 0.5801879
# Achieved with
# ALTERNATIVE transformation method vst
# ALTERNATIVE gene-level ranking metric t-Test
# ALTERNATIVE exponent 1.5




################################################################################
### Phenotype Permutation 9 ####################################################
################################################################################

i <- 9

#########
# step 1: default
#########

# -> adj. p-value = 0.8126698

#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.88228816
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen9 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 9]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen9 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen9, ],
                                    phenotype_labels = phen_bottomly[, 9])

### export
path_filterByExpr_phen9 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                  i,
                                  "/exprdat_filterByExpr_phen_permutation",
                                  i,
                                  ".txt")

write.table(exprdat_prefilt_phen9,
            file = path_filterByExpr_phen9,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value =1
# ->> return to default pre-filtering



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.8018664
# alternative 2: Difference of Classes -> adj. p-value = 0.70637554

# ->> proceed with alternative gene-level ranking metric Difference of Classes




######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.839664
# alternative 2: exponent 1.5 -> adj. p-value = 0.92351174
# alternative 3: exponent 2 -> adj. p-value = 1

# final results: adjusted p-value = 0.70637554
# achieved with
# ALTERNATIVE gene-level ranking metric Difference of classes





################################################################################
### Phenotype Permutation 10 ####################################################
################################################################################

i <- 10


#########
# step 1: default
#########

# -> adj. p-value = 0.65629333


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.77514607
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen10 <-  DGEList(Biobase::exprs(bottomly.eset), group = phen_bottomly[, 10]) %>%
  filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen10 <- voom_trans(Biobase::exprs(bottomly.eset)[prefilt_ind_phen10, ],
                                     phenotype_labels = phen_bottomly[, 10])

### export
path_filterByExpr_phen10 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/Data_task2_Demethylation/Raw/Phen_Permutation",
                                   i,
                                   "/exprdat_filterByExpr_phen_permutation",
                                   i,
                                   ".txt")

write.table(exprdat_prefilt_phen10,
            file = path_filterByExpr_phen10,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.74814177
# -> return to default pre-filtering approach



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.6606892
# alternative 2: Difference of Classes -> adj. p-value =0.8534896

# ->> return to default gene-level ranking metric signal2noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.54858625
# alternative 2: exponent 1.5 -> adj. p-value = 0.30786166
# alternative 3: exponent 2 -> adj. p-value = 0.74843526


# final results: adj. p-value = adj. p-value = 0.30786166
# achieved with
# ALTERNATIVE exponent 1.5










