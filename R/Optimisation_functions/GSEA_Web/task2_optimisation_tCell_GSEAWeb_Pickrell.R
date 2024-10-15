#########################################################################################################################################
### Optimization of GSEA Web Application for Pickrell data set: Mimimize adjusted p-value of gene set "t-cell mediated immunity" (GO:BP)
#########################################################################################################################################

library(edgeR) # for builtin pre-filtering
library(dplyr)


# for this script, the following additional script RNASeq_Transformation.R is required


# load Pickrell data set along with the true and permuted sample conditions
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")

# load the required pre-processing functions
source("./R/Help_functions/PreProcessing_Functions.R")

# load functions that perform a transformation of the RNA-Seq data
source("./R/Help_functions/RNASeq_Transformation.R")


################################################################################
### Generate required ordner structure #########################################
################################################################################

dir.create("./Results/Intermediate_results/GSEA_Web")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Prep")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Prep/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Prep/Phenotypes")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phenotypes")



for(i in 1:10){

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Prep/Phen_Permutation", i))

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i))



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
phen_orig[ pickrell.eset$gender == levels(pickrell.eset$gender)[[1]]] <- 0
phen_orig[ pickrell.eset$gender == levels(pickrell.eset$gender)[[2]]] <- 1

path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phenotypes/Phenotype_original.txt")

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
path_dat_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Original/exprdat_default_phen_original.txt")

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
path_dat_vst_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Original/exprdat_vst_phen_original.txt")


# export
write.table(dat_vst_phenorig,
            file = path_dat_vst_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)








#########
# step 1: default
#########

# -> adj. p-value = 0.8812905

#########
# step 2: vst-transformed gene expression data set
#########

# -> adjusted p-value = 0.45906654
# ->> proceed with vst as RNA-Seq transformation method

#########
# step 3: pre-filtering using edgeR's builtin function filterByExpr (based on default)
#########


# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phenorig <-  DGEList(Biobase::exprs(pickrell.eset), group = pickrell.eset$gender) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phenorig <- variancetransform(Biobase::exprs(pickrell.eset)[prefilt_ind_phenorig, ], phenotype_labels = pickrell.eset$gender)

### export
path_filterByExpr_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Original/exprdat_filterByExpr_phenorig.txt")

write.table(exprdat_prefilt_phenorig,
            file = path_filterByExpr_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# -> adjusted p-value =  0.44699678
# ->> proceed with filtering based on filterByExpr() (based on vst-transformation)



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.44261065
# alternative 2: Difference of Classes -> adj. p-value = 0.45134816

# ->> proceed with alternative gene-level ranking metric t-Test

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.60875833
# alternative 2: exponent 1.5 -> adj. p-value = 0.42853522
# alternative 3: exponent 2 -> adj. p-value = 0.39731848


### final result: adj. p-value = final results: adj. p-value = 0.39731848
# achieved with
# alternative RNA-Seq transformation method vst
# alternative pre-filtering using filterByExpr
# alternative gene-level ranking metric tTest
# alternative exponent 2

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
  phen[phen_pickrell[, i] == "female"] <- 0
  phen[phen_pickrell[, i] == "male"] <- 1

  path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phenotypes/Phenotype_Permutation", i, ".txt")

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
    voom_trans(phenotype_labels = phen_pickrell[, i])

  # generate path
  path_dat <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_default_phen_permutation", i, ".txt")

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
    variancetransform(phenotype_labels = phen_pickrell[, i])

  # generate path
  path_dat_vst <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_vst_phen_permutation", i, ".txt")


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

# -> adj. p-value = 0.8126426



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.68386286
# -> proceed with vst-transformed RNA-Seq data set


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen1 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 1]) %>% filterByExpr()
# perform vst-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen1 <- variancetransform(Biobase::exprs(pickrell.eset)[prefilt_ind_phen1, ], phenotype_labels = phen_pickrell[, 1])

### export
path_filterByExpr_phen1 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen1,
            file = path_filterByExpr_phen1,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.60577023
# -> proceed with pre-filtering using filterByExpr() (based on vst-transformed gene expression data set)



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.6104824
# alternative 2: Difference of Classes  -> adj. p-value = 0.65454173

# ->> return to default ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.98207676
# alternative 2: exponent 1.5 -> adj. p-value = 0.4480666
# alternative 3: exponent 2 -> adj. p-value = 0.38854122



### final results: adj. p-value = 0.38854122
# achieved with
# alternative RNA-Seq transformation vst
# alternative pre-filtering using filterByExpr()
# alternative exponent 2




################################################################################
### Phenotype Permutation 2 ####################################################
################################################################################

i <- 2

#########
# step 1: default
#########

# -> adj. p-value = 0.9996207



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 1
# -> return to voom transformation


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen2 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 2]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen2 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen2, ], phenotype_labels = phen_pickrell[, 2])

### export
path_filterByExpr_phen2 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen2,
            file = path_filterByExpr_phen2,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 1
# -> proceed to default pre-filtering (based on voom-transformation)

########
# step 6: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.9984517 (< 0.9996207 default adjusted p-value)
# alternative 2: Difference of Classes -> adj. p-value = 0.99975663

# -> proceed with alternative gene-level ranking metric t-Test


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.89641523
# alternative 2: exponent 1.5 -> adj- p-value = 0.82541806
# alternative 3: exponent 2 -> adj. p-value = 1


# final results: adj p-value = 0.82541806
# achieved with
# alternative gene-level ranking metric t-Test
# alternative exponent 1.5




################################################################################
### Phenotype Permutation 3 ####################################################
################################################################################

i <- 3


#########
# step 1: default
#########

# -> adj p-value = 0.98444533



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj p-value = 0.9710071
# proceed with RNA-Seq transformation method vst



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr (on vst-transformed gene expression data set)
########


# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen3 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 3]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen3 <- variancetransform(Biobase::exprs(pickrell.eset)[prefilt_ind_phen3, ], phenotype_labels = phen_pickrell[, 3])

### export
path_filterByExpr_phen3 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen3,
            file = path_filterByExpr_phen3,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 1
## -> return to default pre-filtering (based on vst-transformation)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.9711163
# alternative 2: Difference of Classes -> adj. p-value = 0.87568724

# proceed with alternative gene-level ranking metric Difference of Classes

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.82210153
# alternative 2: exponent 1.5 -> adj. p-value = 0.77129024
# alternative 3: exponent 2 -> adj. p-value = 0.841332


### final results: adj p-value = 0.77129024
# achieved with
# ALTERNATIVE ranking metric Difference of Classes
# ALTERNATIVE exponent 1.5




################################################################################
### Phenotype Permutation 4 ####################################################
################################################################################

i <- 4


#########
# step 1: default
#########

# -> adj. p-value = 0.95483476


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.7389655
# -> proceed with alternative RNA-Seq transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen4 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 4]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen4 <- variancetransform(Biobase::exprs(pickrell.eset)[prefilt_ind_phen4, ], phenotype_labels = phen_pickrell[, 4])

### export
path_filterByExpr_phen4 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen4,
            file = path_filterByExpr_phen4,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.7456053
# ->> return to default pre-filtering (based on vst-transformation)




########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.7332202
# alternative 2: Difference of Classes -> adj. p.value = 0.65746826

# -> proceed with alternative gene-level ranking metric Difference of Classes





######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.71072996
# alternative 2: exponent 1.5 -> adj. p-value = 0.8254908
# alternative 3: exponent 2 -> adj. p-value = 0.6730631

# final results: adj. p-value = 0.65746826
# achieved with
# alternative RNA-Seq transformation method vst
# alternative gene-level ranking metric Difference of Classes




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
prefilt_ind_phen5 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 5]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen5 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen5, ], phenotype_labels = phen_pickrell[, 5])

### export
path_filterByExpr_phen5 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

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

# alternative 1: t-Test -> adj. p-value =  1
# alternative 2: Difference of Classes -> adj. p-value = 1

# -> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.94383186
# alternative 2: exponent 1.5 -> adj. p-value = 0.75002885
# alternative 3: exponent 2 -> adj. p-value = 1

# final results: adj. p-value = 0.75002885
# achieved with alternative exponent 1.5


################################################################################
### Phenotype Permutation 6 ####################################################
################################################################################

i <- 6


#########
# step 1: default
#########

# -> adj. p-value = 0.7248434


#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9864163
# ->> return to default RNA-Seq transformation voom



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen6 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 6]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen6 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen6, ], phenotype_labels = phen_pickrell[, 6])

### export
path_filterByExpr_phen6 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen6,
            file = path_filterByExpr_phen6,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# ->> adj. p-value = 1
# ->> return to default pre-filtering method (based on RNA-Seq transformation voom)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.72703534
# alternative 2: Difference of Classes -> adj. p-value = 0.84167004

#->> return to default gene-level Ranking metric Signal2Noise Ratio



######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 0.96868414
# alternative 3: exponent 2 -> adj. p-value = 0.6569779

### final results: adj. p-value = 0.656977
# achieved with alternative exponent 2






################################################################################
### Phenotype Permutation 7 ####################################################
################################################################################

i <- 7

#########
# step 1: default
#########

# -> adj. p-value = 0.95434326



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.9620933
# -> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen7 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 7]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen7 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen7, ], phenotype_labels = phen_pickrell[, 7])

### export
path_filterByExpr_phen7 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen7,
            file = path_filterByExpr_phen7,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# ->  adj. p-value = 1
# ->> return to default pre-filtering (based on voom-transformation)


########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.95436406
# alternative 2: Difference of Classes -> adj. p-value = 0.95512825


# -> return to default gene-level ranking metric Signal2Noise Ratio

######
# step 6: change exponent
######

# alternative 1: exponent 0 ->  adj. p-value = 0.9845177
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.87398493

### end results: adj. p-value = 0.87398493
# Achieved with alternative exponent 2

################################################################################
### Phenotype Permutation 8 ####################################################
################################################################################

i <- 8

#########
# step 1: default
#########

# -> adj. p-value = 0.7738751



#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 0.92897844
# -> return to default RNA-Seq transformation voom



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen8 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 8]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen8 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen8, ], phenotype_labels = phen_pickrell[, 8])

### export
path_filterByExpr_phen8 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen8,
            file = path_filterByExpr_phen8,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 0.92982554
# ->> return to default pre-filtering (based on default RNA-Seq transformation voom)


########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.7694006
# alternative 2: Difference of Classes -> adj. p-value = 0.82668114

# ->> proceed with alternative gene-level ranking metric t-Test

######
# step 5: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 0.9685064
# alternative 2: exponent 1.5 -> adj. p-value = 0.43834278
# alternative 3: exponent 2 -> adj. p-value = 0.68722975

# final results: adj p-value = 0.43834278
# achieved with alternative gene-level ranking metric and alternative exponent 1.5




################################################################################
### Phenotype Permutation 9 ####################################################
################################################################################

i <- 9

#########
# step 1: default
#########

# -> adj. p-value = 0.9857087

#########
# step 2: vst-transformed gene expression data set
#########

# -> adj. p-value = 1
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen9 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 9]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen9 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen9, ], phenotype_labels = phen_pickrell[, 9])

### export
path_filterByExpr_phen9 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen9,
            file = path_filterByExpr_phen9,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 1
# ->> return to default pre-filtering (based on RNA-Seq transformation voom)



########
# step 4: change gene-level ranking metric
########

# alternative 1: t-Test -> adj. p-value = 0.9852834
# alternative 2: Difference of Classes -> adj. p-value = 0.9385907

# ->> proceed with alternative gene-level ranking metric Difference of Classes




######
# step 6: change exponent
######

# alternative 1: exponent 0 -> adj. p-value = 1
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 0.77723676

# final results: adjusted p-value = 0.77723676
# achieved with
# alternative gene-level ranking metric t-Test
# alternative exponent 2





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
prefilt_ind_phen10 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 10]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen10 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen10, ], phenotype_labels = phen_pickrell[, 10])

### export
path_filterByExpr_phen10 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task2_tCell/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen10,
            file = path_filterByExpr_phen10,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> adj. p-value = 1
# -> return to default pre-filtering (based on RNA-Seq transformation method voom)



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
# alternative 2: exponent 1.5 -> adj. p-value = 1
# alternative 3: exponent 2 -> adj. p-value = 1


# final results: adj. p-value = 1
# default results coincide with optimal results










