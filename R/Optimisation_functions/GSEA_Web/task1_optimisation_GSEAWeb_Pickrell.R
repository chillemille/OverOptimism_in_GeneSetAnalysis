################################################################################
### Optimization of GSEA Web Application for Pickrell data set #################
################################################################################

# for this script, the following additional script RNASeq_Transformation.R is required

library(edgeR) # for filterByExpr()

# load Bottomly data set
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")
# load data preprocessing functions
source("./R/Help_functions/RNASeq_Transformation.R")


################################################################################
### Generate folder structure ##################################################
################################################################################

dir.create("./Results/Intermediate_results/GSEA_Web")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Prep")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Prep/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Prep/Phenotypes")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phenotypes")

for(i in 1:10){

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Prep/Phen_Permutation", i))
  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i))
}




#######################################################################################
###pre-filtering function #############################################################
#######################################################################################

pre_filt <- function(expression_data, threshold){

  expression_data_filt<-expression_data[rowSums(expression_data)>=threshold, ]

  # return filtered gene expression data set
  return(expression_data_filt)


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

path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phenotypes/Phenotype_original.txt")

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
path_dat_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Original/exprdat_default_phen_original.txt")

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
path_dat_vst_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Original/exprdat_vst_phen_original.txt")


# export
write.table(dat_vst_phenorig,
            file = path_dat_vst_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)




#########
# step 1: default
#########

# -> 6 DEGS

#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# ->> return to RNA-Seq transformation method voom

#########
# step 3: pre-filtering using edgeR's builtin function filterByExpr
#########

# -> 10 DEGS
# ->> proceed with alternative pre-filtering using filterByExpr

########
# step 4: change gene set database to KEGG
########

# -> 1 DEGS
# -> return to GO with subontology BP

########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 11 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> proceed with alternative gene-level ranking t-Test

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 38 DEGS
# alternative 3: exponent 2 -> 41 DEGS


### final result: 41 DEGS

# achieved with
# alternative pre-filtering using FilterByExpr
# alternative gene-level statistic t-Test
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

  path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phenotypes/Phenotype_Permutation", i, ".txt")

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
  path_dat <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_default_phen_permutation", i, ".txt")

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
  path_dat_vst <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_vst_phen_permutation", i, ".txt")


  # export
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

# -> 1 DEGS

#########
# step 2: vst-transformed gene expression data set
#########

# -> 1 DEGS
# -> return to voom transformation


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen1 <- DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 1]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen1 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen1, ], phenotype_labels = phen_pickrell[, 1])

### export
path_filterByExpr_phen1 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen1,
            file = path_filterByExpr_phen1,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 3 DEGS
# -> proceed with voom transformed and filterByExpr()-filtered gene expression data set


########
# step 4: change gene set database to KEGG
########

# -> 2 DEGS
# ->> return to geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 3 DEGS
# alternative 2: Difference of Classes  -> 0 DEGS

# ->> return to default ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 6 DEGS
# alternative 3: exponent 2 -> 18 DEGS

### final results: 18 DEGS
# - default RNA-Seq transformation voom
# - alternative pre-filtering using filterByExpr()
# - default geneset database GO (BP)
# - default gene-level ranking metric Signal2Noise
# - alternative exponent 2



################################################################################
### Phenotype Permutation 2 ####################################################
################################################################################

i <- 2

#########
# step 1: default
#########

# -> 79 DEGS

#########
# step 2: vst-transformed gene expression data set
#########

# -> 1 DEGS
# -> return to voom transformation


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen2 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 2]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen2 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen2, ], phenotype_labels = phen_pickrell[, 2])

### export
path_filterByExpr_phen2 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen2,
            file = path_filterByExpr_phen2,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 3 DEGS
# -> return to manual pre-filtering with threshold 10


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# ->> return to default gene set database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 80 DEGS
# alternative 2: Difference of Classes -> 61 DEGS

# -> proceed with t-Test (t-statistic) as gene-level ranking metric


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 32 DEGS
# alternative 2: exponent 1.5 -> 8 DEGS
# alternative 3: exponent 2 -> 6 DEGS

# final results: 80 DEGS
# default RNA-Seq transformation voom
# default pre-filtering using manual pre-filtering with threshold 10
# default gene set database GO (BP)
# alternative gene-level ranking metric t-test (t-statistic)
# default exponent 1


################################################################################
### Phenotype Permutation 3 ####################################################
################################################################################

i <- 3


#########
# step 1: default
#########

# -> 0 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 1 DEGS
# ->> proceed with vst-transformed RNA-Seq data set



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr (on vst-transformed gene expression data set)
########


# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen3 <- DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 3]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen3 <- variancetransform(Biobase::exprs(pickrell.eset)[prefilt_ind_phen3, ], phenotype_labels = phen_pickrell[, 3])

### export
path_filterByExpr_phen3 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen3,
            file = path_filterByExpr_phen3,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 0 DEGS
# --> return to manual pre-filtering with threshold 10 (vst-transformed data set)


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# ->> return to geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 1 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# -> return to default gene-level ranking metric Signal2Noise ratio



######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 1 DEGS
# alternative 3: exponent 2 -> 2 DEGS

# final results: 2 DEGS
# alternative RNA-Seq transformation method vst
# default geneset database GO (BP)
# default gene-level ranking metric Signal2Noise ratio
# alternative exponent 2


################################################################################
### Phenotype Permutation 4 ####################################################
################################################################################

i <- 4


#########
# step 1: default
#########

# -> 2 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# ->> return to default transformation voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen4 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 4]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen4 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen4, ], phenotype_labels = phen_pickrell[, 4])

### export
path_filterByExpr_phen4 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen4,
            file = path_filterByExpr_phen4,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 0 DEGS
# ->> return to default pre-filtering using threshold 10 (on voom-transformed gene expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 1 DEGS
# ->> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 3 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> proceed with t-Test (t-statistic) as gene-level ranking metric



######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 1 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

### end results: 3 DEGS
# default RNA-Seq transformation method voom
# default pre-filtering with threshold 10
# default geneset database KEGG
# alternative gene-level ranking metric tTest
# default exponent 1


################################################################################
### Phenotype Permutation 5 ####################################################
################################################################################

i <- 5


#########
# step 1: default
#########

# -> 3 DEGS

#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# ->> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen5 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 5]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen5 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen5, ], phenotype_labels = phen_pickrell[, 5])

### export
path_filterByExpr_phen5 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen5,
            file = path_filterByExpr_phen5,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 2 DEGS
# -> return to default pre-filtering using treshold 10 (voom - transformation)

########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 3 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# -> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

# final results: 3 DEGS
# -> optimal parameter choice coincides with default parameter choice
# -> number of differentially enriched gene sets could not be increased


################################################################################
### Phenotype Permutation 6 ####################################################
################################################################################

i <- 6


#########
# step 1: default
#########

# -> 0 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# ->> return to voom-transformed gene expression data set


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen6 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 6]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen6 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen6, ], phenotype_labels = phen_pickrell[, 6])

### export
path_filterByExpr_phen6 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen6,
            file = path_filterByExpr_phen6,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# ->> 0 DEGS
# -> return to manual pre-filtering (using threshold 10) (voom-transformed expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default gene set database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 0 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# -> return to default gene-level ranking metric Signal2Noise ratio

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 1 DEGS


# final results: 1 DEGS
# default RNA-Seq transformation method voom
# default pre-filtering usign manual filtering with threshold 10
# default geneset database GO (BP)
# default gene-level ranking metric Signal2Noise ratio
# ALTERNATIVE exponent 2


################################################################################
### Phenotype Permutation 7 ####################################################
################################################################################

i <- 7

#########
# step 1: default
#########

# -> 155 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 13 DEGS
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen7 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 7]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen7 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen7, ], phenotype_labels = phen_pickrell[, 7])

### export
path_filterByExpr_phen7 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen7,
            file = path_filterByExpr_phen7,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)
# -> 5 DEGS
# ->> return to manual pre-filtering using threshold 10 (voom- transformed expression data set )


########
# step 4: change gene set database to KEGG
########

# -> 4 DEGS
# ->> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 152 DEGS
# alternative 2: Difference of Classes -> 118 DEGS

# -> return to default gene-level ranking metric signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 70 DEGS

### end results: 155 DEGS
# optimal parameter configuration coincides with default parameter configuration
# number of differentially enriched gene sets could not be increased


################################################################################
### Phenotype Permutation 8 ####################################################
################################################################################

i <- 8

#########
# step 1: default
#########

# -> 0 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# -> return to default transformation method voom



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen8 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 8]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen8 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen8, ], phenotype_labels = phen_pickrell[, 8])

### export
path_filterByExpr_phen8 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen8,
            file = path_filterByExpr_phen8,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 0 DEGS
# ->> return to default pre-filtering with threshold 10


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default gene set database GO (BP)

########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test ->0 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 40 DEGS
# alternative 3: exponent 2 -> 1 DEGS

### final results: 40 DEGS
# default RNA-Seq transformation method voom
# default pre-filtering with threshold 10
# default gene set database GO (BP)
# default gene-level ranking metric Signal2Noise ratio
# alternative exponent 1.5


################################################################################
### Phenotype Permutation 9 ####################################################
################################################################################

i <- 9

#########
# step 1: default
#########

# -> 20 DEGS

#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# -> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen9 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 9]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen9 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen9, ], phenotype_labels = phen_pickrell[, 9])

### export
path_filterByExpr_phen9 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen9,
            file = path_filterByExpr_phen9,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 0 DEGS
# ->> return to default pre-filtering with threshold 10


########
# step 4: change gene set database to KEGG
########

# -> 1 DEGS
# ->> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 17 DEGS
# alternative 2: Difference of Classes -> 6 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio




######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

# ->> end result: 20 DEGS
# optimal parameter configuration coincides with default configuration
# number of differentially enriched gene sets could not be increased




################################################################################
### Phenotype Permutation 10 ####################################################
################################################################################

i <- 10


#########
# step 1: default
#########

# -> 0 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# -> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen10 <-  DGEList(Biobase::exprs(pickrell.eset), group = phen_pickrell[, 10]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen10 <- voom_trans(Biobase::exprs(pickrell.eset)[prefilt_ind_phen10, ], phenotype_labels = phen_pickrell[, 10])

### export
path_filterByExpr_phen10 <- paste0("./Results/Intermediate_results/GSEA_Web/Pickrell/Data_task1/Raw/Phen_Permutation", i, "/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen10,
            file = path_filterByExpr_phen10,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 0 DEGS
# -> return to manual pre-filtering with threshold 10 (voom-transformed data set)


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# ->> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 0 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

# end result: 0 DEGS
# optimal parameter configuration corresponds to default configuration
# number of differentially enriched gene sets could not be increased






