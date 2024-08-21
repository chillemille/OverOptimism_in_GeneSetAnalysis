################################################################################
### Optimization of GSEA Web Application for Bottomly data set #################
################################################################################

library(edgeR) # for filterByExpr()

# load Bottomly data set
source("./R/Prepare_data_and_permutations/Random_Phenotype_Permutations.R")
# load data preprocessing functions
source("./R/Functions/RNASeq_Transformation.R")

################################################################################
### Generate folder structure ##################################################
################################################################################

dir.create("./Results/Intermediate_results/GSEA_Web")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Prep")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Prep/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Prep/Phenotypes")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Original")
dir.create("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phenotypes")

for(i in 1:10){

  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Prep/Phen_Permutation",i))
  dir.create(paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i))
}


################################################################################
###pre-filtering function ######################################################
################################################################################

pre_filt<-function(expression_data, threshold){

  expression_data_filt<-expression_data[rowSums(expression_data)>=threshold,]

  return(expression_data_filt)


}




#######################################################################################
### (I) Full Optimization (pre-processing and internal parameters) ####################
#######################################################################################



#######################################################################################
### Original Phenotype Assignment #####################################################
#######################################################################################


########################
### export required data
########################


##################
### phenotypes ###
##################

# convert levels ("female", "male") to 0 and 1
levels(bottomly.eset$strain)

phen_orig <- c()
phen_orig[bottomly.eset$strain == levels(bottomly.eset$strain)[1]] <- 0
phen_orig[bottomly.eset$strain == levels(bottomly.eset$strain)[2]] <- 1

path_phen_orig <- "./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phenotypes/Phenotype_Original.txt"

write.table(phen_orig,
            file = path_phen_orig,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

####################################
### gene expression data sets ######
####################################

# default pre-process bottomly data
dat_default_phenorig <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
  voom_trans(phenotype_labels = bottomly.eset$strain)

# generate path
path_dat_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Original/exprdat_default_phen_original.txt")

### export
write.table(dat_default_phenorig,
            file = path_dat_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

############################################
### vst transformed gene expression data set
############################################


# default pre-process bottomly data
dat_vst_phenorig <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
  variancetransform(phenotype_labels = bottomly.eset$strain)

# generate path
path_dat_vst_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Original/exprdat_vst_phen_original.txt")


# export
write.table(dat_vst_phenorig,
            file = path_dat_vst_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


####################
### Optimization ###
####################

#########
# step 1: default
#########

# -> 2 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# -> return to default transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phenorig <-  DGEList(bottomly.eset, group = bottomly.eset$strain) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phenorig <- voom_trans(bottomly.eset[prefilt_ind_phenorig,], phenotype_labels = bottomly.eset$strain)

### export
path_filterByExpr_phenorig <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Original/exprdat_filterByExpr_phen_original.txt")

write.table(exprdat_prefilt_phenorig,
            file = path_filterByExpr_phenorig,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 5 DEGS
# proceed with pre-filtering using filterByExpr()


########
# step 4: change gene set database to KEGG
########

# -> 1 DEGS
# ->> return to default gene set database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 5 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 50 +11 = 61 DEGS
# alternative 3: exponent 2 -> 157 + 108 = 265 DEGS

### final results: 265 DEGS
# achieved with
# alternative pre-filtering using filterByExpr()
# alternative exponent 2



#######################################################################################
### Random Permutations of Original Phenotype Assignment ##############################
#######################################################################################


# export phenotype permutations and gene expression data set pre-processed in default manner
for(i in 1:ncol(phen_bottomly)){

  ##############
  ### phenotypes
  ##############

  # convert levels ("female", "male") to 0 and 1
  phen <- c()
  phen[phen_bottomly[,i] == levels(bottomly.eset$strain)[1]] <- 0
  phen[phen_bottomly[,i] == levels(bottomly.eset$strain)[2]] <- 1

  path_phen <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phenotypes/Phenotype_Permutation",i, ".txt")

  write.table(phen,
              file = path_phen,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  ####################################
  ### default gene expression data set
  ####################################

  # default pre-process bottomly data
  dat_default <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
    voom_trans(phenotype_labels = phen_bottomly[,i])

  # generate path
  path_dat <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_default_phen_permutation", i, ".txt")

  ### export
  write.table(dat_default,
              file = path_dat,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)

  ############################################
  ### vst transformed gene expression data set
  ############################################


  # default pre-process bottomly data
  dat_vst <- pre_filt(Biobase::exprs(bottomly.eset), threshold = 10) %>% # default pre-filtering (manual filtering with threshold 10)
    variancetransform(phenotype_labels = phen_bottomly[,i])

  # generate path
  path_dat_vst <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_vst_phen_permutation", i, ".txt")


  # export
  write.table(dat_vst,
              file = path_dat_vst,
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE)


}




#######################################################################################
### Phenotype Permutation 1 ###########################################################
#######################################################################################

i <- 1

#########
# step 1: default
#########

# -> 3 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# -> return to default transformation method voom

########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen1 <-  DGEList(bottomly.eset, group = phen_bottomly[,1]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered pickrell data set
exprdat_prefilt_phen1 <- voom_trans(bottomly.eset[prefilt_ind_phen1,], phenotype_labels = phen_bottomly[,1])

### export
path_filterByExpr_phen1 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen1,
            file = path_filterByExpr_phen1,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# -> 0 DEGS
# ->> return to default pre-filtering using threshold 10 (voom-transformed expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 4 DEGS
# alternative 2: Difference of classes -> 2 DEGS

# ->> proceed with t-Test (t-statistic) as gene-level ranking metric


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

### end result: 4 DEGS
# default RNA-Seq transformation method voom
# default pre-filtering with threshold 10
# default gene set database GO (BP)
# ALTERNATIVE gene-level ranking metric t-Test
# default exponent 1



#######################################################################################
### Phenotype Permutation 2 ###########################################################
#######################################################################################

i <- 2

#########
# step 1: default
#########

# -> 1 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 5 DEGS
# -> proceed with transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen2 <-  DGEList(bottomly.eset, group = phen_bottomly[,2]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen2 <- variancetransform(bottomly.eset[prefilt_ind_phen2,], phenotype_labels = phen_bottomly[,2])

### export
path_filterByExpr_phen2 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen2,
            file = path_filterByExpr_phen2,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 10 DEGS
# proceed with pre-filtering using filterByExpr()


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 10 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 3 DEGS
# alternative 3: exponent 2 -> 6 DEGS

# end results: 10 DEGS
# alternative transformation method vst
# alternative pre-filtering using filterByExpr()



#######################################################################################
### Phenotype Permutation 3 ###########################################################
#######################################################################################

i <- 3


#########
# step 1: default
#########

# -> 0 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# return to default pre-filtering with threshold 10


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen3 <-  DGEList(bottomly.eset, group = phen_bottomly[,3]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen3 <-voom_trans(bottomly.eset[prefilt_ind_phen3,], phenotype_labels = phen_bottomly[,3])

### export
path_filterByExpr_phen3 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen3,
            file = path_filterByExpr_phen3,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# -> 0 DEGS
# ->> return to default pre-filtering with threshold 10


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 0 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# -> return to default gene-level ranking metric signal2noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

# ->> final result: 0 DEGS
# optimal parameter configuration coincides with default configuration
# initual number of differentially enriched gene sets could not be increased


#######################################################################################
### Phenotype Permutation 4 ###########################################################
#######################################################################################

i <- 4


#########
# step 1: default
#########

# -> 1 DEGS



#########
# step 2: vst-transformed gene expression data set
#########

# -> 2 DEGS
# ->> proceed with RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen4 <-  DGEList(bottomly.eset, group = phen_bottomly[,4]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen4 <- variancetransform(bottomly.eset[prefilt_ind_phen4,], phenotype_labels = phen_bottomly[,4])

### export
path_filterByExpr_phen4 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen4,
            file = path_filterByExpr_phen4,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 8 DEGS
# ->> proceed with pre-filtering using filterByExpr() (vst-transformed gene expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# ->> proceed with default gene set database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 8 DEGS
# alternative 2: Difference of Classes -> 10 DEGS

# ->> proceed with alternative gene-level ranking metric Difference of Classes

######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1 -> 13 DEGS
# alternative 3: exponent 2 -> 6 DEGS


# ->> end results: 13 DEGS
# alternative RNA-Seq transformation vst
# alternative pre-filtering using filterByExpr()
# alternative gene-level ranking metric Difference of Classes
# alternative exponent 1.5


#######################################################################################
### Phenotype Permutation 5 ###########################################################
#######################################################################################

i <- 5


#########
# step 1: default
#########

# -> 0 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# ->> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen5 <-  DGEList(bottomly.eset, group = phen_bottomly[,5]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen5 <-voom_trans(bottomly.eset[prefilt_ind_phen5,], phenotype_labels = phen_bottomly[,5])

### export
path_filterByExpr_phen5 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen5,
            file = path_filterByExpr_phen5,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# -> 0 DEGS
# -> return to default pre-filtering with threshold 10 (voom-transformation)


########
# step 4: change gene set database to KEGG
########

# -> 1 DEGS
# ->> proceed with geneset database KEGG


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 1 DEGS
# alternative 2: Difference of Classes -> 1 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 1 DEGS
# alternative 3: exponent 2 -> 1 DEGS

### final results: 1 DEGS
# achieved with alternative geneset database KEGG


#######################################################################################
### Phenotype Permutation 6 ###########################################################
#######################################################################################

i <- 6


#########
# step 1: default
#########

# -> 1 DEGS

#########
# step 2: vst-transformed gene expression data set
#########

# -> 8 DEGS
# -> proceed with RNA-Seq transformation method vst


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen6 <-  DGEList(bottomly.eset, group = phen_bottomly[,6]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen6 <- variancetransform(bottomly.eset[prefilt_ind_phen6,], phenotype_labels = phen_bottomly[,6])

### export
path_filterByExpr_phen6 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen6,
            file = path_filterByExpr_phen6,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

#-> 0 DEGS
# ->> return to default pre-filtering with threshold 10 (vst-transformed gene expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 1 DEGS
# ->> return to default gene set database GO (BP)



########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 8 DEGS
# alternative 2: Difference of Classes -> 2 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 8 DEGS
# alternative 3: exponent 2 -> 1 DEGS

# final result: 8 DEGS
# achieved with alternative RNA-Seq transformation method vst


#######################################################################################
### Phenotype Permutation 7 ###########################################################
#######################################################################################

i <- 7


#########
# step 1: default
#########

# -> 0 DEGS

#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# ->> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen7 <-  DGEList(bottomly.eset, group = phen_bottomly[,7]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen7 <-voom_trans(bottomly.eset[prefilt_ind_phen7,], phenotype_labels = phen_bottomly[,7])

### export
path_filterByExpr_phen7 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen7,
            file = path_filterByExpr_phen7,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 0 DEGS
# ->> return to default pre-filtering with threshold 10 (voom-transformation)


########
# step 4: change gene set database to KEGG
########

# -> 4 DEGS
# -> proceed with geneset database KEGG


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 4 DEGS
# alternatove 2: Difference of Classes -> 1 DEGS

#->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 5 DEGS
# alternative 2: exponent 1.5 -> 1 DEGS
# alternative 3: exponent 2 -> 5 DEGS

# -> end result: 5 DEGS
# achieved with
# ALTERNATIVE gene set database KEGG
# ALTERNATIVE exponents 0 and 2
# -> (tied; in codes R-Coded, in the case of a tie the parameter with the lower index, i.e. exponent 0, is chosen)



#######################################################################################
### Phenotype Permutation 8 ###########################################################
#######################################################################################

i <- 8


#########
# step 1: default
#########

# -> 10 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 33 DEGS
# -> proceed with alternative RNA-Seq transformation method vst



########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen8 <-  DGEList(bottomly.eset, group = phen_bottomly[,8]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen8 <- variancetransform(bottomly.eset[prefilt_ind_phen8,], phenotype_labels = phen_bottomly[,8])

### export
path_filterByExpr_phen8 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen8,
            file = path_filterByExpr_phen8,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 42 DEGS
# -> proceed with filtering using filterByExpr() (vst-transformed gene expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 42 DEGS
# alternative 2: Difference of Classes -> 25 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 9 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 4 DEGS

# -> end result: 42 DEGS
# achieved with
# ALTERNATIVE RNA-Seq transformation method vst
# ALTERNATIVE pre-filtering using filterByExpr()


#######################################################################################
### Phenotype Permutation 9 ###########################################################
#######################################################################################

i <- 9


#########
# step 1: default
#########

# -> 7 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 1 DEGS
# ->> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen9 <-  DGEList(bottomly.eset, group = phen_bottomly[,9]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen9 <-voom_trans(bottomly.eset[prefilt_ind_phen9,], phenotype_labels = phen_bottomly[,9])

### export
path_filterByExpr_phen9 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen9,
            file = path_filterByExpr_phen9,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


# -> 8 DEGS
# ->> proceed with pre-filtering using filterByExpr() (voom-transformed gene expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 2 DEGS
# ->> return to gene set database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 8 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 17 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 1 DEGS

# -> final results: 17 DEGS
# achieved with
# alternative pre-filtering using filterByExpr()
# alternative exponent 0



#######################################################################################
### Phenotype Permutation 9 ###########################################################
#######################################################################################

i <- 10


#########
# step 1: default
#########

# -> 0 DEGS


#########
# step 2: vst-transformed gene expression data set
#########

# -> 0 DEGS
# ->> return to default RNA-Seq transformation method voom


########
# step 3: perform pre-filtering using edgeR's builtin function filterByExpr
########

# generate pre-filtering indicator using filterByExpr()
prefilt_ind_phen10 <-  DGEList(bottomly.eset, group = phen_bottomly[,10]) %>% filterByExpr()
# perform voom-transformation on accordingly pre-filtered bottomly data set
exprdat_prefilt_phen10 <-voom_trans(bottomly.eset[prefilt_ind_phen10,], phenotype_labels = phen_bottomly[,10])

### export
path_filterByExpr_phen10 <- paste0("./Results/Intermediate_results/GSEA_Web/Bottomly/n_DEGS/Data/Raw/Phen_Permutation",i,"/exprdat_filterByExpr_phen_permutation", i, ".txt")

write.table(exprdat_prefilt_phen10,
            file = path_filterByExpr_phen10,
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# -> 0 DEGS
# ->> return to default pre-filtering with threshold 10 (voom-transformed gene expression data set)


########
# step 4: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default geneset database GO (BP)


########
# step 5: change gene-level ranking metric
########

# alternative 1: t-Test -> 0 DEGS
# alternative 2: Difference of Classes -> 0 DEGS


######
# step 6: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

### final results: 0 DEGS
# optimal parameter configuration coincides with default configuration
# number of differentielly enriched gene sets could not be increased




################################################################################
### (II) Optimization of internal parameters of GSEA ###########################
################################################################################

# for each phenotype assignment (i.e. original assignment and random permutations),
# the optimization process is based on the gene expression data set pre-processed
# in default manner
# -> voom-transformation of RNA-Seq data and manual pre-filtering with threshold 10


#######################################################################################
### Original Phenotype Assignment #####################################################
#######################################################################################


#########
# step 1: default
#########

# -> 2 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS
# -> return to default geneset database GO (BP)



########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 2 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> return to default gene-level ranking metric signal2noise ratio


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 52 + 3 = 55 DEGS
# alternative 3: exponent 2 -> 247 + 99 = 346 DEGS


### final results: 246 DEGS
# achieved with ALTERNATIVE exponent 2


#######################################################################################
### Phenotype Permutation 1 ###########################################################
#######################################################################################



#########
# step 1: default
#########

# -> 3 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS
# -> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 4 DEGS
# alternative 2: Difference of Classes -> 2 DEGS

# ->> proceed with gene-level ranking metric t-Test (t-statistic)


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

# -> return to default exponent 1

### final results: 4 DEGS
# achieved with ALTERNATIVE gene-level ranking metric t-Test


#######################################################################################
### Phenotype Permutation 2 ###########################################################
#######################################################################################



#########
# step 1: default
#########

# -> 1 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS
# ->> return to geneset database GO (BP)


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 1 DEGS
# alternative 2: Difference of Classes -> 3 DEGS

# ->> proceed with gene-level ranking metric Difference of Classes


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

## -> return to default exponent 1

### final results: 3 DEGS
# achieved with alternative gene-level ranking metric Difference of Classes


#######################################################################################
### Phenotype Permutation 3 ###########################################################
#######################################################################################



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
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

### final results: 0 DEGS
# optimal parameter choice coincided with default parameter choice
# number of differentially enriched gene sets could not be increased


#######################################################################################
### Phenotype Permutation 3 ###########################################################
#######################################################################################


#########
# step 1: default
#########

# -> 1 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS
# ->> return to default gene set database GO (BP)


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 4 DEGS
# alternative 2: Difference of Classes -> 6 DEGS


# ->> proceed with alernative gene-level ranking metric Difference of Classes


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 5 DEGS
# alternative 3: exponent 2 -> 0 DEGS

### final results: 6 DEGS
# achieved with alternative gene-level ranking metric Difference of Classes

#######################################################################################
### Phenotype Permutation 5 ###########################################################
#######################################################################################


#########
# step 1: default
#########

# -> 0 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS
# ->> proceed with geneset database KEGG


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 1 DEGS
# alternative 2: Difference of Classes -> 1 DEGS

# ->> proceed with default geneset database Signal2Noise ratio


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 1 DEGS
# alternative 3: exponent 2 -> 1 DEGS

### final results: 1 DEGS
# achieved with alternative geneset database KEGG


#######################################################################################
### Phenotype Permutation 6 ###########################################################
#######################################################################################


#########
# step 1: default
#########

# -> 1 DEGS

########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS
# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 1 DEGS
# alternative 2: Difference of Classes -> 0 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio

######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 1 DEGS


### end results: 1 DEGS
# optimal parameter configuration coincides with default parameter configuration
# number of differentially enriched gene sets could not be increased


#######################################################################################
### Phenotype Permutation 7 ###########################################################
#######################################################################################


#########
# step 1: default
#########

# -> 0 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 4 DEGS
# ->> proceed with gene set database KEGG


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 4 DEGS
# alternative 2: Difference of Classes -> 1 DEGS

# ->> return to default gene-level ranking metric Signal2Noise ratio


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 5 DEGS
# alternative 2: exponent 1.5 -> 1 DEGS
# alternative 3: exponent 2 -> 5 DEGS


# final results: 5 DEGS
# achieved with
# alternative gene set database KEGG
# alternative exponents 0 and 2 (tied)


#######################################################################################
### Phenotype Permutation 8 ###########################################################
#######################################################################################


#########
# step 1: default
#########

# -> 10 DEGS


########
# step 2: change gene set database to KEGG
########

# -> 0 DEGS

# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 10 DEGS
# alternative 2: Difference of Classes -> 2 DEGS

# ->> return to defaul gene-level ranking metric Signal2Noise ratio


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 6 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 1 DEGS

### -> return to default exponent 1

### final results: 10 DEGS
# optimal parameter configuration coincides with default parameter configuration
# number of differentially enriched gene sets cannot be increased



#######################################################################################
### Phenotype Permutation 9 ###########################################################
#######################################################################################


#########
# step 1: default
#########

# -> 7 DEGS



########
# step 2: change gene set database to KEGG
########

# -> 1 DEGS
# ->> return to default geneset database GO (BP)


########
# step 3: change gene-level ranking metric
########

# alternative 1: t-Test -> 9 DEGS
# alternative 2: Difference of Classes -> 0 DEGS


# ->> proceed with gene-level ranking metric t-Test (t-statistic)


######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 73 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 4 DEGS

### final results: 73 DEGS
# achieved with
# alternative gene-level ranking metric tTest
# alternative exponent 0



#######################################################################################
### Phenotype Permutation 9 ###########################################################
#######################################################################################


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

## -> return to default gene-level ranking metric Signal2Noise ratio



######
# step 4: change exponent
######

# alternative 1: exponent 0 -> 0 DEGS
# alternative 2: exponent 1.5 -> 0 DEGS
# alternative 3: exponent 2 -> 0 DEGS

### final results: 0 DEGS
# optimal parameter configuration coincides with default configuration
# number of differentially enriched gene sets could not be increased


























