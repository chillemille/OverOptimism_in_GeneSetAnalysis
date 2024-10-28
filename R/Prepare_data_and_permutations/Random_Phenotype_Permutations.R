####################################################################################################
### Upload RNA-Seq gene expression data sets and create random phenotype permutations ##############
####################################################################################################

# source pre-processing functions: especially create_phenpermutations() to generate
# 10 random permutations of the original phenotypes
source("./R/Help_functions/PreProcessing_Functions.R")

############################
### Bottomly et al. data set
############################

load("./GeneExpression_data/bottomly_eset.RData")
# note: we can extract the count data using the command
# Biobase::exprs(bottomly.eset)


# note: this RNA-Seq was downloaded from the following link (Dec 16 2022, 13:54):
# if (!file.exists("bottomly_eset.RData")) download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData",
#                                                       "bottomly_eset.RData")



# create 10 randpm permutations of the true sample conditions
phen_bottomly  <- create_phenpermutations(Biobase::exprs(bottomly.eset), bottomly.eset$strain, 10)
# check dimension
dim(phen_bottomly)
# save in working directory
# save(phen_bottomly, file = "./GeneExpression_Measurements/Save_Phenotype_Permutations_Bottomly.Rdata")

############################
### Pickrell et al. data set
############################

library(tweeDEseqCountData)
data(pickrell)

# create 10 random permutations of the true sample labels
phen_pickrell  <- create_phenpermutations(Biobase::exprs(pickrell.eset), pickrell.eset$gender, 10)

dim(phen_pickrell)
# save in working directory
# save(phen_pickrell, file = "./GeneExpression_Measurements/Phenotype_Permutations_Pickrell.Rdata")
