####################################################################################################
### Upload of RNA-seq gene expression data sets and creation of random phenotype permutations ######
####################################################################################################


############################################################################
### Function to create random permutations of the true sample conditions ###
############################################################################

# arguments: 
# expression_data: gene expression measurements 
# true_sample_conditions: conditions of the samples from the gene expression data 
# nperm: number of permutations to be created 

create_phenpermutations <- function(expression_data, true_sample_conditions, nperm){
  
  # get distribution of true sample conditions
  table(true_sample_conditions) 
  
  # get levels of true sample conditions
  levels(true_sample_conditions)
  
  # create data frame to contain the permuted sample conditions 
  # note that the permuted sample conditions are later filled in columnwisely 
  phen_permutations <- data.frame(to_be_replaced = rep(NA, times = length( true_sample_conditions)))
  
  ### create nperm random phenotype permutations
  for(i in 1:nperm){
    
    # set seed for reproducibility
    set.seed(i)
    
    # randomly generate the permuted positions of the first phenotype 
    sample_index <- sort(sample(x = 1:length(true_sample_conditions), # number of samples from entire gene expression data set 
                                size = table(true_sample_conditions)[1], # sample size of first phenotype 
                                replace = FALSE), # without replacement
                                decreasing = FALSE) # sort indices in ascending manner
    
    # create vector to contain one permutation of the true sample labels 
    phen_vec <- rep(NA, times = length(true_sample_conditions))
    # fill in the first condition with positions according to sample_index 
    phen_vec[sample_index] <- levels(true_sample_conditions)[1]
    # fill in the second condition to the remaining positions 
    phen_vec[is.na(phen_vec)] <- levels(true_sample_conditions)[2]
    
    # check whether frequencies of permuted sample labels coincides with frequency of true sample labels
    if(!all(table(phen_vec) == table(true_sample_conditions))) stop("Error: Frequencies of permuted sample conditions does not 
                                                                    match frequency of true sample conditions")
    
    # add newly created random phenotype permutation
    phen_permutations[, i] <- as.factor(phen_vec)
    
  }
  
  # add identification to all permutations
  colnames(phen_permutations) <- c(paste0("permutation", 1:nperm))
  
  # return data frame of the permuted sample conditions 
  return(phen_permutations)
  
  
  
}

####################################################################################################
### Upload RNA-Seq gene expression data sets and create random phenotype permutations ##############
####################################################################################################

############################
### Bottomly et al. data set
############################

load("./GeneExpression_Measurements/bottomly_eset.RData")
# note: we can extract the count data using the command 
# Biobase::exprs(bottomly.eset)


# note: this RNA-Seq was downloaded from the following link (Dec 16 2022, 13:54): 
# if (!file.exists("bottomly_eset.RData")) download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData", 
#                                                       "bottomly_eset.RData")



# create 10 randpm permutations of the true sample conditions 
phen_bottomly <- create_phenpermutations(Biobase::exprs(bottomly.eset), bottomly.eset$strain, 10)
# check dimension 
dim(phen_bottomly)
# save in working directory
save(phen_bottomly, file = "./GeneExpression_Measurements/Save_Phenotype_Permutations_Bottomly.Rdata")

############################
### Pickrell et al. data set 
############################

library(tweeDEseqCountData)
data(pickrell)

# create 10 random permutations of the true sample labels 
phen_pickrell <- create_phenpermutations(Biobase::exprs(pickrell.eset), pickrell.eset$gender, 10)

dim(phen_pickrell)
# save in working directory 
save(phen_pickrell, file = "./GeneExpression_Measurements/Phenotype_Permutations_Pickrell.Rdata")
