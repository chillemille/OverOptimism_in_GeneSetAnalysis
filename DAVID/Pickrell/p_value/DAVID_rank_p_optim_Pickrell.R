##############################################################################
### Maximization of the adjusted p-value obtained by DAVID for Pickrell data #
##############################################################################

library(DESeq2)
library(limma)
library(edgeR) # for filterByExpr()
library(org.Hs.eg.db)
library(dplyr)

# the link to run DAVID is the following: https://david.ncifcrf.gov/

# set working directory
setwd("./Assessment_OverOptimism")

# load gene expression data set with true phenotype randomly permuted phenotype assignments 
source("./Random_Phenotype_Permutations.R")

# load functions required for data preprocessing 
source("./PreProcessing_Functions.R")



################################################################################
### Original Phenotype Assignment ##############################################
################################################################################


##################
# Step 1: Default (22 differentially expressed genes)
##################

# -> adj. p-value = 1 



#################
# Step 2: Gene List obtained with limma (3 differentially expressed genes)
#################

# -> 



#####################################
# Step 4: Upload alternative universe
#####################################

# -> 0 DEGS 



## -> final results: 


################################################################################
### Random Phenotype Permutations ##############################################
################################################################################


################################################################################
### Phenotype Permutation 1 ####################################################
################################################################################


##################
# Step 1: Default (3 differentially expressed genes)
##################

# -> 



#################
# Step 2: Gene List obtained with limma (0 differentially expressed genes)
#################



#########################################
# Step 3: Change geneset database to KEGG
#########################################

# -> 0 DEGS 



#####################################
# Step 4: Upload alternative universe
#####################################

# -> 0 DEGS 



## -> final results: 


################################################################################
### Phenotype Permutation 2 ####################################################
################################################################################


##################
# Step 1: Default (180 differentially expressed genes)
##################

# -> 


#################
# Step 2: Gene List obtained with limma (107 differentially expressed genes)
#################

# -> 
# ->> 


#########################################
# Step 3: Change geneset database to KEGG
#########################################

# -> 



#####################################
# Step 4: Upload alternative universe
#####################################

# -> 0 DEGS 



## -> final results: 

################################################################################
### Phenotype Permutation 3 ####################################################
################################################################################

# None of the lists of differentially expressed genes, whether generated with 
# DESeq2 nor limma, contained any genes 

# -> DAVID would therefore always return 0 differentially enriched gene sets 



################################################################################
### Phenotype Permutation 4 ####################################################
################################################################################


##################
# Step 1: Default 
##################

# -> 0 DEGS 


#################
# Step 2: Gene List obtained with limma 
#################

# the list of differentially expressed genes contained 0 differentially expressed 
# gene sets when being generated with limma 

# -> skip step 2 and proceed directly to step 3 


#########################################
# Step 3: Change geneset database to KEGG
#########################################

# -> 



#####################################
# Step 4: Upload alternative universe
#####################################

# -> 


################################################################################
### Phenotype Permutation 5 ####################################################
################################################################################



##################
# Step 1: Default 
##################

# -> 


#################
# Step 2: Gene List obtained with limma 
#################

# the list of differentially expressed genes contained 0 differentially expressed 
# gene sets when being generated with limma 

# -> skip step 2 and proceed directly to step 3 


#########################################
# Step 3: Change geneset database to KEGG
#########################################

# -> 



#####################################
# Step 4: Upload alternative universe
#####################################

# -> 0 DEGS 


################################################################################
### Phenotype Permutation 6 ####################################################
################################################################################

# None of the lists of differentially expressed genes, whether generated with 
# DESeq2 nor limma, contained any genes 

# -> DAVID would therefore always return 0 differentially enriched gene sets 



################################################################################
### Phenotype Permutation 7 ####################################################
################################################################################


# None of the lists of differentially expressed genes, whether generated with 
# DESeq2 nor limma, contained any genes 

# -> DAVID would therefore always return 0 differentially enriched gene sets 


################################################################################
### Phenotype Permutation 8 ####################################################
################################################################################


# None of the lists of differentially expressed genes, whether generated with 
# DESeq2 nor limma, contained any genes 

# -> DAVID would therefore always return 0 differentially enriched gene sets 


################################################################################
### Phenotype Permutation 9 ####################################################
################################################################################



##################
# Step 1: Default 
##################

# ->


#################
# Step 2: Gene List obtained with limma 
#################

# the list of differentially expressed genes contained 0 differentially expressed 
# gene sets when being generated with limma 

# -> skip step 2 and proceed directly to step 3 


#########################################
# Step 3: Change geneset database to KEGG
#########################################

# -> 


#####################################
# Step 4: Upload alternative universe
#####################################

# -> 0 DEGS 

################################################################################
### Phenotype Permutation 10 ###################################################
################################################################################


# None of the lists of differentially expressed genes, whether generated with 
# DESeq2 nor limma, contained any genes 

# -> DAVID would therefore always return 0 differentially enriched gene sets 









