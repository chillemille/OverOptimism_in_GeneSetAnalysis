################################################################################
### Run optimizations of the number of differentialla enriched gene sets #######
### for GOSeq, PADOG, clusterProfiler's GSEA and ORA ###########################
################################################################################

# create Results folder in working directory if is does not exist already 
# in this folder, all results are to be stored 

if(!dir.exists("./Results")){dir.create(path = "./Results")}


################################################################################
### GOSeq ######################################################################
################################################################################

# load optimisation functions 
source("./n_DEGS_optim_GOSeq.R")

#############
### Pickrell
#############

# true sample conditions
optim_GOSeq_results_Pickrell_originalphenotype <- GOSeq_joint_optimization(Biobase::exprs(pickrell.eset), 
                                                                           pickrell.eset$gender)

# save results
save(optim_GOSeq_results_Pickrell_originalphenotype, 
     file = "./Results/GOSeq_Results_Pickrell_OriginalPhenotype.RData")

# 10 permutations of sample conditions
optim_GOSeq_results_Pickrell_phenotypepermutation <- apply(FUN = GOSeq_joint_optimization, 
                                                           expression_data = Biobase::exprs(pickrell.eset), 
                                                           X = phen_pickrell, MARGIN = 2)

# save results
save(optim_GOSeq_results_Pickrell_phenotypepermutation, 
     file = "./Results/GOSeq_Results_Pickrell_PhenotypePermutations.RData")

#############
### Bottomly 
#############

# true sample conditions
optim_GOSeq_results_Bottomly_originalphenotype <- GOSeq_joint_optimization(Biobase::exprs(bottomly.eset), 
                                                                           bottomly.eset$strain)

# save results
save(optim_GOSeq_results_Bottomly_originalphenotype, 
     file = "./Results/GOSeq_Results_Bottomly_OriginalPhenotype.RData")


# 10 random permutations of sample conditions
optim_GOSeq_results_Bottomly_phenotypepermutation <- apply(FUN = GOSeq_joint_optimization, 
                                                           expression_data = Biobase::exprs(bottomly.eset), 
                                                           X = phen_bottomly, MARGIN = 2)


# save results
save(optim_GOSeq_results_Bottomly_phenotypepermutation, 
     file = "./Results/GOSeq_Results_Bottomly_PhenotypePermutations.RData")


################################################################################
### PADOG ######################################################################
################################################################################

# load optimisation functions 
source("./n_DEGS_optim_PADOG.R")

# create list of vectors from phen_pickrell 
# -> for use of function lapply() since apply() returns an error 
phen_pickrell_list <- list()

for(i in 1:ncol(phen_pickrell)){
  
  phen_pickrell_list[[i]] <- phen_pickrell[,i]
  
}


# create list of vectors from phen_bottomly 
# -> for use of function lapply() since apply() returns an error 
phen_bottomly_list <- list()

for(i in 1:ncol(phen_bottomly)){
  
  phen_bottomly_list[[i]] <- phen_bottomly[,i]
  
}

#############
### Pickrell 
#############

# original phenotype assignment 
optim_PADOG_results_Pickrell_originalphenotype <- PADOG_optim(Biobase::exprs(pickrell.eset),
                                                              pickrell.eset$gender)

# save results
save(optim_PADOG_results_Pickrell_originalphenotype,
     file = "./Results/PADOG_Results_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optim_PADOG_results_Pickrell_phenotypepermutation <- lapply(FUN = PADOG_optim, expression_data = Biobase::exprs(pickrell.eset), 
                                                            X = phen_pickrell_list)

# save results
save(optim_PADOG_results_Pickrell_phenotypepermutation,
     file = "./Results/PADOG_Results_Pickrell_PhenotypePermutations.RData")

#############
### Bottomly 
#############

# original phenotype assignment 
optim_PADOG_results_Bottomly_originalphenotype <- PADOG_optim(Biobase::exprs(bottomly.eset),
                                                              bottomly.eset$strain)

# save results
save(optim_PADOG_results_Bottomly_originalphenotype,
     file = "./Results/PADOG_Results_Bottomly_OriginalPhenotype.RData")


# 10 random permutations of the sample labels 
optim_PADOG_results_Bottomly_phenotypepermutation <- lapply(FUN = PADOG_optim, expression_data = Biobase::exprs(bottomly.eset), 
                                                            X = phen_bottomly_list)


# save results
save(optim_PADOG_results_Bottomly_phenotypepermutation,
     file = "./Results/PADOG_Results_Bottomly_PhenotypePermutations.RData")


################################################################################
### clusterProfiler's ORA ######################################################
################################################################################


# load optimisation functions 
source("./n_DEGS_optim_cP_ORA.R")


#############
### Pickrell 
#############

# true sample labels
optim_ORA_results_Pickrell_originalphenotype <- cP_ORA_joint_optimization(Biobase::exprs(pickrell.eset), 
                                                                          pickrell.eset$gender)
# save results
save(optim_ORA_results_Pickrell_originalphenotype, 
     file = "./Results/ORA_Results_Pickrell_OriginalPhenotype.RData")



# 10 permuted sample labels
optim_ORA_results_Pickrell_phenotypepermutation <- apply(FUN = cP_ORA_joint_optimization, 
                                                         expression_data = Biobase::exprs(pickrell.eset), 
                                                         X = phen_pickrell, MARGIN = 2)

# save results
save(optim_ORA_results_Pickrell_phenotypepermutation, 
     file = "./Results/ORA_Results_Pickrell_PhenotypePermutations.RData")

#############
### Bottomly 
#############

# true sample labels
optim_ORA_results_Bottomly_originalphenotype <- cP_ORA_joint_optimization(Biobase::exprs(bottomly.eset), 
                                                                          bottomly.eset$strain)

# save results
save(optim_ORA_results_Bottomly_originalphenotype, 
     file = "./Results/ORA_Results_Bottomly_OriginalPhenotype.RData")


# 10 permuted sample labels
optim_ORA_results_Bottomly_phenotypepermutation <- apply(FUN = cP_ORA_joint_optimization, 
                                                         expression_data = Biobase::exprs(bottomly.eset), 
                                                         X = phen_bottomly, 
                                                         MARGIN = 2)


# save results
save(optim_ORA_results_Bottomly_phenotypepermutation, 
     file = "./Results/ORA_Results_Bottomly_PhenotypePermutations.RData")


################################################################################
### clusterProfiler's GSEA #####################################################
#################################################################################

# load optimisation functions 
source("./n_DEGS_optim_cP_GSEA.R")



phen_pickrell_list <- list()

for(i in 1:ncol(phen_pickrell)){
  
  phen_pickrell_list[[i]] <- phen_pickrell[,i]
  
}

phen_bottomly_list <- list()

for(i in 1:ncol(phen_bottomly)){
  
  phen_bottomly_list[[i]] <- phen_bottomly[,i]
  
}


#############
### Pickrell 
#############

# original phenotype assignment 
optim_cP_GSEA_results_Pickrell_originalphenotype <- cP_GSEA_joint_optimization(Biobase::exprs(pickrell.eset), 
                                                                               pickrell.eset$gender)

# save results
save(optim_cP_GSEA_results_Pickrell_originalphenotype, 
     file = "./Results/cP_GSEA_Results_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optim_cP_GSEA_results_Pickrell_phenotypepermutation <- lapply(FUN = cP_GSEA_joint_optimization, 
                                                              expression_data = Biobase::exprs(pickrell.eset), 
                                                              X = phen_pickrell_list)

# save results
save(optim_cP_GSEA_results_Pickrell_phenotypepermutation, 
     file = "./Results/cP_GSEA_Results_Pickrell_PhenotypePermutations.RData")

#############
### Bottomly 
#############

# original phenotype assignment 
optim_cP_GSEA_results_Bottomly_originalphenotype <- cP_GSEA_joint_optimization(Biobase::exprs(bottomly.eset), 
                                                                               bottomly.eset$strain)

# save results
save(optim_cP_GSEA_results_Bottomly_originalphenotype, 
     file = "./Results/cP_GSEA_Results_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optim_cP_GSEA_results_Bottomly_phenotypepermutation <- lapply(FUN = cP_GSEA_joint_optimization, 
                                                              expression_data = Biobase::exprs(bottomly.eset), 
                                                              X = phen_bottomly_list)


# save results
save(optim_cP_GSEA_results_Bottomly_phenotypepermutation, 
     file = "./Results/cP_GSEA_Results_Bottomly_PhenotypePermutations.RData")































