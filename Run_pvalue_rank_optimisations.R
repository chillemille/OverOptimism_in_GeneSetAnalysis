################################################################################
### Run optimizations of the adjusted p-value and rank of two gene sets ########
### for GOSeq, PADOG, clusterProfiler's GSEA and ORA ###########################
################################################################################

# create Results folder in working directory if is does not exist already 
# in this folder, all results are to be stored 

if(!dir.exists("./Results")){dir.create(path = "./Results")}

################################################################################
### GOSeq ######################################################################
################################################################################

# load optimisation functions 
source("./rank_p_OptimisationFunctions_GOSeq.R")


phen_pickrell_list <- list()

for(i in 1:ncol(phen_pickrell)){
  
  phen_pickrell_list[[i]] <- phen_pickrell[,i]
  
}

phen_bottomly_list <- list()

for(i in 1:ncol(phen_bottomly)){
  
  phen_bottomly_list[[i]] <- phen_bottomly[,i]
  
}

### optimization of adjusted p-values 

#############
### Pickrell 
#############

### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimP_GOSeq_tcell_Pickrell_originalphenotype <-  goseq_rank_pvalue_optim(geneset = "GO:0002456", 
                                                                          "GO:BP",
                                                                          metric = "p_adj", 
                                                                          Biobase::exprs(pickrell.eset), 
                                                                          pickrell.eset$gender)

# save results
save(optimP_GOSeq_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimP_GOSeq_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_GOSeq_tcell_Pickrell_phenotypepermutations <-  lapply(FUN = goseq_rank_pvalue_optim, 
                                                             geneset =  "GO:0002456", 
                                                             geneset_database = "GO:BP",
                                                             metric = "p_adj",
                                                             expression_data = Biobase::exprs(pickrell.eset), 
                                                             X = phen_pickrell_list)

# save results
save(optimP_GOSeq_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_GOSeq_tCell_Pickrell_PhenotypePermutations.RData")


### (II) Gene Set "Demethylation" -> "GO:0070988"

# original phenotype assignment 
optimP_GOSeq_Demethylation_Pickrell_originalphenotype <-  goseq_rank_pvalue_optim(geneset= "GO:0070988", "GO:BP",metric = "p_adj", Biobase::exprs(pickrell.eset), pickrell.eset$gender)

# save results
save(optimP_GOSeq_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimP_GOSeq_Demethylation_Pickrell_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations <-  lapply(FUN = goseq_rank_pvalue_optim, 
                                                                     geneset =  "GO:0070988", 
                                                                     metric = "p_adj",
                                                                     expression_data = Biobase::exprs(pickrell.eset), 
                                                                     geneset_database = "GO:BP",
                                                                     X = phen_pickrell_list)

save(optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_GOSeq_Demethylation_Pickrell_PhenotypePermutations.RData")


#############
### Bottomly 
#############

### (I) Gene set "Metabolic Process: GO:0008152"

# original phenotype assignment 
optimP_GOSeq_MetabolicProcess_Bottomly_originalphenotype <- 
  goseq_rank_pvalue_optim(geneset = "GO:0008152", 
                          "GO:BP","p_adj", 
                          Biobase::exprs(bottomly.eset), 
                          bottomly.eset$strain)

# save results
save(optimP_GOSeq_MetabolicProcess_Bottomly_originalphenotype, 
     file = "./Results/optimP_GOSeq_MetabolicProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                       geneset = "GO:0008152", 
                                                                       geneset_database = "GO:BP",
                                                                       metric = "p_adj",
                                                                       expression_data = Biobase::exprs(bottomly.eset), 
                                                                       X = phen_bottomly_list)
# save results 
save(optimP_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations, 
     file = "./Results/optimP_GOSeq_MetabolicProcess_Bottomly_PhenotypePermutations.RData")




### (II) Gene set "Cellular Process: GO:0009987"

# original phenotype assignment 
optimP_GOSeq_CellularProcess_Bottomly_originalphenotype <- 
  goseq_rank_pvalue_optim(geneset = 
                            "GO:0009987", 
                          "GO:BP",
                          metric = "p_adj", 
                          Biobase::exprs(bottomly.eset), 
                          bottomly.eset$strain)

# save results
save(optimP_GOSeq_CellularProcess_Bottomly_originalphenotype, 
     file = "./Results/optimP_GOSeq_CellularProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_GOSeq_CellularProcess_Bottomly_Phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                      geneset = "GO:0009987", 
                                                                      geneset_database = "GO:BP",
                                                                      metric = "p_adj",
                                                                      expression_data = Biobase::exprs(bottomly.eset), 
                                                                      X = phen_bottomly_list)
# save results 
save(optimP_GOSeq_CellularProcess_Bottomly_Phenotypepermutations, 
     file = "./Results/optimP_GOSeq_CellularProcess_Bottomly_PhenotypePermutations.RData")



# optimization of the ranks in the GSA results


#############
### Pickrell 
#############

### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimRank_GOSeq_tcell_Pickrell_originalphenotype <-  goseq_rank_pvalue_optim(geneset= "GO:0002456", 
                                                                             "GO:BP",
                                                                             metric = "rank", 
                                                                             Biobase::exprs(pickrell.eset), 
                                                                             pickrell.eset$gender)

# save results
save(optimRank_GOSeq_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimRank_GOSeq_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_GOSeq_tcell_Pickrell_phenotypepermutations <-  lapply(FUN = goseq_rank_pvalue_optim, 
                                                                geneset =  "GO:0002456", 
                                                                geneset_database = "GO:BP",
                                                                metric = "rank",
                                                                expression_data = Biobase::exprs(pickrell.eset), 
                                                                X = phen_pickrell_list)

# save results
save(optimRank_GOSeq_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_tCell_Pickrell_PhenotypePermutations.RData")


### (II) Gene Set "Demethylation" -> "GO:0070988"

# original phenotype assignment 
optimRank_GOSeq_Demethylation_Pickrell_originalphenotype <-  
  goseq_rank_pvalue_optim(geneset = "GO:0070988", 
                          "GO:BP",
                          metric = "rank", 
                          Biobase::exprs(pickrell.eset), 
                          pickrell.eset$gender)

# save results
save(optimRank_GOSeq_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimRank_GOSeq_Demethylation_Pickrell_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_GOSeq_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                       geneset =  "GO:0070988", 
                                                                       metric = "rank",
                                                                       expression_data = Biobase::exprs(pickrell.eset), 
                                                                       geneset_database = "GO:BP",
                                                                       X = phen_pickrell_list)

save(optimRank_GOSeq_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_Demethylation_Pickrell_PhenotypePermutations.RData")


#############
### Bottomly 
#############

### (I) Gene set "Metabolic Process: GO:0008152"

# original phenotype assignment 
optimRank_GOSeq_MetabolicProcess_Bottomly_originalphenotype <- goseq_rank_pvalue_optim(geneset = "GO:0008152", 
                                                                                       "GO:BP",
                                                                                       "rank", 
                                                                                       Biobase::exprs(bottomly.eset), 
                                                                                       bottomly.eset$strain)

# save results
save(optimRank_GOSeq_MetabolicProcess_Bottomly_originalphenotype, 
     file = "./Results/optimRank_GOSeq_MetabolicProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                          geneset = "GO:0008152", 
                                                                          geneset_database = "GO:BP",
                                                                          metric = "rank",
                                                                          expression_data = Biobase::exprs(bottomly.eset), 
                                                                          X = phen_bottomly_list)
# save results 
save(optimRank_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_MetabolicProcess_Bottomly_PhenotypePermutations.RData")




### (II) Gene set "Cellular Process: GO:0009987"

# original phenotype assignment 
optimRank_GOSeq_CellularProcess_Bottomly_originalphenotype <- goseq_rank_pvalue_optim(geneset = "GO:0009987", 
                                                                                      "GO:BP",
                                                                                      metric = "rank", 
                                                                                      Biobase::exprs(bottomly.eset), 
                                                                                      bottomly.eset$strain)

# save results
save(optimRank_GOSeq_CellularProcess_Bottomly_originalphenotype, 
     file = "./Results/optimRank_GOSeq_CellularProcess_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_GOSeq_CellularProcess_Bottomly_Phenotypepermutations <- lapply(FUN = goseq_rank_pvalue_optim, 
                                                                         geneset = "GO:0009987", 
                                                                         geneset_database = "GO:BP",
                                                                         metric = "rank",
                                                                         expression_data = Biobase::exprs(bottomly.eset), 
                                                                         X = phen_bottomly_list)
# save results 
save(optimRank_GOSeq_CellularProcess_Bottomly_Phenotypepermutations, 
     file = "./Results/optimRank_GOSeq_CellularProcess_Bottomly_PhenotypePermutations.RData")







################################################################################
### PADOG ######################################################################
################################################################################

# load optimisation functions 
source("./rank_p_OptimisationFunctions_PADOG.R")

phen_pickrell_list <- list()

for(i in 1:ncol(phen_pickrell)){
  
  phen_pickrell_list[[i]] <- phen_pickrell[,i]
  
}

phen_bottomly_list <- list()

for(i in 1:ncol(phen_bottomly)){
  
  phen_bottomly_list[[i]] <- phen_bottomly[,i]
  
}

# optimize the adjusted p-values 


#############
### Pickrell 
#############

### (I) Gene Set Primary Immunodeficiency -> "05340"

# original phenotype assignment 
optimP_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype <- 
  PADOG_rankp_optim(geneset = "05340",
                    Biobase::exprs(pickrell.eset), 
                    pickrell.eset$gender, 
                    metric = "p_adj")

# save results
save(optimP_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype, 
     file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations <- 
  lapply(FUN = PADOG_rankp_optim, 
         geneset =  "05340", 
         expression_data = Biobase::exprs(pickrell.eset), 
         metric = "p_adj",
         X = phen_pickrell_list)

# save results
save(optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Pickrell_PhenotypePermutations.RData")


### (II) Gene Set "Graft vs host disease" -> "05332"

# original phenotype assignment 
optimP_PADOG_GraftvsHost_Pickrell_originalphenotype <- PADOG_rankp_optim(geneset = "05332",
                                                                         Biobase::exprs(pickrell.eset), 
                                                                         pickrell.eset$gender, metric = "p_adj")

# save results
save(optimP_PADOG_GraftvsHost_Pickrell_originalphenotype, 
     file = "./Results/optimP_PADOG_GraftvsHost_Pickrell_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations <- lapply(FUN = PADOG_rankp_optim, 
                                                                  geneset =  "05332", 
                                                                  expression_data = Biobase::exprs(pickrell.eset), 
                                                                  metric = "p_adj",
                                                                  X = phen_pickrell_list)

save(optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_PADOG_GraftvsHost_Pickrell_PhenotypePermutations.RData")




#############
### Bottomly 
#############

# We use the same gene sets as for the pickrell data set as the gene sets provided 
# for human and mouse are almost identical 

### (I) Gene Set Primary Immunodeficiency -> " 	05340"

# original phenotype assignment 
optimP_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype <- 
  PADOG_rankp_optim(geneset = "05340",
                    Biobase::exprs(bottomly.eset), 
                    bottomly.eset$strain, 
                    metric = "p_adj")

# save results
save(optimP_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype, 
     file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations <- 
  lapply(FUN = PADOG_rankp_optim, 
         geneset =  "05340", 
         expression_data = Biobase::exprs(bottomly.eset), 
         metric = "p_adj",
         X = phen_bottomly_list)

# save results
save(optimP_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations, 
     file = "./Results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_PhenotypePermutations.RData")


### (II) Gene Set "Graft vs host disease" -> "05332"

# original phenotype assignment 
optimP_PADOG_GraftvsHost_Bottomly_originalphenotype <- 
  PADOG_rankp_optim(geneset = "05332",
                    Biobase::exprs(bottomly.eset), 
                    bottomly.eset$strain, 
                    metric = "p_adj")

# save results
save(optimP_PADOG_GraftvsHost_Bottomly_originalphenotype, 
     file = "./Results/optimP_PADOG_GraftvsHost_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_PADOG_GraftvsHost_Bottomly_phenotypepermutations <- 
  lapply(FUN = PADOG_rankp_optim, 
         geneset =  "05332", 
         expression_data = Biobase::exprs(bottomly.eset), 
         metric = "p_adj",
         X = phen_bottomly_list)

save(optimP_PADOG_GraftvsHost_Bottomly_phenotypepermutations, 
     file = "./Results/optimP_PADOG_GraftvsHost_Bottomly_PhenotypePermutations.RData")



# optimize the ranks 

#############
### Pickrell 
#############

### (I) Gene Set Primary Immunodeficiency -> "05340"

# original phenotype assignment 
optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype <- 
  PADOG_rankp_optim(geneset = "05340",
                    Biobase::exprs(pickrell.eset), 
                    pickrell.eset$gender, 
                    metric = "rank")

# save results
save(optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype, 
     file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations <- 
  lapply(FUN = PADOG_rankp_optim, 
         geneset =  "05340", 
         expression_data = Biobase::exprs(pickrell.eset), 
         metric = "rank",
         X = phen_pickrell_list)

# save results
save(optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_PhenotypePermutations.RData")


### (II) Gene Set "Graft vs host disease" -> "05332"

# original phenotype assignment 
optimRank_PADOG_GraftvsHost_Pickrell_originalphenotype <- PADOG_rankp_optim(geneset = "05332",
                                                                            Biobase::exprs(pickrell.eset), 
                                                                            pickrell.eset$gender, metric = "rank")

# save results
save(optimRank_PADOG_GraftvsHost_Pickrell_originalphenotype, 
     file = "./Results/optimRank_PADOG_GraftvsHost_Pickrell_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_PADOG_GraftvsHost_Pickrell_phenotypepermutations <- lapply(FUN = PADOG_rankp_optim, 
                                                                     geneset =  "05332", 
                                                                     expression_data = Biobase::exprs(pickrell.eset), 
                                                                     metric = "rank",
                                                                     X = phen_pickrell_list)

save(optimRank_PADOG_GraftvsHost_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_PADOG_GraftvsHost_Pickrell_PhenotypePermutations.RData")




#############
### Bottomly 
#############

# We use the same gene sets as for the pickrell data set as the gene sets provided 
# for human and mouse are almost identical 

### (I) Gene Set Primary Immunodeficiency -> "05340"

# original phenotype assignment 
optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype <- 
  PADOG_rankp_optim(geneset = "05340",
                    Biobase::exprs(bottomly.eset), 
                    bottomly.eset$strain, 
                    metric = "rank")

# save results
save(optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype, 
     file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations <- 
  lapply(FUN = PADOG_rankp_optim, 
         geneset =  "05340", 
         expression_data = Biobase::exprs(bottomly.eset), 
         metric = "rank",
         X = phen_bottomly_list)

# save results
save(optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations, 
     file = "./Results/optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_PhenotypePermutations.RData")


### (II) Gene Set "Graft vs host disease" -> "05332"

# original phenotype assignment 
optimRank_PADOG_GraftvsHost_Bottomly_originalphenotype <- 
  PADOG_rankp_optim(geneset = "05332",
                    Biobase::exprs(bottomly.eset), 
                    bottomly.eset$strain, 
                    metric = "rank")

# save results
save(optimRank_PADOG_GraftvsHost_Bottomly_originalphenotype, 
     file = "./Results/optimRank_PADOG_GraftvsHost_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_PADOG_GraftvsHost_Bottomly_phenotypepermutations <- 
  lapply(FUN = PADOG_rankp_optim, 
         geneset =  "05332", 
         expression_data = Biobase::exprs(bottomly.eset), 
         metric = "rank",
         X = phen_bottomly_list)

save(optimRank_PADOG_GraftvsHost_Bottomly_phenotypepermutations, 
     file = "./Results/optimRank_PADOG_GraftvsHost_Bottomly_PhenotypePermutations.RData")









################################################################################
### clusterProfiler's ORA ######################################################
################################################################################


# load optimisation functions 
source("./rank_p_OptimisationFunctions_cP_ORA.R")



phen_pickrell_list <- list()

for(i in 1:ncol(phen_pickrell)){
  
  phen_pickrell_list[[i]] <- phen_pickrell[,i]
  
}

phen_bottomly_list <- list()

for(i in 1:ncol(phen_bottomly)){
  
  phen_bottomly_list[[i]] <- phen_bottomly[,i]
  
}


###### optimization (i.e. minimization) of adjusted p-values 


#############
### Pickrell 
#############



############################################################
### (I) Gene Set "t-cell mediated immunity" -> GO:0002456
############################################################

# original phenotype assignment 
optimP_cP_ORA_tcell_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0002456",
                                                                        metric = "p_adj", 
                                                                        Biobase::exprs(pickrell.eset), 
                                                                        pickrell.eset$gender, 
                                                                        "GO")

# save results
save(optimP_cP_ORA_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimP_ORA_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_cP_ORA_tcell_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                             geneset =  "GO:0002456", 
                                                             metric = "p_adj",
                                                             expression_data = Biobase::exprs(pickrell.eset), 
                                                             geneset_database = "GO",
                                                             X = phen_pickrell_list)

# save results
save(optimP_cP_ORA_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_ORA_tCell_Pickrell_PhenotypePermutations.RData")


############################################################
### (II) Gene Set "Demethylation" -> GO:0070988
############################################################


# original phenotype assignment 
optimP_cP_ORA_Demethylation_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0070988",
                                                                                metric = "p_adj", 
                                                                                Biobase::exprs(pickrell.eset), 
                                                                                pickrell.eset$gender, 
                                                                                "GO")

# save results
save(optimP_cP_ORA_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimP_ORA_Demethylation_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                                     geneset =  "GO:0070988", 
                                                                     metric = "p_adj",
                                                                     expression_data = Biobase::exprs(pickrell.eset), 
                                                                     geneset_database = "GO",
                                                                     X = phen_pickrell_list)

# save results
save(optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_ORA_Demethylation_Pickrell_PhenotypePermutations.RData")



#############
### Bottomly 
#############

### (I) Gene set "t Cell mediated immunity: GO:0002456"

# original phenotype assignment 
optimP_cP_ORA_tCell_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0002456",
                                                                        "p_adj", 
                                                                        Biobase::exprs(bottomly.eset), 
                                                                        bottomly.eset$strain, 
                                                                        "GO")

# save results
save(optimP_cP_ORA_tCell_Bottomly_originalphenotype, 
     file = "./Results/optimP_ORA_tCell_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_cP_ORA_tCell_Bottomly_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                                        geneset = "GO:0002456", 
                                                                        metric = "p_adj",
                                                                        expression_data = Biobase::exprs(bottomly.eset), 
                                                                        geneset_database = "GO",
                                                                        X = phen_bottomly_list)

# save results 
save(optimP_cP_ORA_tCell_Bottomly_phenotypepermutations, 
     file = "./Results/optimP_ORA_tCell_Bottomly_PhenotypePermutations.RData")




### (II) Gene set "Demethylation: "GO:0070988"

# original phenotype assignment 
optimP_cP_ORA_Demethylation_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0070988",
                                                                                metric = "p_adj", 
                                                                                Biobase::exprs(bottomly.eset), 
                                                                                bottomly.eset$strain, 
                                                                                "GO")

# save results
save(optimP_cP_ORA_Demethylation_Bottomly_originalphenotype, 
     file = "./Results/optimP_ORA_Demethylation_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_cP_ORA_Demethylation_Bottomly_Phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                                       geneset = "GO:0070988", 
                                                                       metric = "p_adj",
                                                                       expression_data = Biobase::exprs(bottomly.eset), 
                                                                       geneset_database = "GO",
                                                                       X = phen_bottomly_list)
# save results 
save(optimP_cP_ORA_Demethylation_Bottomly_Phenotypepermutations, 
     file = "./Results/optimP_ORA_Demethylation_Bottomly_PhenotypePermutations.RData")


# optimization (i.e. minimization of the ranks)



#############
### Pickrell 
#############



############################################################
### (I) Gene Set "t-cell mediated immunity" -> GO:0002456
############################################################

# original phenotype assignment 
optimRank_cP_ORA_tcell_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0002456",
                                                                           metric = "rank", 
                                                                           Biobase::exprs(pickrell.eset), 
                                                                           pickrell.eset$gender, 
                                                                           "GO")

# save results
save(optimRank_cP_ORA_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimRank_ORA_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_cP_ORA_tcell_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                                geneset =  "GO:0002456", 
                                                                metric = "rank",
                                                                expression_data = Biobase::exprs(pickrell.eset), 
                                                                geneset_database = "GO",
                                                                X = phen_pickrell_list)

# save results
save(optimRank_cP_ORA_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_ORA_tCell_Pickrell_PhenotypePermutations.RData")


############################################################
### (II) Gene Set "Demethylation" -> GO:0070988
############################################################

# Note: in other GSA tools, demethylation might also be identified by GO:0080111

# original phenotype assignment 
optimRank_cP_ORA_Demethylation_Pickrell_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0070988",
                                                                                   metric = "rank", 
                                                                                   Biobase::exprs(pickrell.eset), 
                                                                                   pickrell.eset$gender, 
                                                                                   "GO")

# save results
save(optimRank_cP_ORA_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimRank_ORA_Demethylation_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_cP_ORA_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                                        geneset =  "GO:0070988", 
                                                                        metric = "rank",
                                                                        expression_data = Biobase::exprs(pickrell.eset), 
                                                                        geneset_database = "GO",
                                                                        X = phen_pickrell_list)

# save results
save(optimRank_cP_ORA_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_ORA_Demethylation_Pickrell_PhenotypePermutations.RData")



#############
### Bottomly 
#############

### (I) Gene set "t Cell mediated immunity: "GO:0002456"

# original phenotype assignment 
optimRank_cP_ORA_tCell_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0002456",
                                                                           "rank", 
                                                                           Biobase::exprs(bottomly.eset), 
                                                                           bottomly.eset$strain, 
                                                                           "GO")

# save results
save(optimRank_cP_ORA_tCell_Bottomly_originalphenotype, 
     file = "./Results/optimRank_ORA_tCell_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_cP_ORA_tCell_Bottomly_phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                                           geneset = "GO:0002456", 
                                                                           metric = "rank",
                                                                           expression_data = Biobase::exprs(bottomly.eset), 
                                                                           geneset_database = "GO",
                                                                           X = phen_bottomly_list)

# save results 
save(optimRank_cP_ORA_tCell_Bottomly_phenotypepermutations, 
     file = "./Results/optimRank_ORA_tCell_Bottomly_PhenotypePermutations.RData")




### (II) Gene set "Demethylation: "GO:0070988"

# original phenotype assignment 
optimRank_cP_ORA_Demethylation_Bottomly_originalphenotype <- ORA_rank_pvalue_optim(geneset = "GO:0070988",
                                                                                   metric = "rank", 
                                                                                   Biobase::exprs(bottomly.eset), 
                                                                                   bottomly.eset$strain, 
                                                                                   "GO")

# save results
save(optimRank_cP_ORA_Demethylation_Bottomly_originalphenotype, 
     file = "./Results/optimRank_ORA_Demethylation_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_cP_ORA_Demethylation_Bottomly_Phenotypepermutations <- lapply(FUN = ORA_rank_pvalue_optim, 
                                                                          geneset = "GO:0070988", 
                                                                          metric = "rank",
                                                                          expression_data = Biobase::exprs(bottomly.eset), 
                                                                          geneset_database = "GO",
                                                                          X = phen_bottomly_list)
# save results 
save(optimRank_cP_ORA_Demethylation_Bottomly_Phenotypepermutations, 
     file = "./Results/optimRank_ORA_Demethylation_Bottomly_PhenotypePermutations.RData")







################################################################################
### clusterProfiler's GSEA #####################################################
#################################################################################

# load optimisation functions 
source("./rank_p_OptimisationFunctions_cP_GSEA.R")


phen_pickrell_list <- list()

for(i in 1:ncol(phen_pickrell)){
  
  phen_pickrell_list[[i]] <- phen_pickrell[,i]
  
}

phen_bottomly_list <- list()

for(i in 1:ncol(phen_bottomly)){
  
  phen_bottomly_list[[i]] <- phen_bottomly[,i]
  
}


##############################################################################################
### Optimization of the adjusted p-values ####################################################
##############################################################################################


#############
### Pickrell 
#############

### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimP_cP_GSEA_tcell_Pickrell_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0002456",
                                                                       "GO", 
                                                                       Biobase::exprs(pickrell.eset), 
                                                                       pickrell.eset$gender, 
                                                                       metric = "p_adj")

# save results
save(optimP_cP_GSEA_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimP_cP_GSEA_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                              geneset = "GO:0002456", 
                                                              geneset_database = "GO",
                                                              expression_data = Biobase::exprs(pickrell.eset), 
                                                              metric = "p_adj",
                                                              X = phen_pickrell_list)
# save results
save(optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_cP_GSEA_tCell_Pickrell_PhenotypePermutations.RData")



### (II) Gene Set "Demethylation" -> GO:0070988

# original phenotype assignment 
optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0070988",
                                                                               "GO", 
                                                                               Biobase::exprs(pickrell.eset), 
                                                                               pickrell.eset$gender, 
                                                                               metric = "p_adj")

# save results
save(optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimP_cP_GSEA_Demethylation_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                                      geneset = "GO:0070988", 
                                                                      geneset_database = "GO",
                                                                      expression_data = Biobase::exprs(pickrell.eset), 
                                                                      metric = "p_adj",
                                                                      X = phen_pickrell_list)
# save results
save(optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimP_cP_GSEA_Demethylation_Pickrell_PhenotypePermutations.RData")






#############
### Bottomly 
#############


### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimP_cP_GSEA_tCell_Bottomly_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0002456",
                                                                       "GO", 
                                                                       Biobase::exprs(bottomly.eset), 
                                                                       bottomly.eset$strain, 
                                                                       metric = "p_adj")

# save results
save(optimP_cP_GSEA_tCell_Bottomly_originalphenotype, 
     file = "./Results/optimP_cP_GSEA_tCell_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_cP_GSEA_tCell_Bottomly_phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                              geneset = "GO:0002456", 
                                                              geneset_database = "GO",
                                                              expression_data = Biobase::exprs(bottomly.eset), 
                                                              metric = "p_adj",
                                                              X = phen_bottomly_list)
# save results 
save(optimP_cP_GSEA_tCell_Bottomly_phenotypepermutations, 
     file = "./Results/optimP_cP_GSEA_tCell_Bottomly_PhenotypePermutations.RData")



### (II) Gene Set "Demethylation" -> GO:0070988

# original phenotype assignment 
optimP_cP_GSEA_Demethylation_Bottomly_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0070988",
                                                                               "GO", 
                                                                               Biobase::exprs(bottomly.eset), 
                                                                               bottomly.eset$strain, 
                                                                               metric = "p_adj")

# save results
save(optimP_cP_GSEA_Demethylation_Bottomly_originalphenotype, 
     file = "./Results/optimP_cP_GSEA_Demethylation_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimP_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                                      geneset = "GO:0070988", 
                                                                      geneset_database = "GO",
                                                                      expression_data = Biobase::exprs(bottomly.eset), 
                                                                      metric = "p_adj",
                                                                      X = phen_bottomly_list)
# save results 
save(optimP_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations, 
     file = "./Results/optimP_cP_GSEA_Demethylation_Bottomly_PhenotypePermutations.RData")





##############################################################################################
### Optimization of the ranks ################################################################
##############################################################################################



#############
### Pickrell 
#############

### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimRank_tcell_Pickrell_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0002456",
                                                                  "GO", 
                                                                  Biobase::exprs(pickrell.eset), 
                                                                  pickrell.eset$gender, metric = "rank")

# save results
save(optimRank_tcell_Pickrell_originalphenotype, 
     file = "./Results/optimRank_cP_GSEA_tCell_Pickrell_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_tcell_Pickrell_phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                         geneset = "GO:0002456", 
                                                         geneset_database = "GO",
                                                         expression_data = Biobase::exprs(pickrell.eset), 
                                                         metric = "rank",
                                                         X = phen_pickrell_list)
# save results
save(optimRank_tcell_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_cP_GSEA_tCell_Pickrell_PhenotypePermutations.RData")



### (II) Gene Set "Demethylation" -> GO:0070988

# original phenotype assignment 
optimRank_Demethylation_Pickrell_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0070988",
                                                                          "GO", 
                                                                          Biobase::exprs(pickrell.eset), 
                                                                          pickrell.eset$gender, 
                                                                          metric = "rank")

# save results
save(optimRank_Demethylation_Pickrell_originalphenotype, 
     file = "./Results/optimRank_cP_GSEA_Demethylation_OriginalPhenotype.RData")

# 10 random phenotype permutations
optimRank_Demethylation_Pickrell_phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                                 geneset = "GO:0070988", 
                                                                 geneset_database = "GO",
                                                                 expression_data = Biobase::exprs(pickrell.eset), 
                                                                 metric = "rank",
                                                                 X = phen_pickrell_list)
# save results
save(optimRank_Demethylation_Pickrell_phenotypepermutations, 
     file = "./Results/optimRank_cP_GSEA_Demethylation_Pickrell_PhenotypePermutations.RData")






#############
### Bottomly 
#############

### (I) Gene Set "t-cell mediated immunity" -> GO:0002456

# original phenotype assignment 
optimRank_tCell_Bottomly_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0002456",
                                                                  "GO", 
                                                                  Biobase::exprs(bottomly.eset), 
                                                                  bottomly.eset$strain, 
                                                                  metric = "rank")

# save results
save(optimRank_tCell_Bottomly_originalphenotype, 
     file = "./Results/optimRank_cP_GSEA_tCell_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_tCell_Bottomly_phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                         geneset = "GO:0002456", 
                                                         geneset_database = "GO",
                                                         expression_data = Biobase::exprs(bottomly.eset), 
                                                         metric = "rank",
                                                         X = phen_bottomly_list)
# save results 
save(optimRank_tCell_Bottomly_phenotypepermutations, 
     file = "./Results/optimRank_cP_GSEA_tCell_Bottomly_PhenotypePermutations.RData")



### (II) Gene Set "Demethylation" -> GO:0070988


# original phenotype assignment 
optimRank_Demethylation_Bottomly_originalphenotype <- cP_GSEA_rankp_optim(geneset = "GO:0070988",
                                                                          "GO", 
                                                                          Biobase::exprs(bottomly.eset), 
                                                                          bottomly.eset$strain, 
                                                                          metric = "rank")

# save results
save(optimRank_Demethylation_Bottomly_originalphenotype, 
     file = "./Results/optimRank_cP_GSEA_Demethylation_Bottomly_OriginalPhenotype.RData")


# 10 random phenotype permutations
optimRank_Demethylation_Bottomly_Phenotypepermutations <- lapply(FUN = cP_GSEA_rankp_optim, 
                                                                 geneset = "GO:0070988", 
                                                                 geneset_database = "GO",
                                                                 expression_data = Biobase::exprs(bottomly.eset), 
                                                                 metric = "rank",
                                                                 X = phen_bottomly_list)
# save results 
save(optimRank_Demethylation_Bottomly_Phenotypepermutations, 
     file = "./Results/optimRank_cP_GSEA_Demethylation_Bottomly_PhenotypePermutations.RData")



###








