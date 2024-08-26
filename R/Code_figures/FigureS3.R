################################################################################
### Results Illustrations for the optimization of the p-values for Bottomly data
### and the random permutations of the true sample labels ######################
################################################################################

library(ggplot2)
library(tidyr)
library(cowplot)

################################################################################
### load help functions to generate results illustrations ######################
################################################################################

source("./R/Figures/helpfunctions_results_illustrations_pvalue_rank.R")


################################################################################
################################################################################
# Load optimization results for the TRUE SAMPLE LABELS #########################
################################################################################
################################################################################

#Load GOSeq results

# Demethylation
load("./Results/Intermediate_results/optimP_GOSeq_MetabolicProcess_Bottomly_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/Intermediate_results/optimP_GOSeq_CellularProcess_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/Intermediate_results/optimP_ORA_Demethylation_Bottomly_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/Intermediate_results/optimP_ORA_tCell_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's GSEA results

# Demethylation
load("./Results/Intermediate_results/optimP_cP_GSEA_Demethylation_Bottomly_OriginalPhenotype.RData")

# t-Cell mediated immunity
load("./Results/Intermediate_results/optimP_cP_GSEA_tCell_Bottomly_OriginalPhenotype.RData")


# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/Intermediate_results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_OriginalPhenotype.RData")

# (ii) Graft vs Host
load("./Results/Intermediate_results/optimP_PADOG_GraftvsHost_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's GSEA results

# Demethylation
load("./Results/Intermediate_results/optimP_cP_GSEA_Demethylation_Bottomly_OriginalPhenotype.RData")

# t-Cell mediated immunity
load("./Results/Intermediate_results/optimP_cP_GSEA_tCell_Bottomly_OriginalPhenotype.RData")

# GSEA & GSEAPreranked results
# Since these are both web-based tools, the optimization of the result was performed manually and outside of R.
# -> We therefore have to preprocess the information on the increase of the differentially enriched gene sets manually
# -> see below





################################################################################
### Prepare results ############################################################
################################################################################

# From each stepwise optimization process, we extract the default and the optimal
# results for the true sample labels


################################################################################
### GOSeq ######################################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_GOSeq_CellularProcess_truephen_default <-
  optimP_GOSeq_CellularProcess_Bottomly_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_GOSeq_CellularProcess_truephen_optim <-
  optimP_GOSeq_CellularProcess_Bottomly_originalphenotype$documentation$p_adj[length(optimP_GOSeq_CellularProcess_Bottomly_originalphenotype$documentation$p_adj)]


# Gene set Demethylation

# default p_adj
padj_GOSeq_MetabolicProcess_truephen_default <-
  optimP_GOSeq_MetabolicProcess_Bottomly_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_GOSeq_MetabolicProcess_truephen_optim <-
  optimP_GOSeq_MetabolicProcess_Bottomly_originalphenotype$documentation$p_adj[length(optimP_GOSeq_MetabolicProcess_Bottomly_originalphenotype$documentation$p_adj)]

################################################################################
### DAVID ######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# default p_adj
padj_DAVID_CellularProcess_truephen_default <- 1
padj_DAVID_CellularProcess_truephen_optim <- 1


# Gene set Demethylation

# default p_adj
padj_DAVID_MetabolicProcess_truephen_default <- 1
# optimal p_adj
padj_DAVID_MetabolicProcess_truephen_optim <- 1

################################################################################
### clusterProfiler's ORA ######################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_cP_ORA_tCell_truephen_default <-
  optimP_cP_ORA_tCell_Bottomly_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_cP_ORA_tCell_truephen_optim <-
  optimP_cP_ORA_tCell_Bottomly_originalphenotype$documentation$p_adj[length(optimP_cP_ORA_tCell_Bottomly_originalphenotype$documentation$p_adj)]


# Gene set Demethylation

# optimal p_adj
padj_cP_ORA_Demethylation_truephen_optim <-
  optimP_cP_ORA_Demethylation_Bottomly_originalphenotype$documentation$p_adj[length(optimP_cP_ORA_Demethylation_Bottomly_originalphenotype$documentation$p_adj)]
# default p_adj
padj_cP_ORA_Demethylation_truephen_default <-
  optimP_cP_ORA_Demethylation_Bottomly_originalphenotype$documentation$p_adj[1]





################################################################################
### PADOG ######################################################################
################################################################################

# Gene set Primary Immunodeficiency

# default p_adj
padj_PADOG_PrimImmun_truephen_default <-
  optimP_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_PADOG_PrimImmun_truephen_optim <-
  optimP_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype$documentation$p_adj[length(optimP_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype$documentation$p_adj)]


# Gene set Graft vs Host

# default p_adj
padj_PADOG_GraftvsHost_truephen_default <-
  optimP_PADOG_GraftvsHost_Bottomly_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_PADOG_GraftvsHost_truephen_optim <-
  optimP_PADOG_GraftvsHost_Bottomly_originalphenotype$documentation$p_adj[length(optimP_PADOG_GraftvsHost_Bottomly_originalphenotype$documentation$p_adj)]


################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_cP_GSEA_tCell_truephen_default <-
  optimP_cP_GSEA_tCell_Bottomly_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_cP_GSEA_tCell_truephen_optim <-
  optimP_cP_GSEA_tCell_Bottomly_originalphenotype$documentation$p_adj[length(optimP_cP_GSEA_tCell_Bottomly_originalphenotype$documentation$p_adj)]


# Gene set Demethylation

# optimal p_adj
padj_cP_GSEA_Demethylation_truephen_optim <-
  optimP_cP_GSEA_Demethylation_Bottomly_originalphenotype$documentation$p_adj[length(optimP_cP_GSEA_Demethylation_Bottomly_originalphenotype$documentation$p_adj)]
# default p_adj
padj_cP_GSEA_Demethylation_truephen_default <-
  optimP_cP_GSEA_Demethylation_Bottomly_originalphenotype$documentation$p_adj[1]



################################################################################
### GSEA #######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# default p_adj
padj_GSEA_tCell_truephen_default <- 0.9645
# optimal p_adj
padj_GSEA_tCell_truephen_optim <- 0.8638


# Gene set Demethylation

# default p_adj
padj_GSEA_Demethylation_truephen_default <- 0.9703
# optimal p_adj
padj_GSEA_Demethylation_truephen_optim <- 0.8922



################################################################################
### GSEAPreranked ##############################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_GSEAPreranked_tCell_truephen_default <- 1
# optimal p_adj
padj_GSEAPreranked_tCell_truephen_optim <- 1


# Gene set Demethylation

# default p_adj
padj_GSEAPreranked_Demethylation_truephen_default <- 1
# optimal p_adj
padj_GSEAPreranked_Demethylation_truephen_optim <- 0.7782


################################################################################
### Illustration plots #########################################################
################################################################################

# In the following section, we plot all default vs. adjusted p-values for all
# optimization processes performed in the Bottomly data set

# We generate ONE plot which contains the default vs. adjusted p-values

# - for all GSA tools (we stratify the plot by the GSA tools)
# - for all gene sets considered (indicated by color)


# At first, we prepare two data sets separately (each containing the adjusted
# p-values for ONE of the two gene sets for each tools, respectively)
# We then combine both data sets later on and generate one plot



# (i) Random permutations of the random permutations of the sample conditions, first data set
padj_GeneSet1_truephen <- data.frame(cP_ORA = c(padj_cP_ORA_tCell_truephen_default, padj_cP_ORA_tCell_truephen_optim),
                                     GOSeq = c(padj_GOSeq_MetabolicProcess_truephen_default, padj_GOSeq_MetabolicProcess_truephen_optim),
                                     DAVID = c(padj_DAVID_MetabolicProcess_truephen_default, padj_DAVID_MetabolicProcess_truephen_optim),
                                     PADOG = c(padj_PADOG_PrimImmun_truephen_default, padj_PADOG_PrimImmun_truephen_optim),
                                     cP_GSEA = c(padj_cP_GSEA_tCell_truephen_default, padj_cP_GSEA_tCell_truephen_optim),
                                     GSEA = c(padj_GSEA_tCell_truephen_default,padj_GSEA_tCell_truephen_optim),
                                     GSEAPreranked = c(padj_GSEAPreranked_tCell_truephen_default, padj_GSEAPreranked_tCell_truephen_optim),
                                     state = c("Default","Maximum")) # state: default vs. optimum




# (ii) Random permutations of the random permutations of the sample conditions, second data set
padj_GeneSet2_truephen <- data.frame(cP_ORA = c(padj_cP_ORA_Demethylation_truephen_default, padj_cP_ORA_Demethylation_truephen_optim),
                                     GOSeq = c(padj_GOSeq_CellularProcess_truephen_default, padj_GOSeq_CellularProcess_truephen_optim),
                                     DAVID = c(padj_DAVID_CellularProcess_truephen_default, padj_DAVID_CellularProcess_truephen_optim),
                                     PADOG = c(padj_PADOG_GraftvsHost_truephen_default, padj_PADOG_GraftvsHost_truephen_optim),
                                     GSEA = c(padj_GSEA_Demethylation_truephen_default,padj_GSEA_Demethylation_truephen_optim),
                                     cP_GSEA = c(padj_cP_GSEA_Demethylation_truephen_default, padj_cP_GSEA_Demethylation_truephen_optim),
                                     GSEAPreranked = c(padj_GSEAPreranked_Demethylation_truephen_default, padj_GSEAPreranked_Demethylation_truephen_optim),
                                     state = c("Default","Maximum")) # state: default vs. optimum




################################################################################
################################################################################
# Load optimization results for the PERMUTED SAMPLE LABELS #####################
################################################################################
################################################################################



#Load GOSeq results

# Demethylation
load("./Results/Intermediate_results/optimP_GOSeq_MetabolicProcess_Bottomly_PhenotypePermutations.RData")

# t Cell mediated immunity
load("./Results/Intermediate_results/optimP_GOSeq_CellularProcess_Bottomly_PhenotypePermutations.RData")


# Load clusterProfiler's ORA results

# t Cell mediated immunity
load("./Results/Intermediate_results/optimP_ORA_tCell_Bottomly_PhenotypePermutations.RData")

# Demethylation
load("./Results/Intermediate_results/optimP_ORA_Demethylation_Bottomly_PhenotypePermutations.RData")


# Load DAVID results : Achtung DAVID fehlt bis jetzt noch


# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/Intermediate_results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_PhenotypePermutations.RData")

# (ii) Graft vs Host
load("./Results/Intermediate_results/optimP_PADOG_GraftvsHost_Bottomly_PhenotypePermutations.RData")


# Load clusterProfiler's GSEA results

# Demethylation
load("./Results/Intermediate_results/optimP_cP_GSEA_Demethylation_Bottomly_PhenotypePermutations.RData")

# t-Cell mediated immunity
load("./Results/Intermediate_results/optimP_cP_GSEA_tCell_Bottomly_PhenotypePermutations.RData")

# GSEA & GSEAPreranked results
# Since these are both web-based tools, the optimization of the result was performed manually and outside of R.
# -> We therefore have to preprocess the information on the increase of the differentially enriched gene sets manually
# -> see below



################################################################################
### Prepare results ############################################################
################################################################################

# From each stepwise optimization process, we extract the default and the optimal
# results for the 10 permutations of the true sample labels


################################################################################
### GOSeq ######################################################################
################################################################################

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_GOSeq_CellularProcess_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_GOSeq_CellularProcess_phenpermutation_optim <- c()


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
padj_GOSeq_MetabolicProcess_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_GOSeq_MetabolicProcess_phenpermutation_optim <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  padj_GOSeq_CellularProcess_phenpermutation_default[i] <-
    optimP_GOSeq_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_GOSeq_CellularProcess_phenpermutation_optim[i] <-
    optimP_GOSeq_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[length(optimP_GOSeq_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  padj_GOSeq_MetabolicProcess_phenpermutation_default[i] <-
    optimP_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_GOSeq_MetabolicProcess_phenpermutation_optim[i] <-
    optimP_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$p_adj)]

}

################################################################################
### DAVID ######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation


# Cellular Process

# store default p_adj for the 10 permutations of the true sample labels
padj_DAVID_CellularProcess_phenpermutation_default <- rep(1, times = 10)

# store optimal p_adj for the 10 permutations of the true sample labels
padj_DAVID_CellularProcess_phenpermutation_optim <- rep(1, times = 10)



# Metabolic Process

# store default p_adj for the 10 permutations of the true sample labels

padj_DAVID_MetabolicProcess_phenpermutation_default <- rep(1, times = 10)


# store optimal p_adj for the 10 permutations of the true sample labels
padj_DAVID_MetabolicProcess_phenpermutation_optim <- rep(1, times = 10)






################################################################################
### clusterProfiler's ORA ######################################################
################################################################################


# Demethylation

# store default p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_Demethylation_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_Demethylation_phenpermutation_optim <- c()


# t Cell mediated immunity
# store default p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_tCell_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_tCell_phenpermutation_optim <- c()


for(i in 1:10){

  # Cellular Process

  # store the default number of differentially enriched gene sets
  padj_cP_ORA_Demethylation_phenpermutation_default[i] <-
    optimP_cP_ORA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_ORA_Demethylation_phenpermutation_optim[i] <-
    optimP_cP_ORA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj)]


  # Metabolic Process

  # store the default number of differentially enriched gene sets
  padj_cP_ORA_tCell_phenpermutation_default[i] <-
    optimP_cP_ORA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_ORA_tCell_phenpermutation_optim[i] <-
    optimP_cP_ORA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$p_adj)]

}

################################################################################
### PADOG ######################################################################
################################################################################

# Primary Immundeficiency

# store default p_adj for the 10 permutations of the true sample labels
padj_PADOG_PrimImmun_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_PADOG_PrimImmun_phenpermutation_optim <- c()


# Graft versus host disease
# store default p_adj for the 10 permutations of the true sample labels
padj_PADOG_GraftvsHost_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_PADOG_GraftvsHost_phenpermutation_optim <- c()


for(i in 1:10){

  # Primary immunodeficiency

  # store the default number of differentially enriched gene sets
  padj_PADOG_PrimImmun_phenpermutation_default[i] <- optimP_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_PADOG_PrimImmun_phenpermutation_optim[i] <- optimP_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations[[i]]$documentation$p_adj)]


  # Graft versus host disease

  # store the default number of differentially enriched gene sets
  padj_PADOG_GraftvsHost_phenpermutation_default[i] <- optimP_PADOG_GraftvsHost_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_PADOG_GraftvsHost_phenpermutation_optim[i] <- optimP_PADOG_GraftvsHost_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_PADOG_GraftvsHost_Bottomly_phenotypepermutations[[i]]$documentation$p_adj)]

}


################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_cP_GSEA_tCell_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_cP_GSEA_tCell_phenpermutation_optim <- c()


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
padj_cP_GSEA_Demethylation_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_cP_GSEA_Demethylation_phenpermutation_optim <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  padj_cP_GSEA_tCell_phenpermutation_default[i] <- optimP_cP_GSEA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_GSEA_tCell_phenpermutation_optim[i] <- optimP_cP_GSEA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_GSEA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$p_adj)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  padj_cP_GSEA_Demethylation_phenpermutation_default[i] <- optimP_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_GSEA_Demethylation_phenpermutation_optim[i] <- optimP_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj)]

}










################################################################################
### GSEA #######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation


# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_GSEA_tCell_phenpermutation_default <- c(0.7845,
                                             0.8578,
                                             1,
                                             0.8031,
                                             0.8984,
                                             1,
                                             1,
                                             0.9167,
                                             1,
                                             0.7764)

# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEA_tCell_phenpermutation_optim <- c(0.5144,
                                           0.6015,
                                           0.6832,
                                           0.609,
                                           0.831,
                                           0.8656,
                                           0.9373,
                                           0.499,
                                           0.6214,
                                           0.4097)



# Demethylation

# store default p_adj for the 10 permutations of the true sample labels

padj_GSEA_Demethylation_phenpermutation_default <- c(0.8478,
                                                     0.9590,
                                                     1,
                                                     0.5624,
                                                     1,
                                                     0.9598,
                                                     0.9158,
                                                     0.7584,
                                                     0.8127,
                                                     0.6563)


# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEA_Demethylation_phenpermutation_optim <- c(0.5555,
                                                   0.9590,
                                                   0.8154,
                                                   0.3003,
                                                   0.8720,
                                                   0.6045,
                                                   0.6694,
                                                   0.5802,
                                                   0.7064,
                                                   0.3079)




################################################################################
### GSEAPreranked ##############################################################
################################################################################

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_tCell_phenpermutation_default <- c(0.6605,
                                                      1,
                                                      0.9798,
                                                      0.5763,
                                                      0.8227,
                                                      0.9372,
                                                      0.8856,
                                                      0.6276,
                                                      0.9904,
                                                      0.8831)

# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_tCell_phenpermutation_optim <- c(0.6605,
                                                    0.7542,
                                                    0.759,
                                                    0.5763,
                                                    0.8227,
                                                    0.8341,
                                                    0.7937,
                                                    0.4338,
                                                    0.0581,
                                                    0.6486)


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_Demethylation_phenpermutation_default <- c(0.9623,
                                                              0.8051,
                                                              1,
                                                              0.5703,
                                                              0.8303,
                                                              0.9600,
                                                              0.9424,
                                                              0.5048,
                                                              1,
                                                              0.2525)
# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_Demethylation_phenpermutation_optim <- c(0.8012,
                                                            0.7094,
                                                            0.8272,
                                                            0.187,
                                                            0.6328,
                                                            0.5234,
                                                            0.8025,
                                                            0.3213,
                                                            0.3808,
                                                            0.2525)

################################################################################
### Illustration plots #########################################################
################################################################################

# In the following section, we plot all default vs. adjusted p-values for all
# optimization processes performed in the Bottomly data set

# We generate ONE plot which contains the default vs. adjusted p-values

# - for all GSA tools (we stratify the plot by the GSA tools)
# - for all gene sets considered (indicated by color)


# At first, we prepare two data sets separately (each containing the adjusted
# p-values for ONE of the two gene sets for each tools, respectively)
# We then combine both data sets later on and generate one plot



# (i) Random permutations of the random permutations of the sample conditions, first data set
padj_GeneSet1_phenpermutation <- data.frame(cP_ORA = c(padj_cP_ORA_tCell_phenpermutation_default, padj_cP_ORA_tCell_phenpermutation_optim),
                                            GOSeq = c(padj_GOSeq_MetabolicProcess_phenpermutation_default, padj_GOSeq_MetabolicProcess_phenpermutation_optim),
                                            DAVID = c(padj_DAVID_CellularProcess_phenpermutation_default, padj_DAVID_CellularProcess_phenpermutation_optim),
                                            PADOG = c(padj_PADOG_PrimImmun_phenpermutation_default, padj_PADOG_PrimImmun_phenpermutation_optim),
                                            cP_GSEA = c(padj_cP_GSEA_tCell_phenpermutation_default, padj_cP_GSEA_tCell_phenpermutation_optim),
                                            GSEA = c(padj_GSEA_tCell_phenpermutation_default,padj_GSEA_tCell_phenpermutation_optim),
                                            GSEAPreranked = c(padj_GSEAPreranked_tCell_phenpermutation_default, padj_GSEAPreranked_tCell_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum
                                            ID = rep(c(1:10),2)) # add the number of the permutation



# (ii) Random permutations of the random permutations of the sample conditions, second data set
padj_GeneSet2_phenpermutation <- data.frame(cP_ORA = c(padj_cP_ORA_Demethylation_phenpermutation_default, padj_cP_ORA_Demethylation_phenpermutation_optim),
                                            GOSeq = c(padj_GOSeq_CellularProcess_phenpermutation_default, padj_GOSeq_CellularProcess_phenpermutation_optim),
                                            DAVID = c(padj_DAVID_MetabolicProcess_phenpermutation_default, padj_DAVID_MetabolicProcess_phenpermutation_optim),
                                            PADOG = c(padj_PADOG_GraftvsHost_phenpermutation_default, padj_PADOG_GraftvsHost_phenpermutation_optim),
                                            GSEA = c(padj_GSEA_Demethylation_phenpermutation_default,padj_GSEA_Demethylation_phenpermutation_optim),
                                            cP_GSEA = c(padj_cP_GSEA_Demethylation_phenpermutation_default, padj_cP_GSEA_Demethylation_phenpermutation_optim),
                                            GSEAPreranked = c(padj_GSEAPreranked_Demethylation_phenpermutation_default, padj_GSEAPreranked_Demethylation_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum
                                            ID = rep(c(1:10),2)) # add the number of the permutation





##########
### ggplot
##########



plot_truelabels <- create_results_illustration_pvalue_rank(padj_GeneSet1_truephen, padj_GeneSet2_truephen,
                                                           "p_adj", "true_labels")

plot_permutedlabels <- create_results_illustration_pvalue_rank(padj_GeneSet1_phenpermutation, padj_GeneSet2_phenpermutation,
                                                               "p_adj", "random_permutations")

plot_grid(plot_permutedlabels, plot_truelabels, labels=c("A", "B"), ncol = 1, nrow = 2)

## uncomment to save
ggsave2("./Results/Figures/FigureS3.eps",
        width = 10,
        height = 12)






