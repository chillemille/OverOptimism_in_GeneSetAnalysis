################################################################################
### Results Illustrations for the optimization of the ranks for Bottomly data
### and the random permutations of the true sample labels ######################
################################################################################

library(ggplot2)
library(tidyr)
library(cowplot)

################################################################################
### load help functions to generate results illustrations ######################
################################################################################

source("./Results_illustrations/helpfunctions_results_illustrations_pvalue_rank.R")

################################################################################
################################################################################
# Load optimization results for the TRUE SAMPLE LABELS #########################
################################################################################
################################################################################

#Load GOSeq results

# Demethylation
load("./Results/optimRank_GOSeq_MetabolicProcess_Bottomly_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/optimRank_GOSeq_CellularProcess_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/optimRank_ORA_tCell_Bottomly_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/optimRank_ORA_Demethylation_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's GSEA results

# Demethylation
load("./Results/optimRank_cP_GSEA_tCell_Bottomly_OriginalPhenotype.RData")

# t-Cell mediated immunity
load("./Results/optimRank_cP_GSEA_Demethylation_Bottomly_OriginalPhenotype.RData")


# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_OriginalPhenotype.RData")

# (ii) Graft vs Host
load("./Results/optimRank_PADOG_GraftvsHost_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's GSEA results



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
rel_rank_GOSeq_CellularProcess_truephen_default <-
  optimRank_GOSeq_CellularProcess_Bottomly_originalphenotype$documentation$rank[1]
# optimal p_adj
rel_rank_GOSeq_CellularProcess_truephen_optim <-
  optimRank_GOSeq_CellularProcess_Bottomly_originalphenotype$documentation$rank[length(optimRank_GOSeq_CellularProcess_Bottomly_originalphenotype$documentation$rank)]


# Gene set Demethylation

# default p_adj
rel_rank_GOSeq_MetabolicProcess_truephen_default <-
  optimRank_GOSeq_MetabolicProcess_Bottomly_originalphenotype$documentation$rank[1]
# optimal p_adj
rel_rank_GOSeq_MetabolicProcess_truephen_optim <-
  optimRank_GOSeq_MetabolicProcess_Bottomly_originalphenotype$documentation$rank[length(optimRank_GOSeq_MetabolicProcess_Bottomly_originalphenotype$documentation$rank)]

################################################################################
### DAVID ######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# default p_adj
rel_rank_DAVID_CellularProcess_truephen_default <- 1
rel_rank_DAVID_CellularProcess_truephen_optim <- 1


# Gene set Demethylation

# default p_adj
rel_rank_DAVID_MetabolicProcess_truephen_default <- 1
# optimal p_adj
rel_rank_DAVID_MetabolicProcess_truephen_optim <- 1

################################################################################
### clusterProfiler's ORA ######################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
rel_rank_cP_ORA_tCell_truephen_default <-
  optimRank_cP_ORA_tCell_Bottomly_originalphenotype$documentation$rank[1]
# optimal p_adj
rel_rank_cP_ORA_tCell_truephen_optim <-
  optimRank_cP_ORA_tCell_Bottomly_originalphenotype$documentation$rank[length(optimRank_cP_ORA_tCell_Bottomly_originalphenotype$documentation$rank)]


# Gene set Demethylation

# optimal p_adj
rel_rank_cP_ORA_Demethylation_truephen_optim <-
  optimRank_cP_ORA_Demethylation_Bottomly_originalphenotype$documentation$rank[length(optimRank_cP_ORA_Demethylation_Bottomly_originalphenotype$documentation$rank)]
# default p_adj
rel_rank_cP_ORA_Demethylation_truephen_default <-
  optimRank_cP_ORA_Demethylation_Bottomly_originalphenotype$documentation$rank[1]





################################################################################
### PADOG ######################################################################
################################################################################

# Gene set Primary Immunodeficiency

# default p_adj
rel_rank_PADOG_PrimImmun_truephen_default <-
  optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype$documentation$rank[1]
# optimal p_adj
rel_rank_PADOG_PrimImmun_truephen_optim <-
  optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype$documentation$rank[length(optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_originalphenotype$documentation$rank)]


# Gene set Graft vs Host

# default p_adj
rel_rank_PADOG_GraftvsHost_truephen_default <-
  optimRank_PADOG_GraftvsHost_Bottomly_originalphenotype$documentation$rank[1]
# optimal p_adj
rel_rank_PADOG_GraftvsHost_truephen_optim <-
  optimRank_PADOG_GraftvsHost_Bottomly_originalphenotype$documentation$rank[length(optimRank_PADOG_GraftvsHost_Bottomly_originalphenotype$documentation$rank)]


################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
rel_rank_cP_GSEA_tCell_truephen_default <-
  optimRank_tCell_Bottomly_originalphenotype$documentation$rank[1]
# optimal p_adj
rel_rank_cP_GSEA_tCell_truephen_optim <-
  optimRank_tCell_Bottomly_originalphenotype$documentation$rank[length(optimRank_tCell_Bottomly_originalphenotype$documentation$rank)]


# Gene set Demethylation

# optimal p_adj
rel_rank_cP_GSEA_Demethylation_truephen_optim <-
  optimRank_Demethylation_Bottomly_originalphenotype$documentation$rank[length(optimRank_Demethylation_Bottomly_originalphenotype$documentation$rank)]
# default p_adj
rel_rank_cP_GSEA_Demethylation_truephen_default <-
  optimRank_Demethylation_Bottomly_originalphenotype$documentation$rank[1]



################################################################################
### GSEA #######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# default p_adj
rel_rank_GSEA_tCell_truephen_default <- 0.9645
# optimal p_adj
rel_rank_GSEA_tCell_truephen_optim <- 0.8638


# Gene set Demethylation

# default p_adj
rel_rank_GSEA_Demethylation_truephen_default <- 0.9703
# optimal p_adj
rel_rank_GSEA_Demethylation_truephen_optim <- 0.8922



################################################################################
### GSEAPreranked ##############################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
rel_rank_GSEAPreranked_tCell_truephen_default <- 1
# optimal p_adj
rel_rank_GSEAPreranked_tCell_truephen_optim <- 1


# Gene set Demethylation

# default p_adj
rel_rank_GSEAPreranked_Demethylation_truephen_default <- 1
# optimal p_adj
rel_rank_GSEAPreranked_Demethylation_truephen_optim <- 0.7782


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
rel_rank_GeneSet1_truephen <- data.frame(cP_ORA = c(rel_rank_cP_ORA_tCell_truephen_default, rel_rank_cP_ORA_tCell_truephen_optim),
                                         GOSeq = c(rel_rank_GOSeq_MetabolicProcess_truephen_default, rel_rank_GOSeq_MetabolicProcess_truephen_optim),
                                         PADOG = c(rel_rank_PADOG_PrimImmun_truephen_default, rel_rank_PADOG_PrimImmun_truephen_optim),
                                         cP_GSEA = c(rel_rank_cP_GSEA_tCell_truephen_default, rel_rank_cP_GSEA_tCell_truephen_optim),
                                         state = c("Default","Maximum")) # state: default vs. optimum






# (ii) Random permutations of the random permutations of the sample conditions, second data set
rel_rank_GeneSet2_truephen <- data.frame(cP_ORA = c(rel_rank_cP_ORA_Demethylation_truephen_default, rel_rank_cP_ORA_Demethylation_truephen_optim),
                                         GOSeq = c(rel_rank_GOSeq_CellularProcess_truephen_default, rel_rank_GOSeq_CellularProcess_truephen_optim),
                                         PADOG = c(rel_rank_PADOG_GraftvsHost_truephen_default, rel_rank_PADOG_GraftvsHost_truephen_optim),
                                         cP_GSEA = c(rel_rank_cP_GSEA_Demethylation_truephen_default, rel_rank_cP_GSEA_Demethylation_truephen_optim),
                                         state = c("Default","Maximum")) # state: default vs. optimum





################################################################################
################################################################################
# Load optimization results for the PERMUTED SAMPLE LABELS #####################
################################################################################
################################################################################

#Load GOSeq results

# Demethylation
load("./Results/optimRank_GOSeq_MetabolicProcess_Bottomly_PhenotypePermutations.RData")

# t Cell mediated immunity
load("./Results/optimRank_GOSeq_CellularProcess_Bottomly_PhenotypePermutations.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/optimRank_ORA_tCell_Bottomly_PhenotypePermutations.RData")

# t Cell mediated immunity
load("./Results/optimRank_ORA_Demethylation_Bottomly_PhenotypePermutations.RData")




# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_PhenotypePermutations.RData")

# (ii) Graft vs Host
load("./Results/optimRank_PADOG_GraftvsHost_Bottomly_PhenotypePermutations.RData")


# Load clusterProfiler's GSEA results

# Demethylation
load("./Results/optimRank_cP_GSEA_Demethylation_Bottomly_PhenotypePermutations.RData")

# t-Cell mediated immunity
load("./Results/optimRank_cP_GSEA_tCell_Bottomly_PhenotypePermutations.RData")





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
rank_GOSeq_CellularProcess_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_GOSeq_CellularProcess_phenpermutation_optim <- c()


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
rank_GOSeq_MetabolicProcess_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_GOSeq_MetabolicProcess_phenpermutation_optim <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  rank_GOSeq_CellularProcess_phenpermutation_default[i] <-
    optimRank_GOSeq_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_GOSeq_CellularProcess_phenpermutation_optim[i] <-
    optimRank_GOSeq_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$rank[length(optimRank_GOSeq_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$rank)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  rank_GOSeq_MetabolicProcess_phenpermutation_default[i] <-
    optimRank_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_GOSeq_MetabolicProcess_phenpermutation_optim[i] <-
    optimRank_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$rank[length(optimRank_GOSeq_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$rank)]

}


################################################################################
### clusterProfiler's ORA ######################################################
################################################################################


# Demethylation

# store default p_adj for the 10 permutations of the true sample labels
rank_cP_ORA_Demethylation_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_cP_ORA_Demethylation_phenpermutation_optim <- c()


# t Cell mediated immunity
# store default p_adj for the 10 permutations of the true sample labels
rank_cP_ORA_tCell_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_cP_ORA_tCell_phenpermutation_optim <- c()


for(i in 1:10){

  # Cellular Process

  # store the default number of differentially enriched gene sets
  rank_cP_ORA_Demethylation_phenpermutation_default[i] <-
    optimRank_cP_ORA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_cP_ORA_Demethylation_phenpermutation_optim[i] <-
    optimRank_cP_ORA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank[length(optimRank_cP_ORA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank)]


  # Metabolic Process

  # store the default number of differentially enriched gene sets
  rank_cP_ORA_tCell_phenpermutation_default[i] <-
    optimRank_cP_ORA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_cP_ORA_tCell_phenpermutation_optim[i] <-
    optimRank_cP_ORA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank[length(optimRank_cP_ORA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank)]

}

################################################################################
### PADOG ######################################################################
################################################################################

# Primary Immundeficiency

# store default p_adj for the 10 permutations of the true sample labels
rank_PADOG_PrimImmun_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_PADOG_PrimImmun_phenpermutation_optim <- c()


# Graft versus host disease
# store default p_adj for the 10 permutations of the true sample labels
rank_PADOG_GraftvsHost_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_PADOG_GraftvsHost_phenpermutation_optim <- c()


for(i in 1:10){

  # Primary immunodeficiency

  # store the default number of differentially enriched gene sets
  rank_PADOG_PrimImmun_phenpermutation_default[i] <- optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_PADOG_PrimImmun_phenpermutation_optim[i] <- optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations[[i]]$documentation$rank[length(optimRank_PADOG_PrimaryImmunodeficiency_Bottomly_phenotypepermutations[[i]]$documentation$rank)]


  # Graft versus host disease

  # store the default number of differentially enriched gene sets
  rank_PADOG_GraftvsHost_phenpermutation_default[i] <- optimRank_PADOG_GraftvsHost_Bottomly_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_PADOG_GraftvsHost_phenpermutation_optim[i] <- optimRank_PADOG_GraftvsHost_Bottomly_phenotypepermutations[[i]]$documentation$rank[length(optimRank_PADOG_GraftvsHost_Bottomly_phenotypepermutations[[i]]$documentation$rank)]

}


################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
rank_cP_GSEA_tCell_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_cP_GSEA_tCell_phenpermutation_optim <- c()


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
rank_cP_GSEA_Demethylation_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
rank_cP_GSEA_Demethylation_phenpermutation_optim <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  rank_cP_GSEA_tCell_phenpermutation_default[i] <- optimRank_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_cP_GSEA_tCell_phenpermutation_optim[i] <- optimRank_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank[length(optimRank_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  rank_cP_GSEA_Demethylation_phenpermutation_default[i] <- optimRank_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rank_cP_GSEA_Demethylation_phenpermutation_optim[i] <- optimRank_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank[length(optimRank_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank)]

}



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
rel_rank_GeneSet1_phenpermutation <- data.frame(cP_ORA = c(rank_cP_ORA_tCell_phenpermutation_default, rank_cP_ORA_tCell_phenpermutation_optim),
                                            GOSeq = c(rank_GOSeq_MetabolicProcess_phenpermutation_default, rank_GOSeq_MetabolicProcess_phenpermutation_optim),
                                            PADOG = c(rank_PADOG_PrimImmun_phenpermutation_default, rank_PADOG_PrimImmun_phenpermutation_optim),
                                            cP_GSEA = c(rank_cP_GSEA_tCell_phenpermutation_default, rank_cP_GSEA_tCell_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum
                                            ID = rep(c(1:10),2)) # add the number of the permutation



# (ii) Random permutations of the random permutations of the sample conditions, second data set
rel_rank_GeneSet2_phenpermutation <- data.frame(cP_ORA = c(rank_cP_ORA_Demethylation_phenpermutation_default, rank_cP_ORA_Demethylation_phenpermutation_optim),
                                            GOSeq = c(rank_GOSeq_CellularProcess_phenpermutation_default, rank_GOSeq_CellularProcess_phenpermutation_optim),
                                            PADOG = c(rank_PADOG_GraftvsHost_phenpermutation_default, rank_PADOG_GraftvsHost_phenpermutation_optim),
                                            cP_GSEA = c(rank_cP_GSEA_Demethylation_phenpermutation_default, rank_cP_GSEA_Demethylation_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum
                                            ID = rep(c(1:10),2)) # add the number of the permutation





################################################################################
### Generate ggplots ###########################################################
################################################################################

plot_truelabels <- create_results_illustration_pvalue_rank(rel_rank_GeneSet1_truephen, rel_rank_GeneSet2_truephen,
                                                           "rel_rank", "true_labels")

plot_permutedlabels <- create_results_illustration_pvalue_rank(rel_rank_GeneSet1_phenpermutation, rel_rank_GeneSet2_phenpermutation,
                                                               "rel_rank", "random_permutations")

plot_relrank <- plot_grid(plot_truelabels, plot_permutedlabels, labels=c("A", "B"), ncol = 1, nrow = 2)

## uncomment to save
ggsave("./Results_illustrations/FigureS4.pdf",
       width = 10,
       height = 12)


