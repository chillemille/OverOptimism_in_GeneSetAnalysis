################################################################################
### Results Illustrations for the optimization of the p-values for Pickrell data
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
load("./Results/optimP_GOSeq_Demethylation_Pickrell_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/optimP_GOSeq_tCell_Pickrell_OriginalPhenotype.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/optimP_ORA_Demethylation_Pickrell_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/optimP_ORA_tCell_Pickrell_OriginalPhenotype.RData")


# Load DAVID results : Achtung DAVID fehlt bis jetzt noch

# Load clusterProfiler's GSEA results

# Demethylation
load("./Results/optimP_cP_GSEA_Demethylation_Pickrell_OriginalPhenotype.RData")

# t-Cell mediated immunity
load("./Results/optimP_cP_GSEA_tCell_Pickrell_OriginalPhenotype.RData")


# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/optimP_PADOG_PrimaryImmunodeficiency_Pickrell_OriginalPhenotype.RData")

# (ii) Graft vs Host
load("./Results/optimP_PADOG_GraftvsHost_Pickrell_OriginalPhenotype.RData")

# GSEA & GSEAPreranked results
# Since these are both web-based tools, the optimization of the result was performed manually and outside of R.
# -> We therefore have to preprocess the information on the increase of the differentially enriched gene sets manually
# -> see below

# load clusterProfiler's GSEA results

# (i) t Cell mediated immunity
load("./Results/optimP_cP_GSEA_tCell_Pickrell_OriginalPhenotype.RData")

# (ii) Demethylation
load("./Results/optimP_cP_GSEA_Demethylation_Pickrell_OriginalPhenotype.RData")

################################################################################
### Prepare results ############################################################
################################################################################

# From each stepwise optimization process, we extract the default and the optimal
# results

################################################################################
### GOSeq ######################################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_GOSeq_tCell_truephen_default <-
  optimP_GOSeq_tcell_Pickrell_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_GOSeq_tCell_truephen_optim <-
  optimP_GOSeq_tcell_Pickrell_originalphenotype$documentation$p_adj[length(optimP_GOSeq_tcell_Pickrell_originalphenotype$documentation$p_adj)]


# Gene set Demethylation

# default p_adj
padj_GOSeq_Demethylation_truephen_default <-
  optimP_GOSeq_Demethylation_Pickrell_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_GOSeq_Demethylation_truephen_optim <-
  optimP_GOSeq_Demethylation_Pickrell_originalphenotype$documentation$p_adj[length(optimP_GOSeq_Demethylation_Pickrell_originalphenotype$documentation$p_adj)]


################################################################################
### DAVID ######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# default p_adj
padj_DAVID_tCell_truephen_default <- 1
padj_DAVID_tCell_truephen_optim <- 1


# Gene set Demethylation

# default p_adj
padj_DAVID_Demethylation_truephen_default <- 1
# optimal p_adj
padj_DAVID_Demethylation_truephen_optim <- 1


################################################################################
### clusterProfiler's ORA ######################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_cP_ORA_tCell_truephen_default <-
  optimP_cP_ORA_tcell_Pickrell_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_cP_ORA_tCell_truephen_optim <-
  optimP_cP_ORA_tcell_Pickrell_originalphenotype$documentation$p_adj[length(optimP_cP_ORA_tcell_Pickrell_originalphenotype$documentation$p_adj)]


# Gene set Demethylation

# optimal p_adj
padj_cP_ORA_Demethylation_truephen_optim <-
  optimP_cP_ORA_Demethylation_Pickrell_originalphenotype$documentation$p_adj[length(optimP_cP_ORA_Demethylation_Pickrell_originalphenotype$documentation$p_adj)]
# default p_adj
padj_cP_ORA_Demethylation_truephen_default <-
  optimP_cP_ORA_Demethylation_Pickrell_originalphenotype$documentation$p_adj[1]



################################################################################
### PADOG ######################################################################
################################################################################

# Gene set Primary Immunodeficiency

# default p_adj
padj_PADOG_PrimImmun_truephen_default <-
  optimP_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_PADOG_PrimImmun_truephen_optim <-
  optimP_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype$documentation$p_adj[length(optimP_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype$documentation$p_adj)]


# Gene set Graft vs Host

# default p_adj
padj_PADOG_GraftvsHost_truephen_default <-
  optimP_PADOG_GraftvsHost_Pickrell_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_PADOG_GraftvsHost_truephen_optim <-
  optimP_PADOG_GraftvsHost_Pickrell_originalphenotype$documentation$p_adj[length(optimP_PADOG_GraftvsHost_Pickrell_originalphenotype$documentation$p_adj)]


################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################


# Gene set t Cell mediated immunity

# default p_adj
padj_cP_GSEA_tCell_truephen_default <-
  optimP_cP_GSEA_tcell_Pickrell_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_cP_GSEA_tCell_truephen_optim <-
  optimP_cP_GSEA_tcell_Pickrell_originalphenotype$documentation$p_adj[length(optimP_cP_GSEA_tcell_Pickrell_originalphenotype$documentation$p_adj)]


# Gene set Demethylation

# optimal p_adj
padj_cP_GSEA_Demethylation_truephen_optim <-
  optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype$documentation$p_adj[length(optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype$documentation$p_adj)]
# default p_adj
padj_cP_GSEA_Demethylation_truephen_default <-
  optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype$documentation$p_adj[1]


################################################################################
### GSEA #######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# default p_adj
padj_GSEA_tCell_truephen_default <- 0.8813
# optimal p_adj
padj_GSEA_tCell_truephen_optim <- 0.3973


# Gene set Demethylation

# default p_adj
padj_GSEA_Demethylation_truephen_default <- 0.1315
# optimal p_adj
padj_GSEA_Demethylation_truephen_optim <- 0.001406



################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_cP_GSEA_tCell_truephen_default <-
  optimP_cP_GSEA_tcell_Pickrell_originalphenotype$documentation$p_adj[1]
# optimal p_adj
padj_cP_GSEA_tCell_truephen_optim <-
  optimP_cP_GSEA_tcell_Pickrell_originalphenotype$documentation$p_adj[length(optimP_cP_GSEA_tcell_Pickrell_originalphenotype$documentation$p_adj)]


# Gene set Demethylation

# optimal p_adj
padj_cP_GSEA_Demethylation_truephen_optim <-
  optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype$documentation$p_adj[length(optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype$documentation$p_adj)]
# default p_adj
padj_cP_GSEA_Demethylation_truephen_default <-
  optimP_cP_GSEA_Demethylation_Pickrell_originalphenotype$documentation$p_adj[1]


################################################################################
### GSEAPreranked ##############################################################
################################################################################

# Gene set t Cell mediated immunity

# default p_adj
padj_GSEAPreranked_tCell_truephen_default <- 0.2928
# optimal p_adj
padj_GSEAPreranked_tCell_truephen_optim <- 0.1993


# Gene set Demethylation

# default p_adj
padj_GSEAPreranked_Demethylation_truephen_default <- 1
# optimal p_adj
padj_GSEAPreranked_Demethylation_truephen_optim <- 0.5797




################################################################################
### Illustration plots #########################################################
################################################################################

# In the following section, we plot all default vs. adjusted p-values for all
# optimization processes performed in the Pickrell data set and the
# TRUE sample labels

# We generate ONE plot which contains the default vs. adjusted p-values

# - for all GSA tools (we stratify the plot by the GSA tools)
# - for all gene sets considered (indicated by color)


# At first, we prepare two data sets separately (each containing the adjusted
# p-values for ONE of the two gene sets for each tools, respectively)
# We then combine both data sets later on and generate one plot



# (i) Random permutations of the random permutations of the sample conditions, first data set
padj_GeneSet1_truephen <- data.frame(cP_ORA = c(padj_cP_ORA_tCell_truephen_default, padj_cP_ORA_tCell_truephen_optim),
                                     GOSeq = c(padj_GOSeq_tCell_truephen_default, padj_GOSeq_tCell_truephen_optim),
                                     DAVID = c(padj_DAVID_tCell_truephen_default, padj_DAVID_tCell_truephen_optim),
                                     PADOG = c(padj_PADOG_PrimImmun_truephen_default, padj_PADOG_PrimImmun_truephen_optim),
                                     cP_GSEA = c(padj_cP_GSEA_tCell_truephen_default, padj_cP_GSEA_tCell_truephen_optim),
                                     GSEA = c(padj_GSEA_tCell_truephen_default,padj_GSEA_tCell_truephen_optim),
                                     GSEAPreranked = c(padj_GSEAPreranked_tCell_truephen_default, padj_GSEAPreranked_tCell_truephen_optim),
                                     state = c("Default","Minimum")) # state: default vs. optimum



# (ii) Random permutations of the random permutations of the sample conditions, second data set
padj_GeneSet2_truephen <- data.frame(cP_ORA = c(padj_cP_ORA_Demethylation_truephen_default, padj_cP_ORA_Demethylation_truephen_optim),
                                     DAVID = c(padj_DAVID_Demethylation_truephen_default, padj_DAVID_Demethylation_truephen_optim),
                                     GOSeq = c(padj_GOSeq_Demethylation_truephen_default, padj_GOSeq_Demethylation_truephen_optim),
                                     PADOG = c(padj_PADOG_GraftvsHost_truephen_default, padj_PADOG_GraftvsHost_truephen_optim),
                                     GSEA = c(padj_GSEA_Demethylation_truephen_default,padj_GSEA_Demethylation_truephen_optim),
                                     cP_GSEA = c(padj_cP_GSEA_Demethylation_truephen_default, padj_cP_GSEA_Demethylation_truephen_optim),
                                     GSEAPreranked = c(padj_GSEAPreranked_Demethylation_truephen_default, padj_GSEAPreranked_Demethylation_truephen_optim),
                                     state = c("Default","Minimum")) # state: default vs. optimum




################################################################################
################################################################################
# Load optimization results for the PERMUTED SAMPLE LABELS #####################
################################################################################
################################################################################

#Load GOSeq results

# Demethylation
load("./Results/optimP_GOSeq_Demethylation_Pickrell_PhenotypePermutations.RData")

# t Cell mediated immunity
load("./Results/optimP_GOSeq_tCell_Pickrell_PhenotypePermutations.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/optimP_ORA_Demethylation_Pickrell_PhenotypePermutations.RData")

# t Cell mediated immunity
load("./Results/optimP_ORA_tCell_Pickrell_PhenotypePermutations.RData")


# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/optimP_PADOG_PrimaryImmunodeficiency_Pickrell_PhenotypePermutations.RData")

# (ii) Graft vs Host
load("./Results/optimP_PADOG_GraftvsHost_Pickrell_PhenotypePermutations.RData")

# GSEA & GSEAPreranked results
# Since these are both web-based tools, the optimization of the result was performed manually and outside of R.
# -> We therefore have to preprocess the information on the increase of the differentially enriched gene sets manually
# -> see below

# load clusterProfiler's GSEA results

# (i) t Cell mediated immunity
load("./Results/optimP_cP_GSEA_tCell_Pickrell_PhenotypePermutations.RData")

# (ii) Demethylation
load("./Results/optimP_cP_GSEA_Demethylation_Pickrell_PhenotypePermutations.RData")


################################################################################
### Prepare results ############################################################
################################################################################

# From each stepwise optimization process, we extract the default and the optimal
# results

################################################################################
### GOSeq ######################################################################
################################################################################



# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_GOSeq_tCell_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_GOSeq_tCell_phenpermutation_optim <- c()


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
padj_GOSeq_Demethylation_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_GOSeq_Demethylation_phenpermutation_optim <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  padj_GOSeq_tCell_phenpermutation_default[i] <-
    optimP_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_GOSeq_tCell_phenpermutation_optim[i] <-
    optimP_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  padj_GOSeq_Demethylation_phenpermutation_default[i] <-
    optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_GOSeq_Demethylation_phenpermutation_optim[i] <-
    optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]

}


################################################################################
### DAVID ######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation


# Cellular Process

# store default p_adj for the 10 permutations of the true sample labels
padj_DAVID_tCell_phenpermutation_default <- rep(1, times = 10)

# store optimal p_adj for the 10 permutations of the true sample labels
padj_DAVID_tCell_phenpermutation_optim <- rep(1, times = 10)



# Metabolic Process

# store default p_adj for the 10 permutations of the true sample labels

padj_DAVID_Demethylation_phenpermutation_default <- rep(1, times = 10)


# store optimal p_adj for the 10 permutations of the true sample labels
padj_DAVID_Demethylation_phenpermutation_optim <- rep(1, times = 10)


################################################################################
### clusterProfiler's ORA ######################################################
################################################################################

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_tCell_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_tCell_phenpermutation_optim <- c()


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_Demethylation_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_cP_ORA_Demethylation_phenpermutation_optim <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  padj_cP_ORA_tCell_phenpermutation_default[i] <-
    optimP_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_ORA_tCell_phenpermutation_optim[i] <-
    optimP_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  padj_cP_ORA_Demethylation_phenpermutation_default[i] <-
    optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_ORA_Demethylation_phenpermutation_optim[i] <-
    optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]

}

################################################################################
### PADOG ######################################################################
################################################################################

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_PADOG_PrimImmun_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_PADOG_PrimImmun_phenpermutation_optim <- c()


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
padj_PADOG_GraftvsHost_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels
padj_PADOG_GraftvsHost_phenpermutation_optim <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  padj_PADOG_PrimImmun_phenpermutation_default[i] <-
    optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_PADOG_PrimImmun_phenpermutation_optim[i] <-
    optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  padj_PADOG_GraftvsHost_phenpermutation_default[i] <-
    optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_PADOG_GraftvsHost_phenpermutation_optim[i] <-
    optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]

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
  padj_cP_GSEA_tCell_phenpermutation_default[i] <-
    optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_GSEA_tCell_phenpermutation_optim[i] <-
    optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  padj_cP_GSEA_Demethylation_phenpermutation_default[i] <-
    optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets
  padj_cP_GSEA_Demethylation_phenpermutation_optim[i] <-
    optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)]

}


################################################################################
### (iv) GSEA ##################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_GSEA_tCell_phenpermutation_default <- c(0.8126,
                                             0.9996,
                                             0.9844,
                                             0.9548,
                                             1,
                                             0.7248,
                                             0.9543,
                                             0.7739,
                                             0.9857,
                                             1)

# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEA_tCell_phenpermutation_optim <- c(0.3885,
                                           0.8254,
                                           0.7713,
                                           0.65,
                                           0.75,
                                           0.657,
                                           0.874,
                                           0.4383,
                                           0.7772,
                                           1)



# Demethylation

# store default p_adj for the 10 permutations of the true sample labels

padj_GSEA_Demethylation_phenpermutation_default <- c(0.9893,
                                                     0.9359,
                                                     0.961,
                                                     1,
                                                     1,
                                                     0.8138,
                                                     0.512,
                                                     1,
                                                     0.7891,
                                                     1)


# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEA_Demethylation_phenpermutation_optim <- c(0.3936,
                                                   0.8779,
                                                   0.6908,
                                                   0.6284,
                                                   0.9679,
                                                   0.5284,
                                                   0.4938,
                                                   0.9203,
                                                   0.6664,
                                                   0.9977)


################################################################################
### (vi) GSEAPreranked #########################################################
################################################################################

# Gene set t Cell mediated immunity

# store default p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_tCell_phenpermutation_default <- c(0.7878,
                                                      0.402,
                                                      0.7536,
                                                      0.6344,
                                                      0.9673,
                                                      0.8538,
                                                      0.6656,
                                                      0.9510,
                                                      0.9398,
                                                      0.8367)
# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_tCell_phenpermutation_optim <- c(0.595,
                                                    0.0677,
                                                    0.3499,
                                                    0.3491,
                                                    0.3898,
                                                    0.8019,
                                                    0.3155,
                                                    0.9510,
                                                    0.4847,
                                                    0.7471)


# Demethylation
# store default p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_Demethylation_phenpermutation_default <- c(0.9516,
                                                              1,
                                                              0.8564,
                                                              0.7073,
                                                              1,
                                                              0.9847,
                                                              0.7614,
                                                              0.99,
                                                              0.9928,
                                                              1)
# store optimal p_adj for the 10 permutations of the true sample labels
padj_GSEAPreranked_Demethylation_phenpermutation_optim <- c(0.7741,
                                                            0.7978,
                                                            0.5727,
                                                            0.2179,
                                                            0.8276,
                                                            0.9058,
                                                            0.4469,
                                                            0.9145,
                                                            0.9904,
                                                            0.1863)

################################################################################
### Illustration plots #########################################################
################################################################################

# In the following section, we plot all default vs. adjusted p-values for all
# optimization processes performed in the Pickrell data set and the PERMUTED
# sample labels

# We generate ONE plot which contains the default vs. adjusted p-values

# - for all GSA tools (we stratify the plot by the GSA tools)
# - for all gene sets considered (indicated by color)


# At first, we prepare two data sets separately (each containing the adjusted
# p-values for ONE of the two gene sets for each tools, respectively)
# We then combine both data sets later on and generate one plot



# (i) Random permutations of the random permutations of the sample conditions, first data set
padj_GeneSet1_phenpermutation <- data.frame(cP_ORA = c(padj_cP_ORA_tCell_phenpermutation_default, padj_cP_ORA_tCell_phenpermutation_optim),
                                            GOSeq = c(padj_GOSeq_tCell_phenpermutation_default, padj_GOSeq_tCell_phenpermutation_optim),
                                            DAVID = c(padj_DAVID_tCell_phenpermutation_default, padj_DAVID_tCell_phenpermutation_optim),
                                            PADOG = c(padj_PADOG_PrimImmun_phenpermutation_default, padj_PADOG_PrimImmun_phenpermutation_optim),
                                            cP_GSEA = c(padj_cP_GSEA_tCell_phenpermutation_default, padj_cP_GSEA_tCell_phenpermutation_optim),
                                            GSEA = c(padj_GSEA_tCell_phenpermutation_default,padj_GSEA_tCell_phenpermutation_optim),
                                            GSEAPreranked = c(padj_GSEAPreranked_tCell_phenpermutation_default, padj_GSEAPreranked_tCell_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum
                                            ID = rep(c(1:10),2)) # add the number of the permutation



# (ii) Random permutations of the random permutations of the sample conditions, second data set
padj_GeneSet2_phenpermutation <- data.frame(cP_ORA = c(padj_cP_ORA_Demethylation_phenpermutation_default, padj_cP_ORA_Demethylation_phenpermutation_optim),
                                            GOSeq = c(padj_GOSeq_Demethylation_phenpermutation_default, padj_GOSeq_Demethylation_phenpermutation_optim),
                                            DAVID = c(padj_DAVID_Demethylation_phenpermutation_default, padj_DAVID_Demethylation_phenpermutation_optim),
                                            PADOG = c(padj_PADOG_GraftvsHost_phenpermutation_default, padj_PADOG_GraftvsHost_phenpermutation_optim),
                                            GSEA = c(padj_GSEA_Demethylation_phenpermutation_default,padj_GSEA_Demethylation_phenpermutation_optim),
                                            cP_GSEA = c(padj_cP_GSEA_Demethylation_phenpermutation_default, padj_cP_GSEA_Demethylation_phenpermutation_optim),
                                            GSEAPreranked = c(padj_GSEAPreranked_Demethylation_phenpermutation_default, padj_GSEAPreranked_Demethylation_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum
                                            ID = rep(c(1:10),2)) # add the number of the permutation




################################################################################
### Generate ggplots ###########################################################
################################################################################

plot_truelabels <- create_results_illustration_pvalue_rank(padj_GeneSet1_truephen, padj_GeneSet2_truephen,
                                                           "p_adj", "true_labels")

plot_permutedlabels <- create_results_illustration_pvalue_rank(padj_GeneSet1_phenpermutation, padj_GeneSet2_phenpermutation,
                                                               "p_adj", "random_permutations")

plot_padj <- plot_grid(plot_truelabels, plot_permutedlabels, labels=c("A", "B"), ncol = 1, nrow = 2)

## uncomment to save
# ggsave2("./Results_illustrations/Figure4.pdf",
#         width = 10,
#         height = 12)
#


