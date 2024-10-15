################################################################################
### Results Illustrations for the optimization of the p-values for Pickrell data
################################################################################

library(ggplot2)
library(tidyr)
library(cowplot)

################################################################################
### load help functions to generate results illustrations ######################
################################################################################

source("./R/Code_figures/helpfunctions_results_illustrations_pvalue_rank.R")

################################################################################
################################################################################
# Load optimization results for the TRUE SAMPLE LABELS #########################
################################################################################
################################################################################

#Load GOSeq results

# Demethylation
load("./Results/Intermediate_results/optimRank_GOSeq_Demethylation_Pickrell_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/Intermediate_results/optimRank_GOSeq_tCell_Pickrell_OriginalPhenotype.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/Intermediate_results/optimRank_ORA_Demethylation_Pickrell_OriginalPhenotype.RData")

# t Cell mediated immunity
load("./Results/Intermediate_results/optimRank_ORA_tCell_Pickrell_OriginalPhenotype.RData")


# Load DAVID results : Achtung DAVID fehlt bis jetzt noch

# Load clusterProfiler's GSEA results

# Demethylation
load("./Results/Intermediate_results/optimRank_cP_GSEA_Demethylation_OriginalPhenotype.RData")

# t-Cell mediated immunity
load("./Results/Intermediate_results/optimRank_cP_GSEA_tCell_Pickrell_OriginalPhenotype.RData")


# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/Intermediate_results/optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_OriginalPhenotype.RData")

# (ii) Graft vs Host
load("./Results/Intermediate_results/optimRank_PADOG_GraftvsHost_Pickrell_OriginalPhenotype.RData")

# GSEA & GSEAPreranked results
# Since these are both web-based tools, the optimization of the result was performed manually and outside of R.
# -> We therefore have to preprocess the information on the increase of the differentially enriched gene sets manually
# -> see below

################################################################################
### Prepare results ############################################################
################################################################################

# From each stepwise optimization process, we extract the default and the optimal
# results

################################################################################
### GOSeq ######################################################################
################################################################################

# Gene set t Cell mediated immunity

# default rel_rank
rel_rank_GOSeq_tCell_truephen_default <-
  optimRank_GOSeq_tcell_Pickrell_originalphenotype$documentation$rank[1]
# optimal rel_rank
rel_rank_GOSeq_tCell_truephen_optim <-
  optimRank_GOSeq_tcell_Pickrell_originalphenotype$documentation$rank[length(optimRank_GOSeq_tcell_Pickrell_originalphenotype$documentation$rank)]


# Gene set Demethylation

# default rel_rank
rel_rank_GOSeq_Demethylation_truephen_default <-
  optimRank_GOSeq_Demethylation_Pickrell_originalphenotype$documentation$rank[1]
# optimal rel_rank
rel_rank_GOSeq_Demethylation_truephen_optim <-
  optimRank_GOSeq_Demethylation_Pickrell_originalphenotype$documentation$rank[length(optimRank_GOSeq_Demethylation_Pickrell_originalphenotype$documentation$rank)]


################################################################################
### DAVID ######################################################################
################################################################################

# extract default and optimal adjusted p-value from the manual documentation

# Gene set t Cell mediated immunity

# default rel_rank
rel_rank_DAVID_tCell_truephen_default  <- 1
rel_rank_DAVID_tCell_truephen_optim  <- 1


# Gene set Demethylation

# default rel_rank
rel_rank_DAVID_Demethylation_truephen_default  <- 1
# optimal rel_rank
rel_rank_DAVID_Demethylation_truephen_optim  <- 1


################################################################################
### clusterProfiler's ORA ######################################################
################################################################################

# Gene set t Cell mediated immunity

# default rel_rank
rel_rank_cP_ORA_tCell_truephen_default <-
  optimRank_cP_ORA_tcell_Pickrell_originalphenotype$documentation$rank[1]
# optimal rel_rank
rel_rank_cP_ORA_tCell_truephen_optim <-
  optimRank_cP_ORA_tcell_Pickrell_originalphenotype$documentation$rank[length(optimRank_cP_ORA_tcell_Pickrell_originalphenotype$documentation$rank)]


# Gene set Demethylation

# optimal rel_rank
rel_rank_cP_ORA_Demethylation_truephen_optim <-
  optimRank_cP_ORA_Demethylation_Pickrell_originalphenotype$documentation$rank[length(optimRank_cP_ORA_Demethylation_Pickrell_originalphenotype$documentation$rank)]
# default rel_rank
rel_rank_cP_ORA_Demethylation_truephen_default <-
  optimRank_cP_ORA_Demethylation_Pickrell_originalphenotype$documentation$rank[1]



################################################################################
### PADOG ######################################################################
################################################################################

# Gene set Primary Immunodeficiency

# default rel_rank
rel_rank_PADOG_PrimImmun_truephen_default <-
  optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype$documentation$rank[1]
# optimal rel_rank
rel_rank_PADOG_PrimImmun_truephen_optim <-
  optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype$documentation$rank[length(optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_originalphenotype$documentation$rank)]


# Gene set Graft vs Host

# default rel_rank
rel_rank_PADOG_GraftvsHost_truephen_default <-
  optimRank_PADOG_GraftvsHost_Pickrell_originalphenotype$documentation$rank[1]
# optimal rel_rank
rel_rank_PADOG_GraftvsHost_truephen_optim <-
  optimRank_PADOG_GraftvsHost_Pickrell_originalphenotype$documentation$rank[length(optimRank_PADOG_GraftvsHost_Pickrell_originalphenotype$documentation$rank)]


################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################


# Gene set t Cell mediated immunity

# default rel_rank
rel_rank_cP_GSEA_tCell_truephen_default <-
  optimRank_tcell_Pickrell_originalphenotype$documentation$rank[1]
# optimal rel_rank
rel_rank_cP_GSEA_tCell_truephen_optim <-
  optimRank_tcell_Pickrell_originalphenotype$documentation$rank[length(optimRank_tcell_Pickrell_originalphenotype$documentation$rank)]


# Gene set Demethylation

# optimal rel_rank
rel_rank_cP_GSEA_Demethylation_truephen_optim <-
  optimRank_Demethylation_Pickrell_originalphenotype$documentation$rank[length(optimRank_Demethylation_Pickrell_originalphenotype$documentation$rank)]
# default rel_rank
rel_rank_cP_GSEA_Demethylation_truephen_default <-
  optimRank_Demethylation_Pickrell_originalphenotype$documentation$rank[1]



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
rel_rank_GeneSet1_truephen  <- data.frame(cP_ORA = c(rel_rank_cP_ORA_tCell_truephen_default, rel_rank_cP_ORA_tCell_truephen_optim),
                                         GOSeq = c(rel_rank_GOSeq_tCell_truephen_default, rel_rank_GOSeq_tCell_truephen_optim),
                                         PADOG = c(rel_rank_PADOG_PrimImmun_truephen_default, rel_rank_PADOG_PrimImmun_truephen_optim),
                                         cP_GSEA = c(rel_rank_cP_GSEA_tCell_truephen_default, rel_rank_cP_GSEA_tCell_truephen_optim),
                                         state = c("Default","Minimum")) # state: default vs. optimum





# (ii) Random permutations of the random permutations of the sample conditions, second data set
rel_rank_GeneSet2_truephen  <- data.frame(cP_ORA = c(rel_rank_cP_ORA_Demethylation_truephen_default, rel_rank_cP_ORA_Demethylation_truephen_optim),
                                         GOSeq = c(rel_rank_GOSeq_Demethylation_truephen_default, rel_rank_GOSeq_Demethylation_truephen_optim),
                                         PADOG = c(rel_rank_PADOG_GraftvsHost_truephen_default, rel_rank_PADOG_GraftvsHost_truephen_optim),
                                         cP_GSEA = c(rel_rank_cP_GSEA_Demethylation_truephen_default, rel_rank_cP_GSEA_Demethylation_truephen_optim),
                                         state = c("Default","Minimum")) # state: default vs. optimum




################################################################################
################################################################################
# Load optimization results for the PERMUTED SAMPLE LABELS #####################
################################################################################
################################################################################

#Load GOSeq results

# Demethylation
load("./Results/Intermediate_results/optimRank_GOSeq_Demethylation_Pickrell_PhenotypePermutations.RData")

# t Cell mediated immunity
load("./Results/Intermediate_results/optimRank_GOSeq_tCell_Pickrell_PhenotypePermutations.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/Intermediate_results/optimRank_ORA_Demethylation_Pickrell_PhenotypePermutations.RData")

# t Cell mediated immunity
load("./Results/Intermediate_results/optimRank_ORA_tCell_Pickrell_PhenotypePermutations.RData")


# Load PADOG results

# (i) Primary Immunodeficiency
load("./Results/Intermediate_results/optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_PhenotypePermutations.RData")

# (ii) Graft vs Host
load("./Results/Intermediate_results/optimRank_PADOG_GraftvsHost_Pickrell_PhenotypePermutations.RData")

# GSEA & GSEAPreranked results
# Since these are both web-based tools, the optimization of the result was performed manually and outside of R.
# -> We therefore have to preprocess the information on the increase of the differentially enriched gene sets manually
# -> see below

# load clusterProfiler's GSEA results

# (i) t Cell mediated immunity
load("./Results/Intermediate_results/optimRank_cP_GSEA_tCell_Pickrell_PhenotypePermutations.RData")

# (ii) Demethylation
load("./Results/Intermediate_results/optimRank_cP_GSEA_Demethylation_Pickrell_PhenotypePermutations.RData")


################################################################################
### Prepare results ############################################################
################################################################################

# From each stepwise optimization process, we extract the default and the optimal
# results

################################################################################
### GOSeq ######################################################################
################################################################################



# Gene set t Cell mediated immunity

# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_GOSeq_tCell_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_GOSeq_tCell_phenpermutation_optim  <- c()


# Demethylation
# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_GOSeq_Demethylation_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_GOSeq_Demethylation_phenpermutation_optim  <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  rel_rank_GOSeq_tCell_phenpermutation_default[i] <-
    optimRank_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_GOSeq_tCell_phenpermutation_optim[i] <-
    optimRank_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  rel_rank_GOSeq_Demethylation_phenpermutation_default[i] <-
    optimRank_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_GOSeq_Demethylation_phenpermutation_optim[i] <-
    optimRank_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank)]

}



################################################################################
### clusterProfiler's ORA ######################################################
################################################################################

# Gene set t Cell mediated immunity

# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_ORA_tCell_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_ORA_tCell_phenpermutation_optim  <- c()


# Demethylation
# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_ORA_Demethylation_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_ORA_Demethylation_phenpermutation_optim  <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  rel_rank_cP_ORA_tCell_phenpermutation_default[i] <-
    optimRank_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_cP_ORA_tCell_phenpermutation_optim[i] <-
    optimRank_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  rel_rank_cP_ORA_Demethylation_phenpermutation_default[i] <-
    optimRank_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_cP_ORA_Demethylation_phenpermutation_optim[i] <-
    optimRank_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank)]

}

################################################################################
### PADOG ######################################################################
################################################################################

# Gene set t Cell mediated immunity

# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_PADOG_PrimImmun_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_PADOG_PrimImmun_phenpermutation_optim  <- c()


# Demethylation
# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_PADOG_GraftvsHost_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_PADOG_GraftvsHost_phenpermutation_optim  <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  rel_rank_PADOG_PrimImmun_phenpermutation_default[i] <-
    optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_PADOG_PrimImmun_phenpermutation_optim[i] <-
    optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$rank)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  rel_rank_PADOG_GraftvsHost_phenpermutation_default[i] <-
    optimRank_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_PADOG_GraftvsHost_phenpermutation_optim[i] <-
    optimRank_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$rank)]

}

################################################################################
### clusterProfiler's GSEA #####################################################
################################################################################

# Gene set t Cell mediated immunity

# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_GSEA_tCell_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_GSEA_tCell_phenpermutation_optim  <- c()


# Demethylation
# store default rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_GSEA_Demethylation_phenpermutation_default  <- c()
# store optimal rel_rank for the 10 permutations of the true sample labels
rel_rank_cP_GSEA_Demethylation_phenpermutation_optim  <- c()


for(i in 1:10){

  # t Cell mediated immunity

  # store the default number of differentially enriched gene sets
  rel_rank_cP_GSEA_tCell_phenpermutation_default[i] <-
    optimRank_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_cP_GSEA_tCell_phenpermutation_optim[i] <-
    optimRank_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_tcell_Pickrell_phenotypepermutations[[i]]$documentation$rank)]


  # Demethylation

  # store the default number of differentially enriched gene sets
  rel_rank_cP_GSEA_Demethylation_phenpermutation_default[i] <-
    optimRank_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets
  rel_rank_cP_GSEA_Demethylation_phenpermutation_optim[i] <-
    optimRank_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank[length(optimRank_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$rank)]

}




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
rel_rank_GeneSet1_phenpermutation  <- data.frame(cP_ORA = c(rel_rank_cP_ORA_tCell_phenpermutation_default, rel_rank_cP_ORA_tCell_phenpermutation_optim),
                                                GOSeq = c(rel_rank_GOSeq_tCell_phenpermutation_default, rel_rank_GOSeq_tCell_phenpermutation_optim),
                                                PADOG = c(rel_rank_PADOG_PrimImmun_phenpermutation_default, rel_rank_PADOG_PrimImmun_phenpermutation_optim),
                                                cP_GSEA = c(rel_rank_cP_GSEA_tCell_phenpermutation_default, rel_rank_cP_GSEA_tCell_phenpermutation_optim),
                                                state = c(rep("Default",10), rep("Minimum",10))) # state: default vs. optimum




# (ii) Random permutations of the random permutations of the sample conditions, second data set
rel_rank_GeneSet2_phenpermutation  <- data.frame(cP_ORA = c(rel_rank_cP_ORA_Demethylation_phenpermutation_default, rel_rank_cP_ORA_Demethylation_phenpermutation_optim),
                                                GOSeq = c(rel_rank_GOSeq_Demethylation_phenpermutation_default, rel_rank_GOSeq_Demethylation_phenpermutation_optim),
                                                PADOG = c(rel_rank_PADOG_GraftvsHost_phenpermutation_default, rel_rank_PADOG_GraftvsHost_phenpermutation_optim),
                                                cP_GSEA = c(rel_rank_cP_GSEA_Demethylation_phenpermutation_default, rel_rank_cP_GSEA_Demethylation_phenpermutation_optim),
                                                state = c(rep("Default",10), rep("Minimum",10))) # state: default vs. optimum


################################################################################
### Generate ggplots ###########################################################
################################################################################

plot_truelabels  <- create_results_illustration_pvalue_rank(rel_rank_GeneSet1_truephen, rel_rank_GeneSet2_truephen,
                                                           "rel_rank", "true_labels")

plot_permutedlabels  <- create_results_illustration_pvalue_rank(rel_rank_GeneSet1_phenpermutation, rel_rank_GeneSet2_phenpermutation,
                                                               "rel_rank", "random_permutations")

plot_grid(plot_permutedlabels, plot_truelabels, labels=c("A", "B"), ncol = 1, nrow = 2)

## uncomment to save
ggsave("./Results/Figures/Figure5.eps",
       width = 10,
       height = 12)

