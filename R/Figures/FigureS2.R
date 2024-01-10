#################################################################################
### Results illustrations: number of differentially enriched gene sets (n_DEGS) #
#################################################################################

library(ggplot2)
library(tidyr)
library(cowplot)

################################################################################
### load help functions to generate results illustrations ######################
################################################################################

source("./R/Figures/helpfunctions_results_illustrations_n_DEGS.R")




################################################################################
#### Illustration of difference in n_DEGS across all tools #####################
################################################################################

#Load GOSeq results

load("./Results/Intermediate_results/GOSeq_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/Intermediate_results/GOSeq_Results_Bottomly_PhenotypePermutations.RData")

# Load clusterProfiler's ORA results

load("./Results/Intermediate_results/ORA_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/Intermediate_results/ORA_Results_Bottomly_PhenotypePermutations.RData")

# Load DAVID results:
# Since DAVID is a web-based application, the results are extracted from the
# corresponding R-script manually

# Load clusterProfiler's GSEA results

load("./Results/Intermediate_results/cP_GSEA_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/Intermediate_results/cP_GSEA_Results_Bottomly_PhenotypePermutations.RData")


# Load PADOG results


load("./Results/Intermediate_results/PADOG_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/Intermediate_results/PADOG_Results_Bottomly_PhenotypePermutations.RData")

# GSEA & GSEAPreranked results
# Since these are both web-based tools, the optimization of the result was performed manually and outside of R.
# -> We therefore have to preprocess the information on the increase of the differentially enriched gene sets manually
# -> see below


########
# GOSeq
########

# true sample labels

# optimal n_DEGS
n_DEGS_GOSeq_bottomly_truephen_optim <-
  optim_GOSeq_results_Bottomly_originalphenotype$documentation$n_DEGS[length(optim_GOSeq_results_Bottomly_originalphenotype$documentation$n_DEGS)]
# default n_DEGS
n_DEGS_GOSeq_bottomly_truephen_default <-
  optim_GOSeq_results_Bottomly_originalphenotype$documentation$n_DEGS[1]


# permuted sample labels

# store optimal n_DEGS for the 10 random permutations
n_DEGS_GOSeq_bottomly_phenpermutation_optim <- c()
# store default n_DEGS for the 10 random permutations
n_DEGS_GOSeq_bottomly_phenpermutation_default <- c()


for(i in 1:10){

  n_DEGS_GOSeq_bottomly_phenpermutation_optim[i] <-
    optim_GOSeq_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[length(optim_GOSeq_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS)]

  n_DEGS_GOSeq_bottomly_phenpermutation_default[i] <-
    optim_GOSeq_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[1]


}

###########
# DAVID ###
###########


# true conditions

# optimal n_DEGS
n_DEGS_DAVID_bottomly_truephen_default <- 2
# default n_DEGS
n_DEGS_DAVID_bottomly_truephen_optim <- 2

# random permutations of conditions

# optimal n_DEGS
n_DEGS_DAVID_bottomly_phenpermutation_default <- rep(0, times = 10)
# true n_DEGS
n_DEGS_DAVID_bottomly_phenpermutation_optim <- rep(0, times = 10)





#######################
# clusterProfiler's ORA
#######################

# true sample labels

# optimal n_DEGS
n_DEGS_cP_ORA_bottomly_truephen_optim <-
  optim_ORA_results_Bottomly_originalphenotype$documentation$n_DEGS[length(optim_ORA_results_Bottomly_originalphenotype$documentation$n_DEGS)]
# default n_DEGS
n_DEGS_cP_ORA_bottomly_truephen_default  <-
  optim_ORA_results_Bottomly_originalphenotype$documentation$n_DEGS[1]


# permuted sample labels

# optimal n_DEGS
n_DEGS_cP_ORA_bottomly_phenpermutation_optim <- c()
# default n_DEGS
n_DEGS_cP_ORA_bottomly_phenpermutation_default <- c()

for(i in 1:10){

  n_DEGS_cP_ORA_bottomly_phenpermutation_optim[i] <-
    optim_ORA_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[length(optim_ORA_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS)]

  n_DEGS_cP_ORA_bottomly_phenpermutation_default[i] <-
    optim_ORA_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[1]


}




########################
# clusterProfiler's GSEA
########################

# true sample labels

# optimal n_DEGS
n_DEGS_cP_GSEA_bottomly_truephen_optim <-
  optim_cP_GSEA_results_Bottomly_originalphenotype$documentation$n_DEGS[length(optim_cP_GSEA_results_Bottomly_originalphenotype$documentation$n_DEGS)]
# default n_DEGS
n_DEGS_cP_GSEA_bottomly_truephen_default <-
  optim_cP_GSEA_results_Bottomly_originalphenotype$documentation$n_DEGS[1]

# permuted sample labels

# optimal n_DEGS
n_DEGS_cP_GSEA_bottomly_phenpermutation_optim <- c()
# default n_DEGS
n_DEGS_cP_GSEA_bottomly_phenpermutation_default <- c()

for(i in 1:10){

  n_DEGS_cP_GSEA_bottomly_phenpermutation_optim[i] <-
    optim_cP_GSEA_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[length(optim_cP_GSEA_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS)]

  n_DEGS_cP_GSEA_bottomly_phenpermutation_default[i] <-
    optim_cP_GSEA_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[1]

}


#######
# PADOG
#######

# true sample labels

# optimal n_DEGS
n_DEGS_PADOG_bottomly_truephen_optim <-
  optim_PADOG_results_Bottomly_originalphenotype$documentation$n_DEGS[length(optim_PADOG_results_Bottomly_originalphenotype$documentation$n_DEGS)]
# default n_DEGS
n_DEGS_PADOG_bottomly_truephen_default <-
  optim_PADOG_results_Bottomly_originalphenotype$documentation$n_DEGS[1]

# permuted sample labels

# optimal n_DEGS
n_DEGS_PADOG_bottomly_phenpermutation_optim <- c()
# default n_DEGS
n_DEGS_PADOG_bottomly_phenpermutation_default <- c()

for(i in 1:10){

  n_DEGS_PADOG_bottomly_phenpermutation_optim[i] <-
    optim_PADOG_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[length(optim_PADOG_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS)]

  n_DEGS_PADOG_bottomly_phenpermutation_default[i] <-
    optim_PADOG_results_Bottomly_phenotypepermutation[[i]]$documentation$n_DEGS[1]


}



#################
# GSEA (web tool)
#################

# true conditions

# optimal n_DEGS
n_DEGS_GSEA_bottomly_truephen_default <- 2
# default n_DEGS
n_DEGS_GSEA_bottomly_truephen_optim <- 265

# random permutations of conditions

# optimal n_DEGS
n_DEGS_GSEA_bottomly_phenpermutation_default <- c(3,1,0,1,0,1,0,10,7,0)
# true n_DEGS
n_DEGS_GSEA_bottomly_phenpermutation_optim <- c(4,10,0,13,1,8,5,42,17,0)



###############
# GSEAPreranked
###############


# optimal n_DEGS
n_DEGS_GSEAPreranked_bottomly_truephen_default <- 0
# default n_DEGS
n_DEGS_GSEAPreranked_bottomly_truephen_optim <- 40

# random permutations of conditions

# optimal n_DEGS
n_DEGS_GSEAPreranked_bottomly_phenpermutation_default <- c(65, 0, 18, 0, 280, 124, 131, 8, 262, 724)
# default n_DEGS
n_DEGS_GSEAPreranked_bottomly_phenpermutation_optim <- c(225, 30, 168, 35, 381, 124, 499, 316, 554, 962)


################################################################################
### Summary ####################################################################
################################################################################


# for each GSA tool, we display the default and optimal number of differentially enriched gene sets in the following
# data frame

# (i) True sample conditions
dat_overview_n_DEGS_bottomly_truephen <- data.frame(cP_ORA = c(n_DEGS_cP_ORA_bottomly_truephen_default, n_DEGS_cP_ORA_bottomly_truephen_optim),
                                                    GOSeq = c(n_DEGS_GOSeq_bottomly_truephen_default, n_DEGS_GOSeq_bottomly_truephen_optim),
                                                    DAVID = c(n_DEGS_DAVID_bottomly_truephen_default, n_DEGS_DAVID_bottomly_truephen_optim),
                                                    PADOG = c(n_DEGS_PADOG_bottomly_truephen_default, n_DEGS_PADOG_bottomly_truephen_optim),
                                                    GSEA = c(n_DEGS_GSEA_bottomly_truephen_default, n_DEGS_GSEA_bottomly_truephen_optim),
                                                    GSEAPreranked = c(n_DEGS_GSEAPreranked_bottomly_truephen_default, n_DEGS_GSEAPreranked_bottomly_truephen_optim),
                                                    cP_GSEA = c(n_DEGS_cP_GSEA_bottomly_truephen_default, n_DEGS_cP_GSEA_bottomly_truephen_optim),
                                                    state = c("Default", "Maximum")) # declare default vs. optimum




# (ii) Random permutations of the true sample conditions (bottomly data set)
dat_overview_n_DEGS_bottomly_phenpermutation <- data.frame(cP_ORA = c(n_DEGS_cP_ORA_bottomly_phenpermutation_default, n_DEGS_cP_ORA_bottomly_phenpermutation_optim),
                                                          GOSeq = c(n_DEGS_GOSeq_bottomly_phenpermutation_default, n_DEGS_GOSeq_bottomly_phenpermutation_optim),
                                                          DAVID = c(n_DEGS_DAVID_bottomly_phenpermutation_default,n_DEGS_DAVID_bottomly_phenpermutation_optim),
                                                          PADOG = c(n_DEGS_PADOG_bottomly_phenpermutation_default, n_DEGS_PADOG_bottomly_phenpermutation_optim),
                                                          GSEA = c(n_DEGS_GSEA_bottomly_phenpermutation_default,n_DEGS_GSEA_bottomly_phenpermutation_optim),
                                                          GSEAPreranked = c(n_DEGS_GSEAPreranked_bottomly_phenpermutation_default, n_DEGS_GSEAPreranked_bottomly_phenpermutation_optim),
                                                          cP_GSEA = c(n_DEGS_cP_GSEA_bottomly_phenpermutation_default, n_DEGS_cP_GSEA_bottomly_phenpermutation_optim),
                                                          state = c(rep("Default",10), rep("Maximum",10)), # state: default vs. optimum
                                                          ID = rep(c(1:10),2)) # add the number of the permutation




################################################################################
### Generate ggplots ###########################################################
################################################################################

plot_truelabels <- create_results_illustration_n_DEGS(dat_overview_n_DEGS_bottomly_truephen, "true_labels", -8)

plot_permutedlabels <- create_results_illustration_n_DEGS(dat_overview_n_DEGS_bottomly_phenpermutation,
                                                          "random_permutations", -14)

plot_grid(plot_permutedlabels, plot_truelabels, labels=c("A", "B"), ncol = 1, nrow = 2)

## uncomment to save
ggsave("./Results/Figures/FigureS2.eps",
       width = 10,
       height = 12)



