#################################################################################
### Results illustrations: number of differentially enriched gene sets (n_DEGS) #
#################################################################################

library(ggplot2)
library(tidyr)

# set working directory 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/Assessment_OverOptimism")


################################################################################
#### Illustration of difference in n_DEGS across all tools #####################
################################################################################

#Load GOSeq results 

load("./Results/GOSeq_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/GOSeq_Results_Bottomly_PhenotypePermutations.RData")

# Load clusterProfiler's ORA results

load("./Results/ORA_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/ORA_Results_Bottomly_PhenotypePermutations.RData")

# Load DAVID results:
# Since DAVID is a web-based application, the results are extracted from the 
# corresponding R-script manually 

# Load clusterProfiler's GSEA results

load("./Results/cP_GSEA_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/cP_GSEA_Results_Bottomly_PhenotypePermutations.RData")


# Load PADOG results 


load("./Results/PADOG_Results_Bottomly_OriginalPhenotype.RData")
load("./Results/PADOG_Results_Bottomly_PhenotypePermutations.RData")

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
n_DEGS_GSEAPreranked_bottomly_truephen_default <- 1 
# default n_DEGS
n_DEGS_GSEAPreranked_bottomly_truephen_optim <- 40

# random permutations of conditions 

# optimal n_DEGS
n_DEGS_GSEAPreranked_bottomly_phenpermutation_default <- c(65, 1, 18, 1, 280, 124, 131, 8, 262, 724)
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

# transform data frame such that each observed n_DEGS value is in a separate column, additionally labelled by the 
# associated GSA tool and state
dat_overview_n_DEGS_bottomly_truephen_long <- pivot_longer(dat_overview_n_DEGS_bottomly_truephen, 
                                                 cols=!c("state"), 
                                                 names_to="GSA_tool", 
                                                 values_to = "n_DEGS")
# transform GSA tools to factors
# -> this way we can fix the order of the tools in the graphic
dat_overview_n_DEGS_bottomly_truephen_long$GSA_tool <- factor(dat_overview_n_DEGS_bottomly_truephen_long$GSA_tool, 
                                                    levels = c("GOSeq", "DAVID", "cP_ORA","PADOG", "GSEA", "GSEAPreranked", "cP_GSEA"))


# In the ggplot, we eventually want to add the names of the GSA tools below the actual x-axis (line annotate())
# for this, we want to insert a line break in "GSEAPreranked" because the term is too long in one word 

# Replace "GSEAPreranked" by "GSEA- \n Preranked"
add_labels_xaxis <- levels(dat_overview_n_DEGS_bottomly_truephen_long$GSA_tool)
add_labels_xaxis[add_labels_xaxis == "GSEAPreranked"] <- "GSEA- \n Preranked"


# Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"


# Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"


# create ggplot 
ggplot(data = dat_overview_n_DEGS_bottomly_truephen_long, 
       aes(x = interaction(GSA_tool, state, lex.order = TRUE), 
           y = n_DEGS, group = 1)) + 
  geom_line(aes(group=GSA_tool), 
            size=0.7, alpha=0.7, col = "#F8766D") + 
  # add scatter for each default and optimal value 
  geom_point(size = 1.2,
             alpha = 0.7, 
             col = "#F8766D") + 
  # add labels "Default" and "Maximum" for each GSA tool on the x-axis 
  scale_x_discrete(labels= rep(c("Default", "Maximum"), 
                               times = 7)) + 
  # rotate labels on the x-axis, remove x-axis title and add space below the plot 
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 9), 
        axis.title.x = element_text(vjust = -13), # move x-axis label downwards to add space for the individual GSA tools 
        plot.margin = margin(t = 1, b = 2, l = 1, r = 1, unit = "cm")) + ## add space below the actual plot (needed for the GSA tool names)
  # Add the tool names to the plot: 
  annotate(geom = "text", 
           x = 1.5 + 2*(0:6), 
           y = -70, 
           label = add_labels_xaxis, size = 3) + 
  # specify range of y-axis
  coord_cartesian(ylim = c(0, 350), 
                  expand = FALSE, 
                  clip = "off") +
  xlab("GSA tools") +
  ylab("Number of differentially enriched gene sets") 



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


# transform data frame such that each observed n_DEGS value is in a separate column, additionally labelled by the 
# associated GSA tool and state
dat_overview_n_DEGS_bottomly_phenpermutation_long <- pivot_longer(dat_overview_n_DEGS_bottomly_phenpermutation, 
                                                                 cols=!c("state", "ID"), 
                                                                 names_to="GSA_tool", 
                                                                 values_to = "n_DEGS")
# transform GSA tools to factors
# -> this way we can fix the order of the tools in the graphic
dat_overview_n_DEGS_bottomly_phenpermutation_long$GSA_tool <- factor(dat_overview_n_DEGS_bottomly_phenpermutation_long$GSA_tool, 
                                                                    levels = c("cP_ORA", "DAVID", "GOSeq", "PADOG", "GSEA", "GSEAPreranked", "cP_GSEA"))


# add the combination of the GSA tool and the permutation ID 
# -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
dat_overview_n_DEGS_bottomly_phenpermutation_long$unique_ID <- paste0(dat_overview_n_DEGS_bottomly_phenpermutation_long$GSA_tool,
                                                                     "_", 
                                                                     dat_overview_n_DEGS_bottomly_phenpermutation_long$ID)

# In the ggplot, we eventually want to add the names of the GSA tools below the actual x-axis (line annotate())
# for this, we want to insert a line break in "GSEAPreranked" because the term is too long in one word 

# # Replace "GSEAPreranked" by "GSEA- \n Preranked"
# add_labels_xaxis <- levels(dat_overview_n_DEGS_bottomly_phenpermutation_long$GSA_tool)
# add_labels_xaxis[add_labels_xaxis == "GSEAPreranked"] <- "GSEA- \n Preranked"
# 
# ## Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
# add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"
# 
# 
# # Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
# add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"


# ggplot 
ggplot(data =dat_overview_n_DEGS_bottomly_phenpermutation_long, 
       aes(x = interaction(GSA_tool, state, lex.order = TRUE), 
           y = n_DEGS, group = 1)) + 
  geom_line(aes(group=unique_ID), size=0.5, alpha=0.7, col = "#F8766D") + 
  geom_point(size = 1, alpha = 1.2, col = "#F8766D") + 
  scale_x_discrete(labels= rep(c("Default", "Maximum"), 
                               times = 7)) + 
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 9), 
        axis.title.x = element_text(vjust = -2), 
        plot.margin = margin(t=1, b =3, l=1, r=1, unit="cm")) + ## add space below the actual plot (needed for the GSA tool names)
  # Add the tool names to the plot: 
  annotate(geom = "text", 
           x = 1.5 + 2*(0:6), 
           y = -350, 
           label = add_labels_xaxis, size =3 ) + 
  coord_cartesian(ylim = c(0,1000),  clip = "off")+ # clip = "off" required to add GSA tool names below the plot
  xlab("GSA tools") +
  ylab("Number of differentially enriched gene sets") +
  scale_y_continuous(breaks = sort(c(0,50, seq(100, 1000, 100))) )



                     