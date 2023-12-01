################################################################################
### Results Illustrations for the optimization of the p-values for Bottomly data
### and the random permutations of the true sample labels ######################
################################################################################

library(ggplot2)
library(tidyr)

################################################################################
### Load corresponding results for Bottomly data set ###########################
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
  rank_cP_GSEA_tCell_phenpermutation_default[i] <- optimRank_cP_GSEA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets 
  rank_cP_GSEA_tCell_phenpermutation_optim[i] <- optimRank_cP_GSEA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank[length(optimRank_cP_GSEA_tCell_Bottomly_phenotypepermutations[[i]]$documentation$rank)] 
  
  
  # Demethylation 
  
  # store the default number of differentially enriched gene sets 
  rank_cP_GSEA_Demethylation_phenpermutation_default[i] <- optimRank_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank[1]
  # store the optimal number of differentially enriched gene sets 
  rank_cP_GSEA_Demethylation_phenpermutation_optim[i] <- optimRank_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank[length(optimRank_cP_GSEA_Demethylation_Bottomly_Phenotypepermutations[[i]]$documentation$rank)] 
  
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
rank_GeneSet1_phenpermutation <- data.frame(cP_ORA = c(rank_cP_ORA_tCell_phenpermutation_default, rank_cP_ORA_tCell_phenpermutation_optim),
                                            GOSeq = c(rank_GOSeq_MetabolicProcess_phenpermutation_default, rank_GOSeq_MetabolicProcess_phenpermutation_optim), 
                                            PADOG = c(rank_PADOG_PrimImmun_phenpermutation_default, rank_PADOG_PrimImmun_phenpermutation_optim),
                                            #cP_GSEA = c(rank_cP_GSEA_tCell_phenpermutation_default, rank_cP_GSEA_tCell_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum 
                                            ID = rep(c(1:10),2)) # add the number of the permutation 


# transform data frame such that each observed adjusted p-value value is in a separate column, additionally labelled by the 
# associated GSA tool and state
rank_GeneSet1_phenpermutation_long <- pivot_longer(rank_GeneSet1_phenpermutation, 
                                                   cols=!c("state", "ID"), 
                                                   names_to="GSA_tool", 
                                                   values_to = "rank")

# Add column that indicates that we consider the "first" gene set for each tool 
rank_GeneSet1_phenpermutation_long$GS <- 1 


# add the combination of the GSA tool, the gene set number for each tool (1 vs. 2) and the permutation ID 
# -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
rank_GeneSet1_phenpermutation_long$unique_ID <- paste0(rank_GeneSet1_phenpermutation_long$GSA_tool,
                                                       "_", 
                                                       rank_GeneSet1_phenpermutation_long$GS, 
                                                       "_",
                                                       rank_GeneSet1_phenpermutation_long$ID)


# (ii) Random permutations of the random permutations of the sample conditions, second data set 
rank_GeneSet2_phenpermutation <- data.frame(cP_ORA = c(rank_cP_ORA_Demethylation_phenpermutation_default, rank_cP_ORA_Demethylation_phenpermutation_optim),
                                            GOSeq = c(rank_GOSeq_CellularProcess_phenpermutation_default, rank_GOSeq_CellularProcess_phenpermutation_optim), 
                                            PADOG = c(rank_PADOG_GraftvsHost_phenpermutation_default, rank_PADOG_GraftvsHost_phenpermutation_optim),
                                            #cP_GSEA = c(rank_cP_GSEA_Demethylation_phenpermutation_default, rank_cP_GSEA_Demethylation_phenpermutation_optim),
                                            state = c(rep("Default",10), rep("Minimum",10)), # state: default vs. optimum 
                                            ID = rep(c(1:10),2)) # add the number of the permutation 


# transform data frame such that each observed adjusted p-value is in a separate column, additionally labelled by the 
# associated GSA tool and state
rank_GeneSet2_phenpermutation_long <- pivot_longer(rank_GeneSet2_phenpermutation, 
                                                   cols=!c("state", "ID"), 
                                                   names_to="GSA_tool", 
                                                   values_to = "rank")

# Indicate that we here consider the "second" gene set for each tool 
rank_GeneSet2_phenpermutation_long$GS <- 2 



# add the combination of the GSA tool, the gene set number for each tool (1 vs. 2) and the permutation ID 
# -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
rank_GeneSet2_phenpermutation_long$unique_ID <- paste0(rank_GeneSet2_phenpermutation_long$GSA_tool,
                                                       "_", 
                                                       rank_GeneSet2_phenpermutation_long$GS, 
                                                       "_", 
                                                       rank_GeneSet2_phenpermutation_long$ID)





# Combine both data sets into one big data set 
rank_allgenesets_phenpermutation_long <- rbind(rank_GeneSet1_phenpermutation_long, 
                                               rank_GeneSet2_phenpermutation_long)

# transform column indicating gene set (gene set 1 vs. gene set 2) to factor
rank_allgenesets_phenpermutation_long$GS <- factor(rank_allgenesets_phenpermutation_long$GS)


# transform GSA tools to factors
# -> this way we can fix the order of the tools in the graphic
rank_allgenesets_phenpermutation_long$GSA_tool <- 
  factor(rank_allgenesets_phenpermutation_long$GSA_tool, 
         levels = c("GOSeq", "cP_ORA",  "PADOG", "cP_GSEA"))



# Replace "GSEAPreranked" by "GSEA- \n Preranked" (line break makes plot easier to understand)
add_labels_xaxis <- levels(rank_allgenesets_phenpermutation_long$GSA_tool)

# Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"


# Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"






##########
### ggplot 
##########


plot_adj_bottomly_phenpermutations <- 
  ggplot(data =rank_allgenesets_phenpermutation_long, 
         aes(x = interaction(GSA_tool, state, lex.order = TRUE), 
             y = rank, group = 1)) + 
  geom_line( aes(group=unique_ID, color = GS), size=0.5, alpha=0.7) + 
  geom_point(aes(color = GS),size = 1.5, alpha = 0.7) + 
  scale_x_discrete(labels= rep(c("Default", "Minimum"), times = 4)) + 
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 9), 
        axis.title.x = element_text(vjust = -2), 
        plot.margin = margin(t=1, b =3, l=1, r=1, unit="cm")) + ## add space below the actual plot (needed for the GSA tool names)
  annotate(geom = "text", x = 1.5 + 2*(0:3), y = -0.2, label = add_labels_xaxis, size =3 )+ 
  coord_cartesian(ylim = c(0,1),  clip = "off")+ # clip = "off" required to add GSA tool names below the plot
  xlab("Computational GSA method") +
  ylab("Relative rank") + 
  labs(color = "Gene set") +
  scale_color_discrete(labels=c("1", 
                                "2")) 










