################################################################################
### Results Illustrations for the optimization of the p-values for Pickrell data
################################################################################

library(ggplot2)
library(tidyr)

# set working directory 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/Assessment_OverOptimism")

################################################################################
# Load optimization results ####################################################
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


# Load DAVID results : Achtung DAVID fehlt bis jetzt noch 

# Load clusterProfiler's GSEA results

# Demethylation 
load("./Results/optimP_cP_GSEA_Demethylation_Pickrell_PhenotypePermutations.RData")

# t-Cell mediated immunity 
load("./Results/optimP_cP_GSEA_tCell_Pickrell_PhenotypePermutations.RData")


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
  padj_GOSeq_tCell_phenpermutation_default[i] <- optimP_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_GOSeq_tCell_phenpermutation_optim[i] <- optimP_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_GOSeq_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
  
  # Demethylation 
  
  # store the default number of differentially enriched gene sets 
  padj_GOSeq_Demethylation_phenpermutation_default[i] <- optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_GOSeq_Demethylation_phenpermutation_optim[i] <- optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_GOSeq_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
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
  padj_cP_ORA_tCell_phenpermutation_default[i] <- optimP_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_cP_ORA_tCell_phenpermutation_optim[i] <- optimP_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
  
  # Demethylation 
  
  # store the default number of differentially enriched gene sets 
  padj_cP_ORA_Demethylation_phenpermutation_default[i] <- optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_cP_ORA_Demethylation_phenpermutation_optim[i] <- optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
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
  padj_PADOG_PrimImmun_phenpermutation_default[i] <- optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_PADOG_PrimImmun_phenpermutation_optim[i] <- optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_PADOG_PrimaryImmunodeficiency_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
  
  # Demethylation 
  
  # store the default number of differentially enriched gene sets 
  padj_PADOG_GraftvsHost_phenpermutation_default[i] <- optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_PADOG_GraftvsHost_phenpermutation_optim[i] <- optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_PADOG_GraftvsHost_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
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
  padj_cP_GSEA_tCell_phenpermutation_default[i] <- optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_cP_GSEA_tCell_phenpermutation_optim[i] <- optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_GSEA_tcell_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
  
  # Demethylation 
  
  # store the default number of differentially enriched gene sets 
  padj_cP_GSEA_Demethylation_phenpermutation_default[i] <- optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_cP_GSEA_Demethylation_phenpermutation_optim[i] <- optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_GSEA_Demethylation_Pickrell_phenotypepermutations[[i]]$documentation$p_adj)] 
  
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


# transform data frame such that each observed adjusted p-value value is in a separate column, additionally labelled by the 
# associated GSA tool and state
padj_GeneSet1_phenpermutation_long <- pivot_longer(padj_GeneSet1_phenpermutation, 
                                                   cols=!c("state", "ID"), 
                                                   names_to="GSA_tool", 
                                                   values_to = "padj")

# Add column that indicates that we consider the "first" gene set for each tool 
padj_GeneSet1_phenpermutation_long$GS <- 1 


# add the combination of the GSA tool, the gene set number for each tool (1 vs. 2) and the permutation ID 
# -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
padj_GeneSet1_phenpermutation_long$unique_ID <- paste0(padj_GeneSet1_phenpermutation_long$GSA_tool,
                                                       "_", 
                                                       padj_GeneSet1_phenpermutation_long$GS, 
                                                       "_",
                                                       padj_GeneSet1_phenpermutation_long$ID)


# Indicate name of the gene set whose adjusted p-value was optimized 

padj_GeneSet1_phenpermutation_long$gene_set <- ifelse(padj_GeneSet1_phenpermutation_long$GSA_tool == "PADOG", 
                                                      "PrimImmun", 
                                                      "tCell") 


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


# transform data frame such that each observed adjusted p-value is in a separate column, additionally labelled by the 
# associated GSA tool and state
padj_GeneSet2_phenpermutation_long <- pivot_longer(padj_GeneSet2_phenpermutation, 
                                                   cols=!c("state", "ID"), 
                                                   names_to="GSA_tool", 
                                                   values_to = "padj")

# Indicate that we here consider the "second" gene set for each tool 
padj_GeneSet2_phenpermutation_long$GS <- 2 



# add the combination of the GSA tool, the gene set number for each tool (1 vs. 2) and the permutation ID 
# -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
padj_GeneSet2_phenpermutation_long$unique_ID <- paste0(padj_GeneSet2_phenpermutation_long$GSA_tool,
                                                       "_", 
                                                       padj_GeneSet2_phenpermutation_long$GS, 
                                                       "_", 
                                                       padj_GeneSet2_phenpermutation_long$ID)



padj_GeneSet2_phenpermutation_long$gene_set <- ifelse(padj_GeneSet1_phenpermutation_long$GSA_tool == "PADOG", 
                                                      "GraftvsHost", 
                                                      "Demethylation")




# Combine both data sets into one big data set 
padj_allgenesets_phenpermutation_long <- rbind(padj_GeneSet1_phenpermutation_long, 
                                               padj_GeneSet2_phenpermutation_long)


# transform GSA tools to factors
# -> this way we can fix the order of the tools in the graphic
padj_allgenesets_phenpermutation_long$GSA_tool <- factor(padj_allgenesets_phenpermutation_long$GSA_tool, 
                                                         levels = c("GOSeq", "DAVID", "cP_ORA",  "PADOG", "cP_GSEA", "GSEA", "GSEAPreranked"))


# transform gene sets to factors 
# -> this way we can fix the order of the gene sets in the legend 
padj_allgenesets_phenpermutation_long$gene_set <- factor(padj_allgenesets_phenpermutation_long$gene_set, 
                                                         levels = c("tCell", "Demethylation","PrimImmun", "GraftvsHost"))


# Replace "GSEAPreranked" by "GSEA- \n Preranked" (line break makes plot easier to understand)
add_labels_xaxis <- levels(padj_allgenesets_phenpermutation_long$GSA_tool)
add_labels_xaxis[add_labels_xaxis == "GSEAPreranked"] <- "GSEA- \n Preranked"


# Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"


# Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"





# add column which indicates default significant threshold used by each tool:
# GSEA (web) and GSEAPreranked: FDR < 0.25, remaining GSA methods: 0.05
padj_allgenesets_phenpermutation_long$sig_threshold <- ifelse(grepl(paste(c("GSEAPreranked","\\bGSEA\\b"), collapse='|'),padj_allgenesets_phenpermutation_long$GSA_tool), 
                                                              yes = 0.25, no = 0.05)


##########
### ggplot 
##########

plot_adj_pickrell_phenpermutations <- 
ggplot(data =padj_allgenesets_phenpermutation_long, 
       aes(x = interaction(GSA_tool, state, lex.order = TRUE), 
           y = padj, group = 1)) + 
geom_line( aes(group=unique_ID, color = gene_set), size=0.5, alpha=0.7) + 
geom_point(aes(color = gene_set),size = 1.5, alpha = 0.7) + 
scale_x_discrete(labels= rep(c("Default", "Minimum"), 
                               times = 6)) + 
theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 9), 
      axis.title.x = element_text(vjust = -2), 
      plot.margin = margin(t=1, b =3, l=1, r=1, unit="cm")) + ## add space below the actual plot (needed for the GSA tool names)
annotate(geom = "text", 
          x = 1.5 + 2*(0:6), 
          y = -0.2, 
          label = add_labels_xaxis, size =3 )+ 
coord_cartesian(ylim = c(0,1),  clip = "off")+ # clip = "off" required to add GSA tool names below the plot
xlab("Computational GSA method") +
ylab("Adjusted p-value") + 
geom_step(data = padj_allgenesets_phenpermutation_long, aes(interaction(GSA_tool, state, lex.order = TRUE), y= sig_threshold), 
          col="gray", linetype = "dashed") + 
  labs(color = "Gene set") +
  scale_color_discrete(labels=c("T cell mediated immunity", 
                                "Demethylation", 
                                "Primary immunodeficiency", 
                                "Graft versus host disease")) 



