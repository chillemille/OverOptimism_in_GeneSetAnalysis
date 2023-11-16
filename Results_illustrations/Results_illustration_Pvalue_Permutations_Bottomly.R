################################################################################
### Results Illustrations for the optimization of the p-values for Bottomly data
### and the random permutations of the true sample labels ######################
################################################################################

library(ggplot2)
library(tidyr)

# set working directory 
setwd("/nfsmb/koll/milena.wuensch/Dokumente/Overoptimism_NEU/NEU/OverOptimism_GSA/Assessment_OverOptimism")

################################################################################
### Load corresponding results for Bottomly data set ###########################
################################################################################



#Load GOSeq results 

# Demethylation
load("./Results/optimP_GOSeq_MetabolicProcess_Bottomly_PhenotypePermutations.RData")

# t Cell mediated immunity 
load("./Results/optimP_GOSeq_CellularProcess_Bottomly_PhenotypePermutations.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/optimP_ORA_MetabolicProcess_Bottomly_PhenotypePermutations.RData")

# t Cell mediated immunity 
load("./Results/optimP_ORA_CellularProcess_Bottomly_PhenotypePermutations.RData")


# Load DAVID results : Achtung DAVID fehlt bis jetzt noch 

# Load clusterProfiler's GSEA results

# Demethylation 
load("./Results/optimP_cP_GSEA_Demethylation_Bottomly_PhenotypePermutations.RData")

# t-Cell mediated immunity 
load("./Results/optimP_cP_GSEA_tCell_Bottomly_PhenotypePermutations.RData")


# Load PADOG results 

# (i) Primary Immunodeficiency 
load("./Results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_PhenotypePermutations.RData")

# (ii) Graft vs Host 
load("./Results/optimP_PADOG_GraftvsHost_Bottomly_PhenotypePermutations.RData")


# Load clusterProfiler's GSEA results

# Demethylation 
load("./Results/optimP_cP_GSEA_Demethylation_Bottomly_PhenotypePermutations.RData")

# t-Cell mediated immunity 
load("./Results/optimP_cP_GSEA_tCell_Bottomly_PhenotypePermutations.RData")

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


# Cellular Process

# store default p_adj for the 10 permutations of the true sample labels  
padj_cP_ORA_CellularProcess_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels 
padj_cP_ORA_CellularProcess_phenpermutation_optim <- c()


# Metabolic Process
# store default p_adj for the 10 permutations of the true sample labels  
padj_cP_ORA_MetabolicProcess_phenpermutation_default <- c()
# store optimal p_adj for the 10 permutations of the true sample labels 
padj_cP_ORA_MetabolicProcess_phenpermutation_optim <- c()


for(i in 1:10){
  
  # Cellular Process
  
  # store the default number of differentially enriched gene sets 
  padj_cP_ORA_CellularProcess_phenpermutation_default[i] <- 
    optimP_cP_ORA_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_cP_ORA_CellularProcess_phenpermutation_optim[i] <- 
    optimP_cP_ORA_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_CellularProcess_Bottomly_Phenotypepermutations[[i]]$documentation$p_adj)] 
  
  
  # Metabolic Process 
  
  # store the default number of differentially enriched gene sets 
  padj_cP_ORA_MetabolicProcess_phenpermutation_default[i] <- 
    optimP_cP_ORA_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[1]
  # store the optimal number of differentially enriched gene sets 
  padj_cP_ORA_MetabolicProcess_phenpermutation_optim[i] <- 
    optimP_cP_ORA_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$p_adj[length(optimP_cP_ORA_MetabolicProcess_Bottomly_phenotypepermutations[[i]]$documentation$p_adj)] 
  
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
padj_GeneSet1_phenpermutation <- data.frame(cP_ORA = c(padj_cP_ORA_MetabolicProcess_phenpermutation_default, padj_cP_ORA_MetabolicProcess_phenpermutation_optim),
                                            GOSeq = c(padj_GOSeq_MetabolicProcess_phenpermutation_default, padj_GOSeq_MetabolicProcess_phenpermutation_optim), 
                                            DAVID = c(padj_DAVID_CellularProcess_phenpermutation_default, padj_DAVID_CellularProcess_phenpermutation_optim),
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

padj_GeneSet1_phenpermutation_long$gene_set <- NA

for(i in 1:nrow(padj_GeneSet1_phenpermutation_long)){
  
  # for GSEA, GSEAPreranked and clusterProfiler's GSEA: the considered gene set 
  # is t Cell mediated immunity 
  if(grepl("GSEA",padj_GeneSet1_phenpermutation_long$GSA_tool[[i]])){
    
    padj_GeneSet1_phenpermutation_long$gene_set[i] <- "tCell"
    
    # for PADOG, the considered gene set is "Primary Immunodeficiency"
  }else if(padj_GeneSet1_phenpermutation_long$GSA_tool[[i]] == "PADOG"){
    
    padj_GeneSet1_phenpermutation_long$gene_set[i] <- "PrimImmun"
    
    # for clusterProfiler's ORA and GOSeq, the considered gene set is 
    # "Metabolic Process"
    }else 
  
      padj_GeneSet1_phenpermutation_long$gene_set[i] <- "MetabolicProcess"
  
  
}



# (ii) Random permutations of the random permutations of the sample conditions, second data set 
padj_GeneSet2_phenpermutation <- data.frame(cP_ORA = c(padj_cP_ORA_CellularProcess_phenpermutation_default, padj_cP_ORA_CellularProcess_phenpermutation_optim),
                                            GOSeq = c(padj_GOSeq_CellularProcess_phenpermutation_default, padj_GOSeq_CellularProcess_phenpermutation_optim), 
                                            DAVID = c(padj_DAVID_MetabolicProcess_phenpermutation_default, padj_DAVID_MetabolicProcess_phenpermutation_optim), 
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



# Indicate name of the gene set whose adjusted p-value was optimized 

padj_GeneSet2_phenpermutation_long$gene_set <- NA

for(i in 1:nrow(padj_GeneSet2_phenpermutation_long)){
  
  # for GSEA, GSEAPreranked and clusterProfiler's GSEA: the considered gene set 
  # is t Cell mediated immunity 
  if(grepl("GSEA",padj_GeneSet2_phenpermutation_long$GSA_tool[[i]])){
    
    padj_GeneSet2_phenpermutation_long$gene_set[i] <- "Demethylation"
    
    # for PADOG, the considered gene set is "Primary Immunodeficiency"
  }else if(padj_GeneSet2_phenpermutation_long$GSA_tool[[i]] == "PADOG"){
    
    padj_GeneSet2_phenpermutation_long$gene_set[i] <- "GraftvsHost"
    
    # for clusterProfiler's ORA and GOSeq, the considered gene set is 
    # "Metabolic Process"
  }else 
    
    padj_GeneSet2_phenpermutation_long$gene_set[i] <- "CellularProcess"
  
  
}




# Combine both data sets into one big data set 
padj_allgenesets_phenpermutation_long <- rbind(padj_GeneSet1_phenpermutation_long, 
                                               padj_GeneSet2_phenpermutation_long)


# transform GSA tools to factors
# -> this way we can fix the order of the tools in the graphic
padj_allgenesets_phenpermutation_long$GSA_tool <- 
  factor(padj_allgenesets_phenpermutation_long$GSA_tool, 
         levels = c("GOSeq", "DAVID", "cP_ORA",  "PADOG", "cP_GSEA", "GSEA", "GSEAPreranked"))


padj_allgenesets_phenpermutation_long$gene_set <- 
  factor(padj_allgenesets_phenpermutation_long$gene_set, 
         levels = c("CellularProcess", "MetabolicProcess", "tCell", "Demethylation", "PrimImmun", "GraftvsHost"))



# Replace "GSEAPreranked" by "GSEA- \n Preranked" (line break makes plot easier to understand)
add_labels_xaxis <- levels(padj_allgenesets_phenpermutation_long$GSA_tool)
add_labels_xaxis[add_labels_xaxis == "GSEAPreranked"] <- "GSEA- \n Preranked"

# Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"


# Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"


# add column which indicates default significant threshold used by each tool:
# GSEA (web) and GSEAPreranked: FDR < 0.25, remaining GSA methods: 0.05
padj_allgenesets_phenpermutation_long$sig_threshold <- 
  ifelse(grepl(paste(c("GSEAPreranked","\\bGSEA\\b"), collapse='|'), padj_allgenesets_phenpermutation_long$GSA_tool), 
         yes = 0.25, 
         no = 0.05)




##########
### ggplot 
##########


plot_adj_bottomly_phenpermutations <- 
  ggplot(data =padj_allgenesets_phenpermutation_long, 
       aes(x = interaction(GSA_tool, state, lex.order = TRUE), 
           y = padj, group = 1)) + 
geom_line( aes(group=unique_ID, color = gene_set), size=0.5, alpha=0.7) + 
geom_point(aes(color = gene_set),size = 1.5, alpha = 0.7) + 
scale_x_discrete(labels= rep(c("Default", "Minimum"), times = 6)) + 
theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 9), 
      axis.title.x = element_text(vjust = -2), 
      plot.margin = margin(t=1, b =3, l=1, r=1, unit="cm")) + ## add space below the actual plot (needed for the GSA tool names)
annotate(geom = "text", x = 1.5 + 2*(0:6), y = -0.2, label = add_labels_xaxis, size =3 )+ 
coord_cartesian(ylim = c(0,1),  clip = "off")+ # clip = "off" required to add GSA tool names below the plot
xlab("Computational GSA method") +
ylab("Adjusted p-value") + 
geom_step(data = padj_allgenesets_phenpermutation_long, aes(interaction(GSA_tool, state, lex.order = TRUE), y= sig_threshold), 
          col="gray", linetype = "dashed") +
labs(color = "Gene set") +
scale_color_discrete(labels=c("Cellular process", 
                              "Metabolic process", 
                              "T cell mediated immunity", 
                              "Demethylation", 
                              "Primary immunodeficiency", 
                              "Graft versus host disease")) 










