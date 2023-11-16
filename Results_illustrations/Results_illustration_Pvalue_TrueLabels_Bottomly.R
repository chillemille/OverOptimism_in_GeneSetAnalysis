################################################################################
### Results Illustrations for the optimization of the p-values for Bottomly data
### and the true sample labels #################################################
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
load("./Results/optimP_GOSeq_MetabolicProcess_Bottomly_OriginalPhenotype.RData")

# t Cell mediated immunity 
load("./Results/optimP_GOSeq_CellularProcess_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's ORA results

# Demethylation
load("./Results/optimP_ORA_MetabolicProcess_Bottomly_OriginalPhenotype.RData")

# t Cell mediated immunity 
load("./Results/optimP_ORA_CellularProcess_Bottomly_OriginalPhenotype.RData")


# Load DAVID results : Achtung DAVID fehlt bis jetzt noch 

# Load clusterProfiler's GSEA results

# Demethylation 
load("./Results/optimP_cP_GSEA_Demethylation_Bottomly_OriginalPhenotype.RData")

# t-Cell mediated immunity 
load("./Results/optimP_cP_GSEA_tCell_Bottomly_OriginalPhenotype.RData")


# Load PADOG results 

# (i) Primary Immunodeficiency 
load("./Results/optimP_PADOG_PrimaryImmunodeficiency_Bottomly_OriginalPhenotype.RData")

# (ii) Graft vs Host 
load("./Results/optimP_PADOG_GraftvsHost_Bottomly_OriginalPhenotype.RData")


# Load clusterProfiler's GSEA results

# Demethylation 
load("./Results/optimP_cP_GSEA_Demethylation_Bottomly_OriginalPhenotype.RData")

# t-Cell mediated immunity 
load("./Results/optimP_cP_GSEA_tCell_Bottomly_OriginalPhenotype.RData")

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
padj_cP_ORA_CellularProcess_truephen_default <- 
  optimP_cP_ORA_CellularProcess_Bottomly_originalphenotype$documentation$p_adj[1]                     
# optimal p_adj
padj_cP_ORA_CellularProcess_truephen_optim <- 
  optimP_cP_ORA_CellularProcess_Bottomly_originalphenotype$documentation$p_adj[length(optimP_cP_ORA_CellularProcess_Bottomly_originalphenotype$documentation$p_adj)] 


# Gene set Demethylation

# optimal p_adj
padj_cP_ORA_MetabolicProcess_truephen_optim <- 
  optimP_cP_ORA_MetabolicProcess_Bottomly_originalphenotype$documentation$p_adj[length(optimP_cP_ORA_MetabolicProcess_Bottomly_originalphenotype$documentation$p_adj)] 
# default p_adj
padj_cP_ORA_MetabolicProcess_truephen_default <- 
  optimP_cP_ORA_MetabolicProcess_Bottomly_originalphenotype$documentation$p_adj[1]                                                                                       





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
padj_GeneSet1_truephen <- data.frame(cP_ORA = c(padj_cP_ORA_MetabolicProcess_truephen_default, padj_cP_ORA_MetabolicProcess_truephen_optim),
                                            GOSeq = c(padj_GOSeq_MetabolicProcess_truephen_default, padj_GOSeq_MetabolicProcess_truephen_optim), 
                                            DAVID = c(padj_DAVID_MetabolicProcess_truephen_default, padj_DAVID_MetabolicProcess_truephen_optim), 
                                            PADOG = c(padj_PADOG_PrimImmun_truephen_default, padj_PADOG_PrimImmun_truephen_optim),
                                            cP_GSEA = c(padj_cP_GSEA_tCell_truephen_default, padj_cP_GSEA_tCell_truephen_optim),
                                            GSEA = c(padj_GSEA_tCell_truephen_default,padj_GSEA_tCell_truephen_optim),
                                            GSEAPreranked = c(padj_GSEAPreranked_tCell_truephen_default, padj_GSEAPreranked_tCell_truephen_optim),
                                            state = c("Default","Maximum")) # state: default vs. optimum 


# transform data frame such that each observed adjusted p-value value is in a separate column, additionally labelled by the 
# associated GSA tool and state
padj_GeneSet1_truephen_long <- pivot_longer(padj_GeneSet1_truephen, 
                                                   cols=!c("state"), 
                                                   names_to="GSA_tool", 
                                                   values_to = "padj")

# Add column that indicates that we consider the "first" gene set for each tool 
padj_GeneSet1_truephen_long$GS <- 1 


# add the combination of the GSA tool, the gene set number for each tool (1 vs. 2) and the permutation ID 
# -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
padj_GeneSet1_truephen_long$unique_ID <- paste0(padj_GeneSet1_truephen_long$GSA_tool,
                                                       "_", 
                                                       padj_GeneSet1_truephen_long$GS)


# Indicate name of the gene set whose adjusted p-value was optimized 

padj_GeneSet1_truephen_long$gene_set <- NA

for(i in 1:nrow(padj_GeneSet1_truephen_long)){
  
  # for GSEA, GSEAPreranked and clusterProfiler's GSEA: the considered gene set 
  # is t Cell mediated immunity 
  if(grepl("GSEA",padj_GeneSet1_truephen_long$GSA_tool[[i]])){
    
    padj_GeneSet1_truephen_long$gene_set[i] <- "tCell"
    
    # for PADOG, the considered gene set is "Primary Immunodeficiency"
  }else if(padj_GeneSet1_truephen_long$GSA_tool[[i]] == "PADOG"){
    
    padj_GeneSet1_truephen_long$gene_set[i] <- "PrimImmun"
    
    # for clusterProfiler's ORA and GOSeq, the considered gene set is 
    # "Metabolic Process"
  }else 
    
    padj_GeneSet1_truephen_long$gene_set[i] <- "MetabolicProcess"
  
  
}



# (ii) Random permutations of the random permutations of the sample conditions, second data set 
padj_GeneSet2_truephen <- data.frame(cP_ORA = c(padj_cP_ORA_CellularProcess_truephen_default, padj_cP_ORA_CellularProcess_truephen_optim),
                                            GOSeq = c(padj_GOSeq_CellularProcess_truephen_default, padj_GOSeq_CellularProcess_truephen_optim), 
                                            PADOG = c(padj_PADOG_GraftvsHost_truephen_default, padj_PADOG_GraftvsHost_truephen_optim),
                                            GSEA = c(padj_GSEA_Demethylation_truephen_default,padj_GSEA_Demethylation_truephen_optim),
                                            cP_GSEA = c(padj_cP_GSEA_Demethylation_truephen_default, padj_cP_GSEA_Demethylation_truephen_optim),
                                            GSEAPreranked = c(padj_GSEAPreranked_Demethylation_truephen_default, padj_GSEAPreranked_Demethylation_truephen_optim),
                                            state = c("Default","Maximum")) # state: default vs. optimum 


# transform data frame such that each observed adjusted p-value is in a separate column, additionally labelled by the 
# associated GSA tool and state
padj_GeneSet2_truephen_long <- pivot_longer(padj_GeneSet2_truephen, 
                                                   cols=!c("state"), 
                                                   names_to="GSA_tool", 
                                                   values_to = "padj")

# Indicate that we here consider the "second" gene set for each tool 
padj_GeneSet2_truephen_long$GS <- 2 



# add the combination of the GSA tool, the gene set number for each tool (1 vs. 2) and the permutation ID 
# -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
padj_GeneSet2_truephen_long$unique_ID <- paste0(padj_GeneSet2_truephen_long$GSA_tool,
                                                       "_", 
                                                       padj_GeneSet2_truephen_long$GS)



# Indicate name of the gene set whose adjusted p-value was optimized 

padj_GeneSet2_truephen_long$gene_set <- NA

for(i in 1:nrow(padj_GeneSet2_truephen_long)){
  
  # for GSEA, GSEAPreranked and clusterProfiler's GSEA: the considered gene set 
  # is t Cell mediated immunity 
  if(grepl("GSEA",padj_GeneSet2_truephen_long$GSA_tool[[i]])){
    
    padj_GeneSet2_truephen_long$gene_set[i] <- "Demethylation"
    
    # for PADOG, the considered gene set is "Primary Immunodeficiency"
  }else if(padj_GeneSet2_truephen_long$GSA_tool[[i]] == "PADOG"){
    
    padj_GeneSet2_truephen_long$gene_set[i] <- "GraftvsHost"
    
    # for clusterProfiler's ORA and GOSeq, the considered gene set is 
    # "Metabolic Process"
  }else 
    
    padj_GeneSet2_truephen_long$gene_set[i] <- "CellularProcess"
  
  
}




# Combine both data sets into one big data set 
padj_allgenesets_truephen_long <- rbind(padj_GeneSet1_truephen_long, 
                                               padj_GeneSet2_truephen_long)


# transform GSA tools to factors
# -> this way we can fix the order of the tools in the graphic
padj_allgenesets_truephen_long$GSA_tool <- 
  factor(padj_allgenesets_truephen_long$GSA_tool, 
         levels = c("GOSeq","cP_ORA", "DAVID",  "PADOG", "cP_GSEA", "GSEA", "GSEAPreranked"))


padj_allgenesets_truephen_long$gene_set <- 
  factor(padj_allgenesets_truephen_long$gene_set, 
         levels = c("CellularProcess", "MetabolicProcess", "tCell", "Demethylation", "PrimImmun", "GraftvsHost"))



# Replace "GSEAPreranked" by "GSEA- \n Preranked" (line break makes plot easier to understand)
add_labels_xaxis <- levels(padj_allgenesets_truephen_long$GSA_tool)
add_labels_xaxis[add_labels_xaxis == "GSEAPreranked"] <- "GSEA- \n Preranked"

# Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"


# Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"


# add column which indicates default significant threshold used by each tool:
# GSEA (web) and GSEAPreranked: FDR < 0.25, remaining GSA methods: 0.05
padj_allgenesets_truephen_long$sig_threshold <- 
  ifelse(grepl(paste(c("GSEAPreranked","\\bGSEA\\b"), collapse='|'), padj_allgenesets_truephen_long$GSA_tool), 
         yes = 0.25, 
         no = 0.05)




##########
### ggplot 
##########


plot_adj_bottomly_truephen <- 
  ggplot(data =padj_allgenesets_truephen_long, 
         aes(x = interaction(GSA_tool, state, lex.order = TRUE), 
             y = padj, group = 1)) + 
  geom_line( aes(group=unique_ID, color = gene_set), size=0.5, alpha=0.7) + 
  geom_point(aes(color = gene_set),size = 1.5, alpha = 0.7) + 
  scale_x_discrete(labels= rep(c("Default", "Minimum"), times = 7)) + 
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 9), 
        axis.title.x = element_text(vjust = -2), 
        plot.margin = margin(t=1, b =3, l=1, r=1, unit="cm")) + ## add space below the actual plot (needed for the GSA tool names)
  annotate(geom = "text", x = 1.5 + 2*(0:6), y = -0.2, label = add_labels_xaxis, size =3 )+ 
  coord_cartesian(ylim = c(0,1),  clip = "off")+ # clip = "off" required to add GSA tool names below the plot
  xlab("Computational GSA method") +
  ylab("Adjusted p-value") + 
  scale_y_continuous(breaks = sort(c(seq(0,1,0.25), 0.05))) +
  geom_step(data = padj_allgenesets_truephen_long, aes(interaction(GSA_tool, state, lex.order = TRUE), y= sig_threshold), 
            col="gray", linetype = "dashed") +
  labs(color = "Gene set") +
  scale_color_discrete(labels=c("Cellular process", 
                                "Metabolic process", 
                                "T cell mediated immunity", 
                                "Demethylation", 
                                "Primary immunodeficiency", 
                                "Graft versus host disease")) 



