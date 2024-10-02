################################################################################
### Functions to generate the results illustration plots goal 2 and 3 ##########
################################################################################

library(ggplot2)
library(dplyr)


################################################################################
### Prepare the data frame to be plotted in the ggplot #########################
################################################################################

# argument goal: c("n_DEGS", "p_adj", "rel_rank")
# sample_labels: c("true_labels", "random_permutations")

prep_data_for_ggplot_pvaluerank <- function(default_to_optim_geneset1, default_to_optim_geneset2, goal, sample_labels){


  # if the data contain the data for the permutations of the sample labels, add column "ID" which contains the
  # which indicates the number of the corresponding permutation
  if(sample_labels == "random_permutations"){

    default_to_optim_geneset1$ID <- rep(c(1:10),2)
    default_to_optim_geneset2$ID <- rep(c(1:10),2)

  }

  # check which of both columns "state" and "ID" are contained in the data sets
  # "state" is always contained, "ID" is contained if the data set contains
  # data on the random permutations
  # we need this indicator "column" to NOT pivot these columns to a longer format
  columns <- c("state", "ID")[c("state", "ID") %in% colnames(default_to_optim_geneset1)]


  # transform data frame such that each observed adjusted p-value value is in a separate column, additionally labelled by the
  # associated GSA tool and state
  default_to_optim_geneset1 <- pivot_longer(default_to_optim_geneset1,
                                            cols=!columns,
                                            names_to="GSA_tool",
                                            values_to = "value")

  # Add column that indicates that we consider the "first" gene set for each tool
  default_to_optim_geneset1$GS <- 1


  # add the combination of
  # - the GSA tool,
  # - the gene set number for each tool (1 vs. 2)
  # - and the permutation ID, IF the data contains the permutations
  # -> this way, each value appears exactly twice and matches the default and optimal value for each tool in each permutation
  default_to_optim_geneset1$unique_ID <- paste0(default_to_optim_geneset1$GSA_tool,
                                                "_",
                                                default_to_optim_geneset1$GS)

  if(sample_labels == "random_permutations"){
    default_to_optim_geneset1$unique_ID <- paste0(default_to_optim_geneset1$unique_ID,
                                                  "_",
                                                  default_to_optim_geneset1$ID)

  }



  # transform data frame such that each observed adjusted p-value is in a separate column, additionally labelled by the
  # associated GSA tool and state
  default_to_optim_geneset2 <- pivot_longer(default_to_optim_geneset2,
                                            cols=!columns,
                                            names_to="GSA_tool",
                                            values_to = "value")

  # Indicate that we here consider the "second" gene set for each tool
  default_to_optim_geneset2$GS <- 2



  # add the combination of the GSA tool, the gene set number for each tool (1 vs. 2) (and the permutation ID
  # -> this way, each value appears exactly twice and matches the default and optimal n_DEGS value for each tool in each permutation
  default_to_optim_geneset2$unique_ID <- paste0(default_to_optim_geneset2$GSA_tool,
                                                "_",
                                                default_to_optim_geneset2$GS)

  if(sample_labels == "random_permutations"){
    default_to_optim_geneset2$unique_ID <- paste0(default_to_optim_geneset2$unique_ID,
                                                  "_",
                                                  default_to_optim_geneset2$ID)

  }


  # Combine both data sets into one big data set
  default_to_optim_allgenesets <- rbind(default_to_optim_geneset1,
                                        default_to_optim_geneset2)


  ############################
  ### fix order of GSA methods
  ############################

  # the order of the methods in the results illustration shall be based on the following
  allmethods <- c("GOSeq", "DAVID", "cP_ORA",  "PADOG", "cP_GSEA", "GSEA", "GSEAPreranked")
  # get the list of tools existent in the current data
  currentmethods <- unique(default_to_optim_allgenesets$GSA_tool)
  # specify the order of the existent methods in the results illustrations
  levels <- currentmethods[order(match(currentmethods,allmethods))]



  # transform GSA tools to factors
  # -> this way we can fix the order of the tools in the graphic
  default_to_optim_allgenesets$GSA_tool <- factor(default_to_optim_allgenesets$GSA_tool,
                                                  levels = levels)


  # transform gene sets to factors
  # -> this way we can fix the order of the gene sets in the legend
  default_to_optim_allgenesets$GS <- factor(default_to_optim_allgenesets$GS)




  if(goal =="p_adj"){
    # for the results illustrations of the adjusted p-value:
    # add column which indicates default significant threshold used by each tool:
    # GSEA (web) and GSEAPreranked: FDR < 0.25, remaining GSA methods: 0.05
    default_to_optim_allgenesets$sig_threshold <- ifelse(grepl(paste(c("GSEAPreranked","\\bGSEA\\b"), collapse='|'),
                                                               default_to_optim_allgenesets$GSA_tool),
                                                         yes = 0.25, no = 0.05)

  }

  return(default_to_optim_allgenesets)


}






################################################################################
### generate ggplot ############################################################
################################################################################

# for the illustration on the adjusted p-value, insert "p_adj" for argument goal


create_results_illustration_pvalue_rank <- function(default_to_optim_geneset1, default_to_optim_geneset2, goal = "p_adj", sample_labels){


  data_prep <- prep_data_for_ggplot_pvaluerank(default_to_optim_geneset1, default_to_optim_geneset2, goal, sample_labels)


  set_ylab <- ifelse(goal == "rel_rank", "(Relative) rank", "Adjusted p-value")


  # Replace "GSEAPreranked" by "GSEA- \n Preranked" (line break makes plot easier to understand)
  add_labels_xaxis <- levels(data_prep$GSA_tool)
  add_labels_xaxis[add_labels_xaxis == "GSEAPreranked"] <- "GSEA- \n Preranked"


  # Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
  add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"


  # Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
  add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"


  plot <-
    ggplot(data = data_prep,
           aes(x = interaction(GSA_tool, state, lex.order = TRUE),
               y = value, group = 1)) +
    geom_line( aes(group=unique_ID, color = GS), size=0.4) +
    geom_point(aes(color = GS),size = 1.3) +
    scale_x_discrete(labels= rep(c("Default", "Minimum"),
                                 times = length(unique(data_prep$GSA_tool)))) +
    theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 9),
          axis.title.x = element_text(vjust = -15, size = 14),
          plot.margin = margin(t=1, b =3, l=1, r=1, unit="cm"), ## add space below the actual plot (needed for the GSA tool names)
          axis.title.y = element_text(size = 14)) +
    annotate(geom = "text",
             x = 1.5 + 2*(0:(length(unique(data_prep$GSA_tool))-1)),
             y = -0.32,
             label = add_labels_xaxis, size =4 )+
    coord_cartesian(ylim = c(0,1),  clip = "off")+ # clip = "off" required to add GSA tool names below the plot
    xlab("GSA method") +
    ylab(set_ylab) +
    labs(color = "Gene set") +
    scale_color_discrete(labels=c("1",
                                  "2"))

  if(goal == "p_adj"){

    plot <- plot +  geom_step( dat = data_prep,
                               aes(interaction(GSA_tool, state, lex.order = TRUE),
                                   y= sig_threshold),
                               col="gray", linetype = "dashed")
  }


  return(plot)

}


