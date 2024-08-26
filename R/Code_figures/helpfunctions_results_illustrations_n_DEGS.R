################################################################################
### Functions to generate the results illustration plots goal 1 ################
################################################################################

library(ggplot2)
library(dplyr)


################################################################################
### Prepare the data frame to be plotted in the ggplot #########################
################################################################################

# sample_labels: c("true_labels", "random_permutations")

prep_data_for_ggplot_n_DEGS <- function(default_to_optim, sample_labels){


  # if the data contain the data for the permutations of the sample labels, add column "ID" which contains the
  # which indicates the number of the corresponding permutation
  if(sample_labels == "random_permutations"){

    default_to_optim$ID <- rep(c(1:10),2)

  }

  # check which of both columns "state" and "ID" are contained in the data sets
  # "state" is always contained, "ID" is contained if the data set contains
  # data on the random permutations
  # we need this indicator "column" to NOT pivot these columns to a longer format
  columns <- c("state", "ID")[c("state", "ID") %in% colnames(default_to_optim)]


  # transform data frame such that each observed adjusted p-value value is in a separate column, additionally labelled by the
  # associated GSA tool and state
  default_to_optim <- pivot_longer(default_to_optim,
                                   cols=!columns,
                                   names_to="GSA_tool",
                                   values_to = "n_DEGS")


  # add the combination of
  # - the GSA tool,
  # - the gene set number for each tool (1 vs. 2)
  # - and the permutation ID, IF the data contains the permutations
  # -> this way, each value appears exactly twice and matches the default and optimal value for each tool in each permutation
  default_to_optim$unique_ID <- default_to_optim$GSA_tool



  if(sample_labels == "random_permutations"){
    default_to_optim$unique_ID <- paste0(default_to_optim$unique_ID,
                                         "_",
                                         default_to_optim$ID)

  }


  ############################
  ### fix order of GSA methods
  ############################

  # the order of the methods in the results illustration shall be based on the following
  allmethods <- c("GOSeq", "DAVID", "cP_ORA",  "PADOG", "cP_GSEA", "GSEA", "GSEAPreranked")
  # get the list of tools existent in the current data
  currentmethods <- unique(default_to_optim$GSA_tool)
  # specify the order of the existent methods in the results illustrations
  levels <- currentmethods[order(match(currentmethods,allmethods))]



  # transform GSA tools to factors
  # -> this way we can fix the order of the tools in the graphic
  default_to_optim$GSA_tool <- factor(default_to_optim$GSA_tool,
                                      levels = levels)



  return(default_to_optim)


}

################################################################################
### generate ggplot ############################################################
################################################################################


# function to square transform the y-axis
# we display the data on the sqrt-scale due to the different magnitudes of the
# points

# note: vpos_methodlabels determines the vertical position of the method labels
# in function annotate() in the plot (this needs to be adjusted to each data since
# of n_DEGS differs considerably between the magnitude settings)
labels_sq <- function(x) {
  paste(x^2)
}


## ggplot
create_results_illustration_n_DEGS <- function(default_to_optim, sample_labels, vpos_methodlabels){


  data_prep <- prep_data_for_ggplot_n_DEGS(default_to_optim, sample_labels)


  # Replace "GSEAPreranked" by "GSEA- \n Preranked" (line break makes plot easier to understand)
  add_labels_xaxis <- levels(data_prep$GSA_tool)
  add_labels_xaxis[add_labels_xaxis == "GSEAPreranked"] <- "GSEA- \n Preranked"


  # Replace "cP_GSEA" by "clusterProfiler's \n GSEA" (line break makes plot easier to understand)
  add_labels_xaxis[add_labels_xaxis == "cP_GSEA"] <- "clusterProfiler's \n GSEA"


  # Replace "cP_ORA" by "clusterProfiler's \n ORA" (line break makes plot easier to understand)
  add_labels_xaxis[add_labels_xaxis == "cP_ORA"] <- "clusterProfiler's \n ORA"



  plot <- ggplot(data = data_prep,
                 aes(x = interaction(GSA_tool, state, lex.order = TRUE),
                     y = sqrt(n_DEGS), group = 1)) +
    geom_line(aes(group=unique_ID), size=0.4,  col ="#F8766D") +
    geom_point(size = 1.3,  col = "#F8766D") +
    scale_x_discrete(labels= rep(c("Default", "Maximum"),
                                 times = length(unique(data_prep$GSA_tool)))) +
    theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 11),
          axis.title.x = element_text(vjust = -15, size = 14),
          plot.margin = margin(t=1, b =3, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
          axis.title.y = element_text(size =14)) +
    # Add the tool names to the plot:
    annotate(geom = "text",
             x = 1.5 + 2*(0:(length(unique(data_prep$GSA_tool))-1)),
             y = vpos_methodlabels,
             label = add_labels_xaxis, size =4 )+
    coord_cartesian(ylim=c(0,round_any(sqrt(max(data_prep$n_DEGS)), 5, ceiling)),clip = "off")+ # clip = "off" required to add GSA tool names below the plot
    xlab("GSA methods") +
    ylab("Number of differentially enriched gene sets")+
    theme(panel.grid.minor = element_blank()) +
    scale_y_continuous(label = labels_sq)


  return(plot)

}


