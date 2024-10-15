library(ggplot2)


uncertainties  <- data.frame(n_uncertainties = c(4,2,2,1,2,3,0,3,2,3,3,2,2,1),
                            GSA_tools = c(rep("GOSeq", 2), rep("DAVID", 2), rep("clusterProfiler's ORA", 2), rep("PADOG",2), rep("clusterProfiler's GSEA",2), rep("GSEA",2), rep("GSEAPreranked",2)),
                            type_uncertainty = rep(c("Parameters", "Data Preprocessing"), 7))

# make column GSA_tools an ordered factor so that ggplot does not change the order automatically
uncertainties$GSA_tools  <- factor(uncertainties$GSA_tools, levels = unique(uncertainties$GSA_tools))
uncertainties$type_uncertainty  <- factor(uncertainties$type_uncertainty, levels = unique(uncertainties$type_uncertainty))


ggplot(uncertainties, aes(n_uncertainties, GSA_tools, fill = type_uncertainty)) +
  geom_bar(stat="identity", position = "stack") +
  labs(fill='Uncertainty') +
  coord_flip() +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  ylab("GSA methods") +
  xlab("Number of exploited uncertainties") +
  scale_fill_manual(values=c('#4AA4DE', 'yellowgreen'))


ggsave("./Results/Figures/Figure2.eps",
       height = 5.5,
       width = 8.5,
       device = "eps")


