#Illustration of results in step diagram

# load clusterProfiler's GSEA results

load("./Results/cP_GSEA_Results_Pickrell_PhenotypePermutations.RData")

# load documentation frame
doc_illustration <- optim_cP_GSEA_results_Pickrell_phenotypepermutation[[6]]$documentation

# insert line break for one optimal parameter
doc_illustration$optimal_parameter[doc_illustration$optimal_parameter == "by cpm>1 in at least 2 samples"] <- "by cpm > 1 in at \n least two samples"
doc_illustration$optimal_parameter[4] <- paste0("method ", doc_illustration$optimal_parameter[4])

# Set optimisation steps to sentence case
doc_illustration$step <- str_to_sentence(doc_illustration$step)
# Revert "id" to "ID"
doc_illustration$step[doc_illustration$step == "Duplicate gene id removal"] <- "Duplicate gene ID removal"
# Add numbering
doc_illustration$step <- paste0(1:length(doc_illustration$step), ". ", doc_illustration$step)


#for facet wrap labels
doc_illustration$phenotype_permutation <-as.factor("Phenotype Permutation 6")
doc_illustration$step <- factor(doc_illustration$step,
                                levels = doc_illustration$step)



#plot step diagram
step_diagram <- ggplot(doc_illustration) +
  geom_step(aes(x=step, y=n_DEGS, group=factor(phenotype_permutation), colour="#F8766D"), direction="hv", linetype=2) +
  geom_point(aes(x=step, y=n_DEGS, colour="#F8766D"), size=1.5) +
  geom_text(aes(x=step, y=n_DEGS), label=doc_illustration$optimal_parameter, nudge_y = 4,  colour="black",size = 4.5)+
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  #xlab("Optimization Step") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ylab("Number of differentially enriched gene sets") +
  xlab("Optimisation step") +
  scale_y_continuous(breaks=seq(0, 50, 10), limits = c(-5,40)) +
  theme(legend.position="none") #+
  facet_wrap(~phenotype_permutation, nrow=1, ncol=1)

# Export to .eps
ggsave(file="./Results_illustrations/FigureS1.eps",
       width = 10,
       height = 7)







