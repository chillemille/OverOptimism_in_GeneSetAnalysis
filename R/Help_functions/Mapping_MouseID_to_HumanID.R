################################################################################
### Mapping Mouse ENSEMBL ID -> Human Ensembl ID -> Human HUGO Gene Symbols ####
################################################################################

library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)



# upload mapping obtained by HCOP tool
# https://www.genenames.org/tools/hcop/
# more precisely: http://ftp.ebi.ac.uk/pub/databases/genenames/hcop/ (downloaded Nov 27 2023)
mapping_human_mouse <- read.csv("./GeneExpression_data/Mapping_ENSEMBL_Human_Mouse.csv",
                                sep = ";")


################################################################################
### Part 1: Framework Mapping Mouse ENSEMBL ID -> Human Ensembl ID #############
################################################################################


# remove irrelevant columns from data set such that only those containing the human gene IDs,
# the mouse ENSEMBL gene IDs as well as the tools by which the respective mapping is supported
# the least will be used as majority vote in the case of ambiguous mappings

ind_colnames <- which(colnames(mapping_human_mouse) %in% c("human_ensembl_gene",   "mouse_ensembl_gene",   "support"))
mapping_human_mouse <- mapping_human_mouse[,   ind_colnames]


# remove all rows from the mapping which contain "-" in either the human ensembl
# or in the mouse ensembl gene ID

table(mapping_human_mouse$human_ensembl_gene != "-") # -> no "-"-values
table(mapping_human_mouse$mouse_ensembl_gene != "-") # -> several "-"-values

# -> remove all rows that contain "-"-values in the mouse ensembl ID
mapping_human_mouse <- mapping_human_mouse[mapping_human_mouse$mouse_ensembl_gene != "-",   ]
dim(mapping_human_mouse)

################################################################################
### (I) Single human ENSEMBL IDs mapped to several mouse ENSEMBL IDs ###########
################################################################################

# search for all human Ensembl IDs that were mapped to several mouse ensembl IDs
mapping_dupl_humanID <- mapping_human_mouse[duplicated(mapping_human_mouse$human_ensembl_gene),  ]
dupl_gene_list_humanID <- unique(mapping_human_mouse$human_ensembl_gene[duplicated(mapping_human_mouse$human_ensembl_gene)])


# subset INITIAL mapping to those genes with a unique mapping from human ensembl IDs
# to mouse ensembl IDs
mapping_unique_humanID <- mapping_human_mouse[!mapping_human_mouse$human_ensembl_gene %in% dupl_gene_list_humanID,  ]


# go through each of the duplicated Ensembl gene IDs and count number of resources that
# support the respective mapping

# create empty data frame which will be filled successively with the final and ambiguous
disambig_mapping_humanID <- data.frame(matrix(NA,   nrow = 0,   ncol = 3))

for(i in 1:length(dupl_gene_list_humanID)){

  # get human ensembl gene ID that is mapped to several
  dupl_gene <- dupl_gene_list_humanID[i]

  # get whole mapping of the human ensembl gene to the several mouse ensembl IDs
  mapping_dupl_gene <- mapping_dupl_humanID[mapping_dupl_humanID$human_ensembl_gene == dupl_gene,  ]

  # count number of supporters for each mappings based on the number of kommas
  # in the column (0 komma -> 1 supporter,   1 komma -> 0 supporters,   ...)
  n_supporters <- str_count(mapping_dupl_gene$support,   pattern = ",  ") + 1

  # get mapping supported by the majority of tools (in case of ties,   choose lower index )
  ind_maxsupport <- min(which(n_supporters == max(n_supporters)))

  # get resulting row of mapping
  final_mapping <- mapping_dupl_gene[ind_maxsupport,  ]

  # bind rows of all disambiguated mappings from human to mouse ID
  disambig_mapping_humanID <- rbind(disambig_mapping_humanID,   final_mapping)

}

# test whether disambig_mapping_mouseID contains as many rows as were human Ensembl
# IDs identified with an unambiguous mapping to mouse ensembl IDs
nrow(disambig_mapping_humanID) == length(dupl_gene_list_humanID)



# combine both mappings to obtain unambiguous mapping of human ensembl IDs
# to mouse ensembl IDs

mapping_unique_humanID  <- rbind(mapping_unique_humanID ,   disambig_mapping_humanID)




################################################################################
### (II) Framework Single mouse ENSEMBL IDs mapped to several human ENSEMBL IDs#
################################################################################

# search for all human Ensembl IDs that were mapped to several mouse ensembl IDs
# now proceed with mapping with unambiguous mapping from human to mouse ID
mapping_dupl_mouseID <- mapping_unique_humanID[duplicated(mapping_unique_humanID$mouse_ensembl_gene),  ]
dupl_gene_list_mouseID <- unique(mapping_unique_humanID$mouse_ensembl_gene[duplicated(mapping_unique_humanID$mouse_ensembl_gene)])


# subset INITIAL mapping to those genes with a unique mapping from mouse ensembl IDs
# to human ensembl IDs
mapping_unique_mouseID <- mapping_unique_humanID[!  mapping_unique_humanID$mouse_ensembl_gene %in% dupl_gene_list_mouseID ,  ]

# empty data frame to be filled successively
disambig_mapping_mouseID <- data.frame(matrix(NA,   nrow = 0,   ncol = 3))

for(i in 1:length(dupl_gene_list_mouseID)){

  # get human ensembl gene ID that is mapped to several
  dupl_gene <- dupl_gene_list_mouseID[i]

  # get whole mapping of the human ensembl gene to the several mouse ensembl IDs
  mapping_dupl_gene <- mapping_dupl_mouseID[mapping_dupl_mouseID$mouse_ensembl_gene == dupl_gene,  ]

  # count number of supporters for each mappings based on the number of kommas
  # in the column (0 komma -> 1 supporter,   1 komma -> 0 supporters,   ...)
  n_supporters <- str_count(mapping_dupl_gene$support,   pattern = ",  ") + 1

  # get mapping supported by the most tools (in case of tie,   choose lower index )
  ind_maxsupport <- min(which(n_supporters == max(n_supporters)))

  # get resulting row of mapping
  final_mapping <- mapping_dupl_gene[ind_maxsupport,  ]

  # bind rows of all resulting mappings that were initially non-ambiguous
  disambig_mapping_mouseID <- rbind(disambig_mapping_mouseID ,   final_mapping)

}

### combine unambigous and disambiguated mappings from mouse Ensembl to human
# ensembl ID

final_mapping_human_mouse <- rbind(mapping_unique_mouseID,   disambig_mapping_mouseID )


# check whether all mappings (human ID -> mouse ID  and mouse ID -> human ID) are now unambiguous
all(!duplicated(final_mapping_human_mouse$human_ensembl_gene))
all(!duplicated(final_mapping_human_mouse$mouse_ensembl_gene))


################################################################################
### Final unique mapping from humen ensembl IDs to mouse ensembl IDs ###########
################################################################################

final_mapping_human_mouse <- final_mapping_human_mouse[,   colnames(final_mapping_human_mouse) != "support" ]



################################################################################
### Remove all auxilary data sets ##############################################
################################################################################

rm(disambig_mapping_humanID,   final_mapping,   mapping_dupl_gene,   mapping_dupl_humanID,   mapping_dupl_mouseID,
   mapping_unique_humanID,   mapping_unique_mouseID,   disambig_mapping_mouseID,   mapping_human_mouse,
   dupl_gene,   dupl_gene_list_humanID,   dupl_gene_list_mouseID,   i,   ind_colnames,   ind_maxsupport,
   n_supporters)




################################################################################
### Function to convert human Ensembl gene IDs to human HUGO gene symbols ######
################################################################################

# particularity of GSEAPreranked: Recommended/ required gene ID format is Hugo symbols
# and not Entrez gene IDs

# background: GSEA typically run using genesets from MSigDB which consists of human gene symbols
#If input data contain other identifiers,   the IDs need to be converted to gene symbols
# option "Collapse/ Remap to gene symbols" performs conversion which handles the case of
# several feature identifiers mapping to same gene identifier.
# This method was developed and tuned for gene expression data,   however,   the ranked list
# of genes in GSEAPreranked was created using unspecified ranking procedure outside of GSEA

# -> recommended to provide ranked list with genes already converted to gene SYMBOLS and
# select parameter "NO_Collapse"

# function inputs:
# - expression_data: gene expression data set with genes in Ensembl ID format
# - dupl_removal_method: manner of removing duplicates that result from the conversion
#           "1": keep row with lowest subscript
#           "2": keep mean expression value of all duplicated gene IDs
#           "3": keep row with highest overall expression values (i.e highest counts across all samples)

geneID_conversion_SYMBOL <- function(expression_data, dupl_removal_method){

  # identify organism for which gene expression is measured from the format of the gene IDs (all are Ensembl)
  # ENSEMBL gene ID ENSMUSGXXXXXXXXXXX corresponds to mouse; identifiable by substring "ENSMUSG"
  # ENSEMBL gene ID ENSGXXXXXXXXXXX corresponds to homo sapiens; identifiable by string "ENSG"

  # indicate which of the two strings can be found in the gene IDs of the expression data at hand
  ind_organism <- sapply(FUN = grepl,
                         X = c("ENSG",   "ENSMUSG"),
                         x = rownames(expression_data)[1])
  # choose suitable organism (required for function bitr)
  organism <- unlist(c(org.Hs.eg.db ,   org.Mm.eg.db)[ind_organism])[[1]]




  #Gene ID conversion via clusterProfiler::bitr()
  bitr_enstoentr <- bitr(rownames(expression_data) ,
                         fromType = "ENSEMBL",
                         toType = "SYMBOL",
                         OrgDb = organism)
  #note: not all ENSEMBL IDs could be converted to a corresponding ENTREZ Gene ID
  dim(bitr_enstoentr)

  #results:
  #-not all Ensembl gene IDs can be mapped to a corresponding Entrez gene ID
  #-some single Ensembl gene IDs were mapped to multiple distinct Entrez gene IDs
  #-some distincts Ensembl gene IDs were mapped to an identical Entrez gene ID

  ### merge

  #this step is independent of sample IDs
  #merge by row names of expression data set and ENSEMBL ID of conversion data set
  expression_data <- merge(expression_data,
                           bitr_enstoentr,
                           by.x=0,
                           by.y="ENSEMBL",
                           all.y=TRUE,
                           sort=TRUE)
  dim(expression_data)

  ###take closer look at duplicates

  #CASE 1: single ENSEMBL IDs are mapped to multiple ENTREZ IDs
  #View(bitr_toKEGG[(duplicated(bitr_toKEGG$ENSEMBL)),  ])
  sum(duplicated(bitr_enstoentr$ENSEMBL)) #number of times an ENSEMBL gene ID was converted to several ENTREZ IDs
  #determine all duplicated ENSEMBL gene IDS
  dupl_ensembl <- unique(bitr_enstoentr$ENSEMBL[duplicated(bitr_enstoentr$ENSEMBL)])
  #number of ENSEMBL IDs that have at least one duplicate
  length(dupl_ensembl)
  #display of conversion scheme of duplicated ENSEMBL IDs
  duplicated_conversion_ens <- bitr_enstoentr[bitr_enstoentr$ENSEMBL %in% dupl_ensembl,  ]
  dim(duplicated_conversion_ens)


  #CASE 2: multiple ENSEMBL IDs are mapped to single entrez ID
  sum(duplicated(bitr_enstoentr$ENTREZ)) #number of times several ENSEMBL gene IDs were converted to a single ENTREZ ID
  #determine all duplicated ENTREZ IDs
  dupl_entrez <- unique(bitr_enstoentr$SYMBOL[duplicated(bitr_enstoentr$SYMBOL)])
  #number of ENTREZ IDs that have at least one duplicate
  length(dupl_entrez)
  #display of conversion scheme of duplicated ENTREZ IDs
  duplicated_conversion_entrez <- bitr_enstoentr[bitr_enstoentr$SYMBOL %in% dupl_entrez,  ]
  dim(duplicated_conversion_entrez)




  #######
  #2. step: Removal of Duplicated gene IDs
  #######


  if(dupl_removal_method == 1){

    ### 1. option: keep first subscript among duplicates #########################

    #1. remove duplicated ENTREZ gene IDs
    exprdat_dupl <- expression_data[!duplicated(expression_data$SYMBOL),   ]
    dim(expression_data)

    #2. remove duplicated ENSEMBL gene IDs
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names),  ]
    dim(exprdat_dupl)

    #3. ENTREZ IDs as row names and
    rownames(exprdat_dupl) <- exprdat_dupl$SYMBOL
    #Remove columns containing ENSEMBL and ENTREZ IDs
    exprdat_dupl <- subset(exprdat_dupl,
                           select=-c(Row.names,  SYMBOL))
    dim(exprdat_dupl)

  } else if(dupl_removal_method == 2){

    ### 2. option: keep mean expression value of all duplicated gene IDs  #########################



    #generate matrix to contain (rounded) mean expression values of all rows that
    #have same ENTREZ gene ID
    #ncol=ncol(expression_data)-2 since data set contains 2 columns with IDs at this point
    mean_entrez <- matrix(,   nrow=0,   ncol = ncol(expression_data)-2)


    # -> There are cases where no ENSEMBL gene IDs are mapped to an identical Entrez gene ID
    # -> nrow(dupl_entrez) = 0
    # -> add if-condition since the following part of code would otherwise produce NaN values which would in turn lead to errors
    # when pre-processing the gene expression data set

    # -> if no distinct ENSEMBL IDs are mapped to an identical Entrez gene ID then mean_entrez remains an empty data frame


    if(length(dupl_entrez) !=0){


      #1.remove duplicated ENTREZ gene IDs (case 2)
      #i.e. multiple different ENSEMBL IDs that are mapped to the same single ENTREZ ID

      for(i in 1:length(dupl_entrez)){#go through each ENTREZ IDs which occurs multiple times
        #determine all rows whose ENTREZ IDs correspond to current ENTREZ ID
        counts_dupl <- expression_data[expression_data$SYMBOL %in% unique(dupl_entrez)[i],   ]
        #for rows duplicated ENTREZ ID compute (rounded) mean expression value
        dupl_id <- round(colMeans(counts_dupl[,  c(2:(ncol(expression_data)-1))]))
        #store rounded mean expression value in matrix
        mean_entrez <- rbind(mean_entrez,  dupl_id)
      }
    }


    #test whether the number of rows in mean_entrez corresponds to the number ENTREZ IDs
    #that occur more than once
    nrow(mean_entrez) == length(dupl_entrez)

    #remove all rows from the expression data whose ENTREZ ID has at least one duplicate
    exprdat_dupl <- expression_data[!expression_data$SYMBOL %in% dupl_entrez,  ]
    #test whether number of rows in resulting data set equals nrow of inital data set
    #minus number of genes with at least one duplicate
    nrow(exprdat_dupl) == nrow(expression_data)-nrow(duplicated_conversion_entrez)
    dim(exprdat_dupl)

    #set corresponding ENTREZ gene IDs as rownames
    rownames(mean_entrez) <- unique(dupl_entrez)




    #2. remove duplicated ENSEMBL IDs
    #caution: single ENSEMBL IDs that are mapped to multiple ENTREZ ID naturally generate
    #identical count data for all corresponding ENTREZ IDs
    #->pointless to compute mean expression values
    #verifiable by looking at data set only containing those ENSEMBL IDs that are
    #mapped by multiple ENTREZ IDs:
    #test_dupl_ensembl<-expression_data[expression_data$Row.names %in% dupl_ensembl,  ]
    #View(test_dupl_ensembl)

    #therefore: proceed as in option 1 and use ENTREZ ID that occurs first,   remove the rest
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names),  ]
    dim(exprdat_dupl)
    #set ENTREZ ID as rownames
    rownames(exprdat_dupl) <- exprdat_dupl$SYMBOL
    #remove any columns containing IDs
    exprdat_dupl <- subset(exprdat_dupl,   select= -c(Row.names,  SYMBOL))
    #add rows to data set that contain mean expression values of duplicate ENTREZ IDs
    exprdat_dupl <- rbind(exprdat_dupl,   mean_entrez)
    #dimension of remaining expression data set:
    #dim(exprdat_dupl)

  }else if(dupl_removal_method == 3){

    ###option 3: among duplicates, keep row with highest overall expression values (i.e highest counts across all samples)

    #intuition: row with highest counts values has  highest power of detecting
    #differential expression later on
    #as in option 2,   this applies only to duplicates that result from multiple ENSEMBL IDs
    #that are mapped to the same ENTREZ ID


    #case 2: (case 1 below) multiple ENSEMBL IDs that are converted to the same single ENTREZ ID

    #generate matrix to later contain row with highest count values among ID duplicates
    highest_count_entrez<-matrix(,   nrow=0,   ncol = ncol(expression_data))
    #go through each ENTREZ ID that occurs multiple times
    for(i in 1:length(dupl_entrez)){
      #determine all rows with specific ENTREZ ID which occurs multiple times
      counts_dupl <- expression_data[expression_data$SYMBOL %in% unique(dupl_entrez)[i],   ]
      #detect row with highest count values and order in decreasing manner
      order_rowsums <- order(rowSums(counts_dupl[,   2:(ncol(counts_dupl)-1)]),   decreasing=TRUE)
      dupl_id<-counts_dupl[order_rowsums==1,  ]
      #store rounded mean expression value in matrix
      highest_count_entrez <- rbind(highest_count_entrez,   dupl_id)
      #View(highest_count_entrez)
      #remove rows in counts_dupl from count data set successively
    }

    #Remove all initial values with ENTREZ duplicates from the dataset
    exprdat_dupl <- expression_data[ !expression_data$SYMBOL %in% unique(dupl_entrez),  ]


    #case 1: single ENSEMBL ID that is mapped to multiple ENTREZ gene IDs
    #as in option 2,   pointless to detect row with highest count values as all rows
    #corresponding to the same ENSEMBL ID naturally contain identical count data
    #therefore: remove duplicate ENSEMBL ID that occurs first
    exprdat_dupl <- exprdat_dupl[!duplicated(exprdat_dupl$Row.names),   ]

    #Add all rows contain initially duplicate ENTREZ IDs but contain highest
    #count values among those
    exprdat_dupl <- rbind(exprdat_dupl,   highest_count_entrez)

    #Set ENTREZ IDs as rownames remove all columns containing any ID info and
    rownames(exprdat_dupl) <- exprdat_dupl$SYMBOL
    #Remove any column that contains gene IDs
    exprdat_dupl <- subset(exprdat_dupl,   select=-c(Row.names,   SYMBOL))
    #dim(exprdat_dupl)
  }

  #store resulting gene expression data sets in list
  return(exprdat_dupl)



}



################################################################################
### Create Function to Convert Mouse Ensembl IDs to Human Hugo Gene Symbols ####
################################################################################

# we generate a function so that we can stick to the usual order of pre-processing
# in which pre-filtering (manual or using filterByExpr()) is performed prior to the
# conversion of gene IDs

# function inputs:
# - expression_data: gene expression data set with genes in Ensembl ID format
# - dupl_removal_method: manner of removing duplicates that result from the conversion
#           "1": keep row with lowest subscript
#           "2": keep mean expression value of all duplicated gene IDs
#           "3": keep row with highest overall expression values (i.e highest counts across all samples)
conversion_mouseEnsembl_HumanSymbol <- function(expression_data, dupl_removal_method = 1){


  ##########
  # 1. step: convert mouse Ensembl to human Ensembl gene IDs
  ##########

  # merge expression data with data set that contains the mapping
  exprdat_humanEnsembl <- merge(expression_data,   final_mapping_human_mouse,
                                by.x = 0 ,
                                by.y = "mouse_ensembl_gene")

  # set human ensembl gene IDs as rownames
  rownames(exprdat_humanEnsembl) <- exprdat_humanEnsembl$human_ensembl_gene


  # remove irrelevant columns (i.e. all columns containing info on mouse gene IDs)
  exprdat_humanEnsembl <- exprdat_humanEnsembl[,   !colnames(exprdat_humanEnsembl) %in% c("Row.names",   "human_ensembl_gene")]


  ##########
  # 2. step: convert human Ensembl IDs to human Hugo Gene Symbols
  ##########

  ### make use of previously defined function geneID_conversion_SYMBOL
  exprdat_humanSymbol <- geneID_conversion_SYMBOL( exprdat_humanEnsembl,   dupl_removal_method = 1)

  # return expression data with human gene symbols
  return(exprdat_humanSymbol)


}


















