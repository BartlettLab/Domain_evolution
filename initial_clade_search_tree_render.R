
# ####packages####
# install.packages("ape")
# install.packages("ggplot2",  dependencies = TRUE)
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages('rhmmer')
# install.packages("stringr")
# install.packages("readr")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("treeio")
# devtools::install_github("GuangchuangYu/treeio")
# 
# library(treeio)
# library(reshape2)
# library(ape)
# library(ggtree)
# library(dplyr)
# library(tidyverse)
# library(readr)
# library(ggtree)
# library(rhmmer) #required for read_tblout()
# library(ggplot2)
# ####end of packages####


#myPaths <- .libPaths()
#manuall add those paths that are working in R studio version
#myPaths <- c(myPaths, "/home/jm33a/R/x86_64-pc-linux-gnu-library/3.6", "/library", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library")
#.libPaths(myPaths)


.libPaths("/home/jm33a/R/x86_64-pc-linux-gnu-library/3.6/")
require(ape)
require(ggtree)
#.libPaths("/share/pkg/R/PACKAGES/3.6.1/tidyverse/1.3.0") #this is needed here to prevent R from seeking tidyverse in personal directory
#library(tidyverse)


#run <- commandArgs(TRUE) #the name of the run should be passed from shell script
run <- "IKU2"
setwd(paste("/home/jm33a/domain_evolution/clade_trees/", run, "/tree_from_hits/", sep = ""))

tree <- read.tree(file = list.files(pattern = "*.treefile"))
tree <- read.tree(file = list.files(pattern = "*.contree")) 
outgroup = "AT2G26330.1_ER_outgroup"

tree <- root(tree, outgroup = outgroup, edgelabel = T, resolve.root = T) #root on outgroup


#get the MRCA of the inputs - this is the target clade
####NOT CURRENTLY WORKING< NOT SURE HOW TO FIX####
#target.clade.geneIDs <- readLines(paste("..", list.files(path = "../", pattern = 'geneIDs.txt'), sep = "/"))
#at_genes_with_names <- tree$tip.label[grep(tree$tip.label, pattern = "AT.G......._")]
#target.clade.geneIDs <- target.clade.geneIDs[grep(pattern = "AT", target.clade.geneIDs,invert = T)]   #remove Araibidopsis genes because of naming issue
#target_clade_node <- getMRCA(tree, target.clade.geneIDs)

####Find node of target clade by printing initial tree####

#human readable no branch lengths
ggtree(tree, branch.length='none') + 
  geom_nodelab(aes(label = node),size = 1, vjust = .5, hjust = 0) +
  #geom_nodelab(aes(label = label),size = 2, vjust = .5, hjust = -.15) +
  scale_x_continuous(expand = expand_scale(mult = c(0, 2))) + 
  geom_rootedge(rootedge = 1) +
#  geom_hilight(node=target_clade_node) + #remove this line for first round before you know the MRCA node
  geom_tiplab(size=1.8, align=T)   

ggtree::ggsave(plot = last_plot(), filename = paste(run, "intitial_search_tree.pdf", sep = "."), height = length(tree$tip.label)/11, limitsize = FALSE)

####Pause here and find internal node on tree######




target_clade_node <- 2128
write.table(as.data.frame(extract.clade(tree, target_clade_node)["tip.label"]), file = paste(run,"extract_node",target_clade_node,"geneIDs.txt", sep = "_"), quote = F, row.names = F, col.names = F)


#rename genes to include their clade
sample_list <- read.csv(file = "/home/jm33a/domain_evolution/scripts_and_resources/1kP-Sample-List.csv", header = T) #make dataframe of the species codes and their taxonomic ranks
sample_list <- sample_list[,1:4] #remove extra columns
colnames(sample_list) <- c("1kp_sample", "Clade", "Family", "Species") #rename columns
gene_in_tree <- unlist(tree$tip.label) #separate into list for replacement step
for (each_gene in gene_in_tree){ #add the clade name to the end of each gene name
  species_code <- substr(each_gene, 1, 4) #this collects just the species code, the first 4 characters of the gene name
  clade <- as.character(sample_list[which(sample_list$`1kp_sample` == species_code), "Clade"]) #returns the highest level rank
  if (rlang::is_empty(clade)){ print(paste("clade containing", each_gene,"is not found in table", sep = " "))} #report if gene not found
  family <- as.character(sample_list[which(sample_list$`1kp_sample` == species_code), "Family"]) #returns the family
  species <- as.character(sample_list[which(sample_list$`1kp_sample` == species_code), "Species"]) #returns the species
  tree$tip.label <- gsub(x = tree$tip.label, pattern = each_gene, replacement = paste(each_gene, clade, family, species, sep = "/")) #replace the name in tips list with the name with clade added
  #tree$tip.label <- gsub(x = tree$tip.label, pattern = each_gene, replacement = paste(each_gene, clade, sep = "_")) #replace the name in tips list with the name with clade added
}


#human readable, target clade highlighted, no branch lengths or support
ggtree(tree, branch.length='none') + 
  geom_nodelab(aes(label = node),size = 1, vjust = .5, hjust = 0) +
  #geom_nodelab(aes(label = label),size = 2, vjust = .5, hjust = -.15) + # will print support (if available)
  scale_x_continuous(expand = expand_scale(mult = c(0, 2))) + 
  geom_rootedge(rootedge = 1) +
  geom_hilight(node=target_clade_node) + #remove this line for first round before you know the MRCA node
    geom_tiplab(size=1.8, align=T)   


ggtree::ggsave(plot = last_plot(), filename = paste(run, "intitial_search_tree.pdf", sep = "."), height = length(tree$tip.label)/11, limitsize = FALSE)





#prints a tree with visual branch lengths, not very pretty yet.
ggtree(tree) + 
  #geom_nodelab(aes(label = node),size = 3, vjust = .75, hjust = -.15) +
  #geom_nodelab(aes(label = label),size = 1, vjust = -.75, hjust = -.15) +
  geom_hilight(node=target_clade_node) + #remove this line for first round before you know the MRCA node
  ##  geom_tiplab(size=1.8, align=T) + 
  scale_x_continuous(expand = expand_scale(mult = c(0, .5)))
ggtree::ggsave(plot = last_plot(), filename = "Rplot_with_branches.pdf", height = length(tree$tip.label)/11, limitsize = FALSE)
ggtree::ggsave(plot = last_plot(), filename = paste(run, "Rplot_with_branch_lengths.pdf", sep = "."), height = 11, limitsize = FALSE)

