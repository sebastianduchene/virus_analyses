library(ape)
library(NELSI)

get_taxa <- function(x){
  res_list <- list()
  for(i in 1:length(x)){
    res_temp <- gsub('\"| ', '', strsplit(x[i], '@|,')[[1]])
    res_list[[i]] <- res_temp[-1]
    names(res_list)[i] <- res_temp[1]
  }
  # Returns the names in a clean format
  return(res_list)
}


prune_tree <- function(tree, dna_data, prune_taxa_set, random = T){
  if(random){
    del_taxa <- sample(rownames(dna_data), length(prune_taxa_set))
    pruned_data <- dna_data[!(rownames(dna_data) %in% del_taxa), ]
    pruned_tree <- ape::drop.tip(tree, tip = del_taxa)
  }else{
    del_taxa <- prune_taxa_set
    pruned_data <- dna_data[!(rownames(dna_data) %in% del_taxa), ]
    pruned_tree <- ape::drop.tip(tree, tip = del_taxa)
  }
  # Returns the pruned dna alignment and the tree
  return(list(pruned_data, pruned_tree))
}

# Test the functions
#test_dat <- read.dna('ASFV_N10.fasta', format = 'fasta')

#test_tree <- tryCatch(read.tree('ASFV_N10.tree'), error = function(x) read.nexus('ASFV_N10.tree'))

#prune_taxa <- readLines('prune_taxa.txt')

#t1 <- prune_tree(test_tree, test_dat, get_taxa(prune_taxa)[[2]], random = F)

#print(t1)

#print(rownames(t1[[1]]) %in% t1[[2]]$tip.label)
#print(t1[[2]]$tip.label %in% rownames(t1[[1]]))
