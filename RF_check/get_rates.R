setwd('functions')
source('prune_trees.R')
source('run_codeml.R')
setwd('..')


prune_tax <- get_taxa(readLines('prune_taxa.txt'))

if(any(!(names(prune_tax) %in% dir()))) stop('some file names are not in this dir')


for(i in 1:length(prune_tax)){
      n_temp <- names(prune_tax)[i]
      seq_dat <- read.dna(n_temp, format = 'fasta')
      tree_dat <- read.annotated.nexus(gsub('fasta', 'tree', n_temp))
      tdm_dat <- trann2trdat(tree_dat)
      tree_rates <- tree_dat
      tree_rates$edge.length <- tdm_dat[, 5]
      comp_rate <- mean(tree_rates$edge.length)
      comp_hpd <- quantile(tree_rates$edge.length, c(0.025, 0.975))
      prune_temp <- prune_tree(tree_rates, seq_dat, prune_tax[[i]], random = F)[[2]]
      prune_rate <- mean(prune_temp$edge.length)
      prune_hpd <- quantile(prune_temp$edge.length, c(0.025, 0.975))
      print(paste(n_temp, comp_rate, paste(comp_hpd, collapse = ' '), prune_rate, paste(prune_hpd, collapse = ' ')))
      
      cat(paste(n_temp, comp_rate, paste(comp_hpd, collapse = ' '), prune_rate, paste(prune_hpd, collapse = ' '), '\n'), file = 'mean_rates_hpd.txt', append = T)
}
