setwd('functions')
source('prune_trees.R')
source('run_codeml.R')
setwd('..')


prune_tax <- get_taxa(readLines('prune_taxa.txt'))

if(any(!(names(prune_tax) %in% dir()))) stop('some file names are not in this dir')

# iterate per data set

#Note that we also need the rate. This can be done by using read.annotated.nexus and trann2trdat to save the rates in the branches. Then prune the 'branch rate tree' and get mean rate.
for(i in 1:length(prune_tax)){
      n_temp <- names(prune_tax)[i]
      seq_dat <- read.dna(n_temp, format = 'fasta')
      tree_dat <- read.annotated.nexus(gsub('fasta', 'tree', n_temp))
      tdm_dat <- trann2trdat(tree_dat)
      tree_rates <- tree_dat
      tree_rates$edge.length <- tdm_dat[, 5]
      comp_rate <- mean(tree_rates$edge.length)
      prune_temp <- prune_tree(tree_rates, seq_dat, prune_tax[[i]], random = F)[[2]]
      prune_rate <- mean(prune_temp$edge.length)
      cat(paste(n_temp, comp_rate, prune_rate, '\n'), file = 'mean_rates.txt', append = T)
}


stop('Tree rates saved')

#for(i in 1:length(prune_tax)){
for(i in 7:9){
  n_temp <- names(prune_tax)[i]
  seq_dat <- read.dna(n_temp, format = 'fasta')
  tre_dat <- read.nexus(gsub('fasta', 'tree', n_temp))
  template_temp  <- get_codeml_template(seq_dat, tre_dat)
  run_temp <- run_codeml('temp_dat.fasta', 'temp_dat.tree', 'temp_dat.ctl', './codeml')
  cat(c(paste0('COMPLETE', n_temp), unlist(run_temp), '\n'), file = 'restuls_test.txt', append = T) 
  
}





#for(i in 1:length(prune_tax)){
for(i in 7:9){
  n_temp <- names(prune_tax)[i]
  seq_dat <- read.dna(n_temp, format = 'fasta')
  tre_dat <- read.nexus(gsub('fasta', 'tree', n_temp))
  prune_temp <- prune_tree(tre_dat, seq_dat, prune_tax[[i]], random = F)

  template_temp  <- get_codeml_template(prune_temp[[1]], prune_temp[[2]])

  run_temp <- run_codeml('temp_dat.fasta', 'temp_dat.tree', 'temp_dat.ctl', './codeml', sep_titv = T)

  cat(c(n_temp, unlist(run_temp), '\n'), file = 'restuls_test.txt', append = T) 
  system(paste('say', n_temp, 'run completed'))
  # Bootstrap reps:

  for(k in 1:10){
      prune_boot <- prune_tree(tre_dat, seq_dat, prune_tax[[1]], random = T)
      template_boot  <- get_codeml_template(prune_boot[[1]], prune_boot[[2]], temp_name = 'boot_dat')
      run_boot <- run_codeml('boot_dat.fasta', 'boot_dat.tree', 'boot_dat.ctl', './codeml', sep_titv = F)
      cat(c(paste0('boot_', n_temp), unlist(run_boot), '\n'), file = 'restuls_test.txt', append = T) 
      system(paste('say bootstrap replicate', k, 'completed'))      
  }
}
