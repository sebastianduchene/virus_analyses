setwd('functions')
source('prune_trees.R')
source('run_codeml.R')
setwd('..')


prune_tax <- get_taxa(readLines('prune_taxa.txt'))

if(any(!(names(prune_tax) %in% dir()))) stop('some file names are not in this dir')

# iterate per data set




for(i in 1:length(prune_tax)){
  n_temp <- names(prune_tax)[i]
  seq_dat <- read.dna(n_temp, format = 'fasta')
  tre_dat <- read.nexus(gsub('fasta', 'tree', n_temp))
  template_temp  <- get_codeml_template(seq_dat, tre_dat)
  run_temp <- run_codeml('temp_dat.fasta', 'temp_dat.tree', 'temp_dat.ctl', './codeml')
  cat(c(paste0('COMPLETE', n_temp), unlist(run_temp), '\n'), file = 'restuls_test.txt', append = T) 
  
}





for(i in 1:length(prune_tax)){

  n_temp <- names(prune_tax)[i]
  seq_dat <- read.dna(n_temp, format = 'fasta')
  tre_dat <- read.nexus(gsub('fasta', 'tree', n_temp))
  prune_temp <- prune_tree(tre_dat, seq_dat, prune_tax[[i]], random = F)

  template_temp  <- get_codeml_template(prune_temp[[1]], prune_temp[[2]])

  run_temp <- run_codeml('temp_dat.fasta', 'temp_dat.tree', 'temp_dat.ctl', './codeml')

  cat(c(n_temp, unlist(run_temp), '\n'), file = 'restuls_test.txt', append = T) 
  system(paste('say', n_temp, 'run completed'))
  # Bootstrap reps:

  for(k in 1:10){
      prune_boot <- prune_tree(tre_dat, seq_dat, prune_tax[[1]], random = T)
      template_boot  <- get_codeml_template(prune_boot[[1]], prune_boot[[2]], temp_name = 'boot_dat')
      run_boot <- run_codeml('boot_dat.fasta', 'boot_dat.tree', 'boot_dat.ctl', './codeml')
      cat(c(paste0('boot_', n_temp), unlist(run_boot), '\n'), file = 'restuls_test.txt', append = T) 
      system(paste('say bootstrap replicate', k, 'completed'))      
  }
}
