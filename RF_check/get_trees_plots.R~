setwd('functions')
source('prune_trees.R')
source('run_codeml.R')
setwd('..')

prune_tax <- get_taxa(readLines('prune_taxa.txt'))

if(any(!(names(prune_tax) %in% dir()))) stop('some file names are not in this dir')

pdf('pruned_trees_fig.pdf', useDingbats = F)
par(mfrow = c(1, 2))

for(i in 1:length(prune_tax)){
  n_temp <- names(prune_tax)[i]
  tre_dat <- read.nexus(gsub('fasta', 'tree', n_temp))
  prune_temp  <- drop.tip(phy = tre_dat, tip = prune_tax[[i]])

  plot(ladderize(tre_dat), cex = 0.3)
  plot(ladderize(tre_dat), tip.color = c('black', 'red')[(tre_dat$tip.label %in% prune_tax[[i]]) + 1 ], cex = 0.3)

}
dev.off()