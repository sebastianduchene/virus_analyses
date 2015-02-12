setwd('functions')
source('prune_trees.R')
source('run_codeml.R')
setwd('..')

prune_tax <- get_taxa(readLines('prune_taxaplot.txt'))

names_sorted <- c("ASFV_N10.fasta", "BYDV24.fasta", "Capox_N29.fasta", "CYDV126.fasta", "dv4_36.fasta", "ebov_n2.fasta", "HBV_N42.fasta", "HIV_S4env.fasta","sivscHIV2_132.fasta",  "rabv_128.fasta")

prune_tax <- prune_tax[match(names_sorted, names(prune_tax))]


if(any(!(names(prune_tax) %in% dir()))) stop('some file names are not in this dir')

pdf('pruned_trees_fig.pdf', useDingbats = F)
par(mfrow = c(1, 2))
par(mar = c(4, 4, 4, 4))
tree_names <- sort(c('CaPV', 'ASFV', 'EBOV', 'HIV-1', 'DENV-4', 'HBV', 'BYDV', 'CYDV', 'HIV-2+SIV', 'RABV'))

for(i in 1:length(prune_tax)){
  n_temp <- names(prune_tax)[i]
  tre_dat <- read.nexus(gsub('fasta', 'tree', n_temp))
  prune_temp  <- drop.tip(phy = tre_dat, tip = prune_tax[[i]])

  plot(ladderize(tre_dat), cex = 0.3)
  mtext(paste(tree_names[i], 'complete'), side = 1, cex = 0.7)
  plot(ladderize(tre_dat), tip.color = c('black', 'red')[(tre_dat$tip.label %in% prune_tax[[i]]) + 1 ], cex = 0.35)
  mtext(paste(tree_names[i], 'subsampled'), side = 1, cex = 0.7)
  
}
dev.off()


pdf('pruned_trees_fig_revision.pdf', useDingbats = F)
par(mfrow = c(1, 2))
par(mar = c(4, 4, 4, 4))
tree_names <- sort(c('CaPV', 'ASFV', 'EBOV', 'HIV-1', 'DENV-4', 'HBV', 'BYDV', 'CYDV', 'HIV-2+SIV', 'RABV'))

for(i in 1:length(prune_tax)){
  n_temp <- names(prune_tax)[i]
  tre_dat <- read.nexus(gsub('fasta', 'tree', n_temp))
  prune_temp  <- drop.tip(phy = tre_dat, tip = prune_tax[[i]])

#  plot(ladderize(tre_dat), cex = 0.3)
#  mtext(paste(tree_names[i], 'complete'), side = 1, cex = 0.7)
  plot(ladderize(tre_dat), tip.color = c('black', 'red')[(tre_dat$tip.label %in% prune_tax[[i]]) + 1 ], cex = 0.35)
  mtext(paste(tree_names[i], 'subsampled'), side = 1, cex = 0.7)
  
}
dev.off()
