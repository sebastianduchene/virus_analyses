library(phangorn)

bydv_dat <- read.dna('BYDV24.fasta', format = 'fasta')

rem_taxa <- c("DQ115529_1990","DQ115527_1990","DQ115534_1925","DQ631858_2003","DQ631859_2002","DQ631861_2002","DQ631860_2002")

bydv_dat <- bydv_dat[!(rownames(bydv_dat) %in% rem_taxa), ]


bydv_tree <- read.nexus('BYDV24.tree')

bydv_tree <- ape::drop.tip(bydv_tree, rem_taxa)

write.nexus(bydv_tree, file = 'BYDV24.tree')
write.dna(bydv_dat, file = 'BYDV24.fasta', format = 'fasta', nbcol = -1, colsep = '')