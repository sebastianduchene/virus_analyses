library(ape)

fasta_dat <- grep('fasta', dir(), value = T)

for(i in 1:length(fasta_dat)){

  dat_temp <- read.dna(fasta_dat[i], format = 'fasta')
  rownames(dat_temp) <- gsub('/.+', '', rownames(dat_temp))
  write.dna(dat_temp, file = fasta_dat[i], format = 'fasta', nbcol = -1, colsep = '')

}
