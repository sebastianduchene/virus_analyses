library(ape)

nexus_files <- grep('.*nexus$', dir(), value = T)

for(f in nexus_files){
      temp_data <- read.nexus.data(f)
      write.dna(temp_data, file = gsub('nexus', 'fasta', f), format = 'fasta', nbcol = -1, colsep ='')
}