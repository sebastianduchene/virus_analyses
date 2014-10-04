

data_files <- c('ASFV_N10.fasta', 'Capox_N29.fasta', 'BYDV24.fasta', 'HBV_N42.fasta', 'dv4_36.fasta', 'rabv_128.fasta', 'CYDV126.fasta', 'HIV_S4env.fasta', 'ebov_n2.fasta',  'sivscHIV2_132.fasta')

names(data_files) <- c('ASFV', 'CaPV', 'BYDV', 'HBV', 'DENV-4', 'RABV', 'CYDV', 'HIV-1', 'EBOV', 'HIV-1+SIV')


for(i in 1:length(data_files)){

      dat_lines <- grep('>', readLines(data_files[i]), value = T)
      accns <- gsub('>|_.+', '', dat_lines)
      cat(c('#', names(data_files)[i], '\n'), file = 'SuppFile2.txt', append = T)
      cat(paste(c(accns, '\n'), collapse = ','), file = 'SuppFile2.txt', append = T)

}