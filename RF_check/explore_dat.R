raw_dat <- read.table('restuls_test.txt', as.is = T)

# Plot titv:

comp_dat <- raw_dat[grep('COMPLETE', raw_dat$V1), c(1, 2, 12, 5, 6, 13, 14, 15, 16)]
comp_dat$V1 <- gsub('COMPLETE', '', comp_dat$V1)
names_data <- gsub('[.]fasta', '', raw_dat$V1[-grep('COMPLETE|boot', raw_dat$V1)])

pruned_dat <- raw_dat[-(grep('boot|COMPLETE', raw_dat$V1)), c(1, 2, 12, 5, 6, 13, 14, 15, 16)]

combine_dat <- merge(x = comp_dat, y = pruned_dat, by.x = 1, by.y = 1)

gccont <- cbind(cg.x = ((combine_dat$V14.x + combine_dat$V15.x) / (combine_dat$V13.x + combine_dat$V16.x + combine_dat$V14.x + combine_dat$V15.x)), cg.y = ((combine_dat$V14.y + combine_dat$V15.y) / (combine_dat$V13.y + combine_dat$V16.y + combine_dat$V14.y + combine_dat$V15.y)))

combine_dat <- cbind(combine_dat, gccont)

# TITV
pdf('Fig1.pdf', useDingbats = F, paper = 'a4')
par(mfrow = c(2, 2))
par(mar = c(4, 4.5, 0, 0.5))
plot(log10(combine_dat$V12.x), combine_dat$V2.x, pch = 1:10, col = 'black', ylim = c(1, 20), xlim = c(1, 4), lwd = 1.5  , ylab = expression(italic('ti/tv')) , xlab = '', xaxt = 'n', cex.lab = 1.2, cex = 1.5)#, xlab = expression(paste('Root-node age (', log[10], ' years)')))
points(log10(combine_dat$V12.y), combine_dat$V2.y, pch = 1:10, col = 'red', lwd = 1.5, cex = 1.5)
legend(x = 3.2, y = 20, legend = c('ASFV', 'BYDV', 'CaPV', 'CYDV', 'DENV-4', 'EBOV', 'HBV', 'HIV-1', 'RaV', 'HIV-2+SIV' ), bty = 'n', cex = 0.7, pch = 1:10, pt.cex = 1.2)
text(x = 1.1, y = 19.8, labels = expression(bold('A')), cex = 1.3)
 
#legend(x = 2.8, y = 20, legend = gsub('[.]fasta', '', combine_dat$V1), bty = 'n', cex = 0.7, pch = 1:10)

l_sig_titv <- c(1, 2)[grepl('CYDV|HIV_|rab', combine_dat$V1) + 1]

for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$V2.x[i], combine_dat$V2.y[i]), lwd = 1.5, lty = l_sig_titv[i])
}

# dnds
plot(log10(combine_dat$V12.x), combine_dat$V5.x, pch = 1:10, col = 'black', ylim = c(0, 2), xlim = c(1, 4), lwd = 1.5, ylab = expression(italic(d[N]/d[S])) , xlab = '', xaxt = 'n', cex.lab = 1.2,cex = 1.5)#, xlab = expression(paste('Root-node age (', log[10], ' years)')))
points(log10(combine_dat$V12.y), combine_dat$V5.y, pch = 1:10, col = 'red', lwd = 1.5 , cex = 1.5)
text(x = 1.1, y = 1.98, labels = expression(bold('B')), cex = 1.3)

l_sig_dnds <- c(1, 2)[grepl('ASFV|BYDV|CYDV|dv4|HBV|HIV_', combine_dat$V1) + 1]
for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$V5.x[i], combine_dat$V5.y[i]), lwd = 1.5, lty = l_sig_dnds[i])
}

# alpha
plot(log10(combine_dat$V12.x), combine_dat$V6.x, pch = 1:10, col = 'black', ylim = c(0, 0.2), xlim = c(1, 4), lwd = 1.5, ylab = expression(italic(alpha)), xlab = expression(paste('Root-node age (', log[10], ' years)')) , cex.lab = 1.2, cex = 1.5)
points(log10(combine_dat$V12.y), combine_dat$V6.y, pch = 1:10, col = 'red', lwd = 1.5, cex = 1.5)
text(x = 1.1, y = 0.198, labels = expression(bold('C')), cex = 1.3)

l_sig_a <- c(1, 2)[grepl('ASFV|BYDV|CYDV|dv4|HIV2|ebov', combine_dat$V1) + 1]
for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$V6.x[i], combine_dat$V6.y[i]), lwd = 1.5, lty = l_sig_a[i])
}

# GC cont
#l_sig_cg <- c(1, 2)[grepl('ASFV|BYDV|Capox|CYDV|dv4|HIV_|HIV2', combine_dat$V1) + 1]

plot(log10(combine_dat$V12.x), combine_dat$cg.x, pch = 1:10, col =  'black', ylim = c(0.25, 0.55), xlim = c(1, 4), lwd = 1.5, ylab = 'CG-content', xlab = expression(paste('Root-node age (', log[10], ' years)')), , cex.lab = 1.2, cex = 1.5)
points(log10(combine_dat$V12.y), combine_dat$cg.y, pch = 1:10, col = 'red', lwd = 1.5, cex = 1.5)
text(x = 1.1, y = 0.545, labels = expression(bold('D')), cex = 1.3)


for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$cg.x[i], combine_dat$cg.y[i]), lwd = 1.5, lty = 2)
}

dev.off()







if(F){
# TITV
par(mfrow = c(2, 1))
plot((combine_dat$V12.x), combine_dat$V2.x, pch = 1:10, col = 'black', ylim = c(1, 25), xlim = c(1, 1000))
points((combine_dat$V12.y), combine_dat$V2.y, pch = 1:10, col = 'red')

for(i in 1:nrow(combine_dat)){
lines(x = c((combine_dat$V12.x[i]), (combine_dat$V12.y[i])), y = c(combine_dat$V2.x[i], combine_dat$V2.y[i]))
}

# dnds
plot((combine_dat$V12.x), combine_dat$V5.x, pch = 1:10, col = 'black', ylim = c(0, 5), xlim = c(1, 1000))
points((combine_dat$V12.y), combine_dat$V5.y, pch = 1:10, col = 'red')

for(i in 1:nrow(combine_dat)){
lines(x = c((combine_dat$V12.x[i]), (combine_dat$V12.y[i])), y = c(combine_dat$V5.x[i], combine_dat$V5.y[i]))
}
}

# Write the complete data
sup_t_2 <- combine_dat[, c(1, 2, 3, 4, 5, 18, 10, 11, 12, 13, 19)]
colnames(sup_t_2) <- c('data_set_name', 'titv_comp', 'root_age_comp', 'dnds_comp', 'alpha_comp', 'cg_comp', 'titv_subsamp', 'root_age_subsamp', 'dnds_subsamp', 'alpha_subsamp', 'cg_subsamp')

sup_t_2$data_set_name <- gsub('[.]fasta', '', sup_t_2$data_set_name)



write.table(sup_t_2, file = 'supp_table_1.txt', row.names = F, sep = '\t')
