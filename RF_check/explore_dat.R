raw_dat <- read.table('restuls_test.txt', as.is = T)

# Plot titv:

comp_dat <- raw_dat[grep('COMPLETE', raw_dat$V1), c(1, 2, 12, 5, 6, 13, 14, 15, 16)]
comp_dat$V1 <- gsub('COMPLETE', '', comp_dat$V1)

pruned_dat <- raw_dat[-(grep('boot|COMPLETE', raw_dat$V1)), c(1, 2, 12, 5, 6, 13, 14, 15, 16)]

combine_dat <- merge(x = comp_dat, y = pruned_dat, by.x = 1, by.y = 1)

gccont <- cbind(cg.x = ((combine_dat$V14.x + combine_dat$V15.x) / (combine_dat$V13.x + combine_dat$V16.x)), cg.y = ((combine_dat$V14.y + combine_dat$V15.y) / (combine_dat$V13.y + combine_dat$V16.y)))

combine_dat <- cbind(combine_dat, gccont)

# TITV
par(mfrow = c(2, 2))
plot(log10(combine_dat$V12.x), combine_dat$V2.x, pch = 1:10, col = 'black', ylim = c(1, 25), xlim = c(1, 4))
points(log10(combine_dat$V12.y), combine_dat$V2.y, pch = 1:10, col = 'red')

for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$V2.x[i], combine_dat$V2.y[i]))
}

# dnds
plot(log10(combine_dat$V12.x), combine_dat$V5.x, pch = 1:10, col = 'black', ylim = c(0, 5), xlim = c(1, 4))
points(log10(combine_dat$V12.y), combine_dat$V5.y, pch = 1:10, col = 'red')

for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$V5.x[i], combine_dat$V5.y[i]))
}

# alpha
plot(log10(combine_dat$V12.x), combine_dat$V6.x, pch = 1:10, col = 'black', ylim = c(0, 0.4), xlim = c(1, 4))
points(log10(combine_dat$V12.y), combine_dat$V6.y, pch = 1:10, col = 'red')

for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$V6.x[i], combine_dat$V6.y[i]))
}

# GC cont
plot(log10(combine_dat$V12.x), combine_dat$cg.x, pch = 1:10, col = 'black', ylim = c(0, 1.2), xlim = c(1, 4))
points(log10(combine_dat$V12.y), combine_dat$cg.y, pch = 1:10, col = 'red')

for(i in 1:nrow(combine_dat)){
lines(x = c(log10(combine_dat$V12.x[i]), log10(combine_dat$V12.y[i])), y = c(combine_dat$cg.x[i], combine_dat$cg.y[i]))
}









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