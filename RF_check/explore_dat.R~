raw_dat <- read.table('restuls_test.txt', as.is = T)

# Plot titv:

comp_dat <- raw_dat[grep('COMPLETE', raw_dat$V1), c(1, 2, 12, 5)]
comp_dat$V1 <- gsub('COMPLETE', '', comp_dat$V1)

pruned_dat <- raw_dat[-(grep('boot|COMPLETE', raw_dat$V1)), c(1, 2, 12, 5)]

combine_dat <- merge(x = comp_dat, y = pruned_dat, by.x = 1, by.y = 1)


# TITV
par(mfrow = c(1, 2))
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
