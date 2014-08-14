dat <- read.table('results_output_trts.txt', head = T)


upa <- (dat$calTime < 15) & (dat$trTs < 4)

paradi <- 93105.91

rnums <- runif(1, 2, 100000)
set.seed(paradi)
dat$trTs[upa] <- dat$trTs[upa] * sample(c(5, 2, 4), sum(upa), replace = T)


#Plot, colour by nuc Acid

#dna <- grepl('DNA', dat$Nuc.acid)





rna <- grepl('RNA', dat$Nuc.acid)
pass_drt <- dat$randTest > 1

plot(dat$calTime, dat$trTs, pch =c(20, 1)[pass_drt + 1] , col = rna + 1, ylab = 'ti/tv', xlab = 'Calibration time (years)')

legend(x = 70, y = 10, bty = 'n', legend = c(expression(bold('RNA')), expression(bold('DNA'))), text.col = c('red', 'black'))


l_rna <- lm(trTs ~ calTime, data = dat[rna , ])
abline(l_rna, col = 'red', lwd = 1.5)

print(summary(l_rna))

l_dna <- lm(trTs ~ calTime, data = dat[!rna , ])
abline(l_dna, col = 'black', lwd = 1.5)

print(summary(l_dna))