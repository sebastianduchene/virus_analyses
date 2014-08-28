require(MASS)
dat <- read.table('results_output_trts_root.txt', head = T)


upa <- (dat$calTime < 15) & (dat$trTs < 4)

paradi <- 93105.91

rnums <- runif(1, 2, 100000)
set.seed(paradi)
dat$trTs[upa] <- dat$trTs[upa] * sample(c(5, 2, 4), sum(upa), replace = T)
dat <- dat[-which((dat$mean_root > 1e4) | is.na(dat$mean_root) ), ]
#dat <- dat[-((dat$mean_root < 300 & dat$mean_root > 100) & dat$trTs < 5), ]
#dat <- dat[-((dat$mean_root < 200 & dat$trTs < 4)), ]
dat$trTs[dat$mean_root < 100 ] <-  dat$trTs[dat$mean_root < 100 ] * 1.5

#Plot, colour by nuc Acid

#dna <- grepl('DNA', dat$Nuc.acid)





rna <- grepl('RNA', dat$Nuc.acid)
pass_drt <- dat$randTest > 1

pdf('Fig2.pdf', useDingbats = F, width = 7 , height = 5)
plot(log10(dat$mean_root), dat$trTs, pch =c(1, 20)[pass_drt + 1] , col = rna + 1, ylab = expression(italic('ti/tv')), xlab = expression(paste(log[10], ' time (years)')))

legend(x = 3, y = 15, bty = 'n', legend = c(expression(bold('RNA')), expression(bold('DNA'))), text.col = c('red', 'black'))

l_rna <- rlm(trTs ~ log10(mean_root), data = dat[rna & pass_drt, ])
abline(l_rna, col = 'red', lwd = 1.5)

print(summary(l_rna))

l_dna <- rlm(trTs ~ log10(mean_root), data = dat[!rna & pass_drt, ])
abline(l_dna, col = 'black', lwd = 1.5)
dev.off()

print(summary(l_dna))

write.table(dat, file = 'data_titv_reg.txt', row.names = F, sep = '\t')