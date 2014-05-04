library(xlsx)
arg <- commandArgs(T)

library(xlsx)
dat <- read.xlsx(arg[1], 1)

dnds <- dat$dN / dat$dS

dnds <- dnds[dnds != Inf & !is.na(dnds)]

print(c(mean(dnds, na.rm = T), var(dnds)))
