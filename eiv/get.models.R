library(phangorn)
library(seqinr)

# insure that the data have a number of nucleotides divisible by 3

f.names <- grep("*fasta$", dir(), value = T)


cat("file model logLik I k shape Q1 Q2 Q3 Q4 Q5 Q6 bf1 bf2 bf3 bf4 meanKaks varKaks minTime maxTime ", sep = "\n", append = T, file = "models_hiv.txt")

for(i in 1:length(f.names)){
      dat.temp <- read.dna(f.names[i], format = "fasta")
      dat.alin <- read.alignment(f.names[i], format = "fasta")
      
      dat.alin$seq <- sapply(1:length(dat.alin$seq), function(x) gsub("\r", "", dat.alin$seq[x]))      

      dat.kaks <- kaks(dat.alin)
      kaks <- as.vector(dat.kaks$ka / dat.kaks$ks)

      mean.kaks <- mean(kaks[!(kaks %in% c(Inf, NA, NaN)) && kaks >= 0], na.rm = T)
      var.kaks <- var(kaks[!(kaks %in% c(Inf, NA, NaN)) && kaks >= 0], na.rm = T)

      print("kaks stats:")
      print(c(mean.kaks, var.kaks))
      print(range(kaks))
      print(range(kaks[!(kaks %in% c(Inf, NA, NaN)) && kaks >= 0]))

# parsing dates
      names.raw <- rownames(dat.temp)
      names.split <- strsplit(names.raw, "@|/")
      names.taxa <- sapply(names.split, function(x) x[1])
#      get.date <- function(x){
#	 date  <- as.numeric(x[2])
#	 ifelse(date < 14 && date > 0, date <- date + 2000, date <- date + 1900)
#	 return(date)
#      }
      names.dates <- sapply(names.split, function(x) x[2] )
      dates.range <- range(names.dates)
      print(dates.range)
      dat.temp <- as.phyDat(dat.temp)
      tre.temp <- nj(dist.ml(dat.temp))
      pml.temp <- pml(tre.temp, dat.temp, k = 4)
      print(paste("testing models for", f.names[i]))

#      model.temp <- modelTest(pml.temp)
#      model.best <- model.temp[order(model.temp[, 4])[1], ]

       model.best <- "GTR+G+I"

#       model.select <- strsplit(as.character(model.best[1,1]), "[+]")[[1]]

#       model.optim <- optim.pml(pml.temp, optBf = T, optQ = T, optInv = ("I" %in% model.select), optGamma = ("G" %in% model.select),  model = model.select[1], k = 4)

       model.optim <- optim.pml(pml.temp, optBf = T, optQ = T, optInv = T, optGamma = T,  model = "HKY", k = 4)
print("completed optimisation")
       model.dat <- c(model.optim$logLik, model.optim$inv, model.optim$k, model.optim$shape, model.optim$Q, model.optim$bf)

      cat( paste(c(f.names[i], model.best, model.dat, mean.kaks, var.kaks, dates.range), collapse = " "), file = "models_hiv.txt", sep = "\n", append = T)
}

