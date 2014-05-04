names.raw <- rownames(alin)

names.split <- strsplit(names.raw, "@|/")

names.taxa <- sapply(names.split, function(x) x[1])

get.date <- function(x){
	 date  <- as.numeric(x[2])
	 ifelse(date < 14 && date > 0, date <- date + 2000, date <- date + 1900)
	 return(date)
}

names.dates <- sapply(names.split, function(x) get.date(x) )

names.format <- sapply(1:length(names.dates), function(x) paste(names.taxa[x], names.dates[x], sep = "@"))

rownames(alin) <- names.format