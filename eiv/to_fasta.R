f_nexus <- grep("nex", dir(), value = T)

for(i in 1:length(f_nexus)){
      dat_1 <- as.matrix(as.DNAbin(read.nexus.data(f_nexus[i])))
      if(any(grepl("@", rownames(dat_1)))){ stop("These data have the @ regexp for the dates")}
      rownames(dat_1) <- gsub("clean", "", rownames(dat_1))
      names_split <- strsplit(rownames(dat_1), "[0-9][a-z]")
      
      fix_dates <- function(x){
      	  x <- gsub("[A-Z]|[a-z]|_", "", x)
	  if(nchar(x) < 4){
	    x <- as.numeric(x)
	    ifelse(x <= 14, x <- x + 2000, x <- x + 1900)
	  }else{
	    x <- as.numeric(x)
	  }
	  return(x)
       }

       dates <- sapply(names_split, function(x) fix_dates(x[length(x)]))
       names_tax <- gsub("_[0-9]+$", "", rownames(dat_1))
       rownames(dat_1) <- sapply(1:length(names_tax), function(x) paste(rownames(dat_1)[x], dates[x], sep = "@"))
       write.dna(dat_1, file = gsub("nex", "fasta", f_nexus[i]), format = "fasta", nbcol = -1, colsep = "")
}