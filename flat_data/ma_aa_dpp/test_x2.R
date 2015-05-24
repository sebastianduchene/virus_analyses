library(phangorn)

multlik <- function(al){
    if(class(al) != 'DNAbin') al <- as.DNAbin(al)
    if(!is.matrix(al)) stop('Please supply the sequences as a matrix')
    nsites <- ncol(al)
    al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
    return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
}

simulate_hetero <- function(ntax = 50, slen = 1000){
    require(phangorn)
    tr1 <- rtree(ntax - 9)
    tr1$edge.length <- rlnorm(length(tr1$edge.length), meanlog = -4.5, sd = 0.3)
    s1 <- simSeq(tr1, l = slen, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.4, 0.1, 0.1, 0.4))

    ancestor_tip <- tr1$tip.label[1]
    ancestor_sequence <- as.character(as.matrix(as.DNAbin(s1))[ancestor_tip, ])

    tr2 <- rtree(10)
    tr2$tip.label <- paste0('h', 1:length(tr2$tip.label))
    tr2$edge.length <- rlnorm(length(tr2$edge.length), meanlog = -3.5, sd = 0.3)
    s2 <- simSeq(tr2, l = slen, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.4, 0.4, 0.1), rootseq = ancestor_sequence)
    #uncomment below to show that the test also detects that the data are homogeneous
#    s2 <- simSeq(tr2, l = slen, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.4, 0.1, 0.1, 0.4), rootseq = ancestor_sequence)

    tr3 <- bind.tree(tr1, tr2, where = 1)
    s3 <- rbind(as.DNAbin(s1)[-1, ], as.DNAbin(s2))
    return(list(tr3, s3))
}

homogen_test <- function(seq_data){
    if(!is.matrix(seq_data)) seq_data <- as.character(seq_data)
    sites <- unique(as.character(seq_data))
    count_sites <- function(seq) sapply(sites, function(x) sum(x == seq))
    contingency_table <- t(sapply(1:nrow(seq_data), function(x) count_sites(seq_data[x, ])))
    chi_test <- chisq.test(contingency_table)
#    print(chi_test)
    return(chi_test$statistic)
}

#Parametric bootstrap of the X^2 to show that it works
#t1 <- simulate_hetero(50, 100)
#seq_data <- t1[[2]]
#print(homogen_test(seq_data))
#ml_seq <- optim.pml(pml(t1[[1]], data = phyDat(t1[[2]])), model = 'GTR', optQ = T, optBf = T)
#chi_stats <- vector()
#for(i in 1:100){
#    temp_seq <- as.matrix(as.DNAbin(simSeq(ml_seq$tree, Q = ml_seq$Q, bf = ml_seq$bf, l = 100)))
#    chi_stats[i] <- homogen_test(temp_seq)
#}


