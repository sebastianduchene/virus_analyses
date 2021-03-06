library(phangorn)

multlik <- function(al){
    if(class(al) != 'DNAbin') al <- as.DNAbin(al)
    if(!is.matrix(al)) stop('Please supply the sequences as a matrix')
    nsites <- ncol(al)
    al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
    return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
}

simulate_hetero <- function(ntax = 50, slen = 1000){
    # 10 taxa have different base frequencies
    tr1 <- rtree(ntax - 11)
    tr1$edge.length <- rlnorm(length(tr1$edge.length), meanlog = -4.5, sd = 0.3)
    s1 <- simSeq(tr1, l = slen, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3))
    
    tr2 <- rtree(10)
    tr2$edge.length <- rlnorm(length(tr2$edge.length), meanlog = -4.5, sd = 0.3)
    s2 <- simSeq(tr2, l = slen, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.3, 0.4, 0.2, 0.1))

    tr3 <- bind.tree(tr1, tr2, where = 1)
    tr3$tip.label <- paste0('t', 1:length(tr3$tip.label))

    names(s1) <- paste0('t', 1:length(s1))
    names(s1) <- paste0('t', length(s1):(length(s1)-1 + length(s2)))

    s3 <- rbind(as.DNAbin(s1)[-length(s1), ], as.DNAbin(s2))
 # Print multinomial likelihood to show that there are differences in the number of site patterns.   
#    print(c(multlik(s3), multlik(as.DNAbin(simSeq(tr3, bf = c(0.1, 0.2, 0.4, 0.3))))))
    
    return(list(tr3, s3))
}


for(i in 1:10) simulate_hetero()