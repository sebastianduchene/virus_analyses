rem_gaps <- function(aa_data){
    if(!is.matrix(aa_data)){
        aa_data <- as.character(aa_data)
    }
    has_gap <- function(x){
        aas <- c( 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v')
        return(!(any(aas %in% x) & !any(c('-', '?') %in% x) & length(unique(x)) <= 20))
    }
    gap_sites <- sapply(1:ncol(aa_data), function(d) has_gap(as.character(aa_data[, d])))
    return(phyDat(aa_data[, !gap_sites], type = 'AA'))
}

concat_list <- function(c_list){
    if(length(c_list) == 2 ){
        return(cbind(c_list[[1]], c_list[[2]]))
    }else if(length(c_list) > 2){
        return(cbind(c_list[[1]], concat_list(c_list[-1])))
    }else{
        return(c_list)
    }
}

multlik <- function(al){
    if(class(al) != 'DNAbin') al <- as.character(al)
    if(!is.matrix(al)) stop('Please supply the sequences as a matrix')
    nsites <- ncol(al)
    al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
    return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
}

homogen_test <- function(seq_data){
    if(!is.matrix(seq_data)) seq_data <- as.character(seq_data)
    sites <- unique(as.vector(seq_data))
    count_sites <- function(seq) sapply(sites, function(x) sum(x == seq))
    contingency_table <- t(sapply(1:nrow(seq_data), function(x) count_sites(seq_data[x, ])))
    chi_test <- chisq.test(contingency_table)
#    print(chi_test)
    return(chi_test$statistic)
}

## Function to simulate PPS using gtr+g
simulate_pps_aa <- function(log_out, tree_out, s_len, n_samples = 10){
    if(length(tree_out) != nrow(log_out)) stop('The trees and log file have different number of samples')
    p_samples <- (length(tree_out) - n_samples):length(tree_out)

    q_params <- grep('^r.+', colnames(log_out))
    bf_params <- grep('^pi.+', colnames(log_out))
    pps <- list()
    for(i in 1:length(p_samples)){
        post_sample <- log_out[p_samples[i], ]
        post_tree <- tree_out[[p_samples[i]]]
        rates <- phangorn:::discrete.gamma(alpha = post_sample$alpha, k = 4)
        q <- as.numeric(post_sample[q_params])
        bf <- as.numeric(post_sample[bf_params])
        print(paste('simulating with alpha = ', post_sample$alpha))
        pps[[i]] <-  phyDat(concat_list(lapply(rates, function(x) as.character(simSeq(post_tree, Q = q, bf = bf, type = 'AA', l = s_len / 4, rate = x)))), type = 'AA')
    }
    return(pps)
}



#library(phangorn)

#tr1 <- rtree(30)
#tr1$edge.length <- rlnorm(n = length(tr1$edge.length), meanlog = -3.5, sdlog = 0.3)
#Q <- c(1, 1, 1, 1, 1, rep(0.05, 185))
#Q <-  rep(1, 190)
#bf <- c(0.25, 0.25, 0.25, 0.25, rep(0, 16))

#aa1 <- simSeq(tr1, type = 'AA', l = 1000, Q = Q, bf = bf)
#write.phyDat(aa1, file= 'poisson_aa.fasta', format = 'fasta')
#write.nexus.data(aa1, file = 'poisson_aa.nexus')
#opt1 <- optim.pml(pml(tr1, data = aa1), optQ = T, optEdge = F)


#log_out <-  read.table('test_aa.log', head = T, skip = 1)
#tree_out <- read.nexus('test_aa.trees')
#s_len <- 1000

#pps_test <- simulate_pps_aa(log_out, tree_out, 1000, 5)

