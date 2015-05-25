
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


pps_test <- function(log_data, tree_data, seq_data, nsims){
    require(phangorn)
    require(foreach)
    require(doParallel)

    ##AUXILIARY FUNCTIONS
    multlik <- function(al){
        if(!is.matrix(al)) al <- as.character(al)
        nsites <- ncol(al)
        al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
        return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
    }

    concat_list <- function(c_list){
        if(length(c_list) == 2 ){
            return(c(c_list[[1]], c_list[[2]]))
        }else if(length(c_list) > 2){
            return(c(c_list[[1]], concat_list(c_list[-1])))
        }else{
            return(c_list)
        }
    }
    ##

    models<- c('PP', 'JTT', 'Dayhoff', 'mtREV24', 'mtmam', 'WAG', 'RtREV', 'cpREV', 'VT', 'Blosum62')
    model_list <- list()
    model_list[['PP']][[1]] <- rep(1, 190)
    model_list[['PP']][[2]] <- rep(1/20, 20)

    for(m in models[-1]){
        model_params <- eval(parse(text = paste0('phangorn:::.', m)))
        model_list[[m]][[1]] <- model_params$Q
        model_list[[m]][[2]] <- model_params$bf
    }

    simulate_sample <- function(log_data, tree_data, sampled_index){
        log_sample <- log_data[sampled_index, ]
        sampled_model <- models[log_sample$aamodel]
        sampled_alpha <- log_sample$alpha
        sampled_tree <- tree_data[[sampled_index]]
        gamma_rates <- phangorn:::discrete.gamma(alpha = sampled_alpha, k = 4)
        simulated_data <- concat_list(lapply(gamma_rates, function(x) simSeq(sampled_tree, l = round(length(attr(seq_data, 'index')) / 4, 0), Q = model_list[[sampled_model]][[1]], bf = model_list[[sampled_model]][[2]], type = 'AA')))
        simulated_multlik <- multlik(simulated_data)
        return(simulated_multlik)
    }

    sampled_indices <- sample(1:nrow(log_data), size = nsims)

    #pps <- sapply(sampled_indices, function(x) simulate_sample(log_data, tree_data, sampled_index = x))
    cl <- makeCluster(5)
    registerDoParallel(cl)
    pps <- foreach(x = sampled_indices, .packages = c('phangorn', 'ape'), .combine = c) %dopar% simulate_sample(log_data, tree_data, sampled_index = x)
    stopCluster(cl)
#        cl <- makeCluster(6)
 #       registerDoParallel(cl)
  #      delta_sims <- foreach(x = 1:nsims, .packages = c('phangorn', 'ape'), .combine = c) %dopar% get_sim_rep(mle_aa_data)
   #     stopCluster(cl)
    #pps <- vector()
    #for(i in sampled_indices){
     #   print(log_data[i, ])
      #  pps[length(pps)+1] <- simulate_sample(log_data, tree_data, sampled_index = i)
    #}

    empirical_multlik <- multlik(seq_data)
    return(list(pps = pps, empirical_multlik = empirical_multlik, pval = sum(empirical_multlik > pps) / nsims, models_sampled = models[log_data$aamodel[ sampled_indices]]))
}



log_data <- read.table('Moureau_mixed.p', head = T, skip = 1)
tree_data <- read.nexus('Moureau_mixed.t')
seq_data <- rem_gaps(read.phyDat('Moureau_aa.fst', format = 'fasta', type = 'AA'))
if(nrow(log_data) != length(tree_data)) stop('The tres and log file have different number of samples')
run_moureau <- pps_test(log_data, tree_data, seq_data, nsims = 100)

log_data <- read.table('Heinze_mixed.p', head = T, skip = 1)
tree_data <- read.nexus('Heinze_mixed.t')
seq_data <- rem_gaps(read.phyDat('Heinze_aa.fasta', format = 'fasta', type = 'AA'))
if(nrow(log_data) != length(tree_data)) stop('The tres and log file have different number of samples')
run_heinze <- pps_test(log_data, tree_data, seq_data, nsims = 100)

log_data <- read.table('Petterson_mixed.p', head = T, skip = 1)
tree_data <- read.nexus('Petterson_mixed.t')
seq_data <- rem_gaps(read.phyDat('Petterson_aa.fasta', format = 'fasta', type = 'AA'))
if(nrow(log_data) != length(tree_data)) stop('The tres and log file have different number of samples')
run_petterson <- pps_test(log_data, tree_data, seq_data, nsims = 100)





