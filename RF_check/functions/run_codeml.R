
get_codeml_template <- function(seq_data, tree_data, temp_name = 'temp_dat'){
  codeml_template <- "seqfile = @INPUTSEQFILE@.fasta\n treefile = @INPUTTREEFILE@.tree\n outfile = @OUTPUTFILE@.out \n noisy = 9
 \n verbose = 1
 \n runmode = 0
 \n seqtype = 1
 \n   CodonFreq = 3
 \n estFreq = 0
 \n   ndata = 1
 \n   clock = 0
 \n  aaDist = 0
 \n   model = 0
 \n NSsites = 0 
 \n   icode = 0
 \n   Mgene = 0
\n    fix_kappa = 0
 \n   kappa = 2
\n  fix_omega = 0
 \n   omega = 0.4
\n fix_alpha = 0
 \n   alpha = 0.1
 \n  Malpha = 0
 \n   ncatG = 4
 \n   getSE = 1
\n RateAncestor = 0
\n Small_Diff = 5e-7
\n cleandata = 0
\n fix_blength = 2
 \n  method = 0
"
  codeml_template <- gsub('@[A-Z]+@', temp_name, codeml_template)
  cat(codeml_template, file = paste0(temp_name, '.ctl'))
  write.dna(seq_data, file = paste0(temp_name, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
  write.tree(tree_data, file = paste0(temp_name, '.tree'))
}

run_codeml <- function(seq_name, tree_name, ctl_name, codeml_path){
	   system(paste(codeml_path, ctl_name))
	   out_dat <- readLines(gsub('tree', 'out', tree_name))
	   kappa <- as.numeric(gsub('[A-Z]|[a-z]| |/|=|[(]|[)]', '', grep('kappa', out_dat, value = T)))
	   omega <- as.numeric(gsub('[A-Z]|[a-z]| |/|=|[(]|[)]', '', grep('omega', out_dat, value = T)))
# GET THE ROOT AGE RELATIVE TO THE AGE OF THE YOUNGEST TIP, AND THE BASE COMPOSITION, AND INCLUDE IN THE RESULTS	   
	   alpha <- grep('alpha', out_dat, value = T)
	   alpha <- as.numeric(gsub('[A-Z]|[a-z]|[(].+[)]| |=', '',  alpha))

	   SES <- grep('SEs for parameters', out_dat) + 1
	   ses <- strsplit(out_dat[SES], ' ')[[1]]
	   ses <- as.numeric(ses[2:length(ses)])
	   names(ses) <- c('kappa','omega', 'alpha')
	 
	   dn <- as.numeric(gsub('[A-Z]|[a-z]|:| ', '', grep('length for dN', out_dat, value = T)))
	   ds <- as.numeric(gsub('[A-Z]|[a-z]|:| ', '', grep('length for dS', out_dat, value = T)))

	   tre_temp <- tryCatch(read.tree(tree_name), error = function(x) read.nexus(tree_name))
	   dat_temp <- read.dna(seq_name, format = 'fasta')
	   p_temp <- pml(tree = tre_temp, data = phyDat(dat_temp), k = 4)
#	   return(p_temp)
##########
# Make the function below to work around the optimisation using the empirical BF
#	   opt_temp <- tryCatch(optim.pml(p_temp, optQ = T, optBf = T, optGamma = T), error = function(x) optim.pml(p_temp, optQ = T, optBf = T, optGamma = F)
#	   Ti_temp <- opt_temp$Q[c(2, 5)]
#	   Tv_temp <- opt_temp$Q[c(1, 3, 4, 6)]
########
	Ti_temp = NA
	Tv_temp = NA

	   root_age <- max(allnode.times(tre_temp))
	   base_comp <- base.freq(dat_temp)
#	   print(kappa)
#	   print(omega)
#	   print(alpha)
#	   print(dn)
#	   print(ds)
#	   print(ses)
	   return(list(kappa = kappa, ti = Ti_temp, tv = Tv_temp, omega = omega, alpha = alpha, dn = dn, ds = ds, SES = ses, root_age = root_age, base_comp = base_comp))
}





#source('functions/prune_trees.R')

#prune_tax <- get_taxa(readLines('prune_taxa.txt'))

#data1 <- read.dna('ASFV_N10.fasta', format = 'fasta')
#tree1 <- read.nexus('ASFV_N10.tree')

#prune_test <- prune_tree(tree1, data1, prune_tax$ASFV_N10.fasta)

#codeml_template <- get_codeml_template(prune_test[[1]], prune_test[[2]])

#run_test <- run_codeml('temp_dat.fasta', 'temp_dat.tree', 'temp_dat.ctl', './codeml')