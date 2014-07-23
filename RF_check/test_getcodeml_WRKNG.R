get_codeml_template <- function(seq_data, tree_data, temp_name = 'temp_dat'){
  codeml_template <- "seqfile = @INPUTSEQFILE@.fasta\n treefile = @INPUTTREEFILE@.tree\n outfile = @OUTPUTFILE@.out \n noisy = 9
 \n verbose = 1
 \n runmode = 0
 \n seqtype = 1
 \n   CodonFreq = 0
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
 \n   ncatG = 5
 \n   getSE = 1
\n RateAncestor = 0
\n Small_Diff = 5e-7
\n cleandata = 0
\n fix_blength = 2
 \n  method = 0
"
  codeml_template <- gsub('@[A-Z]+@', temp_name, codeml_template)
  return(codeml_template)
}

# MAKE FUNCTION TO RUN CODEML
# MAKE FUNCTION TO COLLECT RESUTLS FROM CODEML




source('functions/prune_trees.R')

prune_tax <- get_taxa(readLines('prune_taxa.txt'))

data1 <- read.dna('ASFV_N10.fasta', format = 'fasta')
tree1 <- read.nexus('ASFV_N10.tree')


prune_test <- prune_tree(tree1, data1, prune_tax$ASFV_N10.fasta)

wow <- get_codeml_template(prune_test[[1]], prune_test[[2]])
