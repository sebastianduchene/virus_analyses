source('pps_dpp.R')

trees_log <- read.nexus('test_aa.nxs.run2.t')
params_log <- read.table('test_aa.nxs.run2.p', head = T, skip = 1)

s1 <- simulate_pps_aa(log_out = params_log, tree_out = trees_log, s_len = 1000, n_samples = 20)


