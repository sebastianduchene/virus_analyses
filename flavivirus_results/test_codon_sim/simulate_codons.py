import pyvolve
from Bio import SeqIO
import numpy as np

tr_sim = pyvolve.read_tree(file = 'tr_sim.tree')
freqs = [1/61. for i in range(61)]

#gy_model = pyvolve.Model('GY', {'omega':[0.5, 1, 1.5]  , 'state_freqs':freqs})
#gy_partition = pyvolve.Partition(models = gy_model, size = 30)
#gy_evolver = pyvolve.Evolver(partitions = gy_partition, tree = tr_sim)
#gy_evolver()

print 'gtr'
custom_mu = {"AC":0.5, "AG":0.25, "AT":1.23, "CG":0.55, "CT":1.22, "GT":0.47}
freqs = [0.1, 0.45, 0.3, 0.15]

gtr_model = pyvolve.Model('nucleotide', {'mu': custom_mu, 'state_freqs':freqs}, alpha = 1, num_categories = 4)
#gtr_model =  pyvolve.Model('nucleotide', alpha = 0.2, num_categories = 4, )
gtr_partition = pyvolve.Partition(models = gtr_model, size = 5000)

gtr_evolver = pyvolve.Evolver(partitions = gtr_partition, tree = tr_sim)



gtr_evolver()
print 'gtr'

