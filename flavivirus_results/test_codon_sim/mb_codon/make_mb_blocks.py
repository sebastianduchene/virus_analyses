
# coding: utf-8

# In[17]:

import re, os


# In[18]:

nexus_files = [i for i in os.listdir('.') if len(re.findall('nexus$', i)) > 0]


# In[19]:

print nexus_files
run_block = '''
Begin Mrbayes;
lset nucmodel=codon omegavar=M3;
mcmcp printfreq = 10 ngen = 5000000;
mcmc;
end;
'''


# In[20]:

for f in nexus_files:
    f_temp = open(f, 'r').readlines()
    f_temp.append(run_block)
    open(f, 'w').writelines(f_temp)


# In[8]:




# In[ ]:



