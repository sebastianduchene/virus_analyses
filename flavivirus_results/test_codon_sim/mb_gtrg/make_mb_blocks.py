
# coding: utf-8

# In[5]:

import re, os


# In[9]:

nexus_files = [i for i in os.listdir('.') if len(re.findall('nexus$', i)) > 0]
print os.getcwdu()


# In[10]:

print nexus_files
run_block = '''
begin mrbayes;
lset nst = 6 Rates = Gamma;
mcmcp ngen= 5000000 printfreq = 100;
mcmc;
end;
'''


# In[11]:

for f in nexus_files:
    f_temp = open(f, 'r').readlines()
    f_temp.append(run_block)
    open(f, 'w').writelines(f_temp)


# In[ ]:



