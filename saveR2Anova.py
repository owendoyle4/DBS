#!/usr/bin/env python
# coding: utf-8

# In[42]:


from scipy import io
import numpy as np
import os

def cpickle_load(flname, field=None):
    '''
    Load a pickle storing a dictionary of data
    '''
    import cPickle
    with open(flname, 'r') as fl:
        data = cPickle.load(fl)
        if isinstance(data, list):
            pass
        elif len(data.keys()) == 1:
            data = data[data.keys()[0]]
        elif field is not None:
            data = data[field]
    return data


# In[43]:


#SAVE UNIMANUAL RSQ VALUES FOR MATLAB ANOVA TEST
#MUST RUN BLOCK ABOVE SO rsq_dict IS THE PROPER DICTIONARY

#initialize lists of r-squared. values
fdi_contra     = []
fdi_ipsi       = []
fdi_bi_contra  = []
fdi_bi_ipsi    = []

flex_contra    = []
flex_ipsi      = []
flex_bi_contra = []
flex_bi_ipsi   = []

#hand/condition indicies
c_index = 0
i_index = 1
b_index = 2

for emg_type in ['fdi', 'flex']:
    rsq_dict = r'/Users/owen/Box/DBS_data/results/r-squared/rsq_10hz_{}_move'.format(emg_type)
    rsq_dict = cpickle_load(rsq_dict)
    print rsq_dict.keys()
    for ID in rsq_dict.keys():
        pat = rsq_dict[ID]
        
        if emg_type == 'fdi':
                                                                 #hand    cond   elec(all)
            fdi_contra     = np.concatenate((fdi_contra,      pat[c_index,c_index,:]), axis=None)
            fdi_ipsi       = np.concatenate((fdi_ipsi,        pat[i_index,i_index,:]), axis=None)
            fdi_bi_contra  = np.concatenate((fdi_bi_contra,   pat[c_index,b_index,:]), axis=None)
            fdi_bi_ipsi    = np.concatenate((fdi_bi_ipsi,     pat[i_index,b_index,:]), axis=None)
        
        elif emg_type == 'flex':
                                                                #hand    cond   elec(all)
            flex_contra     = np.concatenate((flex_contra,    pat[c_index,c_index,:]), axis=None)
            flex_ipsi       = np.concatenate((flex_ipsi,      pat[i_index,i_index,:]), axis=None)
            flex_bi_contra  = np.concatenate((flex_bi_contra, pat[c_index,b_index,:]), axis=None)
            flex_bi_ipsi    = np.concatenate((flex_bi_ipsi,   pat[i_index,b_index,:]), axis=None)


# In[44]:


print fdi_contra.size,  fdi_ipsi.size,  fdi_bi_contra.size,  fdi_bi_ipsi.size
print flex_contra.size, flex_ipsi.size, flex_bi_contra.size, flex_bi_ipsi.size


# In[45]:


print 'FDI Unimanual: ',  fdi_contra.mean(),     fdi_ipsi.mean()
print 'Flex Unimanual: ', flex_contra.mean(),    flex_ipsi.mean()
print '--------------------------------------------------------'
print 'FDI Bimanual: ',   fdi_bi_contra.mean(),  fdi_bi_ipsi.mean()
print 'Flex Bimanual: ',  flex_bi_contra.mean(), flex_bi_ipsi.mean()


# In[47]:


#SAVE
os.chdir('/Users/owen/Box/DBS_data/results/r-squared')

io.savemat('anova_uni_rsq', dict(fdi_contra     = fdi_contra,
                                 fdi_ipsi       = fdi_ipsi,
                                 flex_contra    = flex_contra,
                                 flex_ipsi      = flex_ipsi))

io.savemat('anova_bi_rsq',  dict(fdi_bi_contra  = fdi_bi_contra,
                                 fdi_bi_ipsi    = fdi_bi_ipsi,
                                 flex_bi_contra = flex_bi_contra,
                                 flex_bi_ipsi   = flex_bi_ipsi))


# In[ ]:




