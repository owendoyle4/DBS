#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy.io
import numpy as np
import os

def cpickle_load(flname, field=None):
    '''
    Load a pickle storing a dictioanry of data
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

def cpickle_dump(flname, datadict):
    '''
    Dump a dictionary containing some datas into a pickle file
    '''
    import cPickle
    assert isinstance(datadict, dict)
    print 'Saving %s...' % flname
    fl = open(flname, 'wb')
    cPickle.dump(datadict, fl)
    fl.close()
    return


# In[2]:


ID = ['20','25','29','31','33','35','38','39','40']
#Get chan_labels
ch_names = {}
for pos in range(len(ID)):
    pat = ID[pos]
    matcell =  scipy.io.loadmat(   r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/chan_labels_BI0{}'.format(pat,pat))
    chan_labels = matcell['chan_labels'][0:14] #bc some pats have more than 14 brain chans
    
    labels = []
    for pos_label in range(len(chan_labels)):
        labels.append(str(chan_labels[pos_label][0][0]))
    labels = np.array(labels, dtype = np.string_ ) #np array strings of chan names
    
    ch_names['BI0{}'.format(pat)] = labels


# In[3]:


#create list of unique electrode labels
first = True
label_masks = {}

for key, labels in ch_names.iteritems():
    print key
    if first:
        all_labels = labels
        first = False
    else:
        all_labels = np.append(all_labels, labels)

unique_labels = np.unique(all_labels)
unique_labels


# In[4]:


#Create nested dictionary: outer key = electrode labels, inner key = patient ID, stores mask
for label in unique_labels:
    pat_masks = {}
    
    for key, labels in ch_names.iteritems():
        pat_masks[key] = (labels == label)
            
    label_masks[label] = pat_masks


# In[5]:


print label_masks.keys()
label_masks


# In[6]:


#Save masks
os.chdir(r"/Users/owen/Box/DBS_data/sdata_files/elec_masks")
cpickle_dump('elec_label_masks', label_masks)


# In[10]:


#create a new R-squared dictionary organized by electrode type
rsq_dict = {}

#make sure fdi comes second so the flex values for pat39 and 40 are overridden
for name in ["rsq_10hz_flex_move", "rsq_10hz_fdi_move"]:
    temp_dict = r'/Users/owen/Box/DBS_data/results/r-squared/{}'.format(name)
    temp_dict = cpickle_load(temp_dict)
    print temp_dict.keys()
        
    rsq_dict.update(temp_dict)
    


# In[11]:


#check that the patients who have both EMGs use fdi values and not flex
temp_dict = r'/Users/owen/Box/DBS_data/results/r-squared/rsq_10hz_fdi_move'
temp_dict = cpickle_load(temp_dict)

print "***True for both is needed***"
print (temp_dict['BI040'] == rsq_dict['BI040']).all()
print (temp_dict['BI039'] == rsq_dict['BI039']).all()


# In[13]:


rsq_elec = dict()
pat_IDS  = dict()

for elec in label_masks:
    first_pat = True
    pats = []
    for pat in label_masks[elec]:
        pat_mask = label_masks[elec][pat]
        pat_rsq = rsq_dict[pat][:,:,pat_mask]
        pats = pats + [pat]*pat_rsq.shape[2]
        
        if first_pat:
            rsq_elec[elec] = pat_rsq
            first_pat = False
        else:
            rsq_elec[elec] = np.concatenate([rsq_elec[elec], pat_rsq], axis = 2)
    pat_IDS[elec] = pats
    print elec, pats

rsq_elec['pats'] = pat_IDS


# In[14]:


#confirm that positions are preserved
elec = 'M1'
pat  = 'BI029'

print rsq_elec['pats'][elec]
elec_idx = [pat_id == pat for pat_id in rsq_elec['pats'][elec]]
print elec_idx
print '------------------'

pat_idx = label_masks[elec][pat]
print pat_idx
print '------------------'

print rsq_elec[elec][:,:,elec_idx] == rsq_dict[pat][:,:,pat_idx]


# In[ ]:


#Save new rsq dictionary
os.chdir(r"/Users/owen/Box/DBS_data/results/r-squared")
cpickle_dump('rsq_elec', rsq_elec)


# In[ ]:




