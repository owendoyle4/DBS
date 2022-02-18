#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.stats import zscore
import os
from scipy import io
import matplotlib.axes as ax

np.set_printoptions(precision=4, suppress=True)

# defines our correlation "r" value (square this to get R^2)
def columnwise_correlation(ypred, y, zscorea=True, zscoreb=True, axis=0):
    r'''Compute correlations efficiently
    Examples
    --------
    >>> x = np.random.randn(100,2)
    >>> y = np.random.randn(100,2)
    >>> cc = columnwise_correlation(x, y)
    >>> cc.shape
    (2,)
    >>> c1 = np.corrcoef(x[:,0], y[:,0])[0,1]
    >>> c2 = np.corrcoef(x[:,1], y[:,1])[0,1]
    >>> assert np.allclose(cc, np.r_[c1, c2])
    Notes
    -----
    Recall that the correlation cofficient is defined as
    .. math::
       \rho_{x, y} = \frac{cov(X,Y)}{var(x)var(y)}
    Since it is scale invariant, we can zscore and get the same
    .. math::
       \rho_{x, y} = \rho_{zscore(x), zscore(y)} = \frac{cov(X,Y)}{1*1} =
       \frac{1}{N}\frac{\sum_i^n \left(x_i - 0 \right) \left(y_i - 0 \right)}{1*1} =
       \frac{1}{N}\sum_i^n \left(x_i * y_i \right)
    '''
    if zscorea:
        y = zscore(y, axis=axis)
    if zscoreb:
        ypred = zscore(ypred, axis=axis)
    corr = (y * ypred).mean(axis=axis)
    return corr

# CPICKLE loads dict
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

# CPICKLE SAVES dict
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

# # os.chdir('C:\Users\natal\Box\DBS_data\sdata_files\BI020')
# ls #what is in this folder

# # alg2 - 10Hz LPF


# In[6]:


# load the data (emg and neural data)
# preprocessed w/ organized_dbs_data_11_03_20.m 
#######################################################################################
pat = 41 #CHOOSE PATIENT NUMBER
# type_emg = 'flex' #fdi or flex
filename = '10hz' #add
########################################################################################

# SDATA.data
filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/data_{}.mat'.format(pat,filename)
print filepath #add
test = {}

with h5py.File(filepath,'r') as f:
    for k, v in f.items():
        test[k] = np.array(v, dtype = np.float32)
        
dbs_data = test['data']
dbs_data = dbs_data.T

# SDATA event times 
filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/events.mat'.format(pat)
test = {}

with h5py.File(filepath,'r') as f:
    for k, v in f.items():
        test[k] = np.array(v, dtype = np.int)
        
events = test['events']
events = events.T

# load the neural ft 
# load LMP (which was extracted/processed with extracting_DBS_power.m)
filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/dbs_lmp_1k.mat'.format(pat)
test = {}

with h5py.File(filepath,'r') as f:
    for k, v in f.items():
        test[k] = np.array(v, dtype = np.float32)
        
dbs_lmp = test['lmp'] #1khz is dbs_lmp_10
dbs_lmp = dbs_lmp.T

### set neural feat ### 
lmp = dbs_lmp 


# In[7]:


# using move times
# filepath = r'C:\Users\natal\Box\DBS_data\sdata_files\BI0{}\move_BI0{}_temp.mat'.format(pat,pat)
filepath =    r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/move_BI0{}_temp.mat'.format(pat,pat)
print filepath
test = {}

with h5py.File(filepath,'r') as f:
    for k, v in f.items():
        test[k] = np.array(v, dtype = np.int)
        
move_times = test['move_BI0{}_temp'.format(pat)]
move_times = move_times.T

move_indicies = move_times

print move_indicies

move_indicies.shape[0]


# In[8]:


#make bimanual mask

#column indx in move_indicies:
num_blocks = move_indicies.shape[0]
cond_label = 0
start_block = 1 
end_block = 2

data_length = dbs_data.shape[0] #length of all conditons

# initialize mask for both blocks of bimanual condition
biman_mask = np.zeros(data_length, bool) #False

# creates masks based on the label (column 1 of move_indicies)
num_bi = 0
for block in range(num_blocks): #go thru # rows of move_indicies
    start = move_indicies[block,start_block]
    stop = move_indicies[block,end_block]
    
    if move_indicies[block,cond_label] == 2 and [start,stop] != [0,100]: #completely unusable blocks are arbitrarily assigned start times of 0ms and stop times of 100ms, exlude these blocks     
        biman_mask[start:stop] = True
        print(start, stop)
        num_bi += 1

        
# HOW LONG?
print sum(biman_mask)


# In[9]:


#make sure mask works
contra_emg = dbs_data[:,14]
plt.plot(contra_emg)
plt.show()
plt.plot(contra_emg[biman_mask])
plt.show()


# In[10]:


## adding time lags to the emg data
contra_emg = dbs_data[:,14]
ipsi_emg   = dbs_data[:,15]

num_emg_lags = 1000
num_samples = contra_emg.shape[0]

contra_lags  = np.zeros((num_emg_lags,num_samples))
contra_emg_padded = np.hstack([np.zeros([500]),np.squeeze(contra_emg[:-500])])

ipsi_lags  = np.zeros((num_emg_lags,num_samples))
ipsi_emg_padded = np.hstack([np.zeros([500]),np.squeeze(ipsi_emg[:-500])])

for lag in range(num_emg_lags):
    contra_lags[lag,:] = np.hstack([contra_emg_padded[1000-lag:],np.zeros(1000-lag)])
    ipsi_lags[lag,:] = np.hstack([ipsi_emg_padded[1000-lag:],np.zeros(1000-lag)])
    
## to calculate cross correlation
corr_ipsi   = np.zeros((1000, 2)) #1000 time lags, for each hand (contra and ipsi)
corr_contra = np.zeros((1000, 2))


hand = 0
for contant_hand in [contra_emg, ipsi_emg]:
    for lag in range(num_emg_lags):
        corr_contra[lag, hand] = columnwise_correlation(contra_lags[lag,biman_mask], contant_hand[biman_mask])
        corr_ipsi[lag, hand]   = columnwise_correlation(ipsi_lags[lag,biman_mask],   contant_hand[biman_mask])
    hand += 1
            

print corr_ipsi.shape
corr_ipsi = np.squeeze(corr_ipsi)
corr_contra = np.squeeze(corr_contra)


# In[ ]:


if pat == 38 or pat == 39:
    num_bi = 3
else:
    num_bi = 2

## adding time lags to the emg data : 100 lags before
contra_emg = dbs_data[:,14]
ipsi_emg   = dbs_data[:,15]

num_emg_lags = 1000
num_samples = contra_emg.shape[0]

contra_lags  = np.zeros((num_emg_lags,num_samples))
contra_emg_padded = np.hstack([np.zeros([500]),np.squeeze(contra_emg[:-500])])

ipsi_lags  = np.zeros((num_emg_lags,num_samples))
ipsi_emg_padded = np.hstack([np.zeros([500]),np.squeeze(ipsi_emg[:-500])])

for lag in range(num_emg_lags):
    contra_lags[lag,:] = np.hstack([contra_emg_padded[1000-lag:],np.zeros(1000-lag)])
    ipsi_lags[lag,:] = np.hstack([ipsi_emg_padded[1000-lag:],np.zeros(1000-lag)])
    
## to calculate cross correlation
corr_ipsi   = np.zeros((1000, num_bi, 2)) #1000 time lags, for each block, for each hand (contra and ipsi)
corr_contra = np.zeros((1000, num_bi, 2))


block_num = 0
for block in range(len(move_indicies)):
    if move_indicies[block,0] == 2: #if block is bimanual
        block_start = move_indicies[block,1]
        block_end   = move_indicies[block,2]
#         elec_temp = contra_lags[499,block_start:block_end]
#         elec_temp = ipsi_lags[499,block_start:block_end]

        hand = 0
        for elec_temp in [contra_lags, ipsi_lags]:
            elec_temp = elec_temp[499,block_start:block_end]
        
            for lag in range(num_emg_lags):
                corr_ipsi[lag,block_num, hand]   = columnwise_correlation(ipsi_lags[lag,block_start:block_end],elec_temp)
                corr_contra[lag,block_num, hand] = columnwise_correlation(contra_lags[lag,block_start:block_end],elec_temp)
            hand += 1
            
        block_num += 1

print corr_ipsi.shape
corr_ipsi = np.squeeze(corr_ipsi)
corr_contra = np.squeeze(corr_contra)


# In[11]:


plt.plot(corr_contra[:,0])
plt.title("Pat {}: Contra corr with contra zero lag (check)".format(pat))
plt.show()

plt.plot(corr_ipsi[:,0])
plt.title("Pat {}: Ipsi corr with contra zero lag".format(pat))
plt.show()

####

plt.plot(corr_contra[:,1])
plt.title("Pat {}: Contra corr with ipsi zero lag".format(pat))
plt.show()

plt.plot(corr_ipsi[:,1])
plt.title("Pat {}: Ipsi corr with ipsi zero lag (check)".format(pat))
plt.show()


# In[ ]:


## adding time lags to the emg data : 100 lags before
contra_emg = dbs_data[:,14]
ipsi_emg   = dbs_data[:,15]

num_emg_lags = 1000
num_samples = contra_emg.shape[0]

contra_lags  = np.zeros((num_emg_lags,num_samples))
contra_emg_padded = np.hstack([np.zeros([500]),np.squeeze(contra_emg[:-500])])

ipsi_lags  = np.zeros((num_emg_lags,num_samples))
ipsi_emg_padded = np.hstack([np.zeros([500]),np.squeeze(ipsi_emg[:-500])])

for lag in range(num_emg_lags):
    contra_lags[lag,:] = np.hstack([contra_emg_padded[1000-lag:],np.zeros(1000-lag)])
    ipsi_lags[lag,:] = np.hstack([ipsi_emg_padded[1000-lag:],np.zeros(1000-lag)])

    
## to calculate cross correlation
corr_ipsi = np.zeros((1000))
corr_contra = np.zeros((1000))

# elec_temp = contra_lags[499,:]
elec_temp = ipsi_lags[499,:]

block_start = 
block_start = 

for lag in range(num_emg_lags):
    corr_ipsi[lag] = columnwise_correlation(ipsi_lags[lag,1000:17000],elec_temp[1000:17000])
    corr_contra[lag] = columnwise_correlation(contra_lags[lag,1000:17000],elec_temp[1000:17000])

corr_ipsi = np.squeeze(corr_ipsi)
corr_contra = np.squeeze(corr_contra)


# In[ ]:




