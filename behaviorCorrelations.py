# author: Owen Doyle

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats import zscore
from scipy import io
import matplotlib.axes as ax

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


# In[3]:


#load behvaioral data from matlab
times = io.loadmat(r'/Users/owen/Box/DBS_data/behavior/times.mat')

pats =  times['pats'][0]
contra_bl_corrected = times['contra_bl_corrected']
ipsi_bl_corrected   = times['ipsi_bl_corrected']
contra_improvements = [round(imp[0],2) for imp in times['contra_improvements']] #un-nest the array and round to 2 decimal places
ipsi_improvements   = [round(imp[0],2) for imp in times['ipsi_improvements']]   #un-nest the array and round to 2 decimal places


print pats
print contra_bl_corrected
print ipsi_bl_corrected
print contra_improvements
print ipsi_improvements


# In[4]:


#load fdi r2 values
fdi_rsq_dict_cd = r'/Users/owen/Box/DBS_data/results/r-squared/rsq_10hz_fdi_move'
fdi_rsq_dict    = cpickle_load(fdi_rsq_dict_cd)

#load flex r2 values
flex_rsq_dict_cd = r'/Users/owen/Box/DBS_data/results/r-squared/rsq_10hz_flex_move'
flex_rsq_dict    = cpickle_load(flex_rsq_dict_cd)

#combine fdi and flex, where fdi is used for patients in both dictionaries
rsq_dict = flex_rsq_dict.copy()
rsq_dict.update(fdi_rsq_dict)

# print rsq_dict['BI040'] == fdi_rsq_dict['BI040']
# print rsq_dict['BI040'] == flex_rsq_dict['BI040']


# In[5]:


#load in electrode label masks
filepath = '/Users/owen/Box/DBS_data/sdata_files/elec_masks/elec_label_masks'
elec_label_masks = cpickle_load(filepath, field=None)

print elec_label_masks.keys()
print elec_label_masks['M1'].keys()
print elec_label_masks['M1']['BI020']


# In[6]:


## correlate average unimanual contra and unimanual ipsi M1 R2 with behavioral improvement
elec = 'M1'
uniorbi = 'bi'

c_hand = 0
i_hand = 1

if uniorbi == 'uni':
    c_cond, i_cond = 0,1
elif uniorbi == 'bi':
    c_cond, i_cond = 2,2

contra_R2s = []
ipsi_R2s   = []

for pat in pats:
    pat_ID = "BI0" + str(pat)
    R2s = rsq_dict[pat_ID]
    elec_mask = elec_label_masks[elec][pat_ID]
    
    #contra average R2
    elec_contra_R2 = np.mean(R2s[c_hand,c_cond,elec_mask])
    
    #ipsi average R2
    elec_ipsi_R2 = np.mean(R2s[i_hand,i_cond,elec_mask])
    
    contra_R2s += [elec_contra_R2]
    ipsi_R2s   += [elec_ipsi_R2]
    
 


# In[7]:


contra_R2s


# # ECoG Plots

# In[8]:

#plots
plt.figure(0)
plt.plot(pats, zscore(contra_R2s))
plt.plot(pats, zscore(contra_improvements))
plt.title("{} {} Contra R2 and Contra Behavioral Change in Baseline".format(uniorbi,elec))

plt.figure(1)
plt.plot(pats, zscore(ipsi_R2s))
plt.plot(pats, zscore(ipsi_improvements))
plt.title("{} {} Ipsi R2 and Ipsi Behavioral Change in Baseline".format(uniorbi,elec))

print "contra: ", pearsonr(contra_R2s,contra_improvements)
print "ipsi: ",pearsonr(ipsi_R2s,contra_improvements)


# In[9]:


#scatterplots
plt.figure(0)
plt.scatter(contra_R2s,contra_improvements)
plt.title("{} {} Contra R2 and Contra Behavioral Change in Baseline".format(uniorbi,elec))
plt.xlabel("{} {} Contra R2".format(uniorbi,elec))
plt.ylabel("Change from Baseline in last session (s)")

plt.figure(1)
plt.scatter(ipsi_R2s,ipsi_improvements)
plt.title("{} {} Ipsi R2 and Ipsi Behavioral Change in Baseline".format(uniorbi,elec))
plt.xlabel("{} {} Ipsi R2".format(uniorbi,elec))
plt.ylabel("Change from Baseline in last session (s)")


# # STN Plots

# In[47]:


## correlate average unimanual contra and unimanual ipsi M1 R2 with behavioral improvement
elec = 'ring_ventral' 
uniorbi = 'bi'


#'ring_dorsal'
#'medial_dorsal', 'anterior_dorsal', 'lateral_dorsal'
#'medial_ventral','anterior_ventral','lateral_ventral'
#'ring_ventral'

c_hand = 0
i_hand = 1

if uniorbi == 'uni':
    c_cond, i_cond = 0,1
elif uniorbi == 'bi':
    c_cond, i_cond = 2,2



contra_R2s = []
ipsi_R2s   = []

for pat in pats:
    pat_ID = "BI0" + str(pat)
    R2s = rsq_dict[pat_ID]
    elec_mask = elec_label_masks[elec][pat_ID]
    
    #contra average R2
    elec_contra_R2 = np.mean(R2s[c_hand,c_cond,elec_mask])
    
    #ipsi average R2
    elec_ipsi_R2 = np.mean(R2s[i_hand,i_cond,elec_mask])
    
    contra_R2s += [elec_contra_R2]
    ipsi_R2s   += [elec_ipsi_R2]
    


# In[51]:


len(elec_label_masks['M1']['BI020'])


# In[48]:


#scatterplots
plt.figure(0)
plt.scatter(contra_R2s,contra_improvements)
plt.title("{} {} Contra R2 and Contra Behavioral Change in Baseline".format(uniorbi,elec))
plt.xlabel("{} {} Contra R2".format(uniorbi,elec))
plt.ylabel("Change from Baseline in last session (s)")

plt.figure(1)
plt.scatter(ipsi_R2s,ipsi_improvements)
plt.title("{} {} Ipsi R2 and Ipsi Behavioral Change in Baseline".format(uniorbi,elec))
plt.xlabel("{} {} Ipsi R2".format(uniorbi,elec))
plt.ylabel("Change from Baseline in last session (s)")

print "contra R2: ", pearsonr(contra_R2s,contra_improvements)[0]
print "ipsi R2: ",pearsonr(ipsi_R2s,contra_improvements)[0]


# In[55]:


ventral_mask,       dorsal_mask        = [False]*14, [False]*14
ventral_mask[6:10], dorsal_mask[10:14] = [True]*4,   [True]*4


# In[78]:


## correlate average unimanual contra and unimanual ipsi M1 R2 with behavioral improvement
depth = 'ventral' #ventral or dorsal
uniorbi = 'bi'


#'ring_dorsal'
#'medial_dorsal', 'anterior_dorsal', 'lateral_dorsal'
#'medial_ventral','anterior_ventral','lateral_ventral'
#'ring_ventral'

c_hand = 0
i_hand = 1

if uniorbi == 'uni':
    c_cond, i_cond = 0,1
elif uniorbi == 'bi':
    c_cond, i_cond = 2,2

if depth == 'ventral':
    elec_mask = ventral_mask
elif depth == 'dorsal':
    elec_mask = dorsal_mask

contra_R2s = []
ipsi_R2s   = []

for pat in pats:
    pat_ID = "BI0" + str(pat)
    R2s = rsq_dict[pat_ID]
    
    #contra average R2
    elec_contra_R2 = R2s[c_hand,c_cond,elec_mask]
    
    #ipsi average R2
    elec_ipsi_R2 = R2s[i_hand,i_cond,elec_mask]
    
    contra_R2s += [elec_contra_R2]
    ipsi_R2s   += [elec_ipsi_R2]
    


# In[79]:


contra_improvements_vd = [i*np.ones(4) for i in contra_improvements]
ipsi_improvements_vd   = [i*np.ones(4) for i in ipsi_improvements]


# In[80]:

#scatterplots
plt.figure(0)
plt.scatter(contra_R2s,contra_improvements_vd)
plt.title("{} {} Contra R2 and Contra Behavioral Change in Baseline".format(uniorbi,depth))
plt.xlabel("{} {} Contra R2".format(uniorbi,elec))
plt.ylabel("Change from Baseline in last session (s)")

plt.figure(1)
plt.scatter(ipsi_R2s,ipsi_improvements_vd)
plt.title("{} {} Ipsi R2 and Ipsi Behavioral Change in Baseline".format(uniorbi,depth))
plt.xlabel("{} {} Ipsi R2".format(uniorbi,elec))
plt.ylabel("Change from Baseline in last session (s)")

print "contra R2: ", pearsonr(contra_R2s,contra_improvements)[0]
print "ipsi R2: ",pearsonr(ipsi_R2s,contra_improvements)[0]
