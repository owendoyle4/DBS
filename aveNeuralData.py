# author: Owen Doyle

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats import zscore
import os
from scipy import io
import matplotlib.axes as ax

np.set_printoptions(precision=4, suppress=True)

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


# In[2]:


# load in electrode label masks
filepath = '/Users/owen/Box/DBS_data/sdata_files/elec_masks/elec_label_masks'
elec_label_masks = cpickle_load(filepath, field=None)

# check outputs
print elec_label_masks.keys()
print elec_label_masks['M1'].keys()
print elec_label_masks['M1']['BI020']


# In[3]:


#use to check that peak data is properly imported
###############################
pat = 29
filename = 'peak_times_200'
###############################

filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/{}.mat'.format(pat,filename)
print filepath 
data = {}

with h5py.File(filepath,'r') as f:
    for k, v in f.items():
        data[k] = np.array(v, dtype = np.float32)

contra_peaks   = data['contra_peaks']    
ipsi_peaks     = data['ipsi_peaks']
bimanual_peaks = data['bimanual_peaks']
window         = int(data['epoch_width'][0][0])

print contra_peaks.shape
print ipsi_peaks.shape
print bimanual_peaks.shape
print window

print "for each condition: shape = (hands x #peaks). Peaks for both hands are saved for each condition, even unimanual."


# In[9]:


# Load peak times, neural data, and emg data
def load_peaks_data_chans(pats_list, window):
    ##############################################################################################################
    # pats is a list of patient ID numbers                                                                       #
    # ex: pats = [20,25,29,31,33,35,38,39,40,41]                                                                 #
    #                                                                                                            #
    # -returns peaks, a dictionary which contains the times of peaks of both emgs for all conditions of all      #
    # patients along with each patient's neural and emg timeseries, window for epochs, and channel labels        #
    ##############################################################################################################
    
    peaks = {}
    for pat in pats_list:
        #get peak times
        filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/peak_times_{}.mat'.format(pat,str(window))
        matlab_data = {}
        with h5py.File(filepath,'r') as f:
            for k, v in f.items():
                matlab_data[k] = np.array(v, dtype = np.float32)

        #get emg data
        filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/data_10hz.mat'.format(pat)
        with h5py.File(filepath,'r') as f:
            for k, v in f.items():
                matlab_data[k] = np.array(v, dtype = np.float32)
        
        #get neural data
        filepath =    r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/dbs_lmp_1k.mat'.format(pat)
        with h5py.File(filepath,'r') as f:
            for k, v in f.items():
                matlab_data[k] = np.array(v, dtype = np.float32)

        #get channel labels
        filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/chan_labels_BI0{}.mat'.format(pat,pat)
        chan_labels = io.loadmat(filepath)
        
        #add the lpf emg to the lpf neural data set --> 0->5 = ecog, 6->13 = stn, 14->15 = emg
        matlab_data['lmp'] = np.vstack((matlab_data['lmp'], matlab_data['data'][14:16,:]))
        
        #save info to peaks dictionary
        peaks[pat] = {}
        peaks[pat]['data']     = matlab_data['lmp']
        peaks[pat]['contra']   = matlab_data['contra_peaks']    
        peaks[pat]['ipsi']     = matlab_data['ipsi_peaks']
        peaks[pat]['bimanual'] = matlab_data['bimanual_peaks']
        peaks[pat]['window']   = int(data['epoch_width'][0][0])
        peaks[pat]['chans']    = [chan_labels['chan_labels'][num][0][0] for num in range(16)] #6 ecog + 8 stn + 2 emg
    
    return peaks


# In[5]:


#make trials for unimanual data
def make_trials(peaks, uniorbi):
    ################################################################################
    # peaks is a dictionary returned by load_peaks_data_chans()                    #
    # uniorbi is 'uni' or 'bi' to choose between using unimanual or bimanual data  # 
    #                                                                              #                                              #
    # -returns trials, a dictionary which contains the timeseries of the neural    #
    # and emg epochs (trials) for both hands for each patient.                     #
    ################################################################################
    
    #set variables
    contra_idx = 0
    ipsi_idx   = 1
    num_elecs  = 16 #ecog, stn, and emg
    window   = peaks[peaks.keys()[0]]['window'] #width of trial, ms
    trials = {}
    
    #add each patient's data to trials
    for pat in peaks.keys():
        trials[pat] = {}
        
        pat_data = peaks[pat]['data'] #neural and emg timeseries
        
        for hand in ['contra', 'ipsi']:
            hand_idx = eval(hand + '_idx')
            
            #condition = hand for unimanual conditions
            if uniorbi == 'uni':
                cond = hand
            else:
                cond = 'bimanual'

            #get patient's peaks from 'peaks' dictionary for his hand
            pat_peaks = peaks[pat][cond][hand_idx]
            pat_peaks = [peak for peak in pat_peaks if peak] #only use the peaks that are not zero (excludes any buffer zeros at the end of either the contra or ipsi column)
            num_peaks = len(pat_peaks)

            #data shape: #elec x #peaks x time
            trials[pat][hand] = np.zeros((num_elecs, num_peaks, window))
            
            #add each peak's trial to the structure for all electrodes
            for peak_idx in range(num_peaks):
                peak = pat_peaks[peak_idx]
                start = int(peak-(window/2)) 
                stop  = int(peak+(window/2)) 

                trials[pat][hand][:,peak_idx,:] = pat_data[0:num_elecs,start:stop]

                    
    return trials


# In[6]:


def make_trial_avg(trials):
    #############################################################################################################################
    # trials is a dictionary returned by make_trials()                                                                          #
    # note: uniorbi is not an argument, make_trial_avg() will match the uniorbi argument passed to make_trials()                #
    #                                                                                                                           #
    # -returns pat_avgs, a dictionary which contains each electrode's average epoch (trial) for both hands for each patient     #
    # -returns elec_avgs, a dictionary which contains each electrode's average epoch (trial) for both hands across all patients #
    # -returns num_elecs, a dictionary which contains the number of electrodes used in each average in elec_avgs                #
    #############################################################################################################################
    
    pats = trials.keys()
    
    #same for all patients, so just use the first patient
    pat        = pats[0]
    num_elecs  = trials[pat]['contra'].shape[0]  #arbitrarily use contra
    window     = trials[pat]['contra'].shape[2]  #arbitrarily use contra
    
    #initialize dictionaries
    pat_avgs   = {}
    elec_avgs  = {}
    num_elecs  = {} #saves the overall number of each type of electrodes across all patients
    elec_count = {} #only used for indexing into elec_avgs
    
    #initialize electrode averages and retrieve the # of each type of electrode    
    for elec in elec_label_masks.keys():        
        elec_avgs[elec] = {}
        elec_count[elec] = 0
        num_elec = 0
        
        for pat in pats:
            pat_ID = "BI0" + str(pat)
            num_elec += sum(elec_label_masks[elec][pat_ID])
            
        elec_avgs[elec]['contra'] = np.zeros((num_elec, window))
        elec_avgs[elec]['ipsi']   = np.zeros((num_elec, window))
        num_elecs[elec] = num_elec
        
    for pat in pats:
        #average across trials within patient, within electrode channels (ex: if there are two M1 elecs for one pat, they stay separate)
        pat_avgs[pat] = {}
        pat_avgs[pat]['contra'] = np.mean(trials[pat]['contra'],1)
        pat_avgs[pat]['ipsi']   = np.mean(trials[pat]['ipsi'],1)
    
        #group by electrode type across patients (ex: if there are two M1 elecs for pat1 and one M1 elec for pat 2, they are all combined into one group)
        for elec_idx in range(14): #doesn't include 2 emg electrodes
            elec = peaks[pat]['chans'][elec_idx]

            elec_avgs[elec]['contra'][elec_count[elec],:] = pat_avgs[pat]['contra'][elec_idx,:]
            elec_avgs[elec]['ipsi'][elec_count[elec],:]   = pat_avgs[pat]['ipsi'][elec_idx,:]
            
            elec_count[elec] += 1
    
    #average across the group to have two averages (contra and ipsi) per electrode type
    for elec in elec_label_masks.keys():
        elec_avgs[elec]['contra'] = np.mean(elec_avgs[elec]['contra'],0)
        elec_avgs[elec]['ipsi']   = np.mean(elec_avgs[elec]['ipsi'],0)
            
    return pat_avgs, elec_avgs, num_elecs


# In[10]:


uni_pats = [25, 35, 38] #CAN INCLUDE 41 ONCE ELEC LABELS ARE IN
bi_pats  = [20, 25, 33, 40] #CAN INCLUDE 41 ONCE ELEC LABELS ARE IN
all_pats = [20, 25, 29, 31, 33, 35, 38, 40] #not using 39 or 41 yet
#enter parameters
pats    = all_pats
window  = 200

#call functions to get peak for each electrode trials
peaks      = load_peaks_data_chans(pats, window)
uni_trials = make_trials(peaks, 'uni')
bi_trials  = make_trials(peaks, 'bi')

#make averages across trials for visualization
uni_pat_avgs, uni_elec_avgs, num_elecs = make_trial_avg(uni_trials)
bi_pat_avgs, bi_elec_avgs, num_elecs = make_trial_avg(bi_trials)


# In[11]:


#Show each structure's content
print '----------------------------------------------------'
print "'peaks' information"
print "patient list: ", peaks.keys()
print "for each patient: ", peaks[pats[0]].keys()
print "contra: ", peaks[pats[0]]['contra'].shape
print '----------------------------------------------------'
print
print
print '----------------------------------------------------'
print "'trials' information (same for uni and bi)"
print "patient list: ",uni_trials.keys()
print "for each patient: ", uni_trials[pats[0]].keys()
print "contra: ", uni_trials[pats[0]]['contra'].shape
print '----------------------------------------------------'
print 
print
print '----------------------------------------------------'
print "'pat_avgs' information (same for uni and bi)"
print "patient list: ",uni_pat_avgs.keys()
print "for each patient: ", uni_pat_avgs[pats[0]].keys()
print "contra: ", uni_pat_avgs[pats[0]]['contra'].shape
print '----------------------------------------------------'
print
print
print '----------------------------------------------------'
print "'elec_avgs' information (same for uni and bi)"
print uni_elec_avgs.keys()
print "for each electrode: ", uni_elec_avgs['M1'].keys()
print "contra: ", uni_elec_avgs['lateral_ventral']['contra'].shape
print '----------------------------------------------------'


# # Patient average plots

# In[12]:


#CHECK AVERAGE OF EMGS - Peak should be time locked to 150 (center of window)
pat = 25
uniorbi = 'uni'

#choose unimanual or bimanual based on selection
if uniorbi == 'uni':
    pat_avgs = uni_pat_avgs
elif  uniorbi == 'bi':
    pat_avgs = bi_pat_avgs
    
plt.figure()
plt.plot(pat_avgs[pat]['contra'][14,:])
plt.title('Patient {}: Average {}manual Contra EMG Trial'.format(pat,uniorbi))

plt.figure()
plt.plot(pat_avgs[pat]['ipsi'][15,:])
plt.title('Patient {}: Average {}manual Ipsi EMG Trial'.format(pat,uniorbi))
print


# In[13]:


#ECOG PLOTS
elec_idx = 5

for elec in range(6):
    if elec_idx == elec:
        lw = 3
    else:
        lw = 0.5
    
    #During unimanual contra movement
    plt.figure(0)
    plt.plot(pat_avgs[pat]['contra'][elec,:], linewidth = lw)
    plt.title('Patient {}: ECoG during {}manual Contra'.format(pat,uniorbi))
    
    #During unimanual ipsi movement
    plt.figure(1)
    plt.plot(pat_avgs[pat]['ipsi'][elec,:], linewidth = lw)
    plt.title('Patient {}: ECoG during {}manual Ipsi'.format(pat,uniorbi))
    
plt.figure(0)
plt.legend(peaks[pat]['chans'][0:6], bbox_to_anchor=(1.05, 1))
print peaks[pat]['chans'][elec_idx]


# In[14]:


#Plot one ecog elec
elec = 5
lw = 3 #lineweight

contra_sig = pat_avgs[pat]['contra'][elec,:]
ipsi_sig = pat_avgs[pat]['ipsi'][elec,:]

#During unimanual contra movement
plt.figure(0)
plt.plot(contra_sig, linewidth = lw)
plt.title('Patient {}: ECoG during {}manual Contra'.format(pat,uniorbi))

#During unimanual ipsi movement
plt.figure(1)
plt.plot(ipsi_sig, linewidth = lw)
plt.title('Patient {}: ECoG during {}manual Ipsi'.format(pat,uniorbi))

plt.figure(0)
# plt.legend(peaks[pat]['chans'][elec], bbox_to_anchor=(1.05, 1))
print peaks[pat]['chans'][elec]

print pearsonr(contra_sig, ipsi_sig)


# In[15]:


#STN PLOTS
for elec in range(6,14):
    if elec == 6:
        c = 'k'
    elif elec in range(7,10):
        c = 'b'
    elif elec in range(10,13):
        c = 'r'
    else:
        c = 'teal'
    
    #During unimanual contra movement
    plt.figure(0)
    plt.plot(pat_avgs[pat]['contra'][elec,:], linewidth = 0.5, color = c)
    plt.title('Patient {}: STN during {}manual Contra'.format(pat,uniorbi))
    
    #During unimanual ipsi movement
    plt.figure(1)
    plt.plot(pat_avgs[pat]['ipsi'][elec,:], linewidth = 0.5, color = c)
    plt.title('Patient {}: STN during {}manual Ipsi'.format(pat,uniorbi))

plt.figure(0)
plt.legend(peaks[pat]['chans'][6:14], bbox_to_anchor=(1.05, 1))


# # Electrode average plots

# In[17]:


uniorbi = 'bi'

#choose unimanual or bimanual based on selection
if uniorbi == 'uni':
    elec_avgs = uni_elec_avgs
elif  uniorbi == 'bi':
    elec_avgs = bi_elec_avgs
    
#separate by ecog and STN
# elecs = trial_avgs.keys()
elecs = elec_avgs.keys()
ecog_elecs = {}
stn_elecs  = {}

for elec in elecs:
    if 'lfp' in elec:
        stn_elecs[elec] = elec_avgs[elec]
    else:
        ecog_elecs[elec] = elec_avgs[elec]


# In[26]:


#ECoG plots
for elec in ecog_elecs.keys():
    #During contra movement
    plt.figure(0)
    plt.plot(elec_avgs[elec]['contra'], linewidth = 0.5)
    plt.title('ECoG during {}manual Contra'.format(uniorbi))
   
    
    #During ipsi movement
    plt.figure(1)
    plt.plot(elec_avgs[elec]['ipsi'], linewidth = 0.5)
    plt.title('ECoG during {}manual Ipsi'.format(uniorbi))

plt.figure(0)
lgd_txt = ["{} (n={})".format(elec,num_elecs[elec]) for elec in ecog_elecs.keys()]
plt.legend(lgd_txt, bbox_to_anchor=(1.05, 1))


# # Plotting individual electrode averages and held out vs predicted test sets

# across task generalization

# In[27]:


#Compare unimanual c/i to bimanual c/i (descriptive for across task generalization)
pat = 25
elec = 3
cond = 'ipsi'

lw = 3

uni_sig = uni_pat_avgs[pat][cond][elec,:]
bi_sig = bi_pat_avgs[pat][cond][elec,:]

plt.figure
plt.plot(uni_sig, linewidth = lw)
plt.plot(bi_sig, linewidth = lw)
plt.title("patient {}, {}({}) {}".format(pat,peaks[pat]['chans'][elec],elec, cond))
plt.legend(['uni','bi'], loc = 'upper left')
plt.margins(x=0)

print "R: ", pearsonr(uni_sig, bi_sig)[0]
print "R^2: ", pearsonr(uni_sig, bi_sig)[0]**2


# In[28]:


### PREDICTED VS HELD OUT FOR ACROSS TASK GENERALIZATION ###
if pat in [31,33,38]: #flex_pats
    emg_type = 'flex'
else:
    emg_type = 'fdi'
    
if cond == 'contra':
    cond_idx = 0
else:
    cond_idx = 1


#load in train uni test uni data
change_id = r'/Users/owen/Box/DBS_data/results/BI0{}/BI0{}_10hz_{}_move'.format(pat,pat,emg_type)
reg_dict = cpickle_load(change_id)

#load in generalization data
change_id = r'/Users/owen/Box/DBS_data/results/uni_bi_gen/BI0{}_uni_bi_{}_{}'.format(pat,cond,emg_type)
gen_dict = cpickle_load(change_id)


# In[29]:


### PLOT ###
testset_num = 2
#try different test sets (0->4)
gen_quality = 'poorly' #generalization quality: 'well' or 'poorly'
reg_color = 'darkblue'
gen_color = 'darkgreen'
hand_idx = cond_idx #because this is only unimanual

#choose the test sets
if gen_quality == 'well':
    #use the best looking test set (most highly correlated) for both the regular and generalized timeseries
#     reg_test_set =testset_num
    reg_test_set = np.argmax(reg_dict['cc_test'][hand_idx, cond_idx,:, elec])
    gen_test_set = np.argmax(gen_dict['cc_test_gen'][:,elec])

elif gen_quality == 'poorly':
    #use the best looking test set (most highly correlated) for the regular timeseries and 
#     reg_test_set =testset_num
    reg_test_set = np.argmax(reg_dict['cc_test'][hand_idx, cond_idx,:, elec])
    #use the worst test set for generalized timeseries
    gen_test_set = np.argmin(gen_dict['cc_test_gen'][:,elec])

#save test sets and get the time to plot
predicted_reg = zscore(reg_dict['prediction_test'][hand_idx, cond_idx, reg_test_set, elec])
heldout_reg   = zscore(reg_dict['neural_test'][hand_idx, cond_idx, reg_test_set, elec])

predicted_gen = zscore(gen_dict['prediction_test_gen'][gen_test_set,elec])
heldout_gen   = zscore(gen_dict['neural_test_gen'][gen_test_set,elec])

if min(len(predicted_reg), len(predicted_gen)) > 4000:
    time = 4000
else:
    time = min(len(predicted_reg), len(predicted_gen))
    
    
fig, axs = plt.subplots(2, sharey=True)    
#plot regular pred vs held out
axs[0].plot(predicted_reg[0:time],   color = reg_color,   linewidth = 5, alpha = .5)
axs[0].plot(heldout_reg[0:time], color = 'black', linewidth = 4, alpha = .8)
axs[0].set_ylim([-4,4])
axs[0].set_title('BI0{}({}), {} unimanual, testset={}'.format(pat,elec,cond,reg_test_set))
axs[0].margins(x=0)

print "Unimanual R:",pearsonr(predicted_reg[:],heldout_reg)[0]
print "Unimanual R^2:", pearsonr(predicted_reg[:],heldout_reg)[0]**2

#plot generalization pred vs held out
axs[1].plot(predicted_gen[0:time], color = gen_color,   linewidth = 5, alpha = .5)
axs[1].plot(heldout_gen[0:time],   color = 'black', linewidth = 4, alpha = .8)
axs[1].set_ylim([-4,4])
axs[1].set_title('BI0{}({}), {} uni to bi generalization, testset={}'.format(pat,elec,cond,gen_test_set))
axs[1].margins(x=0)

print "Generalized R: ", pearsonr(predicted_gen[:],heldout_gen)[0]
print "Generalized R^2: ", pearsonr(predicted_gen[:],heldout_gen)[0]**2

fig.set_size_inches(12,5)
fig.tight_layout()


# across arm generalization

# In[32]:


#Compare unimanual contra to unimanual ipsi (descriptive for across arm generalization)
pat = 20
elec = 4
train_cond = 'ipsi' #training hand

lw = 3

contra_sig = uni_pat_avgs[pat]['contra'][elec,:]
ipsi_sig = uni_pat_avgs[pat]['ipsi'][elec,:]

plt.figure
plt.plot(contra_sig, linewidth = lw)
plt.plot(ipsi_sig, linewidth = lw)
plt.title("patient {}, {}({})".format(pat,peaks[pat]['chans'][elec],elec))
plt.legend(['contra','ipsi'], loc = 'lower right')
plt.margins(x=0)

print "R: ", pearsonr(contra_sig, ipsi_sig)[0]
print "R^2: ", pearsonr(contra_sig, ipsi_sig)[0]**2


# In[35]:


### PREDICTED VS HELD OUT FOR ACROSS ARM GENERALIZATION ###
if train_cond == 'contra':
    test_cond = 'ipsi'
else:
    test_cond = 'contra'

if pat in [31,33,38]: #flex_pats
    emg_type = 'flex'
else:
    emg_type = 'fdi'
    
if train_cond == 'contra':
    cond_idx = 0
else:
    cond_idx = 1


#load in train uni test uni data
reg_id = r'/Users/owen/Box/DBS_data/results/BI0{}/BI0{}_10hz_{}_move'.format(pat,pat,emg_type)
reg_dict = cpickle_load(reg_id)


#load in generalization data
gen_id = r'/Users/owen/Box/DBS_data/results/uni_uni_gen/train_{}/BI0{}_uni_uni_{}_{}'.format(train_cond,pat,train_cond,emg_type)
gen_dict = cpickle_load(gen_id)


# In[36]:


### PLOT ###

#try different test sets (0->4)
gen_quality = 'poorly' #generalization quality 'well', 'median' or 'poorly'
reg_color = 'darkblue'
gen_color = 'darkgreen'
hand_idx = cond_idx #because this is only unimanual

#choose the test sets
if gen_quality == 'well':
    #use the best looking test set (most highly correlated) for both the regular and generalized timeseries
    reg_test_set = np.argmax(reg_dict['cc_test'][hand_idx, cond_idx,:, elec])
    gen_test_set = np.argmax(gen_dict['cc_test_gen'][:,elec])

elif gen_quality == 'median':
    #use the median test set (sorted by cc) for both the regular and generalized timeseries
    reg_test_set = np.argsort(reg_dict['cc_test'][hand_idx, cond_idx,:, elec])[2] #2 = median index
    gen_test_set = np.argsort(gen_dict['cc_test_gen'][:,elec])[2]
        
elif gen_quality == 'poorly':
    #use the best looking test set (most highly correlated) for the regular timeseries and 
    reg_test_set = np.argmax(reg_dict['cc_test'][hand_idx, cond_idx,:, elec])
    #use the worst test set for generalized timeseries
    gen_test_set = np.argmin(gen_dict['cc_test_gen'][:,elec])

#save test sets and get the time to plot
predicted_reg = zscore(reg_dict['prediction_test'][hand_idx, cond_idx, reg_test_set, elec])
heldout_reg   = zscore(reg_dict['neural_test'][hand_idx, cond_idx, reg_test_set, elec])

predicted_gen = zscore(gen_dict['prediction_test_gen'][gen_test_set,elec])
heldout_gen   = zscore(gen_dict['neural_test_gen'][gen_test_set,elec])

if min(len(predicted_reg), len(predicted_gen)) > 4000:
    time = 4000
else:
    time = min(len(predicted_reg), len(predicted_gen))

    
fig, axs = plt.subplots(2, sharey=True)
#plot regular pred vs held out
axs[0].plot(predicted_reg[0:time],   color = reg_color,   linewidth = 5, alpha = .5)
axs[0].plot(heldout_reg[0:time], color = 'black', linewidth = 4, alpha = .8)
axs[0].set_ylim([-4,4])
axs[0].set_title('BI0{}({}), {} unimanual, testset={}'.format(pat,elec,train_cond,reg_test_set))
axs[0].margins(x=0)

print "Unimanual R:",pearsonr(predicted_reg[:],heldout_reg)[0]
print "Unimanual R^2:", pearsonr(predicted_reg[:],heldout_reg)[0]**2

#plot generalization pred vs held out
axs[1].plot(predicted_gen[0:time], color = gen_color,   linewidth = 5, alpha = .5)
axs[1].plot(heldout_gen[0:time],   color = 'black', linewidth = 4, alpha = .8)
axs[1].set_ylim([-4,4])
axs[1].set_title('BI0{}({}), {} to {} generalization, testset={}'.format(pat,elec,train_cond,test_cond,gen_test_set))
axs[1].margins(x=0)

print "Generalized R: ", pearsonr(predicted_gen[:],heldout_gen)[0]
print "Generalized R^2: ", pearsonr(predicted_gen[:],heldout_gen)[0]**2

fig.set_size_inches(12,5)
fig.tight_layout()


# In[ ]:




