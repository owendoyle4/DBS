#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import h5py
from scipy.stats.stats import pearsonr
from scipy.stats import zscore
import os

np.set_printoptions(precision=4, suppress=True)

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


# #######################################################################################
# pat = 31 #CHOOSE PATIENT NUMBER
# type_emg = 'flex'
# data_filename = '10hz' #add
# name_end = 'move'
# ########################################################################################
# train_cond ='contra'
# #condition######################################################################
# # train_cond = 0 # c:0 or i:1  # CHANGE THIS TO IPSI OR CONTRA 
# ####################################################################################

# # train_hand = train_cond


# In[3]:


def train_uni1_test_uni2(train_cond, pat, data_filename, type_emg, name_end):
    #######################################################################################
    # Train on the ipsilateral data and test on contralateral data, or vice versa.
    # train_cond = 'contra' or 'ipsi'
    # pat = patient ID
    # data_filename = '10hz'
    # type_emg = 'FDI' or 'flex'
    # name_end = further file identification
    #######################################################################################
    
    # init dictionaries
    weights_reg = {}
    cc_test_gen = {}
    prediction_test_gen = {}
    neural_test_gen = {}

    # depending on input 'train_cond' - choose condition
    if train_cond == 'contra':
        train_cond = 0
        test_cond  = 1
    elif train_cond == 'ipsi':
        train_cond = 1
        test_cond  = 0

    train_hand = train_cond
    test_hand  = test_cond

    # LOAD DATA with preprocessed EMG
    filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/data_{}.mat'.format(pat,data_filename)
    test = {}
    with h5py.File(filepath,'r') as f:
        for k, v in f.items():
            test[k] = np.array(v, dtype = np.float32)
    dbs_data = test['data']
    dbs_data = dbs_data.T
#     print filepath

    # load LMP (which was extracted/processed with extracting_DBS_power.m)
    filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/dbs_lmp_1k.mat'.format(pat)
    test = {}
    with h5py.File(filepath,'r') as f:
        for k, v in f.items():
            test[k] = np.array(v, dtype = np.float32)
    dbs_lmp = test['lmp']
    dbs_lmp = dbs_lmp.T

    neural_feat = dbs_lmp #

    # load move times
    filepath = r'/Users/owen/Box/DBS_data/sdata_files/BI0{}/move_BI0{}_temp.mat'.format(pat,pat)
    test = {}
    with h5py.File(filepath,'r') as f:
        for k, v in f.items():
            test[k] = np.array(v, dtype = np.int)
    move_times = test['move_BI0{}_temp'.format(pat)]
    move_times = move_times.T

    # EMG LAGS
    num_emg_lags = 1000
    num_hand = 2
    # separate emg channels from neural channels 
    emg = dbs_data[:,14:16] # fdi for pats with fdi & flex

    num_samples = emg.shape[0] 

    #initialize lags
    emg_lags = np.zeros((num_hand,num_emg_lags,num_samples))
    emg_lags_padded = np.vstack((np.zeros((500,2)),emg[:-500,:]))

    # adding time lags to the emg data : 500 lags before
    for hand in range(num_hand):
         for lag in range(num_emg_lags):
            emg_lags[hand,lag,:] = np.hstack([emg_lags_padded[1000-lag:,hand],np.zeros(1000-lag)])

    emg_lags = emg_lags.T

    #initialize#
    num_elec = dbs_lmp.shape[1]
    neural_mat = []
    emg_mat = []

    # number of mutually exclusive test sets (each comprising of 20% of the data)
    num_test_set = 5
    test_set_run = [[0,1,2,3,4],[1,2,3,4,0],[2,3,4,0,1],[3,4,0,1,2],[4,0,1,2,3]]

    cc_test_gen = np.zeros((num_test_set,num_elec)) # correlation coefficient across folds and test sets


    # MASK FOR TEST HAND/COND
    data_length = dbs_data.shape[0] #length of all conditons
    test_mask = np.zeros(data_length, bool) #intialize as False
    find_test_cond = np.where(move_times[:,0] == test_cond) # find bimanual blocks in move_times
    find_test_cond = find_test_cond[0]


    #loop through bi blocks
    for block in find_test_cond:
        start = move_times[block,1]
        stop = move_times[block,2]
        test_mask[start:stop] = True
    mask = test_mask # TEST MASK (all blocks of the testing condition/hand)

    neural_mat = neural_feat[mask,:]

    # get TRAIN UNI,TEST UNI weights
    change_id = r'/Users/owen/Box/DBS_data/results/BI0{}/BI0{}_{}_{}_{}'.format(pat,pat,data_filename,type_emg,name_end)
    uni_result_dict = cpickle_load(change_id)
    weights_reg = uni_result_dict['weights_reg']

    # USE TRAINING WEIGHTS ON TESTING EMG IN ENCODING MODEL
    neural_feat = dbs_lmp
    emg_hand = emg_lags[mask,:,test_cond]

    #### CREATING TEST SET ####
    for test_set in range(num_test_set):

        test_run = test_set_run[test_set]

        #split neural and emg data into 5 equal test sets
        split_data_neural = np.array_split(neural_mat, num_test_set)
        split_data_emg = np.array_split(emg_hand, num_test_set)

        ### for EMG ###
        test_temp_emg = []
        test_temp_emg.extend(split_data_emg[test_run[4]]) #last fold = index 4 of test_run (held out test set)

        test_emg = np.squeeze(test_temp_emg)

        ### for neural ###
        test_temp_neural = []
        test_temp_neural.extend(split_data_neural[test_run[4]])

        test_neural = np.squeeze(test_temp_neural)

        #using previously saved TRAIN UNI, TEST UNI weights to predict bimanual
        for elec in range(test_neural.shape[1]):

            # loading weights from when we trained w/ UNIMANUAL 
            uni_weights = weights_reg[train_hand,train_cond,test_set,elec]

            # create the prediction for the test set
            # SAME AS LINE BELOW: prediction_test_gen[test_set,elec] = np.squeeze(np.sum([test_emg*uni_weights],axis = 2))
            prediction_test_gen[test_set,elec] = np.sum(test_emg*uni_weights,axis = 1)
            neural_test_gen[test_set,elec] = test_neural[:,elec]


            # record correlation from prediction set
            cc_test_gen[test_set,elec] = columnwise_correlation(prediction_test_gen[test_set,elec], neural_test_gen[test_set,elec], zscorea=True, zscoreb=True, axis=0)

    return cc_test_gen, prediction_test_gen, neural_test_gen


# In[4]:


# ID_list = [20,25,29,35,39,40,41,99] #FDI 
ID_list = [31,32,33,37,38,39,40] #Flex

data_filename0 = '10Hz' #data_filename to load
type_emg0 = 'flex' #MAKE SURE THE LIST YOU INSERT IS THE SAME EMG TYPE
name_end0 = 'move'
print data_filename0 , type_emg0, name_end0


# In[5]:


for train_condition in ['contra', 'ipsi']:
    for curr_pat in range(len(ID_list)):
        pat_id = ID_list[curr_pat]

        [cc_test_gen, prediction_test_gen, neural_test_gen ] = train_uni1_test_uni2(train_cond = train_condition, pat = pat_id, data_filename = data_filename0, type_emg = type_emg0, name_end = name_end0)

        ## Make dictionary to save our vars ##
        name = 'uni_uni'

        os.chdir(r'/Users/owen/Box/DBS_data/results/uni_uni_gen/train_{}'.format(train_condition))
        uni_uni_dict = {}
        uni_uni_dict = {'cc_test_gen':cc_test_gen, 'prediction_test_gen':prediction_test_gen, 'neural_test_gen':neural_test_gen}

        ########################################################################################333
        uni_uni_file_name = 'BI0{}_{}_{}_{}'.format(pat_id, name, train_condition, type_emg0 ) ##########################################################3
        print uni_uni_file_name
        cpickle_dump(uni_uni_file_name, uni_uni_dict)
    

