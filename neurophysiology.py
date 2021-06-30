# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 11:42:17 2019

@author: mahajnal
"""

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

import quantities as pq
import neo
#import elephant as el

import sklearn



# preprocess utils


def medianfromhistogram(x,f):
    # x are points in the value
    # f are the frequencies of those values
    # will iteratively find the median, and returns median and the cdf
    cdf_normalizer = np.sum(f)
    cdf = 0
    found = False
    for i,xi in enumerate(x):
        cdf += f[i]/cdf_normalizer
        if not found and cdf>0.5: m = xi; found=True
    return m,cdf


def meanfromhistogram(x,f):
    # x are points in the value
    # f are the frequencies of those values
    # returns the mean value in x
    m = x @ f / np.sum(f)
    return m,None




def dataindexhistogram(x,bins):
    di = []
    for bx,b in enumerate(bins[:-1]):
        di.append( np.where( np.logical_and(x>=bins[bx],x<bins[bx+1] ) )[0] )
    return di
        



def time2index(t,dt):
    i = int((t/dt).magnitude)
    return i


def downsample(a,r):
    # print(a.shape)
    # print(a.reshape(-1, r).shape)
    # print(a.reshape(-1, r).mean(axis=1).shape)
    s = a.shape
    newshape = [s[0]//r, r]
    newshape.extend(s[1:])
    a_r = a.reshape( newshape ).mean(axis=1)
    return a_r







def smooth(x,kernelwidth=None,mode='valid'):
    if kernelwidth==None:
        kernelwidth=3
    k = create_kernel(kernelwidth,1*pq.dimensionless)
#    k = np.arange(kernelwidth)
#    k = np.hstack((k,kernelwidth,k[::-1]))+1
#    k = k/k.sum()
    s = np.convolve(x.ravel(),k,mode=mode)
    return s


def decay(data,axs=None, smoothing = True, method='linear'):
    # method can be 'loglin' or 'exp'

    if smoothing: s = smooth(data,7)
    else: s = data   # no smoothing
    

    if method=='linear':      #   y = m*x + c
        y = s
        a = np.arange(len(y))
        A = np.vstack([a, np.ones(len(a))]).T
        m, c = np.linalg.lstsq(A, y, rcond=None)[0]
        b = 0

    elif method=='loglin':
        y = np.log(s)  # transform to log
        #y = x
        a = np.arange(len(y))
        A = np.vstack([a, np.ones(len(a))]).T
        m, c = np.linalg.lstsq(A, y, rcond=None)[0]
        b = 0
    
    elif method=='expl2':
        def exp_param(x, t, y):
            nu = 1/len(x)
            return        x[2] + x[1] * np.exp(x[0] * t)    -    y    +   nu * x.T @ x
        
        t = np.arange(len(s))   # parametrize the domain
        x0 = np.array([0.5,1,0.5])      # initialize the parameters
        res_robust = sp.optimize.least_squares(exp_param, x0, loss='soft_l1', f_scale=0.1, args=(t,s))
        m,c,b = res_robust.x

    elif method=='exp':
        def exp_param(x, t, y):
            return        x[2] + x[1] * np.exp(x[0] * t)    -    y
        
        t = np.arange(len(s))   # parametrize the domain
        x0 = np.array([0.5,1,0.5])      # initialize the parameters
        res_robust = sp.optimize.least_squares(exp_param, x0, loss='soft_l1', f_scale=0.1, args=(t,s))
        m,c,b = res_robust.x


    if not axs==None:    
        axs.plot(a,y,alpha=0.2)
        axs.set_ylabel('forward test accuracy')
        axs.set_xticks([0,50,100])
        axs.set_xticklabels(['0','500','1000'])
        axs.set_xlabel('forward timecourse [ms]')

    return m,c,b



def connectedarea(x,i):
    # find connected areas in a signal or image
    x = sp.ndimage.morphology.binary_dilation(x, iterations=i)
    x = sp.ndimage.morphology.binary_erosion(x, iterations=i)
    return x

def removedots(x,i):
    # find isolated points and remove them, while leaving areas unchanged
    x = sp.ndimage.morphology.binary_erosion(x, iterations=i)
    x = sp.ndimage.morphology.binary_dilation(x, iterations=i)
    return x




def interpolate_switchoffs(x,d_th=200,np_th=10, sm_l=5):
    # filter out the signal errors, primarily from eye pupil size data, which is with 2 ms sample rate and around 7000 values
    # d_th differencial threshold, difference in value
    # np_th neighbourhood proximity threshold, number of neighbors in index
    # sm_l width of smoothing in index
    difftrial = np.diff(np.r_[x[0],x])
    bad_start_idc = np.where(difftrial<-d_th)[0]
    bad_end_idc = np.where(difftrial>d_th)[0]
#        bad_flat_idc = np.where(np.abs(difftrial)<1)[0]
    bad_idc = np.union1d( bad_start_idc, bad_end_idc )
#        bad_idc = np.union1d( bad_idc, bad_flat_idc )
    bad_clusters = []
    for tx,_ in enumerate(bad_idc):       # this will only start if there are any bad points
        if tx==0:              # start with the first element of the first cluster
            cluster = [bad_idc[tx]]
            continue
        if bad_idc[tx] - bad_idc[tx-1] < np_th:       # check if it's in the neighbourhood of previous ones
            cluster.append(bad_idc[tx])
        else:
            bad_clusters.append(cluster)
            cluster = [bad_idc[tx]]          # create a new cluster
        if tx==len(bad_idc)-1:                    # if this is the last element of the bad points, put up the cluster
            bad_clusters.append(cluster)

    for cluster in bad_clusters:
        c0 = max(0,cluster[0]-1)
        c1 = min(cluster[-1]+1,len(x)-1)
        x[c0:c1] = x[c0] + (x[c1] - x[c0]) / (c1-c0) * ( np.linspace(0,c1-c0,c1-c0) )


    # smooth
    x_smoothed = smooth(x,kernelwidth=sm_l)
    # x[sm_l*3:-sm_l*3] = x_smoothed[sm_l*3:-sm_l*3]    # no need for this since 'valid' filling method in numpy convolve

    # avoid growing the boundary effect:
    for j in range(sm_l*3):
        x[j] = np.mean(x[:j*2])
        x[-j-1] = np.mean(x[-j*2-1:])

    return x







def localminimamaxima(x,sm=5):
    xs = smooth(x,sm)
    
    mn = np.where(np.diff(np.sign(np.diff(xs)))==+2)[0]
    mx = np.where(np.diff(np.sign(np.diff(xs)))==-2)[0]
    
    return mn,mx






def detectpeaks(x):
    w = np.array([-1.-3,30.,-3,-1.])
    w = w / np.sum(w)
    y = np.convolve(x,w,mode='same')
    return y



def sweep_disjunct_peaks(x,th=1,r=5):
    # x signal
    # th threshold in units of x
    # r range to check refractory after peaks in units of sampling period
    # returns p: 0 1 0 0 1 0    1 at peak times
    p = np.zeros(len(x))
    for i in range(len(x)):
        p[i] = np.abs(x[i])>th
        if np.any(p[i-min(i,r):i-1]==1): p[i] = 0
    return p





def create_kernel(window_width,sampling_period):
    s = 3*window_width/2 # std of gauss; equals the 99% CDF;
    hw = 6*s/2    # half_width for the kernel cutoff
    kt = np.arange(-hw,hw+sampling_period,sampling_period)*sampling_period.units  # the kernel timepoints
    g = 1./(np.sqrt(2*np.pi)*s) * np.exp(  - (kt/s)**2  )
    kernel = g/g.sum().magnitude/sampling_period.magnitude      # norm the kernel with the cutoff and dt
    return kernel


def instantanousrate(spiketimes,window_width,sampling_period):

    kernel = create_kernel(window_width,sampling_period)
    
    # time points in time dimension
    t = np.arange(spiketimes.t_start,spiketimes.t_stop+sampling_period,sampling_period)
    
    
    spikes = np.zeros(int(np.round( (spiketimes.t_stop - spiketimes.t_start)/sampling_period ) + 1))


    for spike in spiketimes.time_slice(spiketimes.t_start, spiketimes.t_stop):
        index = int(np.round((spike - spiketimes.t_start)/sampling_period))
        spikes[index] += 1


    r = np.convolve(spikes, kernel, 'same') * pq.kHz.rescale('Hz')


    return r,t



def countspikes_quantity(spiketimes,t_start,t_end,resolution=1*pq.ms):
    # create a resolution width bin spike count from spike times input; all must be pq.ms or other time quantity
    # this works for any quantity object
    spiketimes = spiketimes.rescale('ms')
    spikecounts = np.zeros( ( ( t_end.rescale('ms')-t_start.rescale('ms') ) / resolution.rescale('ms') ).magnitude.astype(np.int32) )
    # spikecounts[ (spiketimes/resolution).magnitude.astype(np.int64) ] = 1      # cannot account for multiple spikes in a single bin
    for t in spiketimes:           # slower but can account for multiple spikes in a single bin
        spikecounts[   (t/resolution).magnitude.astype(np.int64) ] += 1
    
    return spikecounts




def countspikes(spiketimes,window_width,sampling_period,binsonly=False):
    # this works only for neo spiketimes signal type

    t = np.arange(spiketimes.t_start,spiketimes.t_stop+window_width,window_width)
    c,bins = np.histogram(spiketimes,bins=t,\
                          range=(spiketimes.t_start,(spiketimes.t_stop+window_width)))
    
#    c,bins = np.histogram(spiketimes,bins=int(np.round((spiketimes.t_stop-spiketimes.t_start)/window_width+1)),\
#                          range=(0,spiketimes.t_stop.magnitude))
    # c = (c/window_width).rescale('Hz')
    # c = (c).rescale('Hz')
    c = c * pq.dimensionless
    # print(window_width,c)

    
#    print(int(np.round((spiketimes.t_stop-spiketimes.t_start)/window_width+1)),bins,t)
    if binsonly:
        return c
    else:
        return c,t




def makecontinuoussignal(block,T,continuous_method='count',species='monkey',normalize=True):
    print('recalculating continuous signals for %s: %s...'%(species,continuous_method))

    n_channels = len(block.segments[0].spiketrains)
    if species=='monkey':
        channels_label = ' MUA '
    elif species=='mouse':
        channels_label = ' SUA '
    if continuous_method=='instfr' or continuous_method=='ks2ifr' or continuous_method=='ks3ifr': sampling = T['dt']
    elif continuous_method=='count' or continuous_method=='ks2count' or continuous_method=='ks3count': sampling = T['bin']
        
    chasigs = []
    chasigse = []
    for trx,trial in enumerate(block.segments):
        sigs = []
        for chx in block.channel_indexes[0].index:
            if continuous_method=='instfr' or continuous_method=='ks2ifr' or continuous_method=='ks3ifr':       # instantenous firing rate with kernel smoothing
                sig,t = instantanousrate(  spiketimes=trial.spiketrains[chx],\
                                           window_width=T['bin'],sampling_period=T['dt']  )
                if trx%5==0 and chx==0: print('trial:',trx)
            elif continuous_method=='count' or continuous_method=='ks2count' or continuous_method=='ks3count':      # spike count in disjoint bins
                sig,t = countspikes(  spiketimes=trial.spiketrains[chx],\
                                           window_width=T['bin'],sampling_period=T['dt']  )
#                np.histogram(trial.spiketrains[chx],bins=int(T['endtime']/T['dt']),range=(0,T['endtime'].magnitude))
                if trx%50==0 and chx==0: print('trial:',trx)
            else: return None
            sigs.append(sig)
        asig = neo.AnalogSignal(np.asarray(sigs).T, t_start=T['starttime'],\
             name='%d'%n_channels+channels_label+'continuous %s signals, segments'%continuous_method,\
             sampling_period=sampling, units=pq.Hz)
        asig.channel_index = block.channel_indexes[0]
        trial.analogsignals.append(asig)
        chasigs.append(asig)
        chasigse.extend(asig)
    block.channel_indexes[0].analogsignals.append(neo.AnalogSignal(np.asarray(chasigse),\
                         name='%d'%n_channels+channels_label+'continuous signals, full length',\
                         t_start=T['starttime'], sampling_period=sampling, units=pq.Hz))
    
    
        
    if normalize:
        m,s,_,_ = getfiringratestatistics(block,T)
        # change to zscores, with a stimulus-epoch-exluded baseline signal to calculate the mean
        # change both the single trial and the channel-concatenated representation:
        block.channel_indexes[0].analogsignals[0] = (block.channel_indexes[0].analogsignals[0] - m)/s
        for trial in block.segments:
            trial.analogsignals[0] = (trial.analogsignals[0] - m)/s

    return block







# calculate brain depth, works for mouse silicon probes:
def identifylayers(block,T,depths):
    n_neurons = block.segments[0].analogsignals[0].shape[1]

    edges = np.arange(T['starttime'],T['endtime']+T['dt'],T['dt'])
    n_bins = len(edges)-1
    
    counts = np.zeros((n_neurons,n_bins))
    counts_smooth = np.zeros((n_neurons,n_bins))
    firststimuluspeak = np.zeros((n_neurons)) * pq.ms
    layer = -np.ones((n_neurons),dtype=np.int16) # default unknown values
    
    for n in range(n_neurons):
        for trial in block.segments:
            count, _ = np.histogram(trial.spiketrains[n],bins=edges)
            counts[n,:] += count
        firststimuluspeak[n] = (np.argmax( counts[n,T['stimstart_idx']:] ))*T['dt']
        if depths[n]>250 and depths[n]<=500 and firststimuluspeak[n]<80*pq.ms: layer[n] = 1
        elif depths[n]>500: layer[n] = 2
    
    # counts_smooth[n,:] = np.convolve(counts[n,:],np.ones(5)/5,'same')
    # print(firststimuluspeak)
    
    
    # plot PSTHs
    # plt.figure()
    # plt.plot(counts.T)
    # plt.figure()
    # plt.plot(counts_smooth.T)
    
    # plot first stimulus peak versus cortical depth
    # fig,axs = plt.subplots(1)
    # axs.scatter(firststimuluspeak,depths)
    # axs.invert_yaxis()
    # axs.set_xlim(0,200)
    # axs.set_xlabel('maximum peak after stimulus [ms]')
    # axs.set_ylabel('cortical depth [$\mu$m]')
            
    
    return layer



def identifycelltype(n_neurons,T,templates,returnfullfeatures=False):

    # templates = np.load(datafilepath+'templates.npy')
    # n_neurons = block.segments[0].analogsignals[0].shape[1]

    troughtopeak_time = np.zeros((n_neurons))
    troughtopeak_amplituderatio = np.zeros((n_neurons))
    celltypes = np.zeros((n_neurons),dtype=np.int16)
    for n in range(n_neurons):
        troughtopeak_time[n] = np.abs(np.argmax(templates[n]) - np.argmin(templates[n]))
        troughtopeak_amplituderatio[n] = np.abs(np.max(templates[n])/np.min(templates[n]))
    
    celltypes = (troughtopeak_time>10).astype(np.int16)           #  1 = broad spiking, 0 = narrow spiking

        
    # plot templates
    # for n in range(n_neurons):
    #     plt.figure()
    #     plt.plot(templates[n])
    #     plt.title('neuron %d'%(n+1))


    # plot clusters
    # fig,axs = plt.subplots(1)
    # axs.scatter(troughtopeak_time,troughtopeak_amplituderatio)
    # axs.set_ylim(0,1)
    # axs.set_xlabel('trough to peak time [ms]')
    # axs.set_ylabel('trough to peak amplitude ratio')

    if returnfullfeatures:
        return troughtopeak_time,troughtopeak_amplituderatio
    else:
        return celltypes





def getfiringratestatistics(block,T):
    
    channels = pq.Quantity([ trial.analogsignals[0] for trial in block.segments ], units=block.segments[0].analogsignals[0].units)
#    print(block.segments[0].analogsignals[0].units, channels.units)
#    print(np.hstack( (channels[:,T['epochstart_idx']:T['stimstart_idx'],:], channels[:,T['epochend_idx']+2:-T['trim_idx'],:])).units)
#    print(T['epochstart_idx'],T['stimstart_idx'],T['stimend_idx'],T['epochend_idx'])
    baselines = pq.Quantity(np.hstack( (channels[:,T['epochstart_idx']:T['stimstart_idx'],:], channels[:,T['epochend_idx']+2:-T['trim_idx'],:])).magnitude, units=channels.units)
    evokedlines = channels[:,T['stimstart_idx']:T['stimend_idx'],:]

    m_b = np.mean(baselines,axis=(0,1))
#    print(channels.units)
#    print(baselines.units,evokedlines.units,m_b.units)
#    print( block.channel_indexes[0].analogsignals[0][3,:],block.channel_indexes[0].analogsignals[0][3,:])
#    print(  np.sum(  (block.channel_indexes[0].analogsignals[0]-m_b)**2,axis=0).units)
    s_b = np.sqrt( np.sum(  (block.channel_indexes[0].analogsignals[0]-m_b)**2,axis=0)\
                    /block.channel_indexes[0].analogsignals[0].shape[0] )

    m_e = np.mean(evokedlines,axis=(0,1))
    s_e = np.sqrt( np.sum(  (block.channel_indexes[0].analogsignals[0]-m_e)**2,axis=0)\
                    /block.channel_indexes[0].analogsignals[0].shape[0] )
    
    

    return m_b,s_b,m_e,s_e











# handle bad channels

def removechannelsfromsession(block,channels_to_remove):
    for tx,trial in enumerate(block.segments):
        block.segments[tx].analogsignals[0] = removechannels(block.segments[tx].analogsignals[0], channels_to_remove)
    block.channel_indexes[0].analogsignals[0] = removechannels(block.channel_indexes[0].analogsignals[0], channels_to_remove)
    return block


def keepchannelsinsession(block,channels_to_tokeep):
    for tx,trial in enumerate(block.segments):
        block.segments[tx].analogsignals[0] = keepchannels(block.segments[tx].analogsignals[0], channels_to_tokeep)
    block.channel_indexes[0].analogsignals[0] = keepchannels(block.channel_indexes[0].analogsignals[0], channels_to_tokeep)
    return block



def removechannels(asig,channels_to_remove):
    # asig is a timecourse x channels AnalogSignal array
    # channels_to_remove starts with 1, not with 0
#    if type(asig)!=neo.core.analogsignal.AnalogSignal: raise(TypeError)
    channels = list(np.arange(asig.shape[1],dtype='int16')+1)
    for ch in channels_to_remove:
        channels.remove(ch)
    return asig[:,  np.array(channels,dtype='int16')-1 ]    # convert to 0 start index

def keepchannels(asig,channels_to_keep):
    # asig is a timecourse x channels AnalogSignal array
    # channels_to_keep starts with 1, not with 0
#    if type(asig)!=neo.core.analogsignal.AnalogSignal: raise(TypeError)
    return asig[:,  np.array(channels_to_keep,dtype='int16')-1 ]      # convert to 0 start index









# Fourier spectrum utilities


def removefrequencies(x,frequencies,f_s,filtertype,analog=False):
    w0=np.asarray(frequencies)/(f_s/2.)              # frequency(ies) stats.normalized to the nyquist freq of sampling
    # we will also need an initial condition for the filtered data, to preserve phase information
#    print(w0)
    if filtertype == 'bp':
        b, a = sp.signal.butter(4,w0,'bandpass',analog=analog)                 # bandpass filter, 
        y = sp.signal.filtfilt(b,a,x)#,zi=sp.signal.lfilter_zi(b, a)*x[0] )
    else:
        y = x         # we will iteratively apply the requested filters one by one, if not bandpass
        for k,f in enumerate(frequencies):
    #        print(w0)
            
            if filtertype == 'n':
                b, a = sp.signal.iirnotch(w0,35.,analog=analog)             # single frequency cut
            if filtertype == 'hp':
                b, a = sp.signal.butter(3,w0,'highpass',analog=analog)      # highpass filter for DC remove
            if filtertype == 'lp':
                b, a = sp.signal.butter(3,w0,analog=analog)                 # lowpass filter for noise remove
            
            y = sp.signal.filtfilt(b,a,y)#,zi=sp.signal.lfilter_zi(b, a)*y[0] )
    
    if type(x)==pd.Series:
        y = pd.Series(y, index=x.index)
    return y





def downsamplesignals(signalcollection, divider):
    # signalcollection is a list of classes, each a list of trials each trial analogsignal:
    # divider is how many signals to 

    n_trajectory_orig = signalcollection[0][0].shape[0]
    n_trajectory_downsampled = divider*(n_trajectory_orig//divider)

    signalcollection_downsampled = [    [  \
        neo.AnalogSignal( downsample(trial[:n_trajectory_downsampled,:], divider),\
            name='averaged '+trial.name,t_start=trial.t_start,sampling_period=trial.sampling_period*divider,units=trial.units)\
               for trial in signalcollection[clx]        ] \
                  for clx in range(len(signalcollection))   ]    
    
    return signalcollection_downsampled




def shift_timecourse_analogsignal(signalcollection,shift):
    # shifts the entire signal, 'shift' in units of time
    
    s = int(shift / signalcollection[0][0].sampling_period)
    # print(signalcollection[0][0].sampling_period,shift,s)
    
    signalcollection_shifted = [ [ np.roll(trial,s,axis=0) \
                                       for trial in signal ] \
                                           for signal in signalcollection ]


    return signalcollection_shifted




def convertintervaltofeature(r,w=5,d=5):
    # converts a portion of the next timepoints into extended feature space, each neuron has forward time activity as well
    # e.g. d=1,w=5: 600 timepoints x 4 neurons   --->   600 timepoints x 20 features, or with d=5: 120 timepoints x 20 features

    n_classes = len(r)
    n_timepoints,n_signals = r[0][0].shape


    r_tf = []
    
    for clx in range(n_classes):
        n_trials = len(r[clx])
        
        c = np.zeros((n_trials,(n_timepoints-w)//d+1,n_signals*w))
        
        for tx,t in enumerate(np.arange(0,n_timepoints-w,d)):
            c[:,tx,:] = np.array(r[clx])[:,t:t+w,:].reshape((n_trials,1,w*n_signals)).squeeze()


        # c = np.array(r[clx])[:,:n_timepoints//w*w,:].reshape( (n_trials*w,n_timepoints//w,n_signalchannels), order='F' )
        

        # return as analogsignals
        c = [ neo.AnalogSignal(  trial, name='features and downsample '+r[clx][0].name,t_start=r[clx][0].t_start,sampling_period=r[clx][0].sampling_period*d,units=r[clx][0].units)\
             for trial in c  ]
    
        r_tf.append(c)
    
    return r_tf





def equalizebyindex(r,indexlist):
    # create for each class a series of masked arrays
    # only possible with index masks, not index addresses!
    # r_eq = [  [  np.array(r[clx])[indexlist[clx,tx],tx,:]  for tx in range(r[clx][0].shape[0])  ]  for clx in range(len(r)) ]
    
    n_observations,n_features = r[0][0].shape
    
    # convert indices into masks for masked array for each class
    mask = [ np.ones((len(r[clx]),n_observations,n_features),dtype=bool) for clx in range(len(r)) ] # [class][trials,observ,feat]
    
    
    for clx in range(len(r)):
        for tx in range(n_observations):
            mask[clx][indexlist[clx,tx],tx,:] = np.zeros((len(indexlist[clx,tx]),n_features))
    
    # print('CENSUS')
    # print(mask[0].sum())
    # print(mask[1].sum())
    
    r_eq = [  np.ma.array(r[clx],mask=mask[clx])  for clx in range(len(r))  ]
    
    

    # print('inroutine',len(r_eq),len(r_eq[0]),r_eq[0][0].shape)
    # r_eq_a = [ [ neo.AnalogSignal(  trial, name='equalized '+r[clx][0].name,\
    #             t_start=r[clx][0].t_start,sampling_period=r[clx][0].sampling_period,units=r[clx][0].units)\
    #                for trial in r_eq[clx] ] for clx in range(len(r_eq))  ]

    analogsignal_properties = {'times':r[0][0].times,'t_start':r[0][0].t_start,\
                               'sampling_period':r[0][0].sampling_period,'units':r[0][0].units}
    
    return r_eq,analogsignal_properties






def virtualizetrials_analogsignal(r,w):
    
    # r: analogsignal collections:  [classes][trials][timecourse,neurons]
    n_classes = len(r)
    n_timepoints,n_signalchannels = r[0][0].shape
    
    r_pb = [] # collect and return
    
    for clx in range(n_classes):
        n_trials = len(r[clx])
        # print(n_trials,np.array(r[clx]).shape, np.array(r[clx])[:,:n_timepoints//w,:].shape)
        # reshape into patched timebins; note that the reshape order is to be done with inner fastest
        c = np.array(r[clx])[:,:n_timepoints//w*w,:].reshape( (n_trials*w,n_timepoints//w,n_signalchannels), order='F' )
        # return as analogsignal
        c = [ neo.AnalogSignal(  trial, name='patched %d timebins'%w+r[clx][0].name,t_start=r[clx][0].t_start,sampling_period=r[clx][0].sampling_period*w,units=r[clx][0].units)\
              for trial in c ]
        r_pb.append(c)
    
    
    return r_pb
    
    

def virtualizetrials_maskedarray(r,w,analogsignal_properties):
    # patch together subsequent timebins to increase number of observations
    # new: r: masked array collections [classes][trials,timecourse,neurons]
    # w: number of timebins to patch together
    # a analogsignal properties; times and sample rate will be updated with the virtualization reshape
    n_classes = len(r)
    n_timepoints,n_signalchannels = r[0][0].shape
    # print(n_classes,n_timepoints,n_signalchannels)

    # print(n_classes,len(r[0]),n_timepoints,n_signalchannels )
    r_pb = [] # collect and return
    
    for clx in range(n_classes):
        n_trials = len(r[clx])
        # print(n_trials,np.array(r[clx]).shape, np.array(r[clx])[:,:n_timepoints//w,:].shape)
        # reshape into patched timebins; note that the reshape order is to be done with inner fastest
        c = r[clx][:,:n_timepoints//w*w,:].reshape( (n_trials*w,n_timepoints//w,n_signalchannels), order='F' )
        # return as analogsignal
        # c = [ neo.AnalogSignal(  c[t], name='patched %d timebins'%w+r[clx][0].name,t_start=r[clx][0].t_start,sampling_period=r[clx][0].sampling_period*w,units=r[clx][0].units)\
        #      for t in range(len(c))]
        # return as masked array
        r_pb.append(c)
    
    a = analogsignal_properties.copy()
    a['times'] = a['times'][::w]
    a['sampling_period'] = a['sampling_period']*w
    
    return r_pb, a











# signal ensembles



def trialaveraging(rl):
    # response should be        trials dimensions    x    timepoints    x          channels
#    print('response shape: ', response.shape)
#    response = response[:,T['epochstart_idx']:T['epochend_idx']+1,:]
#    print('response shape, cut to epoch: ', response.shape)

    response = np.array(rl)
    n_trials = response.shape[0]
    rm = neo.AnalogSignal(  np.mean(response,axis=0),    name='trialaveraged mean'+rl[0].name,t_start=rl[0].t_start,sampling_period=rl[0].sampling_period,units=rl[0].units)
    rs = neo.AnalogSignal(np.std(response,axis=0),name='trialaveraged std'+rl[0].name,t_start=rl[0].t_start,sampling_period=rl[0].sampling_period,units=rl[0].units)
    re = neo.AnalogSignal(rs / np.sqrt(n_trials),name='trialaveraged s.e.m.'+rl[0].name,t_start=rl[0].t_start,sampling_period=rl[0].sampling_period,units=rl[0].units)

    return [rm,rs,re]






def get_auctimenormalized(asig,starttime,endtime,cropzero=True,title=None):
    # check if start and end are indexes or times
    if type(starttime)==pq.Quantity:
        start = int(  ( (starttime-asig.t_start)/asig.sampling_period ).magnitude  )
        end   = int(  ( (  endtime-asig.t_start)/asig.sampling_period ).magnitude  )
    else:
        print('times quantity not found, using as index')
        start = starttime
        end = endtime
    # get rid of negative performance due to noise
    cutasig = asig[start:end,:]
    cutasig.t_start = asig.t_start

    if cropzero:
        cutasig[cutasig<0*cutasig.units] = 0*cutasig.units
    
    aucdensity = np.mean(   cutasig, axis=0)

    #integrate out time to get the average
#    interval = (endtime - starttime)
    
    return aucdensity




def get_aucpositivenegative(m1,m2,s1,s2):
    # return the positive and negative parts of the difference integral with confidence intervals
    # with +:   m1-s1-m2-s2>0
    #      -:   m1+s1-m2+s2<0
    l = len(m1)
    positive = m1-s1-m2-s2
    negative = m1+s1-m2+s2
    p = positive[positive>0].sum()/l
    n = negative[negative<0].sum()/l
    return p,n




def showcontinuousrateexamples():
        
    plt.rcParams.update({'font.size': 24})
    plt.rcParams.update({'legend.fontsize': 14})
    plt.rcParams.update({'lines.linewidth': 1})
    
    spikes = [ [100,105,112,136,156,198,202,210,212,250,300], \
               [100, 232, 251]      ]
    f_labels = ['high','low']
    
    for f_levels in range(2):
        spiketrain = neo.SpikeTrain(np.asarray(spikes[f_levels])*pq.ms,t_start=0*pq.ms, t_stop=400*pq.ms)
    
        plt.figure(figsize=(36,24))
        
        j=0
        for w in [5,20,50]:
            for s in [1,25,100]:
                j+=1
                r,c,t,ct = neph.instantenousrate(spiketrain,w*pq.ms,s*pq.ms)
                plt.subplot(3,3,j)
                plt.eventplot(spiketrain,lineoffsets=-7,linelengths=10,color='mediumvioletred')
                plt.plot(ct,c,linewidth=3,color='darkorange')
                plt.plot(t,r,linewidth=3,color='royalblue')
                plt.title('sampling period %d ms'%s)
                if w==50: plt.xlabel('timecourse [ms]')
                if s==1: plt.ylabel('window %d ms\n\nspike count [count/window %d ms]\nrate [Hz]'%(w,w))
                if j==1: plt.legend(['binned spike count','gauss kernel rate','spike times'])
        
        plt.suptitle('transforming spike trains to rate, %s rate'%f_labels[f_levels])
#        plt.savefig('../results/spikes/spiketorate-illustration,%srate.png'%f_labels[f_levels])



        
if __name__ == '__main__':
    print('neurophysiology analysis methods')