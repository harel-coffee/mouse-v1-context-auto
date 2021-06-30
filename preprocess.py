#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 13:31:31 2019

@author: mahajnal
"""


from config import *


import numpy as np
import scipy as sp
import scipy.signal as signal

import pandas as pd
import h5py
import pickle

import matplotlib.pyplot as plt


import quantities as pq
import neo

# import preprocess as preprocess
import neurophysiology as neph



# I/O utils




def getnneurons(datanames):
    numneurons = {'ME103':8,'ME110':48,'ME113':32,'DT008':10,'DT009':4,'DT014':9,'DT017':34,'DT018':43,'DT019':31,'DT020':12,'DT021':8,'DT022':17,'DT030':6,'DT031':8,'DT032':6}
    return numneurons


def loadexperimentdata(fn):
    fullpath = pathdatamouse+trialsfolder+'trials-start-good-C'+fn+'data.csv'

    trialsdata = pd.read_csv(fullpath,sep=',',\
                 usecols=['start','duration','degree','block','freq','water','punish'])
    
    return trialsdata


def exporttrialtimes(fn):
    g = loadexperimentdata(fn)
    
    trialstarts = g['start'] / 1000 * pq.ms

    fullpath = pathdatamouse+trialsfolder+'trialtimes-'+fn+'.csv'
    trialstarts.to_csv(fullpath,sep=',',header=False,index=False)
    return


def exporthdf5(block,dn,path=cacheprefix+'hdf5/'):
    # file saves for julia flux: timepoints, neurons, trials
    D = np.array([trial.analogsignals[0] for trial in block.segments if trial.annotations['block'] in [2,4]  ])        # this was in trials, timepoints, neurons
    D = np.moveaxis(D,[1,2],[2,1])      # row major
    # D = np.moveaxis(D,[0,1,2],[2,0,1])      # column major
    A = pd.DataFrame([ trial.annotations for trial in block.segments if trial.annotations['block'] in [2,4]   ]).values.astype(np.int16)
    
    # D = np.asfortranarray(D)
    # A = np.asfortranarray(A)
    
    hf = h5py.File(path+'ks2ifr-'+dn+'.h5', 'w')
    hf.close()
    hf = h5py.File(path+'ks2ifr-'+dn+'.h5', 'w')
    hf.create_dataset('sua', data=D, dtype='f')
    hf.create_dataset('meta', data=A, dtype='i')
    hf.close()

    


def loaddatamouse(sessionname='',T=None,continuous_method='count',normalize=True,recalculate=False,exportspiketrains=False):

    normalized = ['o','n'][normalize]      # original, vs. normalized
    
    fst=neo.io.PickleIO(cacheprefix+'neo/neuralactivity-%s_%s,%s-%dms%dms.nio'%(sessionname,continuous_method,normalized,T['dt'].magnitude,T['bin'].magnitude))

    if recalculate:
        
        print('%s recalculating %s from raw...'%(sessionname,continuous_method))
        # data = sp.io.loadmat(pathdatamouse + sessionname + 'data.mat')
        
        ES,NR,CT = mat2panda(sessionname)
        if continuous_method[:3]=='ks2' or continuous_method[:3]=='ks3':        # for the Kilosort2,3: this will hold a # SUA list of spike times

            if exportspiketrains=='load':    # reload from simple unit list of list of spikes saved with exportspiketrains = 'save'
                spiketrains = pickle.load(open(cacheprefix+'phys/spiketrains-%s-%s.pck'%(sessionname,continuous_method),'rb'))
                CT = pickle.load(open(cacheprefix+'phys/cellinfo-%s-%s.pck'%(sessionname,continuous_method),'rb'))
                n_neurons = len(spiketrains)

            else:     # generate spiketrains and metadata from from phy2 export files

                spiketrains,depths,templates = ksphy2neospiketrain(dn=sessionname,T=T)
                n_neurons = len(spiketrains)
                    
                #placeholder of cell identities until calculated
                # CT = {'idents':pd.Series(depths),'waveforms':pd.Series(-np.ones(n_neurons)) }
                # neph.identifylayers(block,T,depths)
                CT = {'depths':np.array(depths),\
                    'waveforms':templates,\
                    'layers':np.array(depths)<200,\
                    'celltypes':neph.identifycelltype(n_neurons,T,templates) }

                if exportspiketrains=='save':  # choose this func argument if want to just export list of all spikes without trial segments; it will quit the function afterwards
                    print('exporting raw spiketrains')
                    pickle.dump(spiketrains,open(cacheprefix+'phys/spiketrains-%s-%s.pck'%(sessionname,continuous_method),'wb'))
                    pickle.dump(CT,open(cacheprefix+'phys/cellinfo-%s-%s.pck'%(sessionname,continuous_method),'wb'))
                    return  # don't start making the neo.analogsignals, this is just an export to file in a simple format


        else:      # instfr is JRClust, and NR holds 0s and 1s, so we will need to find timings with np.where later
            n_neurons = NR.shape[1]
    
        g = loadexperimentdata(sessionname)
        trialstarts = g['start'].values * pq.ms
    
        n_trials = trialstarts.shape[0]

        # n_trials = max_tr       # to test for smaller numbers
        
        block = neo.Block(name='mouse '+sessionname,\
                          depths=CT['depths'],\
                          waveforms=CT['waveforms'],\
                          layers=CT['layers'],\
                          celltypes=CT['celltypes'])
        block.channel_indexes.append(neo.ChannelIndex(name='%d channel SUA continuous signal'%n_neurons, index=np.arange(n_neurons)))
        for trx in range(n_trials):

            trial = neo.Segment(name='trial #%d'%(trx+1),index=trx,description=sessionname+' t%d,b%d'%(trx+1,g['block'][trx]+1),\
                                block=g['block'][trx]+1,visual=g['degree'][trx],audio=g['freq'][trx],\
                                punish=g['punish'][trx]    )     # get the stimulus condition ids
            # epoch nomenclature of neo not used for now:
            # trial.epochs=[]
            # trial.epochs.append(neo.Epoch( times=pq.quantity(T['stimstarttime']),\
                      # durations=pq.quantity(T['stimendtime']-T['stimstarttime']),labels=['stimulus']))

            for chx in range(n_neurons):
                if continuous_method[:3]=='ks2' or continuous_method[:3]=='ks3':
                    spiketimes = spiketrains[chx][ np.logical_and( spiketrains[chx]>=trialstarts[trx]+T['starttime'],  \
                                                                   spiketrains[chx]<=trialstarts[trx]+T['endtime'] ) ] \
                                             - trialstarts[trx]
                else:
                    cutstart = int(( trialstarts[trx]+T['starttime'] ).magnitude)
                    cutend = int(( trialstarts[trx]+T['endtime'] ).magnitude)
                    spiketimes = (np.where( NR.iloc[cutstart:cutend,chx] ) + T['starttime'].magnitude )*pq.ms
                    
    #            print(spiketimes)
                trial.spiketrains.append(neo.SpikeTrain(spiketimes,description='neuron %d'%(chx+1),\
                                                        t_start=T['starttime'], t_stop=T['endtime']) )
            block.segments.append(trial)
        
        # count spikes / smooth with gaussian, and add to analogsignal[0]
        block = neph.makecontinuoussignal(block,T,continuous_method=continuous_method,species='mouse',normalize=normalize)

        print(n_trials,n_neurons)
        print(block.segments[0].analogsignals[0].shape, type( block.segments[0].analogsignals[0] ) )
            
        # add running speed:
        block.channel_indexes.append(neo.ChannelIndex(name='running speed', index=np.array([0])))
        L = len(ES['run'])
        sig = signal.resample(ES['run'],  int(L/(T['dt']/(1*pq.ms)).magnitude) )
        block.channel_indexes[1].analogsignals.append(neo.AnalogSignal(sig,\
                                 name='running speed, full length',\
                                 t_start=0*pq.ms, sampling_period=T['dt'], units=pq.cm/pq.s))

        for trx,trial in enumerate(block.segments):
            cutstart = int(( trialstarts[trx]+T['starttime'] ).magnitude)
            cutend = int(( trialstarts[trx]+T['endtime']+T['dt'] ).magnitude)

            L = len(ES['run'][cutstart:cutend])
            sig = signal.resample(ES['run'].iloc[cutstart:cutend],  int(L/(T['dt']/(1*pq.ms)).magnitude) )

            asig = neo.AnalogSignal(np.asarray(sig), t_start=T['starttime'],\
                                    name='running speed',\
                                    sampling_period=T['dt'], units=pq.cm/pq.s)
            asig.channel_index = block.channel_indexes[1]
            trial.analogsignals.append(asig)

        fst.write_block(block)
    else:
        print('%s loading precalculated %s from raw...'%(sessionname,continuous_method))
        block = fst.read_block()


        if 0:
            # actually calculate layer id
            spiketrains,depths,templates = ksphy2neospiketrain(dn=sessionname,T=T)
            n_neurons = len(spiketrains)
            
    
            CT = {'depths':np.array(depths),\
                  'waveforms':templates,\
                  'layers':np.array(depths)<200,\
                  'celltypes':neph.identifycelltype(n_neurons,T,templates) }
            print(CT)
            return
        
    
    return block







    


def mat2panda(fn):
    # if fn in ['ME103','DT013','DT015','DT020']: mstruct = 'estruct'
    # else: mstruct = 'extracellstruct'
    mstruct = 'extracellstruct/'
    
    objects1 = ['vstim','astim','licking','water','run']
    objects2 = ['unitmatrix']
#    if fn[:2]=='AC': objects3 = ['waveforms']
#    else:
    objects3 = ['idents','waveforms']
        
#    layernames = ['superficial','input','deep5','deep6']
    
    f = h5py.File(pathdatamouse+fn+'data.mat','r')

    r=np.asarray(f.get('/%s/vstim'%mstruct)).ravel()
    # times = pd.TimedeltaIndex(start=0,periods=r.shape[0],freq='ms')     # pandas 0.2x
    times = pd.timedelta_range(start=0,periods=r.shape[0],freq='ms')     # pandas 1.x


    ES = pd.DataFrame({},dtype='f2',index=times)
    for k,o in enumerate(objects1):
        aux = np.asarray(f.get('/%s/'%mstruct+o)).ravel()
        # if o=='run': f.close()
        # print(fn,o,aux.shape)
        ES[o] = aux
#        print( o, aux.shape )



    for k,o in enumerate(objects2):
        aux = np.asarray(f.get('/%s/'%mstruct+o))
#        print( o, aux.shape )
    NR = pd.DataFrame(aux,dtype='f2', index=times)

    
    
    CT = pd.DataFrame()
    
    if fn[:2]=='ME' or fn[:2]=='DT':        # idents and waveforms only present in old JRC->matfile formats; ks2 does not have these yet
        ref = f['/%s/'%mstruct+objects3[0]]
        # print( ref )
        aux = np.zeros(ref.shape[1])
        for l, r in enumerate(ref[0,:]):
            aux[l] = np.asarray(f[r])[0][0]
        
        CT['idents'] = aux

        aux = np.asarray(f.get('/%s/'%mstruct+objects3[1]))
        # print( objects3[1], aux.shape )
        CT['waveforms'] = aux[0]

#    CT['layers']=pd.Categorical.from_codes(aux, layernames, True)


    return ES,NR,CT











def ksphy2neospiketrain(dn,T,driftpruninglevel=0):

    from os.path import isfile

    # LP,V1,ACC
    if dn[:2]=='AC':
        pathdatamouse_brainarea = 'ACC/'       
    elif dn[:2]=='LP':
        pathdatamouse_brainarea = 'LP/'
    else:
        pathdatamouse_brainarea = 'V1/'

    suas = []
    cluster_depths = []
    templates = []
    for sh_idx in [1,2]: #[1,2]: # go through shanks
        sua = []
        if continuous_method[:3]=='ks2':
            datafilepath = pathdatamouse_ks2sorted+pathdatamouse_brainarea+'%sshank%d/'%(dn,sh_idx)
        elif continuous_method[:3]=='ks3':
            datafilepath = pathdatamouse_ks3sorted+pathdatamouse_brainarea+'%sshank%d/'%(dn,sh_idx)

        # check if the file exists, and ignore the shank, if there is no sorting found:
        if not isfile(datafilepath+'spike_times.npy'):
            continue
        
        
        spikes_timing = np.load(datafilepath+'spike_times.npy')
        spikes_incluster = np.load(datafilepath+'spike_clusters.npy')
        spike_templates = np.load(datafilepath+'spike_templates.npy')
        all_templates = np.load(datafilepath+'templates.npy')

        
        cluster_info = pd.read_csv(datafilepath+'cluster_info.tsv',sep='\t',usecols=['id','KSLabel','depth','group'])    # KSLabel or the richer cluster_info.npy contains quality measures as well


        
        # find clusters that are 'good' and not 'noise' and not 'drift'-ing: i.e.exclude manually curated noise clusters and drifts:
        if driftpruninglevel==0:
            good_sua_ix = np.logical_and( cluster_info['KSLabel']=='good',   \
                         np.logical_and( cluster_info['group']!='noise', cluster_info['group']!='drift')    )
        elif driftpruninglevel==1: # more strict criteria, with explicit flatness
            good_sua_ix = np.logical_and( cluster_info['KSLabel']=='good',   \
                         np.logical_and( cluster_info['group']!='noise', cluster_info['group']=='flat')    )
                    
        # print(len(cluster_info),len(good_sua_ix))
        
        cluster_id_list = cluster_info['id'].loc[ good_sua_ix ].values
        
        # cluster_channel = cluster_info['channel'].loc[ good_sua_ix ].values         # does not work in newer phy2
        
        # sort with depths:
        cluster_depth = cluster_info['depth'].loc[ good_sua_ix ].values
        
        
        sortorder = np.argsort(cluster_depth)
        
        cluster_id_list = cluster_id_list[sortorder]
        # cluster_channel = cluster_channel[sortorder]
        cluster_depth = cluster_depth[sortorder]

        template = []        
        for n,cluster_id in enumerate(cluster_id_list):
            # extract spike times for a given sua
            sua.append( ( spikes_timing[np.where(spikes_incluster==cluster_id)] / T['ks2samplingrate'] ).rescale('ms') )
            # print(n,cluster_id,len(sua[n]))
            
            # extract the waveform for this sua
            template_id = np.unique(spike_templates[spikes_incluster==cluster_id])[0]

            closestchannel = np.argmax( np.abs(all_templates[template_id,:,:]).sum(axis=0)  )

            # print( 'n%d, t%d, ch%d, mch%d'%(n,template_id,cluster_channel[n], closestchannel)   )

            template.append(  all_templates[template_id,:,closestchannel]   )

        
        suas.extend(sua)
        cluster_depths.extend(cluster_depth)
        templates.extend(np.array(template))
    
        
    return suas,cluster_depths,templates




    




def exportspiketrainsfromneo(dn,block):
    n_trials = len(block.segments)
    n_neurons = len(block.segments[0].spiketrains)
    print('saving %s, neurons %d, trials %d'%(dn,n_neurons,n_trials))
    session = [  [ neuronspiketrain.magnitude   for neuronspiketrain in trial.spiketrains ] for trial in block.segments ]
    pickle.dump(session,open(cacheprefix+'phys/spiketrains,trialcut-%s.pck'%(dn),'wb'))





def licktrajectories(dn,T,binarize=True):
    # get licking data for a mouse; either original analogue, or binarized version with respect to reward delivery start
    
    prematurereset = 0*pq.ms
    

    ES,NR,CT = mat2panda(dn)
    
    g = loadexperimentdata(dn)
    trialstarts = g['start'].values * pq.ms
    
    n_trials = trialstarts.shape[0]


    lick_channel = ES['licking']
    
    L = len(lick_channel)
    lick_channel = signal.resample(lick_channel,  int(L/(T['dt']/(1*pq.ms)).magnitude) )
    
    
    
    if binarize:
        # preprocess
    #    lick_channel = neph.removefrequencies(lick_channel,[50*pq.Hz],1000*pq.Hz,'lp')
    #    lick_channel = neph.smooth(lick_channel,kernelwidth=40).squeeze()
    #    lick_channel = np.convolve(lick_channel,np.r_[np.ones(30),np.zeros(30)]/60,'same')
    
        # check if significant licking occurs
        lick_channel = lick_channel > 3
        # now connect periodic licking intervals into a binary licking activity signal
        lick_channel = neph.connectedarea(lick_channel,15)
        # remove isolated single licks
        lick_channel = neph.removedots(lick_channel,5)
        
        
    
    
        # put back into series with index as timers
        lick_channel = pd.Series(lick_channel)


    
#    print('resampled',lick_channel.shape)




    # go over each trial
    licks = []
    
    for trx in range(n_trials):
        cutstart = int((( trialstarts[trx]+T['starttime'] )/T['dt']).magnitude)
        cutend = int((( trialstarts[trx]+T['endtime'] )/T['dt']).magnitude)

        lick = lick_channel.iloc[cutstart:cutend]

        if binarize:
            
            # remove licking from prestimulus period
            # find the post-reward lick, and follow backwards until first occurance
            w = np.where(   lick[T['reward_idx']:]  )[0]
            if len(w)>0:              # is there a lick in this period at all?
                w += T['reward_idx']    # find the lick
                wb = np.where( lick[:w[0]]==0 )[0]    # find the first point where this lick occurs by finding the last-pre-nonlick
                if len(wb)>0:
    #                if trx==140: print(w[0]*10-1500,wb[-1]*10-1500)
                    lick[:wb[-1]]=0                  # zero any previous licks: these are flagged as _not_ reward related
                lick[:T['stimstart_idx']+neph.time2index(prematurereset,T['dt'])] = 0
            else:    # if there is no lick at all, signify for later NaN detection
                lick[:] = -1
        
        licks.append(lick)
    
    licks = np.array(licks)
    
    return licks











def get_lickonsets(dn,T,trialset):
    licks = licktrajectories(dn,T,binarize=True)
#    print(licks.shape)
    
    
    # use the followings to show individual per trial lick trakectories
    # import figs
    # fig,ax = plt.subplots(3,3,figsize=(36,24))
    # # np.random.seed(15)
    # tl = (np.random.rand(3,3)*160).astype('int16')
    # tl[2,0] = 6
    # for i in [0,1,2]:
    #     for j in [0,1,2]:
    #         axs = ax[i,j]
    #         axs.plot( np.linspace(-1500,4500,600), licks[tl[i,j],:] )
            
    #         axs.set_ylim(-0.1,5)
            
    #         figs.plottoaxis_stimulusoverlay(axs,T)

    #         axs.set_title('%d: %5.3f '%(tl[i,j],np.min(licks[tl[i,j],:])))
        
    
    # get lick times per trial, and their statistics
    # change to comprehension: [xv if c else yv for c, xv, yv in zip(condition, x, y)]
    lickonsets = []
    for lick in licks[trialset]:
        w = np.where(lick>0)[0]
        # print(type(w),w)
        if len(w):
            lickonsets.append(w[0]*T['dt']+T['starttime'])
        else:
            lickonsets.append(np.NaN)
            
    lickonsets = np.array(lickonsets)
    # print(lickonsets)

    lickonset_m = np.nanmedian(lickonsets)
    lickonset_s = np.nanstd(lickonsets)
    lickonset_e = lickonset_s/np.sqrt( (~np.isnan(lickonsets)).sum() )
    
    
    
    
    
    return lickonsets, lickonset_m, lickonset_s, lickonset_e







def get_lickidents(dn,T,trialset):
    
    waterdeliverytime = T['rewardtime']
    
    lickonsets, _,_ = get_lickonsets(dn,T,trialset)
    
    action = lickonsets >= waterdeliverytime
    
    return action









def licksync_trials(dn,T,trialset,responses,padtime=0*pq.ms):
    
    # constrain the lick times to remain between these endpoints
    border_low = 500*pq.ms                       # time after stimulus, allowing 2000 ms pre lick register
    border_high = T['stimendtime']+1000*pq.ms     # this allows 1000 ms post stimulus time
    # padding is how much padding we want before and after the stimulus
    
    # collect lick onset times
    lickonsets, lickonset_m, lickonset_s, lickonset_e = get_lickonsets(dn,T,trialset)
    
    # padtime=0*pq.ms
    # show the lick onsets on a figure:
#    o,m,s = lickonsets, lickonset_m, lickonset_s
#    
#    print(o)
#    print(o[np.logical_not(np.isnan(o))])
#    
#    
#    fig,ax = plt.subplots(1,1,figsize=(12,8))
#    
#    axs = ax
#    axs.hist(o[np.logical_not(np.isnan(o))], 50)
#    axs.set_xlim(-1500,4500)
#    axs.plot([m,m],[0,30],lw=3,color='mediumvioletred')
#    axs.plot([m-s,m-s],[0,30],lw=3,color='fuchsia')
#    axs.plot([m+s,m+s],[0,30],lw=3,color='fuchsia')
    
    responses_licksynced = responses.copy()
    
    for trx,trial in enumerate(responses_licksynced):
        # if there is no lick sync trials to average lick time (i.e. no timeshift):
        if np.isnan(lickonsets[trx]):
            lickonsets[trx] = lickonset_m
        # sync time to lick; border
        # print(border_low,border_high,lickonsets[trx]*pq.ms)
        licktime = max( ( border_low, min( (border_high,lickonsets[trx]*pq.ms) ) ) )
        # if trx in [133,134,138,139]: print(trx,lickonsets[trx]*pq.ms,licktime,trial.times[0],trial.t_start)
        # DT = licktime - T['rewardtime']             # move with delta-t
        DT = lickonset_m*pq.ms - licktime      # if lick onset is later than the average lick onset advance time (it will happen earlier)
        # if trx in [0,2,3,18]:
        # print(trx,DT,'(',lickonset_m*pq.ms,licktime, ')','>',-T['starttime']-padtime,\
        #        '<',-(T['endtime']-T['stimendtime']-padtime))
        if DT>-T['starttime']-padtime:
            DT = -T['starttime']-padtime
        if DT<-(T['endtime']-T['stimendtime']-padtime):
            DT = -(T['endtime']-T['stimendtime']-padtime)
            
        DTix = -neph.time2index(DT,T['dt'])          # index movement is opposite to time value movement
        # padding before and after stimulus - not possible too long, as prev. and next trial stuff interferes
        startpoint = T['stimstart_idx'] - neph.time2index(padtime,T['dt'])
        endpoint = T['stimend_idx'] + neph.time2index(padtime,T['dt'])
        
        trial_ss = responses_licksynced[trx][startpoint:endpoint+1,0]
        responses_licksynced[trx] = responses_licksynced[trx][startpoint+DTix:endpoint+DTix+1,:]
        trial_ls = responses_licksynced[trx][:,0]
        
        # if trx==18:
        #     print('after',DT,'idx',startpoint+DTix,endpoint+DTix+1)
        #     plt.plot(trial_ss)
        #     plt.plot(trial_ls)

        # if trx==6:
        #     fig,ax = plt.subplots(1,2,figsize=(24,8))
        #     ax[0].plot(trial.times,trial[:,0],lw=2,color='navy')
        #     ax[0].set_xlim(-1500,4500)
        #     ax[1].plot(trial2.times,trial2[:,0],lw=2,color='darkred')
        #     ax[1].set_xlim(-1500,4500)
        #     fig.suptitle('%d ms %d ms   %d'%(licktime,DT,DTix))
        #     print(responses_licksynced[trx].shape, type(responses_licksynced[trx]))
            
        # print(trx,responses_licksynced[trx].shape,DTix)

    # print([ len(trial) for trial in responses_licksynced ])
    # print(np.where(np.array([ len(trial) for trial in responses_licksynced ])==0))

    return responses_licksynced









def convert_trialstocsv(dn,ES=None):
    # dn = 'DT008' for audio first; dn = 'DT014' for visual first, 
    # dn = 'GT102' for the double block experiment
    print('Extracting trial events from %s.'%dn)
    ES,_,_ = mat2panda(dn)
    keylist = ['vstim', 'astim', 'block', 'water', 'punish']
    
    kernelwidth = 1000
    triallength = 3500         # this will be used to exclude non-trial-start events; be safe with +500ms
    
    # detect visual event starts; simple: values are degrees 
    cvstim = ES['vstim']>0

    pddvstim_aux = np.where(cvstim==1)[0]
    
    pddvstim = []
    v_eventtype = []
    for tx,t in enumerate(pddvstim_aux):
        if  tx==0 or t-pddvstim[-1]>triallength:
            pddvstim.append(t)
            v_eventtype.append(  ES['vstim'].iloc[t]  )
    
    print('# visual events',len(pddvstim),len(v_eventtype))

    events_visual = pd.DataFrame({'start': pddvstim, 'degree': v_eventtype})

      
    
    
    
    
    # detect audio event starts; need 2nd derivative here
    castim = np.convolve(ES['astim'],np.r_[np.zeros(kernelwidth),np.ones(kernelwidth)]/kernelwidth,'same')
    dcastim = np.diff( castim )
    cdcastim = dcastim>1e-3
    ddastim = np.diff(cdcastim.astype('int16'))
    
    pddastim_aux = np.where(ddastim==1)[0]+1
    
    pddastim = []
    a_eventtype = []
    for tx,t in enumerate(pddastim_aux):
        if  tx==0 or t-pddastim[-1]>triallength:
            pddastim.append(t)
            if ES['astim'].iloc[t:t+triallength].mean()>4: a_eventtype.append(5000)
            else: a_eventtype.append(10000)
    
    print('# audio events',len(pddastim),len(a_eventtype))

    events_audio = pd.DataFrame({'start': pddastim, 'freq': a_eventtype})
    
    
    
    
    
    # merge the two stimuli
    vis = events_visual['start'].iloc[0]
    ais = events_audio['start'].iloc[0]
    if vis<ais:        # if the first block is attend visual, in the cue block just the first modality
        print('visual - audio')
        events_l = events_visual
        events_r = events_audio
        vi = ais-1        # set the other as join point
        lb = True        # extra steps to include the 5th block for receptive field characterization
    else:
        print('audio - visual')
        events_l = events_audio
        events_r = events_visual
        vi = vis-1
        lb = False
    
    print('visual',vis,'audio',ais,'boundary',vi)

    # merge for multimodal events at 1000 seconds nearest:
    events = pd.merge_asof(events_r,events_l,on='start',tolerance=1000, direction='nearest')
    # add the first block:
    events = pd.concat(  (events_l[events_l['start']<vi], events), join='outer')

    if lb:      # add the fifth block in audio first sessions:
        vil = events['start'].iloc[-1]
        events = pd.concat(  (events, events_l[events_l['start']>vil]), join='outer')
    
    
    # normalize the column order
    events['duration'] = 3000*np.ones(len(events),dtype=np.int16)
    events = events[['start','duration','degree','freq']]






    # add block
    events['block'] = np.zeros(len(events),dtype=np.int16)
    # add water
    events['water'] = np.zeros(len(events),dtype=bool)

    blockids = []
    block = 0
    boundaries = []
    tx = -1
    for ix,e in events.iterrows():
        tx += 1
        # print(tx,ix,block,e['degree'],e['freq'])
        if block==0:                    # first singlemodal
            if not ( np.isnan(e['degree']) or np.isnan(e['freq']) ):
                block += 1
                boundaries.append(tx)
            blockids.append(block)
            continue
        if block==1:                     # first multimodal
            if ( np.isnan(e['degree']) or np.isnan(e['freq']) ):
                block += 1
                boundaries.append(tx)
            blockids.append(block)
            continue
        if block==2:                    # second singlemodal
            if not ( np.isnan(e['degree']) or np.isnan(e['freq']) ):
                block += 1
                boundaries.append(tx)
            blockids.append(block)
            continue
        if block==3:                    # second multimodal
            if np.isnan(e['freq']):
                block += 1
                boundaries.append(tx)
            blockids.append(block)
            continue
        blockids.append(block)          # characterization block

    events['block'] = blockids

    print('boundaries ', boundaries)
    
    # add water
    va = (vis>ais).astype(np.int16) # order of context, 0 -> visual first
    go = [ events['degree']==45, events['freq']==5000 ]
    events['water'].iloc[:boundaries[1]] = go[va].iloc[:boundaries[1]]
    events['water'].iloc[boundaries[1]:boundaries[3]] = go[1-va].iloc[boundaries[1]:boundaries[3]]
    
    # add action and reward
    waterdeliverytime = 2000
    licktrajectories = ES['licking'] > 3     # significant licking
    action = []
    punish = []
    print('here1')
    for tx,e in events.iterrows():
        lick = np.max(licktrajectories[e['start']+waterdeliverytime:e['start']+waterdeliverytime+1000+1500]) # check if mouse licked after water delivery before pre 1500 to next trial
        action.append(lick)
        punish.append(   lick ^ e['water']    )

    events['action'] = action
    events['punish'] = np.array(punish, dtype=bool)


    # show and save
    
    print(events)

    # save
    print(events.to_csv(pathdatamouse+'trials-generated/trials-start-good-C%sdata.csv'%dn,index=False))
    
    # debug save
    # events.to_csv(pathdatamouse+'trials-generated/trials-start-good-C%sdata_safe.csv'%dn,index=False)

    
    
    
    
    # diagnostics
    
    
    # events_astim = np.zeros(len(ES['astim']))
    # events_astim[pddastim] = 1
    
    # events_vstim = np.zeros(len(ES['vstim']))
    # events_vstim[pddvstim] = 1
    

    # t0 = 335000; t1 = 367000                  # choose a balanced 5 trials end of audio S start of Av M blocks
    # t0 += 20000; t1 = t0 + 5000                # modulated
    # t0 += 1500; t1 = t0 + 300         # event start magnified
    # t0 += 1950; t1 = t0 + 100         # event end magnified

    # t0 = 356648; t1 = t0 + 3012                # modulated

    # t0 += 14000; t1 = t0 + 5000        # unmodulated
    # t0 += 500; t1 = t0 + 100         # event start magnified
    # t0 += 1950; t1 = t0 + 100         # event end magnified
    
    


    # diagnose visual
    # plt.plot(np.arange(t0,t1),ES['vstim'].iloc[t0:t1]==45,'o',color='darkgreen')
    # plt.plot(np.arange(t0,t1),ES['vstim'].iloc[t0:t1]==135,'o',color='darkred')

    # plt.plot(np.arange(t0,t1),cvstim[t0:t1],'-',color='skyblue')
    # plt.plot(np.arange(t0,t1),dvstim[t0:t1],'-',color='cyan')
    # plt.plot(np.arange(t0,t1),events_vstim[t0:t1],'o',color='navy')

    
    # diagnose audio
    # a = (ES['astim']>0.04).astype(np.int16)
    # plt.plot(np.arange(t0,t1),a[t0:t1],color='grey')
    
    # plt.plot(np.arange(t0,t1),ddastim[t0:t1],'-',color='mediumvioletred')
    # plt.plot(np.arange(t0,t1),events_astim[t0:t1],'o',color='forestgreen')


    
 















def loadtrainingbehaviouraldata(datafileid, useonlymultimodal=False, recalculate=False, plot=False):
    
    
    if recalculate:
        
        # loaded = sp.io.loadmat(pathdatamouse + 'trainingbehaviour/' + datafileid + '.mat')
        loaded = h5py.File(pathdatamouse + 'trainingbehaviour/' + datafileid + '.mat','r')
        data = loaded['BehavStruct'][0][0]      # data contains all the sessions, visual, audio, then mixed combined...
        sessionids = data.dtype.names           # this holds the session id names,   v01,a01..., and mixed av01, va01, av02... etc.
        L = len(sessionids)
        
        
        # indices and labels
        multimodalsessionlabels = [ s for s in sessionids if s[:2] in ['av','va'] ]
        multimodalsessionvisualfirst = [ int(s[:2]=='va')     for s in sessionids if s[:2] in ['av','va'] ]
        if datafileid=='DT008': # erroneous recording sessions:
            del(multimodalsessionlabels[4]); del(multimodalsessionvisualfirst[4])
            del(multimodalsessionlabels[5]); del(multimodalsessionvisualfirst[5])         # delete the pair so that unimodal stays even
        multimodalsessionvisualfirst = np.array(multimodalsessionvisualfirst)
    
        # session data
        multimodalsessions = [ data[sid] for sid in multimodalsessionlabels ]
        
    
        # load the recording session behav as weel:
        bp = preprocess.assessbehaviouralperformancefullarray(datafileid)
    
    
        
    #    print( datafileid, multimodalsessionlabels, multimodalsessionvisualfirst )
        K = len(multimodalsessions)
        
        performances = [] 
        hits = []
        misses = []
        corrrejs = []
        fals = []
        performance = []
        
        cueperiods_initial = []      # collect cue periods by initial,transition
        cueperiods_transition = []
        
        for sx in range(len(multimodalsessions)): # go through each multimodal session
            L = len(multimodalsessions[sx][2][0])
    #        print('session info:   ', sx,L,multimodalsessionlabels[sx])
        
            if multimodalsessionvisualfirst[sx]:     # get visual
                go = (multimodalsessions[sx][0][0].ravel() == 45).astype('int')
            else:   #  or get audio
                go = (multimodalsessions[sx][1][0].ravel() == 5).astype('int')
            
            lick = np.zeros( L, dtype='int16' )
            for l in range(L):
                aux = multimodalsessions[sx][2][0][l][0].ravel()
                if len(aux)>0:
                    lick[l] = aux[-1]>=2.
            
    
            if len(lick)==len(go):      # if there are some issues with the data, then skip
    
                # now find the first multimodal session, if it exists:
                cueonly = False
                if multimodalsessionvisualfirst[sx]:
    #                print((multimodalsessions[sx][1][0].ravel()>0).astype('int16'))
                    mmpoint = np.where(multimodalsessions[sx][1][0].ravel()>0)[0]
                    if len(mmpoint)>0: mmpoint = mmpoint[0]
                    else: cueonly = True
                else:
    #                print((multimodalsessions[sx][0][0].ravel()>0).astype('int16'))
                    mmpoint = np.where(multimodalsessions[sx][0][0].ravel()>0)[0]
                    if len(mmpoint)>0: mmpoint = mmpoint[0]
                    else: cueonly = True
                
                if cueonly:
                    go_cue = go
                    lick_cue = lick
                else:          # only extend the multimodal data, if there is actual multimodal sessions
                    go_cue = go[:mmpoint]
                    lick_cue = lick[:mmpoint]
                    if useonlymultimodal:
                        go = go[mmpoint:]
                        lick = lick[mmpoint:]
                        
                
                    # for multimodal if exists
                    # hit miss correct rej false alarm
                    n_go = np.sum(go)
                    n_nogo = np.sum(1-go)
                    performance = np.array([ np.sum(go & lick) / n_go, np.sum( go & (1-lick) ) / n_go,\
                                             np.sum( (1-go) & (1-lick) ) / n_nogo, np.sum( (1-go) & lick ) / n_nogo ])
            
                    hits.extend( go & lick )
                    misses.extend( go & (1-lick) )
                    corrrejs.extend( (1-go) & (1-lick) )
                    fals.extend( (1-go) & lick )
                    
                # for cue period (always exists)
                # these single modality cues are to be aligned to trials 1-30
                n_go_cue = np.sum(go_cue)
                n_nogo_cue = np.sum(1-go_cue)
    #            performance_cue = np.array([ np.sum(go_cue & lick_cue) / n_go_cue, np.sum( go_cue & (1-lick_cue) ) / n_go_cue,\
    #                                     np.sum( (1-go_cue) & (1-lick_cue) ) / n_nogo_cue, np.sum( (1-go_cue) & lick_cue ) / n_nogo_cue ])
        
    #            hits.extend( go_cue & lick_cue )
    #            misses.extend( go_cue & (1-lick_cue) )
    #            corrrejs.extend( (1-go_cue) & (1-lick_cue) )
    #            fals.extend( (1-go_cue) & lick_cue )
                
                # here we insert all 30 cue trials, and will estimate p of the binomial distribution over sessions (sessionaveraging).
                n_trials_cue = 30
                if len(go_cue)>=n_trials_cue:
                    go_cue = go_cue[:n_trials_cue]
                    lick_cue = lick_cue[:n_trials_cue]
                    if sx%2==0:  # we should separately record initial and transition cue periods
                        cueperiods_initial.append(  [ (go_cue & lick_cue) | ((1-go_cue) & (1-lick_cue)),\
                                                      (go_cue & (1-lick_cue)) |  ((1-go_cue) & lick_cue)  ]   )
                    else:
                        cueperiods_transition.append(  [ (go_cue & lick_cue) | ((1-go_cue) & (1-lick_cue)),\
                                                         (go_cue & (1-lick_cue)) |  ((1-go_cue) & lick_cue)  ]   )
                    # cueperiods will be (sessions,h+c/m+f,trials)=(8,2,30)
                
    #            print(len(go_cue), len(go))
        
            else: print('NaN'); performance = [-1,-1,-1,-1]
    
            performances.append(performance)
        
    
    
        # add the recording session unimodal blocks
        cueperiods_initial.append(    [      [ bp[0][0][b] in [0,2] for b in range(n_trials_cue) ] ,       [ bp[0][0][b] in [1,3]  for b in range(n_trials_cue)     ]  ] )
        cueperiods_transition.append(    [      [bp[2][0][b] in [0,2] for b in range(n_trials_cue) ] ,       [ bp[2][0][b] in [1,3]  for b in range(n_trials_cue)     ]  ] )
            
        
        
        # collect all performance data from all performance:
        alldata = [hits, misses, corrrejs, fals]
        
        # now load the neural recording session into the end of performances array:
    #    print('final recording')
        recordedlist = assessbehaviouralperformance(datafileid, modality='all', multimodalonly=True)
        r = [ len(rr) for rr in recordedlist ]      # note that this loads both unimodal cue and multimodal task blocks
        
        performance = [    r[0]/(r[0]+r[1]), r[1]/(r[0]+r[1]), r[2]/(r[2]+r[3]), r[3]/(r[2]+r[3])     ]
        performances.append(performance)
        
        multimodalsessionlabels.append('record')
        performances = np.array(performances)
        
       
        a_go = np.sum(hits)+np.sum(misses)
        a_nogo = np.sum(corrrejs)+np.sum(fals)
        h = np.sum(hits)
        m = np.sum(misses)
        c = np.sum(corrrejs)
        f = np.sum(fals)
    
        overallperformance = [  h/a_go, m/a_go, c/a_nogo, f/a_nogo, (h+c)/(a_go+a_nogo), (m+f)/(a_go+a_nogo)  ]
        
        o = overallperformance
        
        overallperformance_sem = [  np.sqrt(o[0]*o[1]/a_go), np.sqrt(o[0]*o[1]/a_go), np.sqrt(o[2]*o[3]/a_nogo), np.sqrt(o[2]*o[3]/a_nogo),\
                                    np.sqrt(o[4]*o[5]/(a_go+a_nogo)), np.sqrt(o[4]*o[5]/(a_go+a_nogo))  ]
        
        
        
        # calculate p value and cr.i.
        cue_est = []
    
        cueperiods_initial = np.array(cueperiods_initial)
        cueperiods_transition = np.array(cueperiods_transition)
        for cueperiods in [cueperiods_initial, cueperiods_transition]:
    #        hc_m = cueperiods[:,0,:].mean(axis=0)
    #        mf_m = cueperiods[:,1,:].mean(axis=0)
            # set prior:
            a = 3; b = 1.2
            a = 1; b = 1
            # get maximum likelihood and the 67% and 95% highest density credible interval
            hc_m, hc_cl1, hc_cl2, hc_cu1, hc_cu2, hc_a, hc_b = neba.estimate_binomialp(cueperiods[:,0,:], prior_a=a, prior_b=b, resolution=100)
            mf_m, mf_cl1, mf_cl2, mf_cu1, mf_cu2, mf_a, mf_b = neba.estimate_binomialp(cueperiods[:,1,:], prior_a=b, prior_b=a, resolution=100)
            cue_est.append( [ hc_m, mf_m, hc_cl1, hc_cl2, hc_cu1, hc_cu2, hc_a, hc_b, mf_cl1, mf_cu1  ] )
        cue_est = np.array(cue_est)
        
        cueperiods = [ cueperiods_initial, cueperiods_transition ]
        # cueperiods will be (init/trans,sessions,h+c/m+f,trials)=(14,8,2,30)
        
        
        # now save everything what was calculated
        pickle.dump( (alldata, overallperformance, overallperformance_sem, performances, cueperiods, cue_est), \
                    open(cacheprefix+'phys/behaviour,trainingsessions-%s.pck'%(datafileid),'wb') )

    else:     # just reload
        alldata, overallperformance, overallperformance_sem, performances, cueperiods, cue_est = \
            pickle.load( open(cacheprefix+'phys/behaviour,trainingsessions-%s.pck'%(datafileid), 'rb')  )

        

    if plot:
        fig,ax = plt.subplots(1,2,figsize=(8,16))
        axs = ax[0]
        axs.plot(performances,lw=2)
        axs.set_xticks(np.arange(K+1))
        axs.set_xticklabels(multimodalsessionlabels,rotation=90)
        axs.set_ylim(-0.1,1.1)
        axs.legend(['hit','miss','corrrej','fal'])
        axs.set_title(datafileid + ', cueing unimodal training performance')

        axs = ax[1]
        axs.plot(performances,lw=2)
        axs.set_xticks(np.arange(K+1))
        axs.set_xticklabels(multimodalsessionlabels,rotation=90)
        axs.set_ylim(-0.1,1.1)
        axs.legend(['hit','miss','corrrej','fal'])
        axs.set_title(datafileid + ', multimodal training performance')


    return alldata, overallperformance, overallperformance_sem, performances, cueperiods, cue_est










def getorderattended(datafileid):

    # test block, depending on experiment and attendant stimulus
    if datafileid=='ME103' or datafileid=='ME108' or datafileid=='ME113' or \
        datafileid=='DT014' or datafileid=='DT017' or datafileid=='DT018' or datafileid=='DT019' or datafileid=='DT020' or \
        datafileid=='DT030' or datafileid=='DT031' or \
        datafileid=='LP001' or datafileid=='LP003' or datafileid=='LP005' or\
        datafileid=='AC001' or datafileid=='AC004' or datafileid=='AC006' or \
        datafileid=='AC007' or datafileid=='AC008' or datafileid=='AC009':
        bltv = [1,2]        # for these experiments, order of attended stimuli is: first attend visual, then attend audio
        blta = [3,4]
    elif datafileid=='ME110' or datafileid=='ME112' or \
        datafileid=='DT008' or datafileid=='DT009' or datafileid=='DT021' or datafileid=='DT022' or datafileid=='DT032' or \
        datafileid=='AC003':
        bltv = [3,4]        # these experiments: first attend audio, then attend visual
        blta = [1,2]
    else:                  # give an error
        print('warning: data file %s, order is not specified'%datafileid)
        bltv = None
        blta = None
    return bltv,blta




def getcontextidents(dn,block):
    ixblv = []
    ixbla = []

    blv,bla = getorderattended(dn)
    
    for trx,trial in enumerate(block.segments):
        if trial.annotations['block'] in blv: ixblv.append(trx)
        elif trial.annotations['block'] in bla: ixbla.append(trx)

    return ixblv,ixbla






def getblockidents(dn,block):
    ixs1 = []; ixm2 = []
    ixs3 = []; ixm4 = []
    for trx,trial in enumerate(block.segments):
        if trial.annotations['block'] == 1: ixs1.append(trx)
        elif trial.annotations['block'] == 2: ixm2.append(trx)
        elif trial.annotations['block'] == 3: ixs3.append(trx)
        elif trial.annotations['block'] == 4: ixm4.append(trx)
    return ixs1,ixm2,ixs3,ixm4





def getstimulusidents(dn,block,multimodalonly=False):
    ixv45 = []; ixv135 = []; ixa45 = []; ixa135 = []
    ixa5000 = []; ixa10000 = []; ixv5000 = []; ixv10000 = []

    blv,bla = getorderattended(dn)
    
    if multimodalonly:
        blv = [blv[1]]
        bla = [bla[1]]
    
    for trx,trial in enumerate(block.segments):
        if trial.annotations['block'] in blv and trial.annotations['visual']==45: ixv45.append(trx)
        elif trial.annotations['block'] in blv and trial.annotations['visual']==135: ixv135.append(trx)
        elif trial.annotations['block'] in bla and trial.annotations['visual']==45: ixa45.append(trx)
        elif trial.annotations['block'] in bla and trial.annotations['visual']==135: ixa135.append(trx)
        
        if trial.annotations['block'] in bla and trial.annotations['audio']==5000: ixa5000.append(trx)
        elif trial.annotations['block'] in bla and trial.annotations['audio']==10000: ixa10000.append(trx)
        elif trial.annotations['block'] in blv and trial.annotations['audio']==5000: ixv5000.append(trx)
        elif trial.annotations['block'] in blv and trial.annotations['audio']==10000: ixv10000.append(trx)
    
    
    return ixv45,ixv135,ixa45,ixa135,ixv5000,ixv10000,ixa5000,ixa10000




    
def getcelltypes(block):
    n_neurons = block.segments[0].analogsignals[0].shape[1]
    ixnarrow = []; ixbroad = []
    ixsuperf = []; ixinput = []; ixdeep5 = []; ixdeep6 = []
    for n in range(n_neurons):
        if block.annotations['waveforms'][n]==0: ixnarrow.append(n)
        else: ixbroad.append(n)
        if block.annotations['idents'][n]==1: ixsuperf.append(n)
        elif block.annotations['idents'][n]==2: ixinput.append(n)
        elif block.annotations['idents'][n]==3: ixdeep5.append(n)
        elif block.annotations['idents'][n]==4: ixdeep6.append(n)
    return ixnarrow,ixbroad,ixsuperf,ixinput,ixdeep5,ixdeep6


def getcelltypesinlist(block,withnewline=True):
    ixct = getcelltypes(block)
    ctgroups = []
    for i in range(4):
        ctgroups.append(np.intersect1d(ixct[1],ixct[i+2]).astype('int16'))  # add broad and then narrow for each group
        ctgroups.append(np.intersect1d(ixct[0],ixct[i+2]).astype('int16'))
    ctcolors = ['deepskyblue','lightcoral','darkviolet','deeppink','blue','red','navy','darkred']    
    if withnewline:
        ctlabels = ['superf\nbroad','superf\nnarrow','input\nbroad','input\nnarrow','deep5\nbroad','deep5\nnarrow','deep6\nbroad','deep6\nnarrow']
    else:
        ctlabels = ['superf broad','superf narrow','input broad','input narrow','deep5 broad','deep5 narrow','deep6 broad','deep6 narrow']

    return ctgroups,ctcolors,ctlabels




# physiology
def get_preferredorientation(dn,block,T,absolute=False):
    
    stimIDs = [ [ [1,2,3,4],[45],[] ],  [ [1,2,3,4],[135],[] ] ]    # 45 and 135 from both single and multi trials
    responses = collect_stimulusspecificresponses(block,stimIDs)
    
    orientationresponse = np.zeros((responses[0][0].shape[1],2))
    for oix in range(2):   # 45 and 135
        if absolute:
            orientationresponse[:,oix]= np.abs(\
               np.array(responses[oix])[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1) - \
               np.array(responses[oix])[:,T['start_idx']:T['stimstart_idx'],:].mean(axis=1) ).mean(axis=0)
        else:
            orientationresponse[:,oix]= np.array(responses[oix])[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1).mean(axis=0)
    
    selectivity = (orientationresponse[:,0] - orientationresponse[:,1])
    
    return selectivity















# automatically remove drift and return unremovably drifting cell mask
def removedrift(dn,T,filtercutoff=20*pq.min,sampling_trials=1,analogchannelid=0,continuous_method='instfr'):
    # the purpose of this is to preprocess and remove drifts or remove the to-zero-Hz or from-zero-Hz cells altogether

    doplot = False
    
    # pre and on integration width
    width = 1500*pq.ms        # integration window to check firing rates to reject cell wholly
    windowix = int((width/T['dt']).magnitude)

    # spike count kernel 1 std width
    width_s = 100*pq.ms

    # load the raw firing data


    comparisongroup = [   [ [1,2],   [],    [] ], [ [3,4],  [],     [] ] ]

    # collect firing rates:
    block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=False,recalculate=False)
    across_responses = preprocess.collect_stimulusspecificresponses(block, comparisongroup)
    responses = np.concatenate(across_responses)
    block_orig = preprocess.loaddatamouse(dn,T,continuous_method,normalize=False,recalculate=False)

    fr_pre = responses[:,T['stimstart_idx']-windowix:T['stimstart_idx'],:].mean(axis=1)

    n_trials, n_neurons = fr_pre.shape
    mask = np.zeros(n_neurons,dtype='int16')
    mask_lsq = np.zeros((n_neurons,5))      # return values of the linear least squares slope estimate


    fr_pre_filt = np.zeros((n_trials,n_neurons))
    
    cut_trials = 20  #   high pass filter how many trials above which slow waves to remove
    sampling_trials = 1 #       one trial sampling
    p_crit = 1e-5               # criteria for exclusion
    p_nonfilt = 0.05            # criteria for not filtering at all (i.e. no slow component.)
    
#    average_trial_length = (3+3*0.8+10*0.2)*1000*pq.ms        # estimate by average performance
#    average_trial_length = 6000*pq.ms                          # estimate by actual cut out trial lengths
#    max_t = n_trials*average_trial_length
    
    # now create the filter sampling rate in the time domain (f will= 1/T)
#    filtercutoff = cut_trials * average_trial_length/1000/10 
#    filtercutoff = (1/(49*pq.Hz)).rescale('ms')
    filtercutoff = filtercutoff.rescale('s')
    sampling_period = (sampling_trials * T['dt']).rescale('s')
    
    print(cut_trials,sampling_trials)
#    print(average_trial_length, max_t)
    print(filtercutoff.rescale('ms'),filtercutoff.rescale('min'),'sampling:',sampling_period)
    print((1/filtercutoff).rescale('Hz'),(1/sampling_period).rescale('Hz'))
    
#    filtercutoff = 21*pq.ms
    
    
    
    
    for n in range(n_neurons):

        # load the channel
        achannel = block.channel_indexes[analogchannelid].analogsignals[analogchannelid][:,n].squeeze()

        # filter with window width:
        kernel = neph.create_kernel(width_s,sampling_period.rescale('ms'))

        achannel = np.convolve(achannel, kernel, 'same')

        # remove frequencies
        achannel_filtered = neph.removefrequencies(achannel,[1/filtercutoff],1/sampling_period,filtertype='hp') * pq.Hz
        
        achannel_filtered = achannel_filtered[:,np.newaxis]#/T['dt']*width_s/2

        # plot whole channel        
        if doplot and n in [0,1,3]:
            fig,ax = plt.subplots(2,2,figsize=(36,16))
            axs = ax[0,0]
            axs.plot(achannel,color='black',alpha=0.7)
            axs.plot(achannel_filtered,color='red',alpha=0.7)
            axs = ax[1,0]
            ss = 80000; se = 100000
            axs.plot(np.arange(ss,se),achannel[ss:se],color='black',alpha=0.7)
            axs.plot(np.arange(ss,se),achannel_filtered[ss:se],color='red',alpha=0.7)
        
            u_fx,u_px = sp.signal.welch(x=achannel.squeeze(),fs=1/sampling_period)
            u_fx_filt,u_px_filt = sp.signal.welch(x=achannel_filtered.squeeze(),fs=1/sampling_period)
            
            axs = ax[0,1]
            axs.plot(u_fx,u_px,color='black',alpha=0.7)
            axs.plot(u_fx_filt,u_px_filt,color='red',alpha=0.7)
            
            axs = ax[1,1]
            fig.delaxes(axs)
            
#        u_fx_pre_filt, achannel = sp.signal.welch(x=fr_pre_filt[:,n],fs=1/sampling_trials)
        
        
        # put back the filtered data
        block.channel_indexes[analogchannelid].analogsignals[analogchannelid][:,n] = achannel_filtered
        tl_idx = block.segments[analogchannelid].analogsignals[analogchannelid].shape[0] # trial length
        
        # check if the cell is valid
        fr_pre_filt[:,n] = neph.removefrequencies(fr_pre[:,n],[1/cut_trials],1/sampling_trials,filtertype='hp')
        l = sp.stats.linregress(np.arange(n_trials)+1,np.abs(fr_pre_filt[:,n]))
        mask[n] = l[3]>p_crit   # p value for regression of standard deviations around the moving mean to detect off-on cells
        mask_lsq[n,:] = l
        
        # don't use filter if not needed, replace with original
        if l[3]>p_nonfilt:
            achannel_filtered = achannel[:,np.newaxis]* pq.Hz

        # put back the filtered data into the segments
        for tx,tr in enumerate(block.segments):
            sigs = achannel_filtered[tl_idx*tx:tl_idx*(tx+1)]
            block.segments[tx].analogsignals[analogchannelid][:,n] = sigs

    
    # plot some trials and neurons
    if doplot:
        times = block.segments[0].analogsignals[0].times
        fig,ax = plt.subplots(1,3,figsize=(36,8))
        for ix,n in enumerate([0,1,3]): # which trials
            axs = ax[ix]
#            axs.plot(times,block_orig.segments[55].analogsignals[0][:,n],lw=2,color='black',alpha=0.7)
#            axs.plot(times,block.segments[55].analogsignals[0][:,n],lw=2,color='red',alpha=0.7)
            axs.plot(times,block_orig.segments[55].analogsignals[0][:,n]/np.max(block_orig.segments[55].analogsignals[0][:,n]),lw=2,color='black',alpha=0.7)
            axs.plot(times,block.segments[55].analogsignals[0][:,n]/np.max(block.segments[55].analogsignals[0][:,n]),lw=2,color='red',alpha=0.7)
            axs.set_xticks([0,3000])

    
    
    return block,mask,mask_lsq





















# stimulus id groupings

def collect_stimulusspecificresponses(block,stimulusIDgroups,logical='and',s=0,correctonly=0, erroronly=0):
    # stimulusIDgroups is a list of stimulus identiers grouped as lists in each element of the top list
    # for example to use all stimuli: [[j+1 for j in range(max_st)]]
    # stimulusIDgroups is a triple tuple/list with (blockid,visualid,audioid) for example one elemnet from the list of stimulusIDgroups:
    #  [[1,2],[45,135],[]]   will use both visual, regardless of audio is present and which in block 1 and 2
    # these are for logical 'and'
    # for logical 'or': we have 6 sts elements, first 3 the and second three is combined with a large OR
    # s is signal: 0 neural, 1 runspeed
    # correctonly: 0 (False)  use all trials,  1 (True)  use only not punished trials, only for AND operator
    # erroronly: 0 (False)  use all trials,  1 (True)  use only punished trials, only for AND operator

    
#    n_channels = len(block.segments[0].analogsignals[s])
    n_trials = len(block.segments)

    responses = []
    
    if logical=='and':
        for sts in stimulusIDgroups:
            # print(sts,sts[0],sts[1],sts[2])
            responses.append( [ \
                  trial.analogsignals[s] for trial in block.segments\
                            if ( trial.annotations['block'] in sts[0] and \
                                (trial.annotations['visual'] in sts[1] or sts[1]==[])  and \
                                (trial.annotations['audio'] in sts[2] or sts[2]==[]) and \
                                (1-trial.annotations['punish']>=correctonly )   and \
                                (trial.annotations['punish']>=erroronly )  )    ]  )
    
    elif logical=='or':
        for sts in stimulusIDgroups:
#            print(sts[0],sts[1],sts[2],sts[3],sts[4],sts[5])
            responses.append( [ \
                  trial.analogsignals[s] for trial in block.segments\
                            if ( ( trial.annotations['block'] in sts[0] and \
                                  (trial.annotations['visual'] in sts[1] or sts[1]==[])  and \
                                  (trial.annotations['audio'] in sts[2] or sts[2]==[]) )    \
                                        or \
                                 ( trial.annotations['block'] in sts[3] and \
                                  (trial.annotations['visual'] in sts[4] or sts[4]==[])  and \
                                  (trial.annotations['audio'] in sts[5] or sts[5]==[]) )   
                                ) \
                              ] )
    
    
    else: print('logical operator not applicable for stimulus collection')

    return responses







def collect_stimulusspecificresponses_choice(block,datafileid,onlyinblock=[2,4],s=0,retval=False):
    # choice is determined by lick vs. withhold decision and action
    # lick: hit and false alarm
    # withhold: correct rejection and miss
    # onlyinblock default: only for multimodal blocks
    # s is signal: 0 neural, 1 runspeed
    
    responses = []
    
    ixhit,ixmiss,ixcorrrej,ixfal = assessbehaviouralperformance(datafileid)

    ixlick = np.hstack((ixhit,ixfal))
    ixnolick = np.hstack((ixcorrrej, ixmiss))
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixlick and trial.annotations['block'] in onlyinblock ] )
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixnolick and trial.annotations['block'] in onlyinblock ] )

    if retval:
        return responses, ixlick, ixnolick
    else:
        return responses
    








def collect_stimulusspecificresponses_choice_detailed(block,datafileid,onlyinblock=[2,4],s=0,retval=False):
    # choice is determined by lick vs. withhold decision and action
    # lick: hit and false alarm
    # withhold: correct rejection and miss
    # onlyinblock default: only for multimodal blocks
    # s is signal: 0 neural, 1 runspeed
    
    responses = []
    
    ixhit,ixmiss,ixcorrrej,ixfal = assessbehaviouralperformance(datafileid)

    # ixlick = np.hstack((ixhit,ixfal))
    # ixnolick = np.hstack((ixcorrrej, ixmiss))
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixhit and trial.annotations['block'] in onlyinblock ] )
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixmiss and trial.annotations['block'] in onlyinblock ] )
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixcorrrej and trial.annotations['block'] in onlyinblock ] )
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixfal and trial.annotations['block'] in onlyinblock ] )

    if retval:
        return responses, ixhit,ixmiss,ixcorrrej,ixfal
    else:
        return responses










def collect_rewardspecificresponses(block,datafileid,onlyinblock=[2,4],s=0):
    # choice is determined by lick vs. withhold decision and action
    # lick: hit and false alarm
    # withhold: correct rejection and miss
    # onlyinblock default: only for multimodal blocks
    # s is signal: 0 neural, 1 runspeed
    
    responses = []


    ixv45,ixv135,ixa45,ixa135,ixv5000,ixv10000,ixa5000,ixa10000 = preprocess.getstimulusidents(datafileid,block,multimodalonly=True)
    ixhit,ixmiss,ixcorrrej,ixfal = assessbehaviouralperformance(datafileid,multimodalonly=True)


    
    # get positive water reward, no reward with correct rej, and timing punishment
    
    ixreward = refoldlayerindices(ixhit,np.r_[ixv45,ixa5000])
    ixnoreward = ixempty = refoldlayerindices(ixfal,np.r_[ixv135,ixa10000])
    ixempty = refoldlayerindices(ixcorrrej,np.r_[ixv135,ixa10000])
    ixpunish = refoldlayerindices(np.r_[ixmiss,ixfal],np.r_[ixv135,ixa10000])
    
   
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixreward and trial.annotations['block'] in onlyinblock ] )
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixnoreward and trial.annotations['block'] in onlyinblock ] )
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixempty and trial.annotations['block'] in onlyinblock ] )
    responses.append( [ trial.analogsignals[s] for k,trial in enumerate(block.segments)\
                       if k in ixpunish and trial.annotations['block'] in onlyinblock ] )


    return responses, ixreward, ixnoreward, ixempty, ixpunish
    









def assessbehaviouralperformance(dn,modality='all',multimodalonly=False):
    
    g = loadexperimentdata(dn)
    

    success = g['punish']==False
    
    # get GO signals as True, NOGO signals as False
    truevisual = g['degree']==45
    trueaudio = g['freq']==5000
    

    bltv,blta = getorderattended(dn)
    if multimodalonly:
        bltv = [bltv[1]]
        blta = [blta[1]]
    bltv = np.array(bltv)-1
    blta = np.array(blta)-1
    if modality=='audio': bltv = []
    if modality=='visual': blta = []

    hit                 = g[  (   success &   truevisual & g['block'].isin(bltv)  )    |    \
                              (   success &   trueaudio  & g['block'].isin(blta)  ) ].index
    miss                = g[  ( ~ success &   truevisual & g['block'].isin(bltv)  )    |    \
                              ( ~ success &   trueaudio  & g['block'].isin(blta)  ) ].index
    correctrejection    = g[  (   success & ~ truevisual & g['block'].isin(bltv)  )    |    \
                              (   success & ~ trueaudio  & g['block'].isin(blta)  ) ].index
    falsealarm          = g[  ( ~ success & ~ truevisual & g['block'].isin(bltv)  )    |    \
                              ( ~ success & ~ trueaudio  & g['block'].isin(blta)  ) ].index

    behaviour = (hit,miss,correctrejection,falsealarm)

    return behaviour







def get_conditioned_behaviour(dn,comparisongroups,taskaspects):

    g = loadexperimentdata(dn)

    g['block']+=1
    success = g['punish']==False

    perftask = []
    for cx,comparison in enumerate(taskaspects):
        perf = []
        for c in comparisongroups[cx]:
            # print(g['block'] in  c[0]))
            cond = (  (g['block'].isin(c[0])  )     |   (c[0]==[])   )    & \
                   (  (g['degree'].isin(c[1])  )     |   (c[1]==[])   )    & \
                   (  (g['freq'].isin(c[2])  )     |   (c[2]==[])   )

            perfindiv = success[cond].sum()/cond.sum()
            perf.append(perfindiv)
        perftask.append(perf)

    perftask = pd.DataFrame(np.array(perftask),index=taskaspects,columns=['go','nogo'])
    

    return perftask









def get_conditioned_behaviour_singletrial(dn,comparisongroups,taskaspects):

    # mal = 5

    g = loadexperimentdata(dn)

    g['block']+=1
    success = g['punish']==False
    success.rename('success', inplace=True)




    behma_task = []     # behavior variables collector
    for cx,comparison in enumerate(taskaspects):
        behma = []
        for chx,c in enumerate(comparisongroups[cx]):    # go nogo
            # print(g['block'] in  c[0]))
            cond = (  (g['block'].isin(c[0])  )     |   (c[0]==[])   )    & \
                   (  (g['degree'].isin(c[1])  )     |   (c[1]==[])   )    & \
                   (  (g['freq'].isin(c[2])  )     |   (c[2]==[])   )

            # behma_task_indiv = np.convolve(success[cond],np.ones(mal)/mal)
            # behma_task_indiv = cond.astype(np.int16)
            
            behma_task_indiv = success[cond]
            behma_task_indiv.rename(['go','nogo'][chx],inplace=True)

            # behma_task_indiv = pd.Series(behma_task_indiv).rolling(window=5,  center=True, min_periods=2).mean()


            behma.append(behma_task_indiv)
        
        behma = pd.concat(behma,axis=1)

        behma_task.append(behma)


    # behma_task = pd.DataFrame(np.array(behma_task),index=taskaspects,columns=['go','nogo'])
    # behma_task = np.array(behma_task)

    # returns          taskaspects,classes,trials array
    return behma_task
    # one hot encode trial type
    # put nan-s on invalid columns
    # make a rolling average over n trials of the columns





    # n_trials = len(g)
    # conds = {'degree':45,'freq':5000}
    # smoothedvars = []
    # for v in vars:
    #     smoothedvar = []
    #     for cx,c in enumerate(conds):
    #         mask = g[c['column']]==c['value']
    #         smoothingvar = np.ones((len(mask)))*np.nan         # length of all trials
    #         smoothingvar[mask] = meta[v][mask]                 # where mask is valid, change to valid variable value
    #         if len(movavlens)>1:
    #             movavlen = movavlens[cx]
    #         else:
    #             movavlen = movavlens[0]

    #         if movavlen>0:
    #             smoothedvar.append(  pd.Series(smoothingvar).rolling(window=movavlen,  center=True, min_periods=min_period).mean()  )
    #         else:
    #             smoothedvar.append(  np.array(smoothingvar)   )

    #     smoothedvars.append(smoothedvar)
    
    # smoothedvars = np.array(smoothedvars)       # (variables, conditions, trials)

    # # smoothedvars = np.swapaxes( np.array(smoothedvars), 0,2)

    # return smoothedvars







    return behma_task










def getblocklengths(dn):
    g = loadexperimentdata(dn)
    
    lengths_block = np.zeros(4,dtype='int16')
    for bx in [0,1,2,3]:
        lengths_block[bx] = sum(g['block']==bx)
    
    return lengths_block
    




def assessbehaviouralperformancefullarray(dn):
    bp = list(assessbehaviouralperformance(dn))
    
    # trial index in blocks, start and endpoints
    tlb = np.cumsum(np.hstack((0,getblocklengths(dn))))  #  ; print('trial start index lists',tlb)


    L = -np.ones((tlb[-1]),dtype='int16')
    for cx,c in enumerate(bp):
        L[c] = cx
    
    a = []
    for bx in [0,1,2,3]:     # go through the blocks
        a.append([ L[   tlb[bx]:tlb[bx+1]  ]   ])

    return a





def refoldlayerindices(ix1,ix2):
    aux = np.intersect1d(ix1,ix2)
    return aux



        
if __name__ == '__main__':
    print('preprocess main')
