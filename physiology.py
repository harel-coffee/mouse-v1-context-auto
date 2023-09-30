#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 10:16:45 2020

@author: mahajnal
"""



from operator import index
from numpy import bool8

from quantities.quantity import _reconstruct_quantity
from scipy.fft import next_fast_len
from scipy.stats import bernoulli
from scipy.optimize import minimize, Bounds
from sqlalchemy import column
from config import *

import preprocess

import sklearn.model_selection as skms








def orientationselectivity(dn,block):
    selectivity = preprocess.get_preferredorientation(dn,block,T)
    # plt.bar(np.arange(len(selectivity)),selectivity)
    # plt.title(dn)
    pickle.dump(selectivity,open(cacheprefix+'phys/orientationselectivity_%s'%dn,'wb'))









def cluster_celltypes(datanames):
    print('processing waveform peaks and troughs...')
    times = []
    amplitudes = []
    
    for dn in datanames:
        block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=False)
        n_neurons = block.segments[0].analogsignals[0].shape[1]
        t,u = neph.identifycelltype(n_neurons,T,block.annotations['waveforms'],returnfullfeatures=True)
        times.extend(t)
        amplitudes.extend(u)
    times = np.array(times)
    amplitudes = np.array(amplitudes)
    mask = amplitudes<=1.
    times = times[mask]
    amplitudes = amplitudes[mask]
    
    
    X = np.c_[times,amplitudes]
    y,m,C,model = nedi.fitgaussianmixture( data=X, means_init=np.array([[20,0.2],[7,0.45]]) )

    pickle.dump(model,open(cacheprefix+'phys/broadnarrowspiking,gmm.pck','wb'))
    
    
    fig,ax = plt.subplots(2,1,figsize=(6,12))
    axs = ax[0]
    # axs.scatter(times[y==0],amplitudes[y==0],color='darkgreen')
    # axs.scatter(times[y==1],amplitudes[y==1],color='darkred')

    figs.plottoaxis_gmm(axs,X,y,m,C,colors=['darkgreen','darkred'])

    axs.set_ylim(0,1)
    axs.set_xlabel('trough to peak time [ms]')
    axs.set_ylabel('trough to peak amplitude ratio')

    axs = ax[1]
    axs.hist(times,bins=20)
    axs.set_xlabel('trough to peak time [ms]')

    save = 0 or globalsave
    if save:    
        fig.savefig(resultpath+'broadnarrow,gmm,histogram.png')
    
    
    
    












# *************************************************************************
#                                REWARD, LICKSYNC









def reward(dn,block):
    # compare normal and lick synced responses to reward related activity
    # get trials in a lick-synced manner; use NaN as lick times for no lick trials
    # then licksync no lick trials to water delivery time

    recalculate = 0 or globalrecalculate

    n_trials = len(block.segments)
    
    width = 50*pq.ms        # decoder folding
    padding = 500*pq.ms


    

    # get the action trial-grouping (choice):
    responses_action, ixlick, ixnolick = preprocess.collect_stimulusspecificresponses_choice(block, dn, retval=True)

    responses_action_ls = []
    responses_action_ls.append(preprocess.licksync_trials(dn,T,ixlick,responses_action[0],padding))
    responses_action_ls.append(preprocess.licksync_trials(dn,T,ixnolick,responses_action[1],padding))

    # print(np.array(responses_action_ls[0]).shape)
    # print([ len(responses_action_ls[0][i]) for i in range(len(responses_action_ls[0])) ])
    # print(np.where(np.array([ len(responses_action_ls[0][i]) for i in range(len(responses_action_ls[0])) ])==0))

    if recalculate:
        actiondecoder = nedi.get_responsedecoding(responses_action,width=width)
        pickle.dump(actiondecoder,open(cacheprefix+'continuous/actiondecodes,syncstim,angles-%s_%s.pck'%(dn,continuous_method),'wb'))
    else:
        actiondecoder = pickle.load(open(cacheprefix+'continuous/actiondecodes,syncstim,angles-%s_%s.pck'%(dn,continuous_method),'rb'))

    if recalculate:
        actiondecoder_ls = nedi.get_responsedecoding(responses_action_ls,width=width)
        pickle.dump(actiondecoder_ls,open(cacheprefix+'continuous/actiondecodes,synclick,angles-%s_%s.pck'%(dn,continuous_method),'wb'))
    else:
        actiondecoder_ls = pickle.load(open(cacheprefix+'continuous/actiondecodes,synclick,angles-%s_%s.pck'%(dn,continuous_method),'rb'))






    # get th reward trial-grouping


    # responses will have 3
    responses_cand,ixreward,ixnoreward,ixempty,ixpunish = preprocess.collect_rewardspecificresponses(block,dn)
    # find the two most interesting case: reward and no reward with the same action:
    responses = [responses_cand[0],responses_cand[1]]
    
    # get the same with the lick synced trials
    responses_ls = []
    responses_ls.append(preprocess.licksync_trials(dn,T,ixreward,responses[0],padding))
    responses_ls.append(preprocess.licksync_trials(dn,T,ixpunish,responses[1],padding))
    
    # print(np.concatenate(responses,axis=0).shape)
    # print(np.concatenate(responses_ls,axis=0).shape)
    # return
    
    if recalculate:
        rewarddecoder = nedi.get_responsedecoding(responses,width=width)
        pickle.dump(rewarddecoder,open(cacheprefix+'continuous/rewarddecodes,syncstim,angles-%s_%s.pck'%(dn,continuous_method),'wb'))
    else:
        rewarddecoder = pickle.load(open(cacheprefix+'continuous/rewarddecodes,syncstim,angles-%s_%s.pck'%(dn,continuous_method),'rb'))

    if recalculate:
        rewarddecoder_ls = nedi.get_responsedecoding(responses_ls,width=width)
        pickle.dump(rewarddecoder_ls,open(cacheprefix+'continuous/rewarddecodes,synclick,angles-%s_%s.pck'%(dn,continuous_method),'wb'))
    else:
        rewarddecoder_ls = pickle.load(open(cacheprefix+'continuous/rewarddecodes,synclick,angles-%s_%s.pck'%(dn,continuous_method),'rb'))
    
    
    
    
    
    # ; for now only multimodal blocks; reward expectation should show up in unimodals
    blv,bla = preprocess.getorderattended(dn)
    # determine trial set to get responses from:
    bix = preprocess.getblockidents(dn,block)
    trialset = np.r_[bix[blv[1]-1],bix[bla[1]-1]]
    # 
    o,m,s,e = preprocess.get_lickonsets(dn,T,trialset)

    
    
    
    
    # plot on figures
    fig,ax = plt.subplots(2,2,figsize=(24,16))
    pdf_l = 0.3; pdf_d = 0.1
    
    for dx,decoder in enumerate([actiondecoder,rewarddecoder]):
    
        axs = ax[dx,0]
        
        
        figs.plottoaxis_decoderrocauc(decoder,axs)
        figs.plottoaxis_chancelevel(axs,0.5)
        axs.set_xlim(-1500,4500)
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.add_patch(plt.Rectangle((T['starttime'],0),-T['starttime']+T['endtime'],0.45,fill=True,color='white'))
        
        # water delivery and lick times distribution
        axs.plot([0,0],[pdf_l,pdf_l+pdf_d],color='grey',lw=2)
        axs.plot([3000,3000],[pdf_l,pdf_l+pdf_d],color='grey',lw=2)
        axs.plot([2000,2000],[pdf_l,pdf_l+pdf_d],color='aqua',lw=2)
        figs.plottoaxis_normalpdf(axs,m,s,pdf_l,pdf_d,color='darkred',alpha=0.3)
        figs.plottoaxis_normalpdf(axs,m,e,pdf_l,pdf_d,color='darkred')
        axs.set_ylim(pdf_l,1)

        axs.set_ylabel(['action: lick - nolick','reward: water - false alarm'][dx])

    
        if dx==0: axs.set_title('stimulus sync')

    
    for dx,decoder_ls in enumerate([actiondecoder_ls,rewarddecoder_ls]):
        # we need synchronized times for comparison, at the mean values
        t = np.arange(0*pq.ms,decoder_ls[1].shape[0]*T['dt'],T['dt'])*pq.ms-padding + T['dt']
        
        axs = ax[dx,1]


        # figs.plottoaxis_decoderrocauc(decoder_ls,axs)
        for tt in [0,1]:
            d_m = decoder_ls[tt][:,0].squeeze()
            d_s = decoder_ls[tt][:,1].squeeze()
            d_e = decoder_ls[tt][:,2].squeeze()
            figs.plottoaxis_plottrajectory(axs,t,m=d_m,s=d_e,colorlist=['dodgerblue','darkorange'][tt],alpha=1,linewidth=2,label=None,fill=True)
            figs.plottoaxis_plottrajectory(axs,t,m=d_m,s=d_s,colorlist=['dodgerblue','darkorange'][tt],alpha=0.7,linewidth=2,label=None,fill=True)

        axs.set_xlim(-1500,4500)
        axs.set_ylim(pdf_l,1)
        figs.plottoaxis_chancelevel(axs,0.5)
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.add_patch(plt.Rectangle((T['starttime'],0),-T['starttime']+T['endtime'],0.45,fill=True,color='white'))
    
        # plot the lick sync times 
        axs.plot([m,m],[pdf_l,pdf_l+pdf_d],color='darkred',lw=2) # this is the lick sync
        figs.plottoaxis_normalpdf(axs,2000,s,pdf_l,pdf_d,color='aqua',alpha=0.3)
        figs.plottoaxis_normalpdf(axs,2000,e,pdf_l,pdf_d,color='aqua')
        figs.plottoaxis_normalpdf(axs,0,s,pdf_l,pdf_d,alpha=0.3)
        figs.plottoaxis_normalpdf(axs,0,e,pdf_l,pdf_d)
        figs.plottoaxis_normalpdf(axs,3000,s,pdf_l,pdf_d,alpha=0.3)
        figs.plottoaxis_normalpdf(axs,3000,e,pdf_l,pdf_d)
        
    
        if dx==0: axs.set_title('lick sync ~%d ms'%m)

            
            

    

    
    fig.suptitle(dn+' reward accuracy')

    save = 1
    if save:
        fig.savefig(resultpath+'action,reward_stim+licksync_%s'%dn+ext)






















# **************************************************************************************
#                                RUNSPEED





def joint_differentiabilityrunspeed(dn,block):

#    n_levels = 1000
    max_runspeed = 100.*pq.cm/pq.s
    delta_runspeed = 1.*pq.cm/pq.s
    runspeedlevels = np.linspace(0*pq.cm/pq.s, max_runspeed, int(max_runspeed/delta_runspeed)+1)
    n_runspeedlevels = len(runspeedlevels)
#    runsmoothkernel = np.array(  [1,2,3,4,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1]  )
#    runsmoothkernel = runsmoothkernel/runsmoothkernel.sum()
    runsmoothkernel = neph.create_kernel(window_width=max_runspeed/50,sampling_period=delta_runspeed)
    
    
    n_trials = len(block.segments)
    n_timecourse, n_neurons = block.segments[0].analogsignals[0].shape


    # setup stimulus:
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [ [  [[ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                             [[ [blv[1]], [45],    [5000] ], [ [blv[1]],[135],     [5000] ]  ] ],\
                          [  [[ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                             [[ [bla[1]],  [45],[5000] ], [ [bla[1]],   [45],[10000] ]    ] ],\
                          [  [[  blv,  [],[] ],   [ bla,   [], [] ] ],\
                             [[ [blv[1]],  [45],[5000] ], [ [bla[1]],   [45],[5000] ]   ] ] \
                        ]
        

#    for cx,comparison in enumerate(['visual','audio','choice']):
#        for sx,stimulusspecificity in enumerate(['all','identical']):





    # data structure is a runspeedlevels   x    timecourse    x     neurons   
    H = np.zeros( (n_runspeedlevels,n_timecourse,n_neurons) )
    
    runlevel_counts = np.zeros((n_runspeedlevels,n_timecourse))
    j=0
    for trial in block.segments:
        j+=1
        if dn=='DT031' and j>292: break
        neural = trial.analogsignals[0]
        run = trial.analogsignals[1]
#        print(j,neural.shape,run.shape,n_trials)
        # find runspeed levels for the trial along timecourse axis
        ixrs = np.floor( (run/delta_runspeed).magnitude ).astype('int16')
        ixrs[ixrs>n_runspeedlevels-1]=n_runspeedlevels-1
        for t in range(n_timecourse):
            H[ixrs[t],t,:] += neural[t,:]
            runlevel_counts[ixrs[t],t] += 1

#    print(runlevel_counts.sum(axis=1))
    rlmarginal = runlevel_counts.sum(axis=1)
#    plt.figure(figsize=(18,6))
#    plt.plot(runlevel_counts.sum(axis=0))
#    plt.figure(figsize=(18,6))
#    plt.plot(runlevel_counts.sum(axis=1))
    
    runlevel_counts[runlevel_counts==0.]=1.    # avoid division by zeros - nothing to normalize
    H/=n_trials
    H = H / runlevel_counts[:,:,np.newaxis]
    CH = H.copy()
    for t in range(n_timecourse):
        for n in range(n_neurons):
            CH[:,t,n] = np.convolve(H[:,t,n], runsmoothkernel    ,'same')
        


    if 0:
        fig,ax=plt.subplots(2,6,figsize=(36,16))
    
        for nx in [0,1,2,3,4]:
    #        print(np.min(H),np.min(CH))
    #        print(np.max(H),np.max(CH))
            axs = ax[0][nx]
            figs.plottoaxis_jointtimecourserunspeedneural(H[:,:,nx],axs)
            figs.setxt(axs,[150,450],['0','3000'])
            axs.set_title('neuron %d'%(nx+1))
            axs = ax[1][nx]
            figs.plottoaxis_jointtimecourserunspeedneural(CH[:,:,nx],axs)
            figs.setxt(axs,[150,450],['0','3000'])
            axs.set_xlabel('time within trial [ms]')
        
        ax[0][0].set_ylabel('running speed [cm/s]')
        ax[1][0].set_ylabel('running speed [cm/s]')
        
        axs=ax[0][5]
        figs.plottoaxis_jointtimecourserunspeedneural(H.mean(axis=2),axs)
        figs.setxt(axs,[150,450],['0','3000'])
        axs.set_title('averaged over neurons')
        axs=ax[1][5]
        figs.plottoaxis_jointtimecourserunspeedneural(CH.mean(axis=2),axs)
        figs.setxt(axs,[150,450],['0','3000'])
        axs.set_xlabel('time within trial [ms]')
        
        
        fig.suptitle(dn+' expected rate over p(r|v,t) cond. dist.: E[r(v,t)p(r|v,t)]\nrunspeed (v) vs. firing activity (r) timecourse (t)')
        
    #    print(neural.name,neural.shape,type(neural))
    #    print(run.name,run.shape,type(run))
    
    
        save=0
        if save:
            fig.savefig(resultpath+'runspeed,neural,timecourse_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)


    if 0:
        fig,ax=plt.subplots(4,3,figsize=(24,32))
        for lidx,layer in enumerate(layernames):
            for swx, spikewidth in enumerate(spikewidthnames):
                axs = ax[lidx,swx]
                cond = (block.annotations['idents']==lidx) & (block.annotations['waveforms']==swx)
#                if not cond.any():
#                    continue
                lmask = np.where( cond )[0]
                figs.plottoaxis_jointtimecourserunspeedneural(CH[:,:,lmask].mean(axis=2),axs)
                figs.setxt(axs,[150,450],['0','3000'])
                if lidx==0: axs.set_title(spikewidth+' %d'%len(lmask))
                else: axs.set_title(' %d'%len(lmask))
                if swx==0: axs.set_ylabel(layer+'\nrunspeed [cm/s]')
                

        
        axs = ax[0,2]
        figs.plottoaxis_jointtimecourserunspeedneural(CH.mean(axis=2),axs)
        figs.setxt(axs,[150,450],['0','3000'])
        axs.set_title('all layers and celltypes')
        
        axs = ax[1,2]
        axs.plot(runspeedlevels,rlmarginal)
        axs.set_title('locomotion marginal distribution')
        axs.set_xlabel('run speed [cm/s]')
        
        fig.suptitle(dn+' expected rate over p(r|v,t) cond. dist.: E_p(r|v,t)[r(v,t)]\nrunspeed (v) vs. firing activity (r) timecourse (t)\naveraged over neurons in layers and for celltypes')

                
        save=0
        if save:
            fig.savefig(resultpath+'runspeed,neural,timecourse,layers+types_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)

    return
                
        
        










def runspeed_regression(dn,block):

    # possible choices:
    #   i.i.d. with linear regression
    #   i.i.d. with gaussian process regression
    #   combine over trajectory: rbf kernel with gaussian processes for some autoregression in running and neural activity
    
    recalculate = 0 or globalrecalculate
    method = 'gp'                  # method = {'linreg','gp'}
    cv = 10
    
    width = 500*pq.ms
    stepx = int(width/T['dt'])
    
    n_neurons = block.segments[0].analogsignals[0].shape[1]
            
#    responses = preprocess.collect_stimulusspecificresponses(block,comparisons,logical='and')
    responses = []
    run = []
    for trial in block.segments:
        if trial.annotations['block'] in [2,4]:
            responses.append(trial.analogsignals[0][:,5])        # try only single neurons
            run.append(trial.analogsignals[1])
            run[-1] -= run[-1]
#            run[-1] += ((trial.analogsignals[0][:,5].magnitude*7+35+np.random.randn(trial.analogsignals[1].shape[0],trial.analogsignals[1].shape[1])*2)*trial.analogsignals[1].units)
            run[-1] += ((trial.analogsignals[0][:,4].magnitude*7)*trial.analogsignals[1].units)
#            run[-1] += (\
#               (np.sin(np.linspace(0,4*np.pi,601))*7+\
#                35+\
#                np.random.randn(trial.analogsignals[1].shape[0],1)*2)*trial.analogsignals[1].units\
#                       )



    if recalculate:
        if method=='linreg':
            regressions = nedi.get_linearregression(responses,run,width=width,cv=cv)
        elif method=='gp':
            regressions = nedi.get_gaussianprocessregression(responses,run,width=width,cv=cv)
        else: print('method unknown')
#        pickle.dump(regressions,open(cacheprefix+'phys/runspeed_singleneuralregression,%s,w%dms-%s,%s.pck'%(method,width.magnitude,dn,continuous_method),'wb'))
    else:
        regressions = pickle.load(open(cacheprefix+'phys/runspeed_singleneuralregression,%s,w%dms-%s,%s.pck'%(method,width.magnitude,dn,continuous_method),'rb'))

    # regressions contains 2 signals: regressions[0] is score (R2), regressions[1] is actual trial by trial runspeed
    # runspeed will be:   regressions[1][trials][train-test] AnalogSignal[trajectories, mean-std]
    
    # get the coefficients
#    wx = int((len(regressions)-2)/n_neurons)
#    W = np.reshape(np.array(regressions[7:]), (wx,n_neurons,regressions[2].shape[0],regressions[2].shape[1]) ).mean(axis=0)

    
    if 1:      # single trial runspeed estimates
        whichtrials = np.random.randint(0,len(responses),size=(4,6))

        fig,ax = plt.subplots(whichtrials.shape[0],whichtrials.shape[1],figsize=(np.array(whichtrials.shape)[::-1]*6))
        fig.suptitle( '%s, runspeed* ~ GP(neural activity)\n$R^2$: train** %5.3f, test*** %5.3f\n*green, **blue, ***orange'%(dn,regressions[0][0][:,0].mean(), regressions[0][1][:,0].mean()) )
        
        times = regressions[1][0][0].times

        for j in range(whichtrials.shape[0]):
            for k in range(whichtrials.shape[1]):
                trx = whichtrials[j,k]
                axs = ax[j,k]
                
                
                # plot runspeed estimates:
                # regressions[1] is the predictions:   ( trials, trajectories, train-test, mean-std )
                d = regressions[1]

                # test estimate
                te_ix = 1
                m = d[trx][te_ix][:,0].squeeze()
                s = (d[trx][te_ix][:,2]).squeeze() #*2/np.sqrt(50)).squeeze()
                axs.plot(times,m,linewidth=2,color='maroon',alpha=0.8)
                axs.fill_between(times,m-s,m+s, color='maroon',alpha=0.3)

                # train estimate
                tr_ix = 0
                m = d[trx][tr_ix][:,0].squeeze()
                s = (d[trx][tr_ix][:,2]).squeeze()#*2/np.sqrt(100)).squeeze()
                axs.plot(times,m,linewidth=2,color='dodgerblue',alpha=0.8)
                axs.fill_between(times,m-s,m+s, color='dodgerblue',alpha=0.3)
                
                # true data:
                axs.plot(run[0].times[::stepx],run[trx][::stepx],linewidth=2,color='seagreen',alpha=0.6)
                axs.plot(run[0].times,run[trx],linewidth=1,color='darkgreen',alpha=0.2)
            
                axs.set_xlim(axs.get_xlim())
                axs.set_ylim(axs.get_ylim())
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.setxt(axs)
                
                axs.set_title('trial # %d'%trx)
    
    
    
    
    if 1:         # trial averages of runspeed estimates
        te_ix = 1
        times = regressions[1][0][0].times
        # true runspeed
        d  = np.mean( np.array(run), axis=0 )[:-1:stepx]
        # estimated runspeed
        r_ = np.array([ regressions[1][trx][te_ix] for trx in range(len(regressions[1])) ])
        d_ = np.mean( r_, axis=0 )

        fig,ax = plt.subplots(1,3,figsize=(32,12))
        
        axs = ax[0]
        axs.plot(times,d,linewidth=2,color='darkgreen',alpha=0.8)
        m = d_[:,0].squeeze()
        s = (d_[:,2]).squeeze()  # *2/np.sqrt(len(regressions[1]))).squeeze()
        axs.plot(times,m,linewidth=2,color='maroon',alpha=0.8)
        axs.fill_between(times,m-s,m+s, color='maroon',alpha=0.3)
        
        axs.set_xlim(axs.get_xlim())
        axs.set_ylim(axs.get_ylim())
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.legend(['run true','run estimate'])
        axs.set_ylabel('run speed [cm/s]')
        axs.set_xlabel('time from stimulus onset [ms]')

        
#        axs.set_title(  'Mean runspeed estimates.\nR2 train: %4.2f+/-%4.2f, test: %4.2f+/-%4.2f'\
#                      %(regressions[0][0][:,1].mean(),regressions[0][0][:,1].mean(),regressions[0][1][:,0].mean(),regressions[0][1][:,2].mean()) )
        



        axs = ax[1]          # now explore the metrics

        figs.plottoaxis_decoderrocauc(regressions[0][:2],axs,colorlist=['dodgerblue','maroon'],plottrain=True)
        figs.setxt(axs)
        axs.set_ylim([-1.05,1.05])
        figs.plottoaxis_chancelevel(axs,0)
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.legend(['train','test'])

        axs.set_ylabel('explained variance R$^2$')
        axs.set_xlabel('time from stimulus onset [ms]')
        axs.set_title('decoding runspeed')
        
#
#        axs = ax[2]    # PC angles           # only good for linear regr. fits (not caluculated yet for GP)
#        figs.plottoaxis_decoderangles(regressions[0][2:7],axs)       # plot 5 of the subspace directions
#        figs.plottoaxis_chancelevel(axs,np.pi/2,'orthogonal')
#        figs.plottoaxis_chancelevel(axs,0.0,'parallel')
#        axs.set_ylim([0-np.pi/2*0.05,np.pi/2*1.05])
#        figs.plottoaxis_stimulusoverlay(axs,T)
#        axs.set_yticks([0,np.pi/2]); axs.set_yticklabels(['$0$','$\pi/2$'])
#        axs.legend(['weights angle to PC%d'%(d+1) for d in range(5) ]);
#        
#        axs.set_ylabel('angles between\nneural coefficients and principal components')
#        axs.set_xlabel('time from stimulus onset [ms]')
#        axs.set_title('projection angles to sliding principal components')
        
        
        
        fig.suptitle(dn+' timecourse of runspeed decoder explanation from neural data\nsliding window %d ms, %d-fold CV'%(width,cv))
    
        save=0
        if save:
            fig.savefig(resultpath+'runspeed,neural,decoding,pcangles_%s_w%dms-%dms%dms_%s'%(continuous_method,width.magnitude,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)










def runspeed_taskconditions(dn,block):
# trial averaging runspeed profiles stratified over task conditions
    
    recalculate = 0 or globalrecalculate
    
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [   [ [ [blv[1]],   [45],    [] ],  [ [blv[1]],[135],     [] ]  ], \
                            [ [ [bla[1]],   [],[5000] ],    [ [bla[1]],   [],[10000] ]    ],    \
                            [ [ [blv[1]],   [],[] ],        [ [bla[1]],   [],[] ]    ],     \
                            [[],[]] \
                        ]

    taskaspects = ['visual','audio','context','choice']
    classlabels = [ ['45°','5kHz','attend visual','lick'],['135°','10kHz','attend audio','withhold lick'] ]
    
    width = 50*pq.ms
#    n_neurons = block.segments[0].analogsignals[0].shape[1]
    
    runspeedslist = []
    decoders = []
    
    for cx,stimulusIDs in enumerate(comparisongroups):
        if cx<3:
            runspeeds = preprocess.collect_stimulusspecificresponses(block,stimulusIDs,s=1)    # s=1: analogsignals[1] is runspeed
#            responses = preprocess.collect_stimulusspecificresponses(block,stimulusIDs)
        else:
            runspeeds = preprocess.collect_stimulusspecificresponses_choice(block,dn,s=1)
#            responses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
        runspeedslist.append(runspeeds)


        if recalculate:
            acrossdecoder = nedi.get_responsedecoding(runspeeds,width=width)
            # now save
            pickle.dump(acrossdecoder,open(cacheprefix+'phys/rundecodes-%s-%s-%s.pck'%(dn,continuous_method,taskaspects[cx]),'wb'))
        else:
            acrossdecoder = pickle.load(open(cacheprefix+'phys/rundecodes-%s-%s-%s.pck'%(dn,continuous_method,taskaspects[cx]),'rb'))
        decoders.append(acrossdecoder)







    colors = ['seagreen','firebrick']
    fig,ax = plt.subplots(2,4,figsize=(32,14))


    for cx,stimulusIDs in enumerate(comparisongroups):
        times = runspeeds[0][0].times

        runspeeds = runspeedslist[cx]
        axs = ax[0,cx]
        for clx in range(2):
            m = np.array(runspeeds[clx]).mean(axis=0).squeeze()
            s = np.array(runspeeds[clx]).std(axis=0).squeeze()
            
            figs.plottoaxis_plottrajectory(axs,times,m=m,s=s/2,colorlist=colors[clx],alpha=0.8,linewidth=2,label=classlabels[clx][cx])
            axs.fill_between(times, m-2*s/np.sqrt(len(runspeeds[clx])), m+2*s/np.sqrt(len(runspeeds[clx])), color=colors[clx],alpha=0.4)

#            figs.plottoaxis_plottrajectory(axs,times,m=m,s=s/2,colorlist='seagreen',alpha=0.8,linewidth=1,label=classlabels[1-clx][cx])
#            axs.fill_between(times, m-2*s/np.sqrt(len(runspeeds[clx])), m+2*s/np.sqrt(len(runspeeds[clx])), color='seagreen',alpha=0.4)


#            axs.legend([classlabels[clx][cx]])
            axs.set_ylim(-5,105)
            axs.set_xlim(-1500,4500)
            figs.setxt(axs)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0)
            
            if cx==0: axs.set_ylabel('runspeed [cm/s]')
            axs.set_title('%s'%taskaspects[cx])

        axs.legend()


        # decoding task variables from runspeeds
        decoder = decoders[cx]
        axs = ax[1,cx]
        figs.plottoaxis_decoderrocauc(decoder,axs)
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,0.5)
        axs.set_yticks([0.5,1.0])
        if cx==0: axs.set_ylabel('decoder classification accuracy')

        axs.set_xlabel('time from stimulus onset [ms]')


        

    fig.suptitle(dn+' trial distribution of running speed in task conditions')

    save=0
    if save:
        fig.savefig(resultpath+'runspeed,acrosstaskconditionclasses_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)











def corr_runspeed_dbnv(dn,block):

    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [   [ [ [2,4],   [45],    [] ],  [ [2,4],[135],     [] ]  ], \
                            [ [ [2,4],   [],[5000] ],    [ [2,4],   [],[10000] ]    ],    \
                            [ [ [blv[1]],   [],[] ],        [ [bla[1]],   [],[] ]    ],     \
                            [[],[]], \
                            [ [ [2,4],   [],[] ]  ],    \
                        ]

    # a special consideration, that the 5th row will be used for all 4 as a 3rd "class" to use all trials together
    # this will be the actual correlation that takes into account the variability related to the task variable

    taskaspects = ['visual','audio','context','choice']
    classlabels = [ ['45°','5kHz','attend visual','lick'],['135°','10kHz','attend audio','withhold lick'],['all','all','all','all'] ]
    taskcolors = [ ['darkgreen','darkred'],['mediumseagreen','orangered'] ]
    variablecolors = ['navy','darkgreen','mediumvioletred','darkorange']
    correlationcolors = ['dodgerblue','lightgreen','fuchsia','gold']

    width = 50*pq.ms
    n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
    depth_idx = 150
    
    runspeedslist = []
    acrossdecoders = []


    # collect runspeeds in trials
    for cx,comparison in enumerate(taskaspects):
        stimulusIDs = comparisongroups[cx]
        if cx!=3:
            runspeeds = preprocess.collect_stimulusspecificresponses(block,stimulusIDs,s=1)    # s=1: analogsignals[1] is runspeed
            # responses = preprocess.collect_stimulusspecificresponses(block,stimulusIDs)
        else:
            runspeeds = preprocess.collect_stimulusspecificresponses_choice(block,dn,s=1)
            # responses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
        
        # now add the alltrial runspeeds as a 3rd class (not usual comparisongroup wording, but with choice it's not possible to do other way)
        runspeeds.append(  preprocess.collect_stimulusspecificresponses(block,comparisongroups[4],s=1)[0]   )

        runspeedslist.append(runspeeds)



    # collect decision vector bases    
    c_db_vectors = np.zeros((n_neurons,len(taskaspects)))
    for cx,comparison in enumerate(taskaspects):
        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
        acrossdecoders.append( acrossdecoder )
        wx = int((len(acrossdecoder)-7)/n_neurons)
        coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
        # mean of on stimulus
        c_db_vectors[:,cx] = coeff[:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,0].mean(axis=1)



    # c_db_matrix projects to the 3 dbnvs each row:   c_db_matrix @ neural_input
    # Q,R = np.linalg.qr(c_db_matrix)     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
    # V = Q                          # activity transformer matrix
    # print('c_db',c_db_matrix.shape,'Q',Q.shape,'R',R.shape)
    # B = Q.T @ c_db_matrix / np.linalg.norm(c_db_matrix,axis=0,keepdims=True)        # transformed dbnvs (in columns) in orthonormal coordinates (in rows)
    # print( 'V',V.shape, 'B', B.shape )




        
    projected_dynamics = []    # this will be (tasks,classes,[trials,trajectory])
    for cx,comparison in enumerate(taskaspects):
        # now project each activity onto the dbnv basis vectorspace

        # collect neural responses
        if cx!=3: acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx])
        else: acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)

        # now project all trials regardless of task variables, to the taskaspect DV (only a single list in acrossresponses)
        acrossresponses.append(preprocess.collect_stimulusspecificresponses(block,comparisongroups[4])[0])

        # extract projector vector relevant to the task        
        V = c_db_vectors[:,cx]

        projected_dynamic = []
        for aix,classresponse in enumerate(acrossresponses):
            # project the entire dynamics trajectory onto the DVs axes for each trial
            if len(classresponse)>0:
                X = np.dot(np.array(classresponse),V)  #    (trials,trajectory) = (trials,trajectory,neurons) * (,,neurons)
            else:
                X = np.empty((0,n_trajectory,n_ba))
            projected_dynamic.append(X)
        
        
        projected_dynamics.append(projected_dynamic)







    # now we have both runspeed and projected neuronal activity onto dbnvs
    # we can compare the two trial by trial, with correlations:
    
    corr = np.zeros((len(taskaspects),3,n_trajectory,2,2))   # (aspects,classes+all,trajectory,{normal,shuffled},{r,p})

    for cx,comparison in enumerate(taskaspects):
        for clx in range(3):      # note that for the non-class separated all-trials correlation is the 3rd "class"
            # print(cx,clx,len(runspeedslist[cx][clx]),len(projected_dynamics[cx][clx]))
            shuffle_mask = np.random.permutation(len(projected_dynamics[cx][clx]))
            for t in range(n_trajectory):
                # print(np.array(runspeedslist[cx][clx]).shape, np.array(projected_dynamics[cx][clx]).shape )
                corr[cx,clx,t,0,:] = sp.stats.pearsonr(  np.array(runspeedslist[cx][clx])[:,t,0],\
                                                       np.array(projected_dynamics[cx][clx])[:,t]   )
                corr[cx,clx,t,1,:] = sp.stats.pearsonr(  np.array(runspeedslist[cx][clx])[:,t,0],\
                                                       np.array(projected_dynamics[cx][clx])[shuffle_mask,t]   )
                










    
    # FIGURE
    
    
    fig,ax = plt.subplots(8,len(taskaspects),figsize=(len(taskaspects)*8,8*6))
    
    
    
    
    # for visual control for runspeed measurement errors

    for cx in [0,1,2,3]:
        ax[0,cx].remove()
        axs = fig.add_subplot(8,1,1)
        times = block.channel_indexes[1].analogsignals[0].times.rescale('min')
        axs.plot(times,block.channel_indexes[1].analogsignals[0],color='r',alpha=0.6)
        axs.set_xlim(times[0],times[-1])
        # axs.set_xlabel('in session time [min]                        ',labelpad=-20)
        axs.set_xlabel('in session time [min]',labelpad=-5)
        axs.set_ylabel('runspeed [cm/s]')
    
    
    
    
    for cx,comparison in enumerate(taskaspects):

        axs = ax[1,cx]
        
        times = runspeeds[0][0].times

        runspeeds = runspeedslist[cx]
        axs = ax[1,cx]
        for clx in range(2):             # we only need to show the two task-specific trial groupings, not the all (3rd class)
            m = np.array(runspeeds[clx]).mean(axis=0).squeeze()
            s = np.array(runspeeds[clx]).std(axis=0).squeeze()
            
            figs.plottoaxis_plottrajectory(axs,times,m=m,s=s/2,colorlist=variablecolors[cx],alpha=0.8-0.4*clx,linewidth=2,label=classlabels[clx][cx])
            axs.fill_between(times, m-2*s/np.sqrt(len(runspeeds[clx])), m+2*s/np.sqrt(len(runspeeds[clx])), color=variablecolors[cx],alpha=0.4-0.15*clx)

#            figs.plottoaxis_plottrajectory(axs,times,m=m,s=s/2,colorlist='seagreen',alpha=0.8,linewidth=1,label=classlabels[1-clx][cx])
#            axs.fill_between(times, m-2*s/np.sqrt(len(runspeeds[clx])), m+2*s/np.sqrt(len(runspeeds[clx])), color='seagreen',alpha=0.4)


        # axs.set_ylim(-5,105)
        axs.set_xlim(-1400,4400)
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        # figs.plottoaxis_chancelevel(axs,0)
        
        if cx==0: axs.set_ylabel('runspeed [cm/s]')
        axs.set_title('%s'%taskaspects[cx])

        axs.legend()
    
    
    
    

        # plot neural projection firing rate averages
        axs = ax[2,cx]
        
        times = acrossresponses[0][0].times

        for clx in range(2):      # we only need to show the two task-specific trial groupings, not the all (3rd class)

            trials = projected_dynamics[cx][clx]            # averaging over trials
            trial_av = trials.mean(axis=0)        # trial averages
            trial_e = trials.std(axis=0)*2/np.sqrt(len(trials))
    
            axs.plot( times, trial_av, lw=3,color=variablecolors[cx],alpha=0.8-0.4*clx,label=classlabels[clx][cx])
            axs.fill_between( times,trial_av-trial_e,trial_av+trial_e,\
                     color=variablecolors[cx], alpha=0.4-0.15*clx)

        axs.set_xlim(-1400,4400)

        axs.legend(frameon=False)

        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,0)

        figs.setxt(axs)
        
        if cx==0: axs.set_ylabel('neural activity\nprojected onto dbnv')


        axs = ax[3,cx]

        figs.plottoaxis_decoderrocauc(acrossdecoders[cx][:2],axs,colorlist=['black',variablecolors[cx]],plottrain=False)       # plot the performance
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,0.5)
        figs.setxt(axs)
        if cx==0: axs.set_ylabel('decoder accuracy')


        # correlation between runspeed and firing activity projected onto relevant dbnv
        for rx in [0,1]:     # origina and shuffled neural
            for sx in [0,1]:        #   r  and  p
                for clx in range(3):   # classes
                    #r
                    axs = ax[4+2*rx+sx,cx]
                    
                    corr_all = corr[cx,clx,:,rx,sx]
                    notsignificant = corr[cx,clx,:,rx,1]>0.05
                    corr_sig = corr[cx,clx,:,rx,sx].copy()
                        
                    
                    corr_sig[notsignificant] = np.nan
                    if clx<2:
                        axs.plot(times,corr_all,color=variablecolors[cx],alpha=0.1-0.05*clx)
                        axs.plot(times,corr_sig,color=variablecolors[cx],lw=2+(clx==2)*2,alpha=1-0.4*clx,label=classlabels[clx][cx])
                    else: 
                        axs.plot(times,corr_all,color=correlationcolors[cx],lw=2,alpha=0.1)
                        axs.plot(times,corr_sig,color=correlationcolors[cx],lw=4,alpha=1,label=classlabels[clx][cx])
            
                if sx==0: axs.set_ylim(-1,1)
                else: axs.set_ylim(0,1)
                axs.set_xlim(-1400,4400)
    
                if cx==0: axs.set_ylabel(['ordered\n','shuffled\n'][rx]+['correlation r','correlation p'][sx])
    
                figs.setxt(axs)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.05)
        
                axs.legend(frameon=False)







    fig.suptitle(dn+'    %d neurons'%n_neurons+'\n(1st row) runspeed over the session (incorrect waiting times cut)\n (2-3rd row) runspeed and dbnv projected neural activity trial averages,\n(4-5th row) correlation of (2) and (3) over trials by timepoints (faint lines: non-significant')

    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'runspeed,dbnv,correlation_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
        # fig.savefig(resultpath+'runspeed,dbnv,correlation-shuffle,iid_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)










def decoder_equalizedrunspeeddistribution(dn,block):
    # load a trial-index indexed list for each timecourse point, where the values are runspeed distributions-equalized trial indices: include or not include


    recalculate = 0 or globalrecalculate
    dostats = 0
    doplot = 0 or globaldoplot

    
    decoderwidth = 5         # number of consecutive timepoints to average runspeed over and feed to decoder as features
    patchwidth = 3             # number of consecutive timepoints to patch together, original sampling rate 100 Hz -> 10 ms bins
    n_bootstrap = 10             # number of multiple subsamplings to average over
    
    runshifts_times = np.arange(-200,201,100)*pq.ms
    # runshifts_times = np.arange(-100,101,200)*pq.ms
    # runshifts_times = np.arange(0,100,100)*pq.ms      # use this for runspeed histograms
    # runshifts = (runshifts_times / T['dt']).astype(np.int16)
    
    max_runspeed = 100.*pq.cm/pq.s
    delta_runspeed = 5.*pq.cm/pq.s
    runspeedlevels = np.linspace(0*pq.cm/pq.s, max_runspeed, int(max_runspeed/delta_runspeed)+1)
    n_runspeedlevels = len(runspeedlevels)-1
    
    
    n_trials = len(block.segments)
    n_trajectory_full, n_neurons = block.segments[0].analogsignals[0].shape
    # trajectory for the final accuracy timecouse: decoderwidth for features and runspeedlevels, and hump together patchwidth number of them
    n_trajectory = n_trajectory_full//decoderwidth//patchwidth

    # setup stimulus:
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [  [ [ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                           [ [ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ] ],\
                           [[  [blv[1]],  [],[] ],   [ [bla[1]],   [], [] ] ],\
                           [ [], [] ] \
                        ]
        

    taskaspects = ['visual','audio','context','choice']
    
    
    equalized_trialindex = np.empty( (len(taskaspects),2,n_trajectory), dtype=np.object)   # this will be a 4 by 2 by 601/decoderwidth list, with the aspects and their two conditions, holding lists of trial indices
    all_trialindex = np.empty( (len(taskaspects),2,n_trajectory), dtype=np.object)
    lens_trialindex = np.zeros( (2,len(taskaspects),2,n_trajectory) )  # ({eq,noeq},taskaspects,class,trajectory)
    stats_trialindex = np.zeros((2,len(taskaspects),2)) # ({eq,noeq},taskaspects,{mean,std})
    
    
    acrossdecoders_runshifts = []
    acrossdecoders_noeq_runshifts = []
    
    for rshx,rsh in enumerate(runshifts_times):
        print('shifting runspeed',rsh)
        
        acrossdecoders = []
        acrossdecoders_noeq = []
        
        
        for cx,comparison in enumerate(taskaspects):
            print('calculating for '+comparison)
            if cx!=2: continue
            stimulusIDs = comparisongroups[cx]
            
            if cx!=3:
                runspeeds_full = preprocess.collect_stimulusspecificresponses(block,stimulusIDs,s=1)    # s=1: analogsignals[1] is runspeed
                responses_full = preprocess.collect_stimulusspecificresponses(block,stimulusIDs)
            else:
                runspeeds_full = preprocess.collect_stimulusspecificresponses_choice(block,dn,s=1)
                responses_full = preprocess.collect_stimulusspecificresponses_choice(block,dn)
            
            # patch da batch
            
            # create runspeed average distributions             
            runspeeds_full = neph.downsamplesignals(runspeeds_full,decoderwidth)
            runspeeds_full = neph.shift_timecourse_analogsignal(runspeeds_full,rsh)
            
            # create width features for decoders, and then create virtual trials below
            responses_full = neph.convertintervaltofeature(responses_full,decoderwidth,decoderwidth)
        
            
    
            if recalculate:
                
                # create a structure that will collect and average multiple runs of subsampling
                acrossdecoder_m =  [ neo.AnalogSignal(  np.zeros( (responses_full[0][0].shape[0]//patchwidth,3) ),\
                       name='subsample average roc', t_start=responses_full[0][0].t_start,\
                       sampling_period=responses_full[0][0].sampling_period*patchwidth, units=responses_full[0][0].units   )  for clx in range(len(responses_full))  ]
                acrossdecoder_noeq_m = acrossdecoder_m.copy()
                
                
                for bx in range(n_bootstrap):        # do the multiple bootstrap runs
                    
                    runspeeds = runspeeds_full       # start each bootstrap clean
                    responses = responses_full
                
                    # save control states
                    runspeeds_noeq = runspeeds
                    responses_noeq = responses



                    # patch together a couple of forward times as observations:
                    runspeeds = neph.virtualizetrials_analogsignal(runspeeds,patchwidth)


                    # go through the trial selection process: equalizing with runspeed distribution
                    # patch batches together to increase number of trials i.e. observations
                    
                    r_dists = np.zeros((n_trajectory,n_runspeedlevels,2))
                    for t in range(n_trajectory): # do separately each timecourse point but use all trials
                        bins_data_indices = []
                        for clx in [0,1]: # go through the two classes
                            r_dists[t,:,clx],_ = np.histogram(np.array(runspeeds[clx])[:,t],bins=runspeedlevels)
                            bins_data_indices.append( neph.dataindexhistogram(np.array(runspeeds[clx])[:,t],bins=runspeedlevels) )
                        
                        # r_dists_eq will hold the difference in the number of samples at each runspeed level
                        # positive numbers mean that we need to undersample the second class, negative to undersample the first class
                        r_dists_eq_t = r_dists[t,:,0]-r_dists[t,:,1]
                        equalized_trialindex[cx, 0, t] = [];equalized_trialindex[cx, 1, t] = []
                        # print(cx,t,r_dists_eq_t)
                        for rx,rl in enumerate(r_dists_eq_t): # do the subsampling for the more numerous class
                            if len(bins_data_indices[0][rx])>0 and len(bins_data_indices[1][rx])>0:
                                if rl==0:
                                    for clx in [0,1]:
                                        equalized_trialindex[cx, clx, t].extend(bins_data_indices[clx][rx])
                                else:
                                    i_more = int(rl<0)
                                    i_less = int(rl>0)
                                    
                                    subsamples = np.random.permutation(len(bins_data_indices[i_more][rx]))[:len(bins_data_indices[i_less][rx])]
                                    # print(rx,rl,'><',i_more,i_less)
                                    # print(bins_data_indices[i_more][rx],bins_data_indices[i_less][rx])
                                    # print(subsamples)
                                    equalized_trialindex[cx, i_less, t].extend(bins_data_indices[i_less][rx])
                                    equalized_trialindex[cx, i_more, t].extend(bins_data_indices[i_more][rx][subsamples])
                        
                        # create a control list, with the same number of trials in each class
                        for clx in [0,1]:
                            all_trialindex[cx, clx, t] = np.random.permutation(len(responses[clx])*patchwidth)[:len(equalized_trialindex[cx,clx,t])]
                        
                        for clx in [0,1]: # collect for statistics
                            lens_trialindex[0,cx,clx,t] = len(equalized_trialindex[cx, clx, t])
                            lens_trialindex[1,cx,clx,t] = len(all_trialindex[cx, clx, t])
                            
                    
                    if 0:
                        
                        if 0:
                            pickle.dump(r_dists,open(cacheprefix+'locomotion/runspeeddistributions_%s-%s-%s.pck'%(dn,continuous_method,comparison),'wb'))
                        
                        r_dists_eq = r_dists[:,:,0]-r_dists[:,:,1]
                        
                        ts0 = 0
                        skip = 4
                        lump = 1
                        n_ts = n_trajectory//skip
                        # ts_m = 10//decoderwidth 
                        
                        
                        fig,ax = plt.subplots(2,n_ts,figsize=(n_ts*7,2*7))
                        
                        for k,t in enumerate(np.arange(0,n_trajectory,skip)):
                            
                            axs = ax[0,k]
                            axs.plot(runspeedlevels[:-1],r_dists[t:t+lump,:,0].sum(axis=0),color='rebeccapurple',label='av')
                            axs.plot(runspeedlevels[:-1],r_dists[t:t+lump,:,1].sum(axis=0),color='deeppink',label='aa')
                            if k==0:
                                axs.legend(frameon=False)
                                axs.set_ylabel('runspeed histogram\n# trials')
                            axs.set_title('timecourse # %d ms'%(runspeeds[0][0].times[t]))
                            
            
                            axs = ax[1,k]
                            axs.plot(runspeedlevels[:-1],r_dists_eq[t:t+lump,:].sum(axis=0),color='mediumvioletred')
                            figs.plottoaxis_chancelevel(axs)
                            if k==0: axs.set_ylabel('runspeed histogram difference\n# trials')
                            axs.set_xlabel('runspeed cm/s')
                        
                        fig.suptitle(dn+' runspeed histogram, %s differences $\Delta$r=%d cm/s, averaging=%dms, patch=%dms'%(comparison,delta_runspeed.magnitude,decoderwidth*T['dt'],decoderwidth*patchwidth*T['dt'].magnitude))
                        if 0:
                            fig.savefig(resultpath+'locomotionequalization_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                        
                        return
        
                    
                    if dostats:
                        # collect stats for loco-equalized:
                        stats_trialindex[0,cx,0] = np.mean(lens_trialindex[0,cx,0,:])
                        stats_trialindex[0,cx,1] = 2*np.std(lens_trialindex[0,cx,0,:])/np.sqrt(n_trajectory)
                        # collect control (this is originally the same at each timepoint):
                        stats_trialindex[1,cx,0] = min((lens_trialindex[1,cx,0,0],lens_trialindex[1,cx,1,0])) # find the class with the smaller number of trials that was available, as only fair
                        if cx==2:
                            print('stats',dn,comparison,' %d +- %4.2f  /  %d = %4.2f %%'%\
                              (stats_trialindex[0,cx,0],stats_trialindex[0,cx,1],\
                               stats_trialindex[1,cx,0],\
                               stats_trialindex[0,cx,0]/stats_trialindex[1,cx,0]*100)  )
                        continue
        
        
        
                    
                    # patch together a couple of forward times as observations:
                    responses = neph.virtualizetrials_analogsignal(responses,patchwidth)
                    responses_noeq = neph.virtualizetrials_analogsignal(responses_noeq,patchwidth)

        
                    # now runspeed-equalize by the collected indices:
                    responses,sample_analogsignal = neph.equalizebyindex(responses,equalized_trialindex[cx,:,:])
                    # patch together a couple of forward times as observations:
                    # responses,sample_analogsignal = neph.virtualizetrials_maskedarray(responses,patchwidth,sample_analogsignal)
                
                    # bootstrap for the controls
                    responses_noeq,_ = neph.equalizebyindex(responses_noeq,all_trialindex[cx,:,:])
                    # responses_noeq,_ = neph.virtualizetrials_maskedarray(responses_noeq,patchwidth,sample_analogsignal)
                    
                    
                    
                    
                    
        
                    acrossdecoder = nedi.get_pointwisedecoder_tdistribution(responses,sample_analogsignal)
                    acrossdecoder_noeq = nedi.get_pointwisedecoder_tdistribution(responses_noeq,sample_analogsignal)
                    
                    
    
                    acrossdecoder_m = [ acrossdecoder_m[clx] + acrossdecoder[clx]/n_bootstrap \
                                       for clx in range(len(responses)) ]
        
                    acrossdecoder_noeq_m = [ acrossdecoder_noeq_m[clx] + acrossdecoder_noeq[clx]/n_bootstrap \
                                       for clx in range(len(responses_noeq)) ]
    
                pickle.dump(acrossdecoder_m,open(cacheprefix+'locomotion/responsedecodes,locomotionequalized,shift%+dms_angles-%s-%s-%s.pck'%(rsh.magnitude,dn,continuous_method,comparison),'wb'))
    
    
                pickle.dump(acrossdecoder_noeq_m,open(cacheprefix+'locomotion/responsedecodes,locomotionequalized-control,shift%+dms_angles-%s-%s-%s.pck'%(rsh.magnitude,dn,continuous_method,comparison),'wb'))
    
            else:
                acrossdecoder_m = pickle.load(open(cacheprefix+'locomotion/responsedecodes,locomotionequalized,shift%+dms_angles-%s-%s-%s.pck'%(rsh.magnitude,dn,continuous_method,comparison),'rb'))
                acrossdecoder_noeq_m = pickle.load(open(cacheprefix+'locomotion/responsedecodes,locomotionequalized-control,shift%+dms_angles-%s-%s-%s.pck'%(rsh.magnitude,dn,continuous_method,comparison),'rb'))
    
    
    
            
            acrossdecoders.append(acrossdecoder_m)
            acrossdecoders_noeq.append(acrossdecoder_noeq_m)
            
        acrossdecoders_runshifts.append(acrossdecoders)
        acrossdecoders_noeq_runshifts.append(acrossdecoders_noeq)




    if doplot:
        taskcolors = ['navy','darkgreen','mediumvioletred','darkorange']
        fig,ax = plt.subplots(4,len(taskaspects),figsize=(len(taskaspects)*6,4*6))
        
        s = len(runshifts_times)//2+1-1  # find the zero shift position
        
        acrossdecoders = acrossdecoders_runshifts[s]
        acrossdecoders_noeq = acrossdecoders_noeq_runshifts[s]
        
        for cx,comparison in enumerate(taskaspects):
            
            axs = ax[0,cx]

            figs.plottoaxis_decoderrocauc(acrossdecoders_noeq[cx][:2],axs,colorlist=['white',taskcolors[cx]])       # plot the performance
            figs.setxt(axs)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)
            if cx==0: axs.legend(['train','test'],frameon=False); axs.set_ylabel('original\nroc accuracy')
            axs.set_title(comparison)



            axs = ax[1,cx]

            figs.plottoaxis_decoderrocauc(acrossdecoders[cx][:2],axs,colorlist=['white',taskcolors[cx]])       # plot the performance
            figs.setxt(axs)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)
            if cx==0: axs.legend(['train','test'],frameon=False); axs.set_ylabel('locomotion dist. equalized\nroc accuracy')


            axs = ax[2,cx]

            figs.plottoaxis_decoderrocauc(acrossdecoders_noeq[cx][:2],axs,colorlist=['white','black'],plottrain=False,onlysem=True,label='original')       # plot the performance
            figs.plottoaxis_decoderrocauc(acrossdecoders[cx][:2],axs,colorlist=['white',taskcolors[cx]],plottrain=False,onlysem=True,label='locomotion dist. equalized')       # plot the performance
            figs.setxt(axs)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)
            if cx==0: axs.legend(frameon=False); axs.set_ylabel('both\nroc accuracy')


            

            axs = ax[3,cx]
            times = acrossdecoders_noeq_runshifts[rshx][cx][1].times
            for rshx,rsh in enumerate(runshifts_times):
                tc = np.array(mcs.to_rgb(taskcolors[cx]))
                tc /= max(tc)
                colors = np.array([[1,1,1],tc]) * (rshx+1) /  (len(runshifts_times)+2)
                axs.plot(times,acrossdecoders_noeq_runshifts[rshx][cx][1][:,0],lw=2,color=colors[0],label='original %+dms'%rsh)       # plot the performance
                axs.plot(times,acrossdecoders_runshifts[rshx][cx][1][:,0],lw=2,color=colors[1],label='loco. eq. %+dms'%rsh)       # plot the performance
            figs.setxt(axs)
            axs.set_ylim([0.45,1.01])
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)
            if cx==0: axs.legend(frameon=False); axs.set_ylabel('both, runspeed shifts\nroc accuracy')



        fig.suptitle(dn+' decoders: original and runspeed distribution-equalized at each timepoint; last row: runspeed shifted')
        
        
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'decoders,locomotionequalized-p%dms_%s-%dms%dms_%s'%(decoderwidth*patchwidth*T['dt'],continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
            # fig.savefig(resultpath+'runspeed,dbnv,correlation-shuffle,iid_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)

            

    
    
    return








def getmovementpctimecourses(dn, multimodalonly=True, calculuslevel='motion'):

    starttime = T['starttime'].magnitude
    stimtime = (T['stimendtime'] - T['stimstarttime']).magnitude
    endtime = T['endtime'].magnitude
    fps = T['videofps']

    movementpcs = preprocess.loadvideopca(dn, calculuslevel=calculuslevel)
    triallist = preprocess.loadexperimentdata(dn, multimodalonly=multimodalonly)
    # triallist['start'] /= 1000
    n_trials = triallist.values.shape[0]
    n_pcs = movementpcs.values.shape[1]-1      # the first column is the timestamp
    movementpcs['time'] *= 1000       # convert to milliseconds
    timestampfps = np.arange( starttime, endtime, 1000./fps )   # milliseconds fps
    ntimestamps = len(timestampfps)      # timestamps cap

    # print(triallist)
    # print(movementpcs)

    # dimensions will be (trials,timecourse,pcs)
    movementpcslist = np.zeros((n_trials,len(timestampfps),n_pcs))

    # find each trial still in the video
    for trx,trialstart in enumerate(triallist['start'][triallist['start']+stimtime<movementpcs['time'].iloc[-1]]):
        slicemask = (movementpcs['time']>=trialstart+starttime) & (movementpcs['time']<trialstart+endtime)              # check if trial is in the video
        if sum(slicemask)<ntimestamps: continue                           # if a trial is not fully within the video, skip
        movementpcslist[trx,:,:] = (movementpcs[slicemask].iloc[:,1:]).iloc[:ntimestamps,:]     # skip the first column (time), and only 6 sec to enter

    
    return movementpcslist, timestampfps











def display_movementpca(dn):

    doplot = 0 or globaldoplot

    triallist = preprocess.loadexperimentdata(dn, multimodalonly=True)
    movementpcslist, timestampfps = getmovementpctimecourses(dn)          # movementpcslist dimensions are (trials,timecourse,pcs)

    # return
    print('movement video derived motion differential PCs')
    print(movementpcslist.shape)
    print(movementpcslist.mean(axis=0))


    criteria = [[c,v,a] for c in [1,3] for v in [45,135] for a in [5000,10000]]
    print(criteria)

    maxpc = 2

    X = []
    ttx = []
    for cx,c in enumerate(criteria):
        mask = (triallist['block']==c[0]) & (triallist['degree']==c[1]) & (triallist['freq']==c[2])
        X.append(movementpcslist[mask,:,:maxpc])
        print('mask criteria sum', c, sum(mask))
        ttx.append(np.where(mask)[0][0])
    
    
    if doplot:
        colors = ['navy','seagreen','pink','purple','darkgreen','turquoise','orange','gold']
        fig, ax = plt.subplots(2,4,figsize=(4*8,2*8))

        for hx in range(2):
            for wx in range(4):
                ix = hx*4+wx
                c = criteria[ix]
                axs = ax[hx,wx]
                for k in range(maxpc):
                    m = X[ix][:,:,k].mean(axis=0)
                    e = X[ix][:,:,k].std(axis=0)/np.sqrt(X[ix][:,:,k].shape[0])
                    axs.plot(timestampfps, m,color=colors[ix],alpha=1.-k/(maxpc+1), label='pc%d'%k)
                    axs.fill_between(timestampfps, m-e, m+e, color=colors[ix], alpha=1./3-k/(3*maxpc+1))
                    # axs.plot(timestampfps, movementpcslist[ttx[ix],:,k],color=colors[ix],alpha=1.-k*0.2, label='pc%d'%k)

                if wx==0: axs.set_ylabel('trialaveraged PC 1-%d +/-sem'%maxpc)
                axs.set_title('%s %d %d'%(['visual','audio'][c[0]==3],c[1],c[2]))
                figs.setxt(axs)
                figs.plottoaxis_stimulusoverlay(axs,T)
                # print(hx,wx, X[ix][:,:,:maxpc].mean(axis=(0,1)) )


        fig.suptitle('%s, PC1..%d of movement video, multimodal trials'%(dn,maxpc))


        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'movementvideo-PCs-visualcolors_%s'%(dn)+ext)     #PCs, all






def decode_movementpca_singleconcat(dn,block,calculuslevel='motion'):
    # create bins of 20 Hz (fps) for neural data
    # use 20 Hz video PC timecourses
    # decode one from the other with single cells and single PCs
    # test shared covariances
    # calculuslevel = {'posture','motion'}   : posture is absolute image intensity, motion is frame to frame intensity difference

    recalculate = 0 or globalrecalculate

    doplot = 1 or globaldoplot

    th = 0.05   # significance threshold


    triallist = preprocess.loadexperimentdata(dn, multimodalonly=True)
    movementpcslist, timestampfps = getmovementpctimecourses(dn, calculuslevel=calculuslevel)          # movementpcslist dimensions are (trials,timecourse,pcs)

    allcomplextrials = [ [[2,4], [], []] ]
    neuralactivity = preprocess.collect_stimulusspecificresponses(block, allcomplextrials)
    downsamplerate = 5   # length * binsize / fps
    neuralactivity = np.array(neph.downsamplesignals(neuralactivity, downsamplerate)[0])

    n_trials = neuralactivity.shape[0]
    n_timecourse = neuralactivity.shape[1]
    n_neurons = neuralactivity.shape[2]
    n_pcs = movementpcslist.shape[2]

    print('downsampled neural activity: ', neuralactivity.shape)
    print('movement PCs: ', movementpcslist.shape)

    movementpcslistconcat = np.reshape(movementpcslist, (-1,n_pcs))
    neuralactivityconcat = np.reshape(neuralactivity, (-1,n_neurons ))
    n_observations = movementpcslistconcat.shape[0]

    # # calculate in 1 sec bins
    # broadsamplerate = 20
    # neuralactivitylowres = neph.downsamplearray(neuralactivity, broadsamplerate, axis=1)    # gather spike counts in 1 sec bins
    # movementpcslistlowres = neph.downsamplearray(movementpcslist, broadsamplerate, axis=1)    # average movement PCs in 1 sec bins
    # n_timecourse_lowres = neuralactivitylowres.shape[1]
    # timestampfpslowres = timestampfps[::broadsamplerate]
    # timestampfpslowres += ((timestampfpslowres[1] - timestampfpslowres[0])//2)




    # calculate per neuron per pc correlations along the trials
    if recalculate:
        # two values:  neural->behaviour, behaviour->neural
        rhos = np.zeros((2,n_neurons,n_pcs))
        ps = np.zeros((2,n_neurons,n_pcs))
        accs = np.zeros((2,n_neurons,n_pcs,2,3))          # (train/test, stats)
        coefs = np.zeros((2,n_neurons,n_pcs,3))           # (stats)
        ys_te = np.zeros((2,n_observations,n_neurons,n_pcs))               # cv perdiction for each observation
        for nx in range(n_neurons):
            print('neuron:', nx)
            for px in range(n_pcs):
                # if nx>10 or px>6: break

                # correlation significance
                rho,p = sp.stats.linregress(neuralactivityconcat[:,nx], movementpcslistconcat[:,px])[2:4]     # rho and p
                rhos[0,nx,px] = rho
                ps[0,nx,px] = p
                rho,p = sp.stats.linregress(movementpcslistconcat[:,px], neuralactivityconcat[:,nx])[2:4]     # rho and p
                rhos[1,nx,px] = rho
                ps[1,nx,px] = p

                # CV regressions
                accs[0,nx,px,:,:],coefs[0,nx,px,:],ys_te[0,:,nx,px] = \
                       nedi.get_singleregression(neuralactivityconcat[:,nx][:,np.newaxis], movementpcslistconcat[:,px][:,np.newaxis])
                accs[1,nx,px,:,:],coefs[1,nx,px,:],ys_te[1,:,nx,px] = \
                       nedi.get_singleregression(movementpcslistconcat[:,px][:,np.newaxis], neuralactivityconcat[:,nx][:,np.newaxis])



        pickle.dump((rhos,ps,accs,coefs,ys_te),\
                open(cacheprefix+'locomotion/%spcs,neurons-single,concatenated_%s.pck'%(calculuslevel, dn), 'wb'))
    else:
        rhos,ps,accs,coefs,ys_te = pickle.load(open(cacheprefix+'locomotion/%spcs,neurons-single,concatenated_%s.pck'%(calculuslevel, dn), 'rb'))





    if doplot:

        # neural -> behaviour
        n_pcs = 6
        n_neurons = 10
        fig, ax = plt.subplots(n_pcs,n_neurons,figsize=(n_neurons*6,n_pcs*6))
        fig.subplots_adjust(wspace=0.,hspace=0.)

        ht = 0.07 # text height
        for nx in range(n_neurons):
            for px in range(n_pcs):
                axs = ax[px,nx]

                axs.text(0.01,0.7+3*ht,'$\\rho$=%6.4f, w=%6.4f$\\pm$%6.4f'%(rhos[0,nx,px],coefs[0,nx,px,0],coefs[0,nx,px,2]), fontsize='xx-small', transform=axs.transAxes)
                axs.text(0.01,0.7+2*ht,'$p$=%12.10f'%ps[0,nx,px], fontsize='xx-small', transform=axs.transAxes)
                axs.text(0.01,0.7+1*ht,'$R^2_{train}$=%6.4f$\\pm$%6.4f'%(accs[0,nx,px,0,0],accs[0,nx,px,0,2]),fontsize='xx-small', transform=axs.transAxes)
                axs.text(0.01,0.7+0*ht,'$R^2_{test}$=%6.4f$\\pm$%6.4f'%(accs[0,nx,px,1,0],accs[0,nx,px,1,2]),fontsize='xx-small', transform=axs.transAxes)
                
                axs.scatter(neuralactivityconcat[:,nx], movementpcslistconcat[:,px],color='dodgerblue',s=1,alpha=0.5)
                axs.scatter(neuralactivityconcat[:,nx], ys_te[0,:,nx,px],color='orange',s=1,alpha=0.5)


                axs.set_xlim(-0.5,3.5)
                axs.set_ylim(-1,1)
                if px<n_pcs-1: axs.set_xticklabels([])
                if nx>0: axs.set_yticklabels([])



                if px==0: axs.set_title('neuron %d'%(nx+1),fontsize='x-small')
                if nx==0: axs.set_ylabel('PC %d [S.U.]'%(px+1))
                if px==n_pcs-1: axs.set_xlabel('firing rate [S.U.]',fontsize='x-small')


        
        plt.suptitle('%s, neural->behaviour'%calculuslevel)


        save = 1 or globalsave
        if save:
            fig.savefig(resultpath+'neural-%s-single,p,cv-n,b_%s'%(calculuslevel, dn)+ext)

    

        # behaviour -> neural
        n_pcs = 10
        n_neurons = 6
        fig, ax = plt.subplots(n_neurons,n_pcs,figsize=(n_pcs*6,n_neurons*6))
        fig.subplots_adjust(wspace=0.,hspace=0.)

        ht = 0.07 # text height
        for nx in range(n_neurons):
            for px in range(n_pcs):
                axs = ax[nx,px]

                axs.text(0.01,0.7+3*ht,'$\\rho$=%6.4f, w=%6.4f$\\pm$%6.4f'%(rhos[1,nx,px],coefs[1,nx,px,0],coefs[1,nx,px,2]), fontsize='xx-small', transform=axs.transAxes)
                axs.text(0.01,0.7+2*ht,'$p$=%12.10f'%ps[1,nx,px], fontsize='xx-small', transform=axs.transAxes)
                axs.text(0.01,0.7+1*ht,'$R^2_{train}$=%6.4f$\\pm$%6.4f'%(accs[1,nx,px,0,0],accs[1,nx,px,0,2]),fontsize='xx-small', transform=axs.transAxes)
                axs.text(0.01,0.7+0*ht,'$R^2_{test}$=%6.4f$\\pm$%6.4f'%(accs[1,nx,px,1,0],accs[1,nx,px,1,2]),fontsize='xx-small', transform=axs.transAxes)
                
                axs.scatter(movementpcslistconcat[:,px], neuralactivityconcat[:,nx],color='dodgerblue',s=1,alpha=0.5)
                axs.scatter(movementpcslistconcat[:,px], ys_te[1,:,nx,px],color='orange',s=1,alpha=0.5)


                axs.set_xlim(-1,1)
                axs.set_ylim(-0.5,3.5)
                if nx<n_neurons-1: axs.set_xticklabels([])
                if px>0: axs.set_yticklabels([])



                if nx==0: axs.set_title('PC %d'%(px+1),fontsize='x-small')
                if px==0: axs.set_ylabel('neuron %d\nfiring rate [S.U.]'%(nx+1))
                if nx==n_neurons-1: axs.set_xlabel('[S.U.]',fontsize='x-small')


        
        plt.suptitle('%s, behaviour->neural'%calculuslevel)


        save = 1 or globalsave
        if save:
            fig.savefig(resultpath+'neural-%s-single,p,cv-b,n_%s'%(calculuslevel, dn)+ext)
























def decode_movementpca(dn,block,calculuslevel='motion'):
    # create bins of 20 Hz (fps) for neural data
    # use 20 Hz video PC timecourses
    # decode one from the other
    # test shared covariances
    # calculuslevel = {'posture','motion'}   : posture is absolute image intensity, motion is frame to frame intensity difference


    recalculatecorrelations = 0 or globalrecalculate
    recalculateneuraltobehaviourdecode = 0 or globalrecalculate
    recalculatebehaviourtoneuraldecode = 0 or globalrecalculate

    doplot = 0 or globaldoplot

    th = 0.05   # significance threshold


    triallist = preprocess.loadexperimentdata(dn, multimodalonly=True)
    movementpcslist, timestampfps = getmovementpctimecourses(dn, calculuslevel=calculuslevel)          # movementpcslist dimensions are (trials,timecourse,pcs)

    allcomplextrials = [ [[2,4], [], []] ]
    neuralactivity = preprocess.collect_stimulusspecificresponses(block, allcomplextrials)
    downsamplerate = 5   # length * binsize / fps
    neuralactivity = np.array(neph.downsamplesignals(neuralactivity, downsamplerate)[0])

    n_trials = neuralactivity.shape[0]
    n_timecourse = neuralactivity.shape[1]
    n_neurons = neuralactivity.shape[2]
    n_pcs = movementpcslist.shape[2]

    print('downsampled neural activity: ', neuralactivity.shape)
    print('movement PCs: ', movementpcslist.shape)




    # calculate per neuron per pc correlations along the trials
    if recalculatecorrelations and globalrecalculate:
        rhos = np.zeros((n_timecourse,n_pcs,n_neurons))
        ps = np.zeros((n_timecourse,n_pcs,n_neurons))
        for nx in range(n_neurons):
            for px in range(n_pcs):
                # if nx>10 or px>6: break
                for tx in range(n_timecourse):
                    rho,p = sp.stats.linregress(neuralactivity[:,tx,nx], movementpcslist[:,tx,px])[2:4]     # rho and p
                    rhos[tx,px,nx] = rho
                    ps[tx,px,nx] = p
                # print('neuron %d, pc %d: mean(rho_t)=%4.2f, mean(p_t)=%6.5f'%(nx+1, px+1, rhos[:,nx,px].mean(), ps[:,nx,px].mean()) )
        pickle.dump((rhos,ps), open(cacheprefix+'locomotion/%spcs,neurons-correlations,timecourse_%s.pck'%(calculuslevel, dn), 'wb'))
    else:
        rhos,ps = pickle.load(open(cacheprefix+'locomotion/%spcs,neurons-correlations,timecourse_%s.pck'%(calculuslevel, dn), 'rb'))
    print('num significant:', np.sum(ps<=th), '/', n_timecourse*n_neurons*n_pcs)



    reguralization = 'ridge'     # 'ridge' or 'lasso'
    alpha = 500.0


    # decode motion from neurons
    if recalculateneuraltobehaviourdecode and globalrecalculate:
        # accuracies will be (time,pcs,traintest,stats)
        accuraciesnm,coefsnm = nedi.get_linearregressionmultivariate(neuralactivity,movementpcslist,timestampfps,'regression', reguralization=reguralization,alpha=alpha)
        pickle.dump((accuraciesnm,coefsnm), open(cacheprefix+'locomotion/%spcs,neurons-decoder,%s,%d,neuraltobehaviour,timecourse_%s.pck'%(calculuslevel, reguralization,np.round(alpha),dn), 'wb'))
    else:
        accuraciesnm,coefsnm = pickle.load(open(cacheprefix+'locomotion/%spcs,neurons-decoder,%s,%d,neuraltobehaviour,timecourse_%s.pck'%(calculuslevel, reguralization,np.round(alpha),dn), 'rb'))



    # decode neurons from motion (glm)
    if recalculatebehaviourtoneuraldecode and globalrecalculate:
        # accuracies will be (time,neurons,traintest,stats)
        accuraciesmn,coefsmn = nedi.get_linearregressionmultivariate(movementpcslist,neuralactivity,timestampfps,'regression')
        pickle.dump((accuraciesmn, coefsmn), open(cacheprefix+'locomotion/%spcs,neurons-decoder,%s,%d,behaviourtoneural,timecourse_%s.pck'%(calculuslevel, reguralization,np.round(alpha),dn), 'wb'))
    else:
        accuraciesmn,coefsmn = pickle.load(open(cacheprefix+'locomotion/%spcs,neurons-decoder,%s,%d,behaviourtoneural,timecourse_%s.pck'%(calculuslevel, reguralization,np.round(alpha),dn), 'rb'))












    
    if doplot:
        n_pcs = 6
        n_neurons = 10
        fig, ax = plt.subplots(n_pcs,n_neurons,figsize=(n_neurons*12,n_pcs*12))
        fig.subplots_adjust(wspace=0.,hspace=0.)


        for nx in range(n_neurons):
            for px in range(n_pcs):
                axs = ax[px,nx]
                snf = ps[:,px,nx]<=th # significance mask

                axs.plot(timestampfps, rhos[:,px,nx],color='orange',lw=2,alpha=0.8,label='r')
                axs.plot(timestampfps[snf], rhos[snf,px,nx],'.',color='orange',lw=2,alpha=1)
                axs.set_ylim(-1,1)
                if px<n_pcs-1: axs.set_xticklabels([])
                if nx>0: axs.set_yticklabels([])

                axst = axs.twinx()
                axst.semilogy(timestampfps, ps[:,px,nx],color='slategrey', lw=2,alpha=0.8,label='p')
                axst.semilogy(timestampfps[snf], ps[snf,px,nx],'.',color='slategrey', lw=2,alpha=1)
                axst.semilogy([timestampfps[0],timestampfps[-1]],[th,th],'--',lw=1,color='slategrey',alpha=0.7)
                axst.set_ylim(1e-5,1)
                if nx<n_neurons-1: axst.set_yticklabels([])

                # if px==0 and nx==0: axs.legend(frameon=False)
        
        
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'neural-%s-correlations_%s'%(calculuslevel, dn)+ext)     #PCs, all

    




    if doplot:

        colors = ['dodgerblue','orange']
        labels = ['train','test']
        fig, ax = plt.subplots(5,5,figsize=(5*8,5*8))

        for hx in range(5):
            for wx in range(5):
                px = hx*5+wx
                axs = ax[hx,wx]

                for k in [0,1]:
                    m = accuraciesmn[:,px,k,0]
                    e = accuraciesmn[:,px,k,2]
                    axs.plot(timestampfps, m, color=colors[k], lw=3,label=labels[k])
                    axs.fill_between(timestampfps, m-e, m+e, color=colors[k], alpha=0.3)

                if wx==0 and hx==0: axs.legend(frameon=False)
                figs.setxt(axs)
                axs.set_ylim(-2,1)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs)
    

                axs.set_title('PC %d'%(px+1))
                if wx==0: axs.set_ylabel('$R^2$')



        fig.suptitle('%s   neurons -> %s PCs linear regression '%(dn, calculuslevel))

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'neural-%s-decoder,%s,%d,n,m_%s'%(calculuslevel, reguralization,np.round(alpha),dn)+ext)     #PCs, all





    if doplot:
    
        colors = ['dodgerblue','orange']
        labels = ['train','test']
        fig, ax = plt.subplots(5,9,figsize=(9*8,5*8))

        for hx in range(5):
            for wx in range(9):
                nx = hx*9+wx
                if nx>43: break
                axs = ax[hx,wx]
                for k in [0,1]:
                    m = accuraciesmn[:,nx,k,0]
                    e = accuraciesmn[:,nx,k,2]
                    axs.plot(timestampfps, m, color=colors[k], lw=3,label=labels[k])
                    axs.fill_between(timestampfps, m-e, m+e, color=colors[k], alpha=0.3)

                if wx==0 and hx==0: axs.legend(frameon=False)
                figs.setxt(axs)
                axs.set_ylim(-2,1)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs)
    

                axs.set_title('neuron %d'%(nx+1))
                if wx==0: axs.set_ylabel('$R^2$')



        fig.suptitle('%s   %s PCs -> neurons linear regression '%(dn, calculuslevel))

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'neural-%s-glm,%s,%d,m,n_%s'%(calculuslevel, reguralization,np.round(alpha),dn)+ext)     #neurons, all



    return






def decode_movementbodyparts(dn,block):
    # create bins of 20 Hz (fps) for neural data
    # use 20 Hz video body parts mean intensity timecourses
    # decode one from the other
    # test shared covariances
    # calculuslevel = {'posture','motion'}   : posture is absolute image intensity, motion is frame to frame intensity difference


    recalculatecorrelations = 0 or globalrecalculate
    recalculateneuraltobehaviourdecode = 0 or globalrecalculate
    recalculatebehaviourtoneuraldecode = 0 or globalrecalculate

    doplot = 1 or globaldoplot

    th = 0.05   # significance threshold

    triallist = preprocess.loadexperimentdata(dn, multimodalonly=True)
    timestampfps = np.arange( T['starttime'].magnitude, T['endtime'].magnitude, 1000./T['videofps'] )   # milliseconds fps
    movementbodypartslist = preprocess.loadmovingtrialsbodyparts(dn)         # dimensions: (trials,timecourse,bodyparts)
    movementbodypartslist = movementbodypartslist[:,:-1,:]
    bodypartnames = ['nose','mouth','eye','ear','forepaw','back']


    allcomplextrials = [ [[2,4], [], []] ]
    neuralactivity = preprocess.collect_stimulusspecificresponses(block, allcomplextrials)
    downsamplerate = 5   # length * binsize / fps
    neuralactivity = np.array(neph.downsamplesignals(neuralactivity, downsamplerate)[0])

    n_trials = neuralactivity.shape[0]
    n_timecourse = neuralactivity.shape[1]
    n_neurons = neuralactivity.shape[2]
    n_bodyparts = movementbodypartslist.shape[2]

    print('downsampled neural activity: ', neuralactivity.shape)
    print('movement body  partss: ', movementbodypartslist.shape)




    # calculate per neuron per pc correlations along the trials
    if recalculatecorrelations and globalrecalculate:
        rhos = np.zeros((n_timecourse,n_bodyparts,n_neurons))
        ps = np.zeros((n_timecourse,n_bodyparts,n_neurons))
        for nx in range(n_neurons):
            for bx in range(n_bodyparts):
                # if nx>10 or px>6: break
                for tx in range(n_timecourse):
                    rho,p = sp.stats.linregress(neuralactivity[:,tx,nx], movementbodypartslist[:,tx,bx])[2:4]     # rho and p
                    rhos[tx,bx,nx] = rho
                    ps[tx,bx,nx] = p
                # print('neuron %d, pc %d: mean(rho_t)=%4.2f, mean(p_t)=%6.5f'%(nx+1, px+1, rhos[:,nx,px].mean(), ps[:,nx,px].mean()) )
        pickle.dump((rhos,ps), open(cacheprefix+'locomotion/motionbodyparts,neurons-correlations,timecourse_%s.pck'%(dn), 'wb'))
    else:
        rhos,ps = pickle.load(open(cacheprefix+'locomotion/motionbodyparts,neurons-correlations,timecourse_%s.pck'%(dn), 'rb'))
    print('num significant:', np.sum(ps<=th), '/', n_timecourse*n_neurons*n_bodyparts)



    reguralization = 'ridge'     # 'ridge' or 'lasso'
    alpha = 1.0


    # decode motion from neurons
    if recalculateneuraltobehaviourdecode and globalrecalculate:
        # accuracies will be (time,pcs,traintest,stats)
        accuraciesnm,coefsnm = nedi.get_linearregressionmultivariate(neuralactivity,movementbodypartslist,timestampfps,'regression', reguralization=reguralization,alpha=alpha)
        pickle.dump((accuraciesnm,coefsnm), open(cacheprefix+'locomotion/motionbodyparts,neurons-decoder,%s,%d,neuraltobehaviour,timecourse_%s.pck'%(reguralization,np.round(alpha),dn), 'wb'))
    else:
        accuraciesnm,coefsnm = pickle.load(open(cacheprefix+'locomotion/motionbodyparts,neurons-decoder,%s,%d,neuraltobehaviour,timecourse_%s.pck'%(reguralization,np.round(alpha),dn), 'rb'))



    # decode neurons from motion (glm)
    if recalculatebehaviourtoneuraldecode and globalrecalculate:
        # accuracies will be (time,neurons,traintest,stats)
        accuraciesmn,coefsmn = nedi.get_linearregressionmultivariate(movementbodypartslist,neuralactivity,timestampfps,'regression')
        pickle.dump((accuraciesmn, coefsmn), open(cacheprefix+'locomotion/motionbodyparts,neurons-decoder,%s,%d,behaviourtoneural,timecourse_%s.pck'%(reguralization,np.round(alpha),dn), 'wb'))
    else:
        accuraciesmn,coefsmn = pickle.load(open(cacheprefix+'locomotion/motionbodyparts,neurons-decoder,%s,%d,behaviourtoneural,timecourse_%s.pck'%(reguralization,np.round(alpha),dn), 'rb'))












    
    if doplot:
        n_bodyparts = 6
        n_neurons = 10
        fig, ax = plt.subplots(n_bodyparts,n_neurons,figsize=(n_neurons*12,n_bodyparts*12))
        fig.subplots_adjust(wspace=0.,hspace=0.)


        for nx in range(n_neurons):
            for bx in range(n_bodyparts):
                axs = ax[bx,nx]
                snf = ps[:,bx,nx]<=th # significance mask

                axs.plot(timestampfps, rhos[:,bx,nx],color='orange',lw=2,alpha=0.8,label='r')
                axs.plot(timestampfps[snf], rhos[snf,bx,nx],'.',color='orange',lw=2,alpha=1)
                axs.set_ylim(-1,1)
                if bx<n_bodyparts-1: axs.set_xticklabels([])
                if nx>0: axs.set_yticklabels([])

                axst = axs.twinx()
                axst.semilogy(timestampfps, ps[:,bx,nx],color='slategrey', lw=2,alpha=0.8,label='p')
                axst.semilogy(timestampfps[snf], ps[snf,bx,nx],'.',color='slategrey', lw=2,alpha=1)
                axst.semilogy([timestampfps[0],timestampfps[-1]],[th,th],'--',lw=1,color='slategrey',alpha=0.7)
                axst.set_ylim(1e-5,1)
                if nx<n_neurons-1: axst.set_yticklabels([])

                # if bx==0 and nx==0: axs.legend(frameon=False)
        
        
        save = 1 or globalsave
        if save:
            fig.savefig(resultpath+'neural-motion,bodyparts-correlations_%s'%(dn)+ext)     #PCs, all

    


    print(accuraciesmn.shape, accuraciesnm.shape)

    if doplot:

        colors = ['dodgerblue','orange']
        labels = ['train','test']
        fig, ax = plt.subplots(2,3,figsize=(3*8,2*8))

        for hx in range(2):
            for wx in range(3):
                bx = hx*3+wx
                axs = ax[hx,wx]

                for k in [0,1]:
                    m = accuraciesnm[:,bx,k,0]
                    e = accuraciesnm[:,bx,k,2]
                    axs.plot(timestampfps, m, color=colors[k], lw=3,label=labels[k])
                    axs.fill_between(timestampfps, m-e, m+e, color=colors[k], alpha=0.3)

                if wx==0 and hx==0: axs.legend(frameon=False)
                figs.setxt(axs)
                axs.set_ylim(-2,1)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs)
    

                axs.set_title(bodypartnames[bx])
                if wx==0: axs.set_ylabel('$R^2$')



        fig.suptitle('%s   neurons -> motion body parts linear regression '%(dn))

        save = 1 or globalsave
        if save:
            fig.savefig(resultpath+'neural-motion,bodyparts-decoder,%s,%d,n,m_%s'%(reguralization,np.round(alpha),dn)+ext)     #PCs, all





    if doplot:
    
        colors = ['dodgerblue','orange']
        labels = ['train','test']
        fig, ax = plt.subplots(5,9,figsize=(9*8,5*8))

        for hx in range(5):
            for wx in range(9):
                nx = hx*9+wx
                if nx>43: break
                axs = ax[hx,wx]
                for k in [0,1]:
                    m = accuraciesmn[:,nx,k,0]
                    e = accuraciesmn[:,nx,k,2]
                    axs.plot(timestampfps, m, color=colors[k], lw=3,label=labels[k])
                    axs.fill_between(timestampfps, m-e, m+e, color=colors[k], alpha=0.3)

                if wx==0 and hx==0: axs.legend(frameon=False)
                figs.setxt(axs)
                axs.set_ylim(-2,1)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs)
    

                axs.set_title('neuron %d'%(nx+1))
                if wx==0: axs.set_ylabel('$R^2$')



        fig.suptitle('%s   motion bodyparts -> neurons linear regression '%(dn))

        save = 1 or globalsave
        if save:
            fig.savefig(resultpath+'neural-motion,bodyparts-glm,%s,%d,m,n_%s'%(reguralization,np.round(alpha),dn)+ext)     #neurons, all



    return












def decode_movementpca_tasks(dn):
    # use 20 Hz video PC timecourses
    # decode task variables (visual,audio,context,choice) from it

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot

    # setup stimulus:
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [   [ [ [], [45],    [] ],        [    [],[135],     [] ]  ], \
                            [ [ [], [],    [5000] ],        [    [],[],     [10000] ]  ], \
                            [ [ [blv[1]],  [],[] ], [ [bla[1]],   [],[] ]    ],     \
                            [ [True], [False] ], \
                            #[ [ [bla[0]],  [],[5000] ], [ [bla[0]],   [],[10000] ]    ],     \
                            [ [], [] ],     \
                        ]
    taskaspects = ['visual','audio','context','choice','shuffled']   #'simpleaudio',
    columnnames = ['degree','freq','block','action','']    # 'freq',
    class1value = [45,5000,blv[1],True,0]    # 5000

    triallist = preprocess.loadexperimentdata(dn, multimodalonly=True)
    triallist['block'] += 1
    movementpcslist, timestampfps = getmovementpctimecourses(dn)          # movementpcslist dimensions are (trials,timecourse,pcs)
    # excludepc_indices = np.array([])
    # includepc_indices = np.array([ px for px in range(movementpcslist.shape[2]) if px not in excludepc_indices ])
    # movementpcslist = movementpcslist[:,:,includepc_indices]
    n_trials,n_timestamps,n_features = movementpcslist.shape
    n_targets = 1



    accuracylist = []
    coeflist = []
    if recalculate:
        accuracies = np.zeros((n_timestamps,n_targets,2,3,2,n_features,len(taskaspects)))
        coefs = np.zeros((n_timestamps,n_targets,n_features,3,2,n_features,len(taskaspects)))
        for cx,task in enumerate(taskaspects):
            print(task)
            comparisonmask = []
            if cx<4:
                for comparison in comparisongroups[cx]:
                    if task != 'choice':
                        comparisonmask.append( \
                        ( (triallist['block'].isin(comparison[0])) | (len(comparison[0])==0) ) & \
                        ( (triallist['degree'].isin(comparison[1])) | (len(comparison[1])==0) ) & \
                        ( (triallist['freq'].isin(comparison[2])) | (len(comparison[2])==0) ) )
                    else:
                        comparisonmask.append( triallist['action'].isin(comparison) )

                targets = np.hstack([ triallist[columnnames[cx]][comparisonmask[i]]==class1value[cx] for i in [0,1]  ])
                predictors = np.vstack([movementpcslist[comparisonmask[i],:,:] for i in [0,1]])
            # elif cx==4: # simple audio in single modality audio; this only works once single trial PCAs are calculated as well
            #     triallist_alltrials = preprocess.loadexperimentdata(dn, multimodalonly=False)
            #     triallist_alltrials['block'] += 1
            #     mask_simpleaudio = triallist_alltrials['block']==bla[0]
            #     movementpcslist_alltrials, _ = getmovementpctimecourses(dn, multimodalonly=False)
            #     predictors = movementpcslist_alltrials[mask_simpleaudio,:,:]
            #     targets = triallist_alltrials[mask_simpleaudio][columnnames[cx]]==class1value[cx]
            elif cx==4: # shuffle
                L = movementpcslist.shape[0]
                targets = np.hstack([ np.ones(L//2), np.zeros(L//2) ])[np.random.permutation(L)]
                predictors = movementpcslist

            # fill out the targets axes, with dimensions = (observations,timepoints,features)
            targets = np.tile(targets, [1,len(timestampfps),1])
            targets = np.swapaxes(targets, 0, 2)

            # for px in np.arange(n_features-1,n_features-3,-1):
            # for px in np.arange(3-1,-1,-1):
            for px in np.arange(n_features):
                pcslice = slice(px,n_features)   # cut the last px number of PCs
                    
                # only use one pc from back to front
                accuracyo,coefo = nedi.get_linearregressionmultivariate(predictors[:,:,px][:,:,np.newaxis],targets,timestampfps,'classification')
                # cumulatively use more and more pcs from back to front
                if px==n_features-1:
                    accuracyc,coefc = accuracyo,coefo
                else:
                    accuracyc,coefc = nedi.get_linearregressionmultivariate(predictors[:,:,pcslice],targets,timestampfps,'classification')

                # (n_timestamps,n_targets,train/test,stats,one/cumulative,pcused,task)
                accuracies[:,:,:,:,0,px,cx] = accuracyo
                accuracies[:,:,:,:,1,px,cx] = accuracyc
                # (n_timestamps,n_targets,n_features,stats,one/cumulative,pcused,task)
                coefs[:,:,px:n_features,:,0,px,cx] = coefo
                coefs[:,:,px:n_features,:,1,px,cx] = coefc

            

        pickle.dump((accuracies,coefs), open(cacheprefix+'locomotion/motionpcs,variables-decreasing-decoder,timecourse_%s.pck'%(dn), 'wb'))
    else:
        accuracies,coefs = pickle.load(open(cacheprefix+'locomotion/motionpcs,variables-decreasing-decoder,timecourse_%s.pck'%(dn), 'rb'))

    


    print(accuracies.shape)
    print(coefs.shape)

    if doplot:
    
        n_displaypc = 6
        labels = ['train','test']
        taskcolors = ['navy','darkgreen','mediumvioletred','darkorange','green','red']
        n_tasks = len(taskaspects)
        
        fig, ax = plt.subplots(1+n_displaypc,n_tasks,figsize=(n_tasks*12,(1+n_displaypc)*8))

        for cx,task in enumerate(taskaspects):
            colors = ['grey',taskcolors[cx]]
            axs = ax[0,cx]
            for k in [0,1]:
                m = accuracies[:,0,k,0,1,0,cx]         # just one dim target on 2nd dimension
                e = accuracies[:,0,k,2,1,0,cx]
                axs.plot(timestampfps, m, color=colors[k], lw=3)
                axs.fill_between(timestampfps, m-e, m+e, color=colors[k], alpha=0.3)

            figs.setxt(axs)
            axs.set_ylim(0.45,1.05)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)


            axs.set_title(task)
            if cx==0: axs.set_ylabel('accuracy')



            for nx in range(min(n_features,n_displaypc)):
                axs = ax[1+nx,cx]
                m = coefs[:,0,nx,0,1,0,cx]
                e = coefs[:,0,nx,2,1,0,cx]
                axs.plot(timestampfps, m, color=taskcolors[cx])
                axs.fill_between(timestampfps, m-e, m+e, color=taskcolors[cx], alpha=0.3)
                figs.setxt(axs)
                axs.set_ylim(-1.5,1.5)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs)
                if cx==0: axs.set_ylabel('weight PC %d'%(nx+1))


        fig.suptitle('%s   diff. motion PCs -> task variables logistic regression '%dn)

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'motionpcs,variables-decoder,timecourse_%s'%(dn)+ext)     #neurons, all




        # show the 
        fig, ax = plt.subplots(6,7,figsize=(7*12,6*8))
        cx = 2   # context only
        for hx in range(7):
            for wx in range(7):
                px = hx*7+wx
                if px>=n_features: break


                axs = ax[hx,wx]
                for k in [0,1]:
                    for ox in [1,0]:   # only and up until
                        colors = [['grey','darkgrey'][ox],['rebeccapurple','mediumvioletred'][ox]]
                        labels = ['only PC %d'%(px+1), 'PC %d-%d'%(px+1,n_features)]
                        m = accuracies[:,0,k,0,ox,px,cx]         # just one dim target on 2nd dimension
                        e = accuracies[:,0,k,2,ox,px,cx]
                        axs.plot(timestampfps, m, color=colors[k], lw=3, label=labels[ox]+[' train',' test'][k])
                        axs.fill_between(timestampfps, m-e, m+e, color=colors[k], alpha=0.3)

                figs.setxt(axs)
                axs.set_ylim(0.45,1.05)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.5)

                axs.legend(frameon=False)
                if cx==0: axs.set_ylabel('accuracy')




        fig.suptitle('%s   diff. motion reverse cumulative and single PCs -> context logistic regression '%dn)

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'motionpcs,context,multipcs-decoder,timecourse_%s'%(dn)+ext)     #neurons, all













def subspace_decode_motiontoneuron_reducedrank(dn, block, calculuslevel='motion'):
    # create bins of 20 Hz (fps) for neural data
    # use 20 Hz video PC timecourses
    # reduced rank linear regression of neurons from motion
    # test for saturating rank

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot

    

    movementpcslist, timestampfps = getmovementpctimecourses(dn, calculuslevel=calculuslevel)          # movementpcslist dimensions are (trials,timecourse,pcs)


    allcomplextrials = [ [[2,4], [], []] ]
    neuralactivity = preprocess.collect_stimulusspecificresponses(block, allcomplextrials)
    downsamplerate = 5   # length * binsize / fps
    neuralactivity = np.array(neph.downsamplesignals(neuralactivity, downsamplerate)[0])

    n_trials = neuralactivity.shape[0]
    n_timecourse = neuralactivity.shape[1]
    n_neurons = neuralactivity.shape[2]
    n_pcs = movementpcslist.shape[2]
    maxrank = min(n_neurons, n_pcs)

    print('downsampled neural activity: ', neuralactivity.shape)
    print('movement PCs: ', movementpcslist.shape)


    # trial timepoints concatenated, keep a singleton time dimension
    movementpcslistconcat = np.expand_dims(np.reshape(movementpcslist, (-1,n_pcs)), axis=1)
    neuralactivityconcat = np.expand_dims(np.reshape(neuralactivity, (-1,n_neurons )), axis=1)


    
    # calculate in 1 sec bins
    broadsamplerate = 20
    neuralactivitylowres = neph.downsamplearray(neuralactivity, broadsamplerate, axis=1)    # gather spike counts in 1 sec bins
    movementpcslistlowres = neph.downsamplearray(movementpcslist, broadsamplerate, axis=1)    # average movement PCs in 1 sec bins
    n_timecourse_lowres = neuralactivitylowres.shape[1]
    timestampfpslowres = timestampfps[::broadsamplerate]
    timestampfpslowres += ((timestampfpslowres[1] - timestampfpslowres[0])//2)
    
    
    
    # decode neurons from motion (glm)
    if recalculate:
        # accuracies will be (time,traintest,stats,ranks)
        # coefs will be (time,n_targets,n_predictors,stats,ranks)

        # motion -> neural
        # timecourse
        accuraciesmnranks = np.zeros((n_timecourse,2,3,maxrank))
        coefsmnranks = np.zeros((n_timecourse,n_neurons,n_pcs,3,maxrank))
        # concatenate
        accuraciesmnranksconcat = np.zeros((1,2,3,maxrank))
        coefsmnranksconcat = np.zeros((1,n_neurons,n_pcs,3,maxrank))
        # low resolution
        accuraciesmnrankslowres = np.zeros((n_timecourse_lowres,2,3,maxrank))
        coefsmnrankslowres = np.zeros((n_timecourse_lowres,n_neurons,n_pcs,3,maxrank))


        # neural -> motion
        # timecourse
        accuraciesnmranks = np.zeros((n_timecourse,2,3,maxrank))
        coefsnmranks = np.zeros((n_timecourse,n_pcs,n_neurons,3,maxrank))
        # concatenate
        accuraciesnmranksconcat = np.zeros((1,2,3,maxrank))
        coefsnmranksconcat = np.zeros((1,n_pcs,n_neurons,3,maxrank))
        # low resolution
        accuraciesnmrankslowres = np.zeros((n_timecourse_lowres,2,3,maxrank))
        coefsnmrankslowres = np.zeros((n_timecourse_lowres,n_pcs,n_neurons,3,maxrank))


        for r in range(maxrank):
            print('rank %d, '%r)

            # timecourse
            accuraciesmn,coefsmn = nedi.get_linearregressionmultivariate(movementpcslist,neuralactivity,timestampfps,'reducedrankregression', rank=r)
            accuraciesmnranks[:,:,:,r] = accuraciesmn[:,0,:,:]
            coefsmnranks[:,:,:,:,r] = coefsmn
            accuraciesnm,coefsnm = nedi.get_linearregressionmultivariate(neuralactivity,movementpcslist,timestampfps,'reducedrankregression', rank=r)
            accuraciesnmranks[:,:,:,r] = accuraciesnm[:,0,:,:]
            coefsnmranks[:,:,:,:,r] = coefsnm
            accuraciesmnlowres,coefsmnlowres = nedi.get_linearregressionmultivariate(movementpcslistlowres,neuralactivitylowres,timestampfpslowres,'reducedrankregression', rank=r)
            accuraciesmnrankslowres[:,:,:,r] = accuraciesmnlowres[:,0,:,:]
            coefsmnrankslowres[:,:,:,:,r] = coefsmnlowres
            accuraciesnmlowres,coefsnmlowres = nedi.get_linearregressionmultivariate(neuralactivitylowres,movementpcslistlowres,timestampfpslowres,'reducedrankregression', rank=r)
            accuraciesnmrankslowres[:,:,:,r] = accuraciesnmlowres[:,0,:,:]
            coefsnmrankslowres[:,:,:,:,r] = coefsnmlowres


            # concatenated
            accuraciesmnconcat,coefsmnconcat = nedi.get_linearregressionmultivariate(movementpcslistconcat,neuralactivityconcat,timestampfps,'reducedrankregression', rank=r)
            accuraciesmnranksconcat[:,:,:,r] = accuraciesmnconcat[:,0,:,:]
            coefsmnranksconcat[:,:,:,:,r] = coefsmnconcat
            accuraciesnmconcat,coefsnmconcat = nedi.get_linearregressionmultivariate(neuralactivityconcat,movementpcslistconcat,timestampfps,'reducedrankregression', rank=r)
            accuraciesnmranksconcat[:,:,:,r] = accuraciesnmconcat[:,0,:,:]
            coefsnmranksconcat[:,:,:,:,r] = coefsnmconcat


        pickle.dump((accuraciesmnranks, coefsmnranks, accuraciesnmranks, coefsnmranks, accuraciesmnrankslowres, coefsmnrankslowres),\
                open(cacheprefix+'locomotion/%spcs,neurons-rrr,behaviourtoneural,timecourse_%s.pck'%(calculuslevel, dn), 'wb'))
        pickle.dump((accuraciesmnranksconcat, coefsmnranksconcat, accuraciesnmranksconcat, coefsnmranksconcat, accuraciesnmrankslowres, coefsnmrankslowres),\
                open(cacheprefix+'locomotion/%spcs,neurons-rrr,behaviourtoneural,concatenated_%s.pck'%(calculuslevel, dn), 'wb'))
    else:
        accuraciesmnranks, coefsmnranks, accuraciesnmranks, coefsnmranks, accuraciesmnrankslowres, coefsmnrankslowres = pickle.load(open(cacheprefix+'locomotion/%spcs,neurons-rrr,behaviourtoneural,timecourse_%s.pck'%(calculuslevel, dn), 'rb'))
        accuraciesmnranksconcat, coefsmnranksconcat, accuraciesnmranksconcat, coefsnmranksconcat, accuraciesnmrankslowres, coefsnmrankslowres = pickle.load(open(cacheprefix+'locomotion/%spcs,neurons-rrr,behaviourtoneural,concatenated_%s.pck'%(calculuslevel, dn), 'rb'))



    # plot

    if doplot:
        colors = ['dodgerblue','orange']
        colorlist = plt.cm.viridis( np.linspace(1.0, 0.33, maxrank) )
        labels = ['train','test']
        titles = ['%s->neural','neural->%s']%(calculuslevel,calculuslevel)
        concattitles = [' timecourse',' concatenated','lowres']
        print('plotting ',calculuslevel     )
        fig, ax = plt.subplots(2,3,figsize=(3*12,2*8))


        for ct,(accuraciesranks,accuraciesranksconcat,accuraciesrankslowres) in enumerate(\
                  zip((accuraciesmnranks,accuraciesnmranks),\
                      (accuraciesmnranksconcat,accuraciesnmranksconcat),\
                      (accuraciesmnrankslowres,accuraciesnmrankslowres))):
            
            
            axs = ax[ct,0]
            for r in range(maxrank):
                colors = ['grey',colorlist[r]]
                labels = [None, 'rank %d'%r]
                for k in [0,1]:
                    m = accuraciesranks[:,k,0,r]
                    e = accuraciesranks[:,k,2,r]
                    axs.plot(timestampfps, m, color=colors[1], lw=1, label=labels[k])
                    axs.fill_between(timestampfps, m-e, m+e, color=colors[k], alpha=0.3)

            axs.legend(frameon=False, loc='upper left', ncol=4, fontsize=8)
            figs.setxt(axs)
            axs.set_ylim(-4,1)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs)
            axs.set_title(titles[ct]+concattitles[ct])


            axs.set_ylabel('$R^2$')




            axs = ax[ct,1]

            colors = ['grey','rebeccapurple']
            labels = [None, 'rank %d'%r]
            ranklist = np.arange(maxrank)+1
            for k in [0,1]:
                m = accuraciesranksconcat[0,k,0,:]
                e = accuraciesranksconcat[0,k,2,:]
                axs.plot(ranklist, m, color=colors[k], lw=1, label=labels[k])
                axs.fill_between(ranklist, m-e, m+e, color=colors[k], alpha=0.3)


            axs.set_ylim(-0.2,+0.2)
            axs.set_title(titles[ct]+concattitles[ct])






            axs = ax[ct,2]
            
            for r in range(maxrank):
                colors = ['grey',colorlist[r]]
                labels = [None, 'rank %d'%r]
                for k in [0,1]:
                    m = accuraciesrankslowres[:,k,0,r]
                    e = accuraciesrankslowres[:,k,2,r]
                    axs.plot(timestampfpslowres, m, color=colors[1], lw=1, label=labels[k])
                    axs.fill_between(timestampfpslowres, m-e, m+e, color=colors[k], alpha=0.3)

            axs.legend(frameon=False, loc='upper left', ncol=4, fontsize=8)
            figs.setxt(axs)
            axs.set_ylim(-4,1)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs)
            axs.set_title(titles[ct]+concattitles[ct])


            axs.set_ylabel('$R^2$')



        figs.plottoaxis_chancelevel(axs)



        fig.suptitle('%s   %s PCs -> neurons   linear reduced rank regressions '%(calculuslevel, dn))

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'neural-%s-rrr,bn,nb,timecourse+concat_%s'%(calculuslevel, dn)+ext)     #neurons, all






    return











def comparestationarycontext(dn, block):
    # separate motion and stationary by video motion energy thresholds
    # compare context representation to all trials

    print('decoding context at various motion threshold selected trials')

    recalculate = 0 or globalrecalculate
    doplot = 1 or globaldoplot


    # collect motion criteria at different thresholds and allowed proportions
    thresholds, proportions, movingtrials = preprocess.loadmovingthresholdtrials(dn)
    n_thresholds = len(thresholds)
    n_proportions = len(proportions)
    n_trials = movingtrials.shape[0]//len(proportions)

    movingtrials = np.reshape(movingtrials,(n_trials,n_proportions,-1))
    stationarytrials = (1-movingtrials).astype(np.bool8)
    # print(movingtrials.shape)
    # print(n_trials-movingtrials.sum(axis=0))


    # collect neural activity
    allcomplextrials = [ [[2,4], [], []] ]
    targets_all = np.hstack((2*np.ones(70),4*np.ones(70)))     # reestablish block ids
    neuralactivity_all = np.array(preprocess.collect_stimulusspecificresponses(block, allcomplextrials)[0])
    wx = int((50*pq.ms/T['dt']).magnitude)  # width of the feature space in time for each neuron
    times = block.segments[0].analogsignals[0].times
    if wx>1:   # create wider feature space of several subsequent timepoints
        times = times[:-wx]
        neuralactivity_all = neph.slidingwindow(neuralactivity_all,wx)


    n_neurons = block.segments[0].analogsignals[0].shape[1]
    n_timestamps = len(times)

    n_trials_contextminimum = 5 # the number of trials required per trial as a minimum to do a decoding

    # calculate the decoder grids
    if recalculate:
        accuracies = np.nan*np.ones((n_timestamps,2,3,n_proportions,n_thresholds))        # (timepoints, train/test, stats, proportions, thresholds)
        coefs = np.nan*np.ones((n_timestamps,n_neurons,3,n_proportions,n_thresholds))     # (timepoints, neurons, stats, proportions, thresholds)

        for px,pr in enumerate(proportions):
            print('at proportions', pr)
            for mx,th in enumerate(thresholds):
                # check if we have enough trials
                if stationarytrials[:,px,mx][:n_trials//2].sum()<n_trials_contextminimum or  \
                   stationarytrials[:,px,mx][n_trials//2:].sum()<n_trials_contextminimum: continue

                # fill out the targets axes, with dimensions = (observations,timepoints,features)
                targets = targets_all[stationarytrials[:,px,mx]]
                predictors = neuralactivity_all[stationarytrials[:,px,mx],:,:]

                targets = np.tile(targets, [1,n_timestamps,1])
                targets = np.swapaxes(targets, 0, 2)

                accuracy,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'classification')
                
                # reshape to average over feature timewidth
                if wx>1:
                    coef = np.reshape(coef, [coef.shape[0], 1, wx, n_neurons, 3]).mean(axis=2)

                # (n_timestamps,train/test,stats,task,symmetry)
                accuracies[:,:,:,px,mx] = accuracy.squeeze()
                # (n_timestamps,n_features,stats,task,symmetry)
                coefs[:,:,:,px,mx] = coef.squeeze()

    
        # do for all trials as control    
        targets = targets_all
        predictors = neuralactivity_all
        targets = np.tile(targets, [1,n_timestamps,1])
        targets = np.swapaxes(targets, 0, 2)
        accuracy,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'classification')
        if wx>1: coef = np.reshape(coef, [coef.shape[0], 1, wx, n_neurons, 3]).mean(axis=2)
        accuracyall = accuracy.squeeze()
        coefsall = coef.squeeze()
            

        pickle.dump((accuracies,coefs,accuracyall,coefsall), open(cacheprefix+'locomotion/stationarytrials,threshold-context-decoder,timecourse_%s.pck'%(dn), 'wb'))
    else:
        accuracies,coefs,accuracyall,coefsall = pickle.load(open(cacheprefix+'locomotion/stationarytrials,threshold-context-decoder,timecourse_%s.pck'%(dn), 'rb'))

    
    
    
    
    
    
    if doplot:
        # show the 
        fig, ax = plt.subplots(1,2,figsize=(2*12,1*8))

        axs = ax[0]
        prop_disp,th_disp = 1,3     # display indices for detailed timecourse
        for mx in range(2):
            for k in [0,1]:
                colors = [['darkgrey','lightgrey'][mx],['rebeccapurple','mediumvioletred'][mx]]
                labels = ['stationary','all']
                if mx==0:
                    m = accuracies[:,k,0,prop_disp,th_disp]         # just one dim target on 2nd dimension
                    e = accuracies[:,k,2,prop_disp,th_disp]
                else:
                    m = accuracyall[:,k,0]
                    e = accuracyall[:,k,2]
                axs.plot(times, m, color=colors[k], lw=3, label=labels[mx]+[' train',' test'][k])
                axs.fill_between(times, m-e, m+e, color=colors[k], alpha=0.3)

            figs.setxt(axs)
            axs.set_ylim(0.45,1.05)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)

            axs.legend(frameon=False)
            axs.set_ylabel('accuracy')



        axs = ax[1]
        colors = plt.cm.viridis( np.linspace(0, 0.8, n_proportions) )
        for px,pr in enumerate(proportions):

            labels = ['stationary','all']
            m = accuracies[:,1,0,px,:].mean(axis=0)         # just one dim target on 2nd dimension
            e = accuracies[:,1,2,px,:].mean(axis=0)
            axs.plot(thresholds, m, color=colors[px], lw=1, label='prop=%4.2f'%proportions[px])
            axs.fill_between(thresholds, m-e, m+e, color=colors[px], alpha=0.3)


        # show the point in the grid on the first subplot
        axs.plot(thresholds[th_disp],accuracies[:,1,0,prop_disp,th_disp].mean(axis=0),'o',color='rebeccapurple')

        axs.set_ylim(0.45,1.05)
        figs.plottoaxis_chancelevel(axs,0.5)

        axs.set_xlabel('motion energy threshold [AU]')
        axs.legend(frameon=False)
        axs.set_ylabel('accuracy')
        
    


        fig.suptitle('%s   decoding context from threshold-stationary trials '%dn)

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'stationarytrials,threshold-context-decoder,timecourse_%s'%(dn)+ext)     #neurons, all

    










def comparemovementdistributioncontext(dn, block, roi='total'):

    # separate motion and stationary by video motion energy thresholds
    # compare context representation to all trials

    print('decoding context at various motion intensity levels')

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot


    # collect motion criteria at different thresholds and allowed proportions

    if roi=='total':
        loader = preprocess.loadmovinglevelstrials
        bodypartnames = ['total']
    elif roi=='bodyparts':
        loader = preprocess.loadmovinglevelstrialsbodyparts
        bodypartnames = ['nose','mouth','eye','ear','forepaw','back']

    levels, movingtrials = loader(dn, bodypartnames)
    print(levels)

    # n_bodyparts = len(bodypartnames)
    n_levels = len(levels)
    n_trials, n_movement_timestamps = movingtrials.shape[:2]
    if len(movingtrials.shape)>2: n_bodyparts = movingtrials.shape[2]
    else: n_bodyparts = 1; movingtrials = movingtrials[:,:,np.newaxis]
    downsamplerate = 5   # length * binsize / fps; video is downsampled compared to neural activity


    # collect neural activity
    allcomplextrials = [ [[2,4], [], []] ]
    neuralactivity_all = np.array(preprocess.collect_stimulusspecificresponses(block, allcomplextrials)[0])
    times = block.segments[0].analogsignals[0].times
    wx = int((50*pq.ms/T['dt']).magnitude)  # width of the feature space in time for each neuron
    if wx>1:   # create wider feature space of several subsequent timepoints
        times = times[:-wx]
        neuralactivity_all = neph.slidingwindow(neuralactivity_all,wx)
    indices = np.arange(len(times))
    print(neuralactivity_all.shape, movingtrials.shape, levels.shape)

    targets_all = np.hstack((2*np.ones(n_trials//2),4*np.ones(n_trials//2)))     # reestablish block ids


    n_neurons = block.segments[0].analogsignals[0].shape[1]
    n_timestamps = len(times)

    n_trials_contextminimum = 5 # the number of trials required per trial as a minimum to do a decoding

    # calculate the decoder grids
    if recalculate:
        accuracies = np.nan*np.ones((n_timestamps,2,3,n_levels,n_bodyparts))        # (timepoints, train/test, stats, levels)
        coefs = np.nan*np.ones((n_timestamps,n_neurons,3,n_levels,n_bodyparts))     # (timepoints, neurons, stats, levels)

        accuraciesrandom = np.nan*np.ones((n_timestamps,2,3,n_levels,n_bodyparts))        # (timepoints, train/test, stats, levels)
        coefsrandom = np.nan*np.ones((n_timestamps,n_neurons,3,n_levels,n_bodyparts))     # (timepoints, neurons, stats, levels)

        # go over body parts
        for bx,bn in enumerate(bodypartnames):
            print('Calculating:',bn)

            # go over timestamps
            for bt,nt in enumerate(indices[:-2:downsamplerate]): # behaviour time, neural time


                for mx,lv in enumerate(levels):      # go over the midpoints

                    mask_levels = movingtrials[:,bt,bx]==mx
                    # check if we have levels at all:
                    if np.sum(mask_levels)==0: continue


                    num_vi = mask_levels[:n_trials//2].sum()
                    num_au = mask_levels[n_trials//2:].sum()
                    # check if we have enough trials in both contexts
                    if num_vi<n_trials_contextminimum or  \
                        num_au<n_trials_contextminimum: continue
                    
                    # equalize the number of trials in both contexts
                    di = num_vi-num_au
                    if di>0:
                        inds = np.argwhere(mask_levels[:n_trials//2])
                        mask_levels[inds[0:np.abs(di)]] = False
                    elif di<0:
                        aux = np.copy(mask_levels)
                        aux[:n_trials//2] = False
                        inds = np.argwhere(aux)         # this will result in the second half checking only
                        mask_levels[inds[0:np.abs(di)]] = False
                    
                    num_vi = mask_levels[:n_trials//2].sum()
                    num_au = mask_levels[n_trials//2:].sum()
                    print(nt,bt,mx,lv, 'hit', num_vi, num_au)

                    # fill out the targets axes, with dimensions = (observations,timepoints,features)
                    predictors = neuralactivity_all[mask_levels,nt:nt+downsamplerate,:]
                    targets = targets_all[mask_levels]

                    targets = np.tile(targets, [1,downsamplerate,1]) # extend dimensions
                    targets = np.swapaxes(targets, 0, 2)

                    accuracy,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'classification')
                    
                    # reshape to average over feature timewidth
                    if wx>1:
                        coef = np.reshape(coef, [coef.shape[0], 1, wx, n_neurons, 3]).mean(axis=2)

                    # (n_timestamps,train/test,stats,task,symmetry)
                    accuracies[nt:nt+downsamplerate,:,:,mx,bx] = accuracy.squeeze()
                    # (n_timestamps,n_features,stats,task,symmetry)
                    coefs[nt:nt+downsamplerate,:,:,mx,bx] = coef.squeeze()


                    # random control
                    chosen_indices_vi = np.random.permutation(n_trials//2)[:num_vi]
                    chosen_indices_au = np.random.permutation(n_trials//2)[:num_au]  + n_trials//2
                    chosen_indices = np.hstack((chosen_indices_vi,chosen_indices_au))
                    predictors = neuralactivity_all[chosen_indices,nt:nt+downsamplerate,:]
                    targets = targets_all[chosen_indices]
                    
                    targets = np.tile(targets, [1,downsamplerate,1]) # extend dimensions
                    targets = np.swapaxes(targets, 0, 2)

                    accuracy_,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'classification')
                    
                    # reshape to average over feature timewidth
                    if wx>1:
                        coef = np.reshape(coef, [coef.shape[0], 1, wx, n_neurons, 3]).mean(axis=2)

                    # (n_timestamps,train/test,stats,task,symmetry)
                    accuraciesrandom[nt:nt+downsamplerate,:,:,mx,bx] = accuracy.squeeze()
                    # (n_timestamps,n_features,stats,task,symmetry)
                    coefsrandom[nt:nt+downsamplerate,:,:,mx,bx] = coef.squeeze()



        # aggregate random
        accuraciesrandom = np.nanmean(accuraciesrandom,axis=(0,3,4))
        coefsrandom = np.nanmean(coefsrandom,axis=(0,3,4))

        # do for all trials as control    
        predictors = neuralactivity_all
        targets = targets_all
        targets = np.tile(targets, [1,n_timestamps,1])
        targets = np.swapaxes(targets, 0, 2)
        accuracy,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'classification')
        if wx>1: coef = np.reshape(coef, [coef.shape[0], 1, wx, n_neurons, 3]).mean(axis=2)
        accuracyall = accuracy.squeeze()
        coefsall = coef.squeeze()
            

        pickle.dump((accuracies,coefs,accuracyall,coefsall,accuraciesrandom,coefsrandom), open(cacheprefix+'locomotion/movementdistribution,levels,%s-context-decoder,timecourse_%s.pck'%(roi,dn), 'wb'))
    else:
        accuracies,coefs,accuracyall,coefsall,accuraciesrandom,coefsrandom = pickle.load(open(cacheprefix+'locomotion/movementdistribution,levels,%s-context-decoder,timecourse_%s.pck'%(roi,dn), 'rb'))

    
    
    
    print(accuracies.shape, times.shape)
    if doplot:

        bodypartsmotions = preprocess.loadmovingtrialsbodyparts(dn)

        fig, ax = plt.subplots(3,n_bodyparts,figsize=(n_bodyparts*12,3*8))
        # accuracies (timepoints, train/test, stats, levels)
        # coefs (timepoints, neurons, stats, levels)

        for bx,bn in enumerate(bodypartnames):

            # plot time course of decoding accuracy
            axs = ax[0,bx]
            colors = ['black','mediumvioletred']
            for mx,lv in enumerate(levels):
                color = np.array([mx/n_levels,0,mx/n_levels])
                for k in [1]:
                    m = accuracies[:,k,0,mx,bx]
                    e = accuracies[:,k,2,mx,bx]
                    axs.plot(times, m, color=['black',color][k], lw=1, label=[None,'%4.2f'%lv][k])
                    axs.fill_between(times, m-e, m+e, color=['black',color][k],alpha=0.3)

            figs.setxt(axs)
            axs.set_ylim(0.45,1.05)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)
            axs.set_xlim(times[0],times[-1])

            # axs.legend(frameon=False)
            axs.set_ylabel('accuracy')
            axs.set_title(bn)




            # cntext decoders mean over timepoints
            axs = ax[1,bx]

            m = np.nanmean(accuracies[:,1,0,:,bx], axis=0)
            e = np.nanmean(accuracies[:,1,2,:,bx], axis=0)
            axs.plot(levels, m, color='mediumvioletred', lw=1)
            axs.fill_between(levels, m-e, m+e, color='mediumvioletred', alpha=0.3)


            axs.set_xlim(0,0.05)
            axs.set_xticks(np.arange(0,0.051,0.01))
            axs.set_ylim(0.45,1.05)
            figs.plottoaxis_chancelevel(axs,0.5)

            axs.set_xlabel('motion energy level [AU]')
            axs.legend(frameon=False)
            axs.set_ylabel('accuracy')
            axs.set_title(bn)
            
    

            # context histograms
            axs = ax[2,bx]
            bins = np.arange(0,0.05,0.0005)
            xv = bodypartsmotions[:n_trials//2,:,bx].flatten()
            xa = bodypartsmotions[n_trials//2:,:,bx].flatten()
            kst = sp.stats.ks_2samp(xv, xa, alternative='two-sided')
            
            for cx,contextlabel in enumerate(['visual','audio']):
                h,_ = np.histogram([xv,xa][cx], bins=bins)
                axs.plot(bins[:-1]+0.0001/2,h,color=['navy','darkgreen'][cx], label=contextlabel)
            axs.text(0.5,0.9, 'KS p=%6.4f'%kst.pvalue, transform=axs.transAxes, ha='left',color='black')
            axs.set_xlim(bins[0],bins[-1])
            axs.set_xticks(np.arange(0,0.051,0.01))
            axs.legend(frameon=False)
            axs.set_xlabel('motion energy level [AU]')
            axs.set_title(bn)


        fig.suptitle('%s   decoding context at movement-levels, %s '%(dn,roi))

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'movementdistribution,levels,%s-context-decoder,timecourse_%s'%(roi,dn)+ext)     #neurons, all

    



    return











def get_mask_cleverness(dn, ma=20, ma_threshold=0.5, visualfirst=True, go='nogo', congruency='incongruent'):
    # selects a subset of trials, and returns this selection mask
    # seperatively for the two contexts in a 2 elements list

    blv,bla = preprocess.getorderattended(dn) # get order of context

    thetas, index_contextchangepoint = behaviour_likelihood_idealobserver(dn, ma=ma, full=False, returnbehaviourma=True) # currently only works for the first changepoint
    # thetas dimensions are ( trials, {all,go,nogo}, {all,congruent,incongruent} )

    cix = 2-(congruency=='congruent')    # index for congruency: 1 or 2  for cong. incongr.       (0 is both together)
    gix = 2-(go=='go')                   # index for goness:     1 or 2  for go nogo              (0 is both together)
    cci = thetas[:,gix,cix]            # clever/clueless index: view on incongruent nogo trials fraction correct, as a proxy to indicate clever/clueless state
    mask_clevers = cci > ma_threshold       # designate periods where they perform incongruent nogo consistently (moving average) better than chance

    mask_clevers = [ mask_clevers[:index_contextchangepoint[0]], mask_clevers[index_contextchangepoint[0]:] ]        # could change to multiple changepoints
    if visualfirst:
        mask_clevers = [ mask_clevers[blv[1]==4], mask_clevers[blv[1]==2] ]       # reorder for visual context first, audio context second

    return mask_clevers





















def behaviour_likelihood_simplemodels(dn):

    doplot = 1 or globaldoplot


    lowprob = 1e-3     # prevent exact 1s and 0s in probabilities

    theta_data_conditional,triallist_conditional, n_trials, labels = \
        preprocess.get_behaviour_conditional_hypercube(dn,multimodalonly=True, lowprob=lowprob)

    # also get: trialindices_conditional

    # print(labels)
    # print(data_likelihood_conditional)
    # print('%s, data MLL: '%dn, marginal_loglikelihood, 'ML: ', np.exp(marginal_loglikelihood))

    # models; cube: (context, visual, audio) 2x2x2 matrix
    # m functions will give a probability mass cube from parameters
    # [[['visual,45,5000' 'visual,45,10000']   
    # ['visual,135,5000' 'visual,135,10000']]
    # [['audio,45,5000' 'audio,45,10000']     
    # ['audio,135,5000' 'audio,135,10000']]]

    model_labels = ['contextual','visual only','audio only','lick only','global bias','chance','local bias']
    model_colors = ['purple','navy','darkgreen','gold','darkorange','red','black']
    # ideal observer always perfect choice given by the model
    # for the lick bias model calculate the random bias based on the data (so not exactly chance):
    # random 0.5 chance for the last, then data dummy
    # bias  maximum likelihood estimate from the data (p = mean of x)
    theta_pre_lickbias = np.concatenate(np.concatenate(np.concatenate(triallist_conditional))).mean()
    theta_pre = [1-lowprob,1-lowprob,1-lowprob,1-lowprob,theta_pre_lickbias,0.5,np.nan]
    n_models = len(theta_pre)

    m_context = lambda c: np.array(    [[[c,   c],  [1-c, 1-c]],       \
                                        [[c, 1-c],  [c, 1-c]]]    )
    m_visual = lambda v: np.array(    [[[v,   v],  [1-v, 1-v]],       \
                                        [[v,   v],  [1-v, 1-v]]]    )
    m_audio = lambda a: np.array(     [[[a, 1-a],  [a, 1-a]],       \
                                        [[a, 1-a],  [a, 1-a]]]    )
    m_random = lambda r: np.array(     [[[r, r],  [r, r]],       \
                                        [[r, r],  [r, r]]]    )


    # solution 1
    #
    # fit the models to find the parameter at the maximum likelihood explaining the data
    # print marginal likelihood of the data with the model
    # thetas = np.zeros((5))
    # model_loglikelihood_conditional = np.zeros((5,2,2,2))
    # optimalization
    # for mx,m in enumerate([m_context,m_visual,m_audio,m_random]):
    #     costfunction = lambda theta: (  ( data_likelihood_conditional - m(theta) )  **2  ).mean()
    #     res = sp.optimize.minimize(costfunction, np.random.rand())
    #     thetas[mx] = res.x
    #     # print(model_labels[mx], res.x)
    #     model_loglikelihood_conditional[mx] = m(thetas[mx])
    # add the data likelihood
    # model_loglikelihood_conditional[4] = data_likelihood_conditional




    # solution 2
    #
    # use the models as ideal observer models, i.e. hit and correct rejection probabilities = 1.0
    # use random model to generate lick only strategy, as well as a theta = 0.5, exact chance model
    #
    # calculate likelihood for each conditioned trial list with the models

    

    model_loglikelihood_conditional = np.zeros((n_models,2,2,2))
    model_loglikelihood_marginal = np.zeros((n_models,3))      # attend visual, attend audio, total sum



    # get overall conditional and marginal likelihoods
    for mx,m in enumerate([m_context,m_visual,m_audio,m_random,m_random,m_random,m_random]):
        # print('\n%s model'%model_labels[mx])
        if mx<n_models-1: # calculate with model thetas
            theta = m( theta_pre[mx] )
        elif mx==n_models-1: # calculate data likelihood with the data conditional thetas
            theta = theta_data_conditional

        for cx in [0,1]:
            for vx in [0,1]:
                for ax in [0,1]:
                    # summary for the whole condition
                    model_loglikelihood_conditional[mx,cx,vx,ax] = neba.f_loglikelihood_bernoulli( \
                        triallist_conditional[cx,vx,ax], theta[cx,vx,ax] )
                    
                    # print(mx,'(',cx,vx,ax,')','L:',model_loglikelihood_conditional[mx,cx,vx,ax], 'p: %4.2f'%theta[cx,vx,ax], '   n=%d'%len(triallist_conditional[cx,vx,ax]) )
                    
            model_loglikelihood_marginal[mx,cx] = model_loglikelihood_conditional[mx,cx].sum()
        model_loglikelihood_marginal[mx,2] = model_loglikelihood_conditional[mx].sum()
        # print('MLL:',model_loglikelihood_marginal[mx])








    # moving average of Ls in windows of trials
    ma = 16 # moving average width for the sliding windows in number of trials
    model_loglikelihood_ma = np.empty((n_models,2),dtype=object)          # for the two contexts
    g = preprocess.loadexperimentdata(dn)
    g['block']+=1
    g['success'] = g['punish']==False
    g['action'] = np.logical_not(np.logical_xor(g['water'], g['success']))
    blv,bla = preprocess.getorderattended(dn)

    for cx in [0,1]:
        h = g[g['block'].isin([[2,4][cx]])]      # use only the multimodal blocks
        n_trials_context = len(h)

        # for model estimated bernoulli theta parameters:
        theta_ma = np.zeros((n_models,n_trials_context))
        
        cubeindex_trials = pd.DataFrame()          
        cubeindex_trials['context'] = h['block']==2       #  first context, later we calculate if it is visual or audio
        cubeindex_trials['degree'] = (h['degree']==45).astype(np.int16)
        cubeindex_trials['freq'] = (h['freq']==5000).astype(np.int16)


        # for performance displays:
        theta_all_ma = np.zeros((n_trials_context))
        theta_congruent_ma = np.zeros((n_trials_context))
        theta_incongruent_ma = np.zeros((n_trials_context))


        for mx,m in enumerate([m_context,m_visual,m_audio,m_random,m_random,m_random,m_random]):
            model_loglikelihood_ma[mx,cx] = np.zeros((n_trials_context))
            for t in range( n_trials_context ):
                sl = slice( max(0,t-ma//2), min(t+ma//2,n_trials_context-1) )

                if mx<n_models-1:             # get the probability from the model based choices
                    response_cube = m( theta_pre[mx] )
                    # create actions for the models; 0==(bla[1]==2) means if the context is visual (0) or audio (1), but order in the session is kept
                    #  bla[1]).astype(np.int16)
                    actions = response_cube[(cubeindex_trials['context']==(blv[1]==2)).astype(np.int16)[sl], cubeindex_trials['degree'][sl], cubeindex_trials['freq'][sl]]
                    theta_ma[mx][t] = actions.mean()

                elif mx==n_models-1:          # local bias model: get the probability of action from the data
                    theta_ma[mx][t] = h['action'][ sl ].mean()
                    
                
                # the above generated simulation creates sliding theta values use this for LL:
                # calculate the log likelihood of the data with the models for this sliding window
                # + add some compensatory numbers for the moving average windows at the edges of the data (   log p^x1+x2 = log p^ N xaverage = N xaverage log p )
                model_loglikelihood_ma[mx,cx][t] = neba.f_loglikelihood_bernoulli( \
                        h['action'][sl].values, theta_ma[mx][t] )  *  ma / len(g[sl])
    
    

        


    # for performance displays:
    theta_all_ma = np.zeros((n_trials,3))
    theta_congruent_ma = np.zeros((n_trials,3))
    theta_incongruent_ma = np.zeros((n_trials,3))
    h = g[g['block'].isin([2,4])]
    n_trials = len(h)
    index_contextchangepoint = h['block'].ne(2).values.argmax()          # get the splitpoint index between contexts relative within the selected mm blocks

    for mx,m in enumerate([m_context,m_visual,m_audio,m_random,m_random,m_random,m_random]):
        for t in range( n_trials ):
            sl = slice( max(0,t-ma//2), min(t+ma//2,n_trials-1) )

            # and also gather the performance moving averages:
            theta_all_ma[t,0] = h['success'][ sl ].mean()
            theta_all_ma[t,1] = h[ sl ][h['water']==True]['success'].mean()
            theta_all_ma[t,2] = h[ sl ][h['water']==False]['success'].mean()

            congruent_trials = h[ sl ]
            mask = ((congruent_trials['degree']==45) & (congruent_trials['freq']==5000)) |  \
                ((congruent_trials['degree']==135) & (congruent_trials['freq']==10000))
            theta_congruent_ma[t,0] = congruent_trials['success'][  mask  ].mean()
            theta_congruent_ma[t,1] = congruent_trials[h['water']==True]['success'][  mask  ].mean()
            theta_congruent_ma[t,2] = congruent_trials[h['water']==False]['success'][  mask  ].mean()

            incongruent_trials = h[ sl ]
            mask = ((incongruent_trials['degree']==45) & (incongruent_trials['freq']==10000)) |  \
                ((incongruent_trials['degree']==135) & (incongruent_trials['freq']==5000))
            theta_incongruent_ma[t,0] = incongruent_trials['success'][  mask  ].mean()
            theta_incongruent_ma[t,1] = incongruent_trials[h['water']==True]['success'][  mask  ].mean()
            theta_incongruent_ma[t,2] = incongruent_trials[h['water']==False]['success'][  mask  ].mean()












    # both indices will have two contexts visual and audio
    indices_congruent = np.array([ [[0,0],[1,1]], [[0,0],[1,0]] ],dtype=np.int16)
    indices_conflicting = np.array([ [[0,1],[1,0]], [[0,1],[1,1]] ],dtype=np.int16)


    if doplot:

        congruency_labels = ['congruent go','congruent nogo','conflicting go','conflicting nogo']
        xticklist = -np.array([0,1,2, 4,5,6, 8])    # corresponding to easy interpretation of the model groups and data
        
        
        
        if 1:
            fig,ax = plt.subplots(5,4,figsize=(4*8,5*8))

            for cx in [0,1]:
                
                p_congruent_go = model_loglikelihood_conditional[:,cx,indices_congruent[cx,0,0],indices_congruent[cx,0,1]]
                p_congruent_nogo = model_loglikelihood_conditional[:,cx,indices_congruent[cx,1,0],indices_congruent[cx,1,1]]
                p_conflicting_go = model_loglikelihood_conditional[:,cx,indices_conflicting[cx,0,0],indices_conflicting[cx,0,1]]
                p_conflicting_nogo = model_loglikelihood_conditional[:,cx,indices_conflicting[cx,1,0],indices_conflicting[cx,1,1]]


                for gfx,p_gfs in enumerate([ [p_congruent_go, p_congruent_nogo], [p_conflicting_go, p_conflicting_nogo]]):  # congruent conflicting
                    for gx in [0,1]:     # go, nogo
                        axs = ax[gfx,cx*2+gx]

                        for mx in range(n_models):
                            # axs.bar(x=xticklist[mx],height=p_gfs[gx][mx],color=model_colors[mx])
                            axs.plot(p_gfs[gx][mx],xticklist[mx], ['o','d'][mx==n_models-1], markersize=24,color=model_colors[mx])
                            axs.text(p_gfs[gx][mx],xticklist[mx]+0.2, '%4.2f'%p_gfs[gx][mx], color=model_colors[mx])
                        
                        if gfx==0:
                            axs.set_title('%s\n%s'%(['visual set','audio set'][cx],congruency_labels[gfx*2+gx]))
                        else:
                            axs.set_title('%s'%(congruency_labels[gfx*2+gx]))

                        # if cx*2+gx+gfx==0: axs.legend(model_labels, frameon=False, loc='lower left')

                        # axs.set_ylim(-0.55,0.55)

                        axs.set_yticks([])
                        if gx+cx==0:
                            axs.set_yticks(xticklist)
                            axs.set_yticklabels(model_labels)
                            axs.set_ylabel('conditional log likelihood\n%s'%['congruent','conflicting'][gfx])
                        axs.set_ylim(xticklist[-1]-0.5,0.5)
                        
                        # figs.plottoaxis_chancelevel(axs,p_gfs[gx][-1])




            ax_center1 = fig.add_subplot(5,3,8)
            ax_center2 = fig.add_subplot(5,1,4)
            ax_center3 = fig.add_subplot(5,1,5)

            for cx in [0,1,2]:
                if cx<2: axs = ax[2,cx*3]
                else:    axs = ax_center1

                for mx in range(n_models):
                    # axs.bar(x=xticklist[mx],height=p_gfs[gx][mx],color=model_colors[mx])
                    axs.plot(model_loglikelihood_marginal[mx,cx],xticklist[mx], ['o','d'][mx==n_models-1], markersize=24,color=model_colors[mx])
                    axs.text(model_loglikelihood_marginal[mx,cx]*0.9,xticklist[mx], '%4.2f'%model_loglikelihood_marginal[mx,cx], color=model_colors[mx],verticalalignment='center')

                    axs.set_yticks([])
                    if cx==0:
                        axs.set_yticks(xticklist)
                        axs.set_yticklabels(model_labels)
                        axs.set_ylabel('marginal log likelihood')
                    axs.set_ylim(xticklist[-1]-0.5,0.5)

                    axs.set_title(['visual set','audio set','total'][cx])





            # log likelihood moving averages

            axs = ax_center2
            ylim_aux = []
            for mx in range(n_models):
                if mx==0: ls='.-'
                else: ls='-'
                model_loglikelihood_ma_bothcontexts = np.concatenate(model_loglikelihood_ma[mx,:])
                axs.plot(model_loglikelihood_ma_bothcontexts, ls, lw=2, color=model_colors[mx], label=model_labels[mx])
                if mx!=3: ylim_aux.append(model_loglikelihood_ma_bothcontexts)       # collect for automatic optimal y_lim placement

            axs.legend(frameon=False)

            ylim_aux = np.concatenate(ylim_aux)#[0,1,2,4,5,6]])#.flatten()
            axs.set_ylim(ylim_aux.mean()-ylim_aux.std()*4,ylim_aux.mean()+ylim_aux.std()*4)
            axs.vlines(index_contextchangepoint-0.5,axs.get_ylim()[0],axs.get_ylim()[1],color='black')

            axs.set_ylabel('log likelihood')




            # performance moving average
            axs = ax_center3
            gng = ['all','go','nogo']
            ls = ['-','-','--']
            for k in [0,1,2]:       # all, go, nogo
                lws = 1/((k>0)*2+1)
                axs.plot(theta_all_ma[:,k], ls[k], lw=5*lws, color='black', label='data performance, all, %s'%gng[k])
                axs.plot(theta_congruent_ma[:,k], ls[k], lw=3*lws, color='darkturquoise', label='data performance, congruent only, %s'%gng[k])
                axs.plot(theta_incongruent_ma[:,k], ls[k], lw=3*lws, color='deeppink', label='data performance, incongruent only, %s'%gng[k])

            axs.vlines(index_contextchangepoint-0.5,axs.get_ylim()[0],axs.get_ylim()[1],color='black')

            axs.legend(frameon=False)

            axs.set_xlim(ax_center2.get_xlim())
            axs.set_ylim(-0.05,1.05)
            figs.plottoaxis_chancelevel(axs,0.5)
            
            axs.set_ylabel('data fraction correct')

            cls = [ ['visual','audio'],['audio','visual'] ][blv[1]==2]
            axs.set_xlabel('attend %s                                               attend %s\ntrial number'%(cls[0],cls[1])+\
                               ', ma=%d'%ma)





            for axs in [ax[2,1],ax[2,2],  ax[3,0],ax[3,1],ax[3,2],ax[3,3],  ax[4,0],ax[4,1],ax[4,2],ax[4,3]]:
                figs.invisibleaxes(axs,which=['right','top','left','bottom'])
    












            fig.suptitle('%s, choice models and data, conditional and marginal log likelihoods'%(dn))
            # fig.tight_layout()



            save = 0 or globalsave
            if save:
                fig.savefig(resultpath+'strategymodel_conditional,marginal,movingaverage-likelihoods_%s'%(dn)+ext)






                    

    return










def behaviour_likelihood_idealobserver(dn, ma=20, full=True, multimodalonly=True,
                returnlikelihoodratio=False, returnincongruentlikelihood=False,
                returnbehaviourma=False, onlybehaviourma=False):
    # print(dn,'behaviour moving average')
    doplot = 0 or globaldoplot


    lowprob = 5e-2     # prevent exact 1s and 0s in probabilities

    g = preprocess.loadexperimentdata(dn,full=full,multimodalonly=multimodalonly)
    g['block']+=1
    g['success'] = g['punish']==False
    g['action'] = np.logical_not(np.logical_xor(g['water'], g['success']))
    blv,bla = preprocess.getorderattended(dn)
    # labels_contextorder = [ ['visual','audio'],['audio','visual'] ][bla[1]==2]        # only works for two context sets (one shift)

    g = g[g['block'].isin([2,4])]      # use only the multimodal blocks
    n_trials = len(g)
    g['conditionedindex'] = np.arange(n_trials)
    # index_contextchangepoint = g['block'].ne(2).values.argmax()          # get the splitpoint index between contexts relative within the selected mm blocks
    index_contextchangepoint = np.where(np.abs(np.diff(g['block']))>0)[0]+1    # splitpoints list for more than a single context set shift
    labels_contextorder = np.array(['visual','audio'])[ (g['block'].iloc[np.r_[0,index_contextchangepoint]]==bla[1]).values.astype(np.int16) ]
    colors_contextorder = np.array(['navy','darkgreen'])[ (g['block'].iloc[np.r_[0,index_contextchangepoint]]==bla[1]).values.astype(np.int16) ]

    congruent_mask = ((g['degree']==45) & (g['freq']==5000)) |  \
           ((g['degree']==135) & (g['freq']==10000))
    idx_congruent_go = g['conditionedindex'].loc[g[congruent_mask & (g['water']==True)].index].values
    idx_congruent_nogo = g['conditionedindex'].loc[g[congruent_mask & (g['water']==False)].index].values
    incongruent_mask = ((g['degree']==45) & (g['freq']==10000)) |  \
           ((g['degree']==135) & (g['freq']==5000))
    idx_incongruent_go = g['conditionedindex'].loc[g[incongruent_mask & (g['water']==True)].index].values
    idx_incongruent_nogo = g['conditionedindex'].loc[g[incongruent_mask & (g['water']==False)].index].values






    key = [['freq','degree'],['degree','freq']][blv[1]==2]
    valuego = [[5000,45],[45,5000]][blv[1]==2]
    signals_contextual = np.concatenate([ g[g['block']==k][key[kx]]==valuego[kx] for kx,k in enumerate([2,4]) ]).astype(np.int16)
    signals_visual = (g['degree']==45).values.astype(np.int16)
    signals_audio = (g['freq']==5000).values.astype(np.int16)
    signals = [signals_contextual,signals_visual,signals_audio]
    n_models = len(signals)+2
    choices = g['action'].values.astype(np.int16)



    labels_models = ['contextual','visual only','audio only','lick bias','chance']
    colors_models = ['purple','navy','darkgreen','goldenrod','red']


    LLs = np.zeros((n_trials,n_models))
    for mx,_ in enumerate(labels_models):
        if mx<3:
            p = signals[mx] * (1-2*lowprob) + lowprob
        elif mx==3:
            p = choices.mean()
        elif mx==4:
            p = 0.5
        LLs[:,mx] = (     choices * np.log(  p  )   +   (1-choices) * np.log( (1-p) )     )





    # for performance displays:
    theta_all_ma = np.zeros((n_trials,3))            # both, go, nogo trials
    theta_congruent_ma = np.zeros((n_trials,3))
    theta_incongruent_ma = np.zeros((n_trials,3))
    h = g[g['block'].isin([2,4])]
    n_trials = len(h)
    # index_contextchangepoint = h['block'].ne(2).values.argmax()          # get the splitpoint index between contexts relative within the selected mm blocks
    index_contextchangepoint = np.where(np.abs(np.diff(g['block']))>0)[0]+1    # splitpoints list for more than a single context set shift

    # start = [0,index_contextchangepoint]
    # stop = [index_contextchangepoint,n_trials]
    start = np.r_[0,index_contextchangepoint]
    stop = np.r_[index_contextchangepoint, n_trials]

    for cx in range(len(start)):
        for t in np.arange(start[cx],stop[cx]):
            sl = slice( max(start[cx],t-ma//2), min(t+ma//2,stop[cx]) )
            watersl = h[sl]['water']==True

            # and also gather the performance moving averages:
            theta_all_ma[t,0] = h['success'][ sl ].mean()
            theta_all_ma[t,1] = h['success'][ sl ][  watersl ].mean()
            theta_all_ma[t,2] = h['success'][ sl ][ ~watersl ].mean()

            congruent_trials = h[ sl ]
            mask = ((congruent_trials['degree']==45) & (congruent_trials['freq']==5000)) |  \
                ((congruent_trials['degree']==135) & (congruent_trials['freq']==10000))
            theta_congruent_ma[t,0] = congruent_trials['success'][  mask ].mean()
            theta_congruent_ma[t,1] = congruent_trials['success'][  mask &  watersl ].mean()
            theta_congruent_ma[t,2] = congruent_trials['success'][  mask & ~watersl ].mean()

            incongruent_trials = h[ sl ]
            mask = ((incongruent_trials['degree']==45) & (incongruent_trials['freq']==10000)) |  \
                ((incongruent_trials['degree']==135) & (incongruent_trials['freq']==5000))
            theta_incongruent_ma[t,0] = incongruent_trials['success'][  mask  ].mean()
            theta_incongruent_ma[t,1] = incongruent_trials['success'][  mask &  watersl ].mean()
            theta_incongruent_ma[t,2] = incongruent_trials['success'][  mask & ~watersl ].mean()




    # followings is for individual mice calculation returns for function calls
    if returnlikelihoodratio or returnbehaviourma:
        R = []

        # strategy ideal observer likelihoods ratios to chance (LL differences)
        # LLs are log likelihoods per trial for each model
        # we return three log likelihood differences: contextual, visual only and audio only to chance ratio
        # return all trials and return
        if returnlikelihoodratio:
            LLDs = LLs[:,0:3] - LLs[:,4][:,np.newaxis]
            R.extend([ LLDs.mean(axis=0), LLDs[idx_incongruent_nogo,:].mean(axis=0) ])


        # moving average behaviours, 
        if returnbehaviourma:
            theta_collect = np.stack([theta_all_ma, theta_congruent_ma, theta_incongruent_ma],axis=2)
            R.extend([theta_collect,index_contextchangepoint])
            # pickle.dump((theta_collect,index_contextchangepoint),open(cacheprefix+'behaviour/fractioncorrect-both,congruent,incongruent-all,go,nogo_%s'%(dn),'wb'))
        
        return R



    if returnincongruentlikelihood:



        return




    # display figure


    if doplot:



        if onlybehaviourma:
            fig,ax = plt.subplots(1,3,figsize=(3*12,1*8))
        else:
            fig,ax = plt.subplots(2,1,figsize=(1*20,2*8))
        
        
        if not onlybehaviourma:
            axs = ax[0]
            for mx in range(n_models):
                if mx==0: l='o-'
                else: l='-'
                axs.plot(LLs[:,mx],l,color=colors_models[mx],lw=2,label='%s  mLL=%4.2f'%(labels_models[mx],LLs[:,mx].mean()))

            mx = 0
            axs.plot(idx_congruent_go, 0.1*np.ones(len(idx_congruent_go)),'^',color='lightgrey',label='congruent, go  %s mLL=%4.2f'%(labels_models[mx],LLs[idx_congruent_go,mx].mean()))
            axs.plot(idx_congruent_nogo, 0.1*np.ones(len(idx_congruent_nogo)),'x',color='lightgrey',label='congruent, nogo  %s mLL=%4.2f'%(labels_models[mx],LLs[idx_congruent_nogo,mx].mean()))

            axs.plot(idx_incongruent_go, 0.1*np.ones(len(idx_incongruent_go)),'^',color='cyan',label='incongruent, go  %s mLL=%4.2f'%(labels_models[mx],LLs[idx_incongruent_go,mx].mean()))
            axs.plot(idx_incongruent_nogo, 0.1*np.ones(len(idx_incongruent_nogo)),'x',color='red',label='incongruent, nogo  %s mLL=%4.2f'%(labels_models[mx],LLs[idx_incongruent_nogo,mx].mean()))

            axs.set_ylim(np.min(LLs.flatten())*1.05,0.12)

            for ichp in index_contextchangepoint:
                axs.vlines(ichp-0.5,axs.get_ylim()[0],axs.get_ylim()[1],ls='--',lw=2,color='black',alpha=0.2)

            axs.set_ylabel('single trial log likelihood')
            # axs.set_xlabel('trial number\n%s              %s'%(labels_contextorder[0],labels_contextorder[1]))
            axs.set_xlabel('trial number')


            axs.legend(bbox_to_anchor=(1,0),loc='lower left',frameon=False)



        if not onlybehaviourma:
            axs = ax[1]
        # performance moving average
        gng = ['all','go','nogo']
        acolors = ['grey','black','red']     # both congruent incongruent
        # dark/light: relevant/irrelevant
        # solid/dash: congruent/incongruent
        # darkturquoise/deeppink: go/nogo         ->        black/red: go/nogo        
        # darkorchid/fuchsia  and    goldenrod/gold:               clever,relevant/clever,irrelevant   and    clueless,relevant/clueless,irrelevant
        for cx in range(len(start)):
            for k in [0,1,2]:       # all, go, nogo      [ 0,1,2 ]
                axs = ax[k]
                lws = 1/((k>0)*2+1)
                label = ['fraction correct, %s'%gng[k], 'fraction correct, %s'%gng[k], 'fraction correct, %s'%gng[k] ][k]
                # all, congruent, incongruent:
                for ex,(theta_,lwsmul,thls,labelpostfix) in enumerate(zip([theta_all_ma, theta_congruent_ma, theta_incongruent_ma],[6,3,3],\
                                            ['-','-','--'],[', both',', congruent only',', incongruent only'])):
                    if k==0 and ex>0: continue
                    axs.plot(np.arange(start[cx],stop[cx]), theta_[start[cx]:stop[cx],k], ls=thls, lw=lwsmul*lws, color=acolors[k], label=[label+labelpostfix,None][cx>0])
                if k>0: axs.legend(['all','congruent','incongruent'],frameon=False,loc='lower right')
                axs.set_title(['performance','performance go','performance nogo'][k])

                ee_icp = np.r_[0,index_contextchangepoint,n_trials]
                for ix in range(len(ee_icp)-1):
                    if ix>0:   axs.vlines(ee_icp[ix]-0.5,-0.05,1.05,ls='--',lw=2,color='black',alpha=0.2)
                    axs.text((ee_icp[ix]+ee_icp[ix+1])/2, 0.01, labels_contextorder[ix], color=colors_contextorder[ix],alpha=0.8,horizontalalignment='center')


                if not onlybehaviourma: axs.set_xlim(ax[0].get_xlim())
                axs.set_ylim(-0.05,1.05)
                figs.plottoaxis_chancelevel(axs,0.5)
                if k==0: axs.set_ylabel('animal %s\nfraction correct'%dn)
                # axs.set_xlabel('trial number, ma=%d'%ma)
                axs.set_xlabel('trial number')
        
        # axs.legend(bbox_to_anchor=(1,0),loc='lower left',frameon=False)

        # cls = [ ['visual','audio'],['audio','visual'] ][blv[1]==2]



        if onlybehaviourma:
            # print()
            fig.suptitle('%s behaviour moving averages'%(dn))
        else:
            fig.suptitle('%s, ideal (%d%%) observer models, single trial and mean log likelihoods'%(dn,100-lowprob*100))
        fig.tight_layout()




        save = 0 or globalsave
        if save:
            if onlybehaviourma:
                fig.savefig(resultpath+'movingaverage%d,performance_%s'%(ma,dn)+ext)
            else:
                fig.savefig(resultpath+'strategymodel_idealobserversingletrial,movingaverageperformance-likelihoods_%s'%(dn)+ext)





    return















def behaviour_likelihood_sigmoidlinearmodels(dn):

    doplot = 1 or globaldoplot


    # load data

    g = preprocess.loadexperimentdata(dn)
    g['block']+=1
    g['success'] = g['punish']==False
    g['action'] = np.logical_not(np.logical_xor(g['water'], g['success']))
    blv,bla = preprocess.getorderattended(dn)

    g = g[g['block'].isin([2,4])]      # use only the multimodal blocks
    n_trials = len(g)
    g['conditionedindex'] = np.arange(n_trials)
    index_contextchangepoint = g['block'].ne(2).values.argmax()          # get the splitpoint index between contexts relative within the selected mm blocks


    congruent_mask = ((g['degree']==45) & (g['freq']==5000)) |  \
           ((g['degree']==135) & (g['freq']==10000))
    idx_congruent_go = g['conditionedindex'].loc[g[congruent_mask & (g['water']==True)].index].values
    idx_congruent_nogo = g['conditionedindex'].loc[g[congruent_mask & (g['water']==False)].index].values
    incongruent_mask = ((g['degree']==45) & (g['freq']==10000)) |  \
           ((g['degree']==135) & (g['freq']==5000))
    idx_incongruent_go = g['conditionedindex'].loc[g[incongruent_mask & (g['water']==True)].index].values
    idx_incongruent_nogo = g['conditionedindex'].loc[g[incongruent_mask & (g['water']==False)].index].values





    key = [['freq','degree'],['degree','freq']][blv[1]==2]
    valuego = [[5000,45],[45,5000]][blv[1]==2]
    signals_contextual = np.concatenate([ g[g['block']==k][key[kx]]==valuego[kx] for kx,k in enumerate([2,4]) ]).astype(np.int16)
    signals_visual = (g['degree']==45).values.astype(np.int16)
    signals_audio = (g['freq']==5000).values.astype(np.int16)
    signals_bias = np.zeros(n_trials)
    signals = [signals_contextual,signals_visual,signals_audio,signals_bias]
    n_models = len(signals)
    choices = g['action'].values.astype(np.int16)




    sigmoid = lambda x: 1/(1+np.exp(-x))
    def negativeloglikelihood(theta,S,d):
        # theta = (w1,w2,w3,...,b)   parameters to find where return is extreme
        # S stimulus conditions for trials, N by S matrix to drive the parametrized models
        # d animal choices for trials, N length vector
        w = theta[:-1]
        b = theta[-1]
        # p = sigmoid(np.dot(S,w)+b)     # vector stimulus
        p = sigmoid(S*w+b)               # single stimulus
        LL = (     d * np.log(  p  )     +    (1-d) * np.log( (1-p) )     ).sum()
        return -LL


    # cross validation schene
    k_fold = 10
    kf = skms.KFold(n_splits=k_fold)




    
    # fit the models to find the parameter at the maximum likelihood explaining the data
    W_hat = np.zeros((n_models,k_fold))
    B_hat = np.zeros((n_models,k_fold))
    LLs = np.zeros((n_models,n_trials,2))          # train test
    # optimization
    for (fold, (idx_tr,idx_test)) in enumerate( kf.split(np.arange(n_trials)) ):
        for mx,signal in enumerate(signals):
            # costfunction = lambda theta: (  ( data_likelihood_conditional - theta )  **2  ).mean()
            n_s = 1
            res = sp.optimize.minimize(fun=negativeloglikelihood, x0=np.random.randn(n_s+1), args=(signal[idx_tr],choices[idx_tr]))
            W_hat[mx,fold] = res.x[:-1][0]
            if mx==3: W_hat[mx,fold] = 0 # bias model does not fit w (signals are zeroed out)
            B_hat[mx,fold] = res.x[-1]

            # now we have estimates
            # calculate the sigmoids and the data likelihood for all folds and models
            for ix,idx in enumerate([idx_tr,idx_test]):
                p = sigmoid(  signals[mx]*W_hat[mx,ix] + B_hat[mx,ix] )[idx]
                aux = (     choices[idx] * np.log(  p  )     +    (1-choices[idx]) * np.log( (1-p) )     )
                if fold==0: # train, mean of all folds
                    LLs[mx,idx,fold] += aux / (k_fold-1)
                elif fold==1:   # test, single held out
                    LLs[mx,idx,fold] = aux




    # plot

    labels_models = ['contextual','visual only','audio only','bias']
    colors_models = ['purple','navy','darkgreen','goldenrod']
    context_order = [ ['visual','audio'],['audio','visual'] ][bla[1]==2]

    fig,ax = plt.subplots(2,4,figsize=(4*8,2*8))

    
    for mx in range(n_models):


        # choice probability
        axs = ax[0,mx]
        t = np.arange(-2,2.,0.01)
        for k in range(k_fold):
            p = sigmoid(  t*W_hat[mx,k] + B_hat[mx,k]  )
            axs.plot(t,p,lw=2,color=colors_models[mx])
        axs.set_ylim(0,1)
        axs.set_xlim(-1,2)
        figs.plottoaxis_crosshair(axs,0.0,0.5)
        figs.plottoaxis_crosshair(axs,1.0,0.5)

        if mx==0: axs.set_ylabel('choice probability, P(lick)')
        # axs.set_yticks([0,0.5,1])
        # axs.set_yticklabels(['0 withhold','0.5 chance', '1 lick'])
        axs.set_xticks([0,1])
        axs.set_xticklabels(['0\nnogo','1\ngo'])
        axs.set_title(labels_models[mx])




        # log likelihood
        axs = ax[1,mx]

        axs.plot(LLs[mx,:,0],'-',color=colors_models[mx],alpha=0.3,label='%s, train, %4.2f'%(labels_models[mx],LLs[mx,:,0].mean()))
        axs.plot(LLs[mx,:,1],'-',color=colors_models[mx],label='%s, test %4.2f'%(labels_models[mx],LLs[mx,:,1].mean()))

        axs.plot(idx_congruent_go, LLs[mx,idx_congruent_go,1],'^',color=colors_models[mx],label='congruent, go %4.2f'%LLs[mx,idx_congruent_go,1].mean())
        axs.plot(idx_congruent_nogo, LLs[mx,idx_congruent_nogo,1],'x',color=colors_models[mx],label='congruent, nogo %4.2f'%LLs[mx,idx_congruent_nogo,1].mean())

        axs.plot(idx_incongruent_go, LLs[mx,idx_incongruent_go,1],'^',color='cyan',label='incongruent, go %4.2f'%LLs[mx,idx_incongruent_go,1].mean())
        axs.plot(idx_incongruent_nogo, LLs[mx,idx_incongruent_nogo,1],'x',color='red',label='incongruent, nogo %4.2f'%LLs[mx,idx_incongruent_nogo,1].mean())

        axs.set_ylim(np.min(LLs.flatten()),0)

        axs.vlines(index_contextchangepoint-0.5,axs.get_ylim()[0],axs.get_ylim()[1],ls='--',lw=2,color='black',alpha=0.2)

        axs.legend(frameon=False)

        if mx==0: axs.set_ylabel('single trial log likelihood')
        axs.set_xlabel('trial number\n%s              %s'%(context_order[0],context_order[1]))





    fig.suptitle('%s, cross validated linear sigmoid model fit, numbers in legend: mean trial LL'%dn)
    fig.tight_layout()


    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'strategymodel_sigmoidlinear,loglikelihood_%s'%(dn)+ext)











def behaviour_symmetry_highperformance(dn):


    # get high performance and low performanc trial indices

    g = preprocess.loadexperimentdata(dn)
    g['block']+=1
    g['success'] = g['punish']==False
    blv,bla = preprocess.getorderattended(dn)

    # apply behaviour masks
    action = ['go','nogo']
    congruency = ['congruent','incongruent']
    mask_clevers_list = [[],[]] # holds context dependent list
    mask_clevers = []           # holds indexed by original trial order
    for c in congruency:
        for a in action:
            mask_clever_contexts = get_mask_cleverness(dn, ma_threshold=0.5, visualfirst=True, go=a, congruency=c)
            for cx,mask_clever_context in enumerate(mask_clever_contexts):
                mask_clevers_list[cx].append(mask_clever_context)
            mask_clever = np.hstack( get_mask_cleverness(dn, ma_threshold=0.5, visualfirst=False, go=a, congruency=c) )
            mask_clevers.append(mask_clever)
    
    mask_clevers_list = [ np.vstack(mask_clevers_list[cx]).T for cx in [0,1]]
    mask_clevers = np.vstack(mask_clevers).T


    mask_contextuals = [ np.prod(mask_clevers_list[cx],axis=1) for cx in [0,1] ]
    # display the number of trials that has above threshold movingaverage for all 4 combinations of congruenct and action
    print(dn, 'V:%d/%d, A:%d/%d'%(np.sum(mask_contextuals[0]), len(mask_contextuals[0]), np.sum(mask_contextuals[1]), len(mask_contextuals[1])))
    mask_highperf = np.bool8(np.prod(mask_clevers, axis=1))









    
    # get congruency indices

    g = preprocess.loadexperimentdata(dn,full=True,multimodalonly=True)
    g['block']+=1
    g['success'] = g['punish']==False
    g['action'] = np.logical_not(np.logical_xor(g['water'], g['success']))
    blv,bla = preprocess.getorderattended(dn)
    # labels_contextorder = [ ['visual','audio'],['audio','visual'] ][bla[1]==2]        # only works for two context sets (one shift)

    g = g[g['block'].isin([2,4])]      # use only the multimodal blocks
    n_trials = len(g)
    g['conditionedindex'] = np.arange(n_trials)
    # index_contextchangepoint = g['block'].ne(2).values.argmax()          # get the splitpoint index between contexts relative within the selected mm blocks
    index_contextchangepoint = np.where(np.abs(np.diff(g['block']))>0)[0]+1    # splitpoints list for more than a single context set shift
    labels_contextorder = np.array(['visual','audio'])[ (g['block'].iloc[np.r_[0,index_contextchangepoint]]==bla[1]).values.astype(np.int16) ]
    colors_contextorder = np.array(['navy','darkgreen'])[ (g['block'].iloc[np.r_[0,index_contextchangepoint]]==bla[1]).values.astype(np.int16) ]

    congruent_mask = ((g['degree']==45) & (g['freq']==5000)) |  \
           ((g['degree']==135) & (g['freq']==10000))
    idx_congruent_go = g['conditionedindex'].loc[g[congruent_mask & (g['water']==True)].index].values
    idx_congruent_nogo = g['conditionedindex'].loc[g[congruent_mask & (g['water']==False)].index].values
    incongruent_mask = ((g['degree']==45) & (g['freq']==10000)) |  \
           ((g['degree']==135) & (g['freq']==5000))
    idx_incongruent_go = g['conditionedindex'].loc[g[incongruent_mask & (g['water']==True)].index].values
    idx_incongruent_nogo = g['conditionedindex'].loc[g[incongruent_mask & (g['water']==False)].index].values




    # assemble signals

    key = [['freq','degree'],['degree','freq']][blv[1]==2]
    valuego = [[5000,45],[45,5000]][blv[1]==2]
    signals_contextual = np.concatenate([ g[g['block']==k][key[kx]]==valuego[kx] for kx,k in enumerate([2,4]) ]).astype(np.int16)
    signals_visual = (g['degree']==45).values.astype(np.int16)
    signals_audio = (g['freq']==5000).values.astype(np.int16)
    
    
    # print(mask_clevers.shape, len(mask_contextuals), mask_highperf.shape, sum(mask_highperf))
    # print(signals_contextual[mask_highperf].shape, signals_contextual[congruent_mask].shape)

    # create mask for context (i.e. attend visual and attend audio, in this order)
    mask_context = g['block'] == blv[1]


    # assemble the data vector for two contexts using masks for context, congruency
    # but only use high performance trials using "mask_highperf" variable
    # and incongruency variations will be for cycled later

    signal_crossparts = [signals_audio, signals_visual]
    signals = [[[],[]],[[],[]]] # signals with two contexts and congruency variations
    action = [[[],[]],[[],[]]]
    meanresponses = np.zeros((2,2,3))
    for cx in [0,1]:
        for ix in [0,1]:
            mask = (mask_context==1-cx) & (congruent_mask==1-ix) & mask_highperf
            signals[cx][ix].append( signals_contextual[mask] )          # contextually correct
            signals[cx][ix].append( signal_crossparts[cx][mask] )       # contextually opposite
            action[cx][ix] = g['action'][mask]
    




    # calculate log likelihoods for each model, in each context and congruency variation

    lowprob = 1e-3     # prevent exact 1s and 0s in probabilities

    # relevancy in models means contextually correct or opposite modality
    # log likelihood mean over all trials
    LLs = np.zeros((2,2,2,6))      # relevancy, context, congruency, model
    # p or parameter:
    mp = np.zeros((2,2,2,6,2))         # relevancy, context, congruency, model, parameter

    for cx in [0,1]:
        for ix in [0,1]:



            for gx in [0,1,2]:
                if gx<2: continue # don't do go and nogo trials separately
                # create signal masks
                if gx<2: # go and nogo signals
                    sl = signals[cx][ix][0]==1-gx
                else: # both together
                    sl = np.tile(True,len(signals[cx][ix][0]))

                choices = action[cx][ix][sl]
                meanresponses[cx,ix,gx] = choices.mean()


                # lick bias: # no oppsosite for this, so rx=1 is omitted
                modelnumber = 0
                relevancy = 0
                p = meanresponses[cx,ix,gx]
                if gx==2: mp[relevancy,cx,ix,modelnumber,0] = p-0.5     # save overall bias for display
                p = np.clip(p,lowprob,1-lowprob)
                LLs[relevancy,cx,ix,modelnumber] = (     choices * np.log(  p  )   +   (1-choices) * np.log( (1-p) )     ).mean()


                # contextually correct:
                modelnumber = 1
                relevancy = 0
                data = signals[cx][ix][0][sl]
                p = data * (1-2*lowprob) + lowprob
                LLs[relevancy,cx,ix,modelnumber] = (     choices * np.log(  p  )   +   (1-choices) * np.log( (1-p) )     ).mean()
                mp[relevancy,cx,ix,modelnumber,:] = p.mean()

                # contextually opposite:
                modelnumber = 1
                relevancy = 1
                data = signals[cx][ix][1][sl]
                p = data * (1-2*lowprob) + lowprob
                LLs[relevancy,cx,ix,modelnumber] = (     choices * np.log(  p  )   +   (1-choices) * np.log( (1-p) )     ).mean()
                mp[relevancy,cx,ix,modelnumber,:] = p.mean()


                # contextually correct with lick bias
                modelnumber = 2
                relevancy = 0
                data = signals[cx][ix][0][sl]      # contextually opposite, but it is already loaded
                # p = data * (1-2*lowprob) + lowprob     + meanresponses[cx,ix,gx]-0.5
                p = data + meanresponses[cx,ix,gx]-0.5
                p = np.clip(p,lowprob,1-lowprob)
                LLs[relevancy,cx,ix,modelnumber] = (     choices * np.log(  p  )   +   (1-choices) * np.log( (1-p) )     ).mean()
                mp[relevancy,cx,ix,modelnumber,:] = p.mean()

                # opposite with lick bias
                modelnumber = 2
                relevancy = 1
                data = signals[cx][ix][1][sl]      # contextually opposite, but it is already loaded
                p = data * (1-2*lowprob) + lowprob     + meanresponses[cx,ix,gx]-0.5
                p = np.clip(p,lowprob,1-lowprob)
                LLs[relevancy,cx,ix,modelnumber] = (     choices * np.log(  p  )   +   (1-choices) * np.log( (1-p) )     ).mean()
                mp[relevancy,cx,ix,modelnumber,:] = p.mean()



    
    
    
    
    # parameteric maximum log likelihood models

    # caculate bernoulli log likelihood (LL) for choices at beta bias, lambda lapse
    # maximise LL, and store the best LL, and the parameter values that gave it



    def p_bias_data(x, b=0, v=[1,0]):
        # x is a vector of 1s and 0s, the target
        # b bias parameter
        # v is p=1 and p=0 classes in this order, default [1,0]
        p = np.zeros(x.shape)
        p[x==v[0]] = 1 + b
        p[x==v[1]] = b
        p = np.clip(p,lowprob,1-lowprob)
        return p

    def p_lapse_data(x, l=0, v=[1,0]):
        # x is a vector of 1s and 0s, the target
        # l lapse parameter
        # v is p=1 and p=0 classes in this order, default [1,0]
        l = np.clip(l,0.0,0.5)       # avoid lapse to reverse to opposite modality, max is fully random
        p = np.zeros(x.shape)
        p[x==v[0]] = 1 - l
        p[x==v[1]] = l
        p = np.clip(p,lowprob,1-lowprob)
        return p
    
    def p_biaslapse_data(x, b=0, l=0, v=[1,0]):
        # x is a vector of 1s and 0s, the target
        # b bias parameter
        # l lapse parameter
        # v is p=1 and p=0 classes in this order, default [1,0]
        l = np.clip(l,0.0,0.5)       # avoid lapse to reverse to opposite modality, max is fully random
        p = np.zeros(x.shape)
        p[x==v[0]] = 1 + b - l
        p[x==v[1]] = b + l
        p = np.clip(p,lowprob,1-lowprob)
        return p



    def modelLL_bias(k,x,b):      # context aware bias model
        # k is choices, x is target
        # positive b means that p is shifted towards 1
        # negative b means that p is shifted towards 0
        p = p_bias_data(x, b=b)
        ll = nedi.ll_bernoulli(k, p).mean()
        return ll

    def modelLL_lapse(k,x,l):      # context aware lapse model
        p = p_lapse_data(x, l)
        ll = nedi.ll_bernoulli(k, p).mean()
        return ll
    
    def modelLL_biaslapse(k,x,b,l):      # context aware bias+lapse model
        p = p_biaslapse_data(x, b, l)
        ll = nedi.ll_bernoulli(k, p).mean()
        return ll






    # parameter sweeps
    biasrange = np.arange(-0.50,0.51,0.01)
    lapserange = np.arange(0.00,1.01,0.01)





    # context aware parametric models

    # bias model
    modelnumber = 3
    LLb = np.zeros((2,2,2,len(biasrange))) # (relevancy, context, congruency, parameterrange)
    for rx in [0,1]:             # relevancy
        for cx in [0,1]:         # context
            for ix in [0,1]:     # congruency
                for bx,bias in enumerate(biasrange):
                    choices = action[cx][ix]
                    target = signals[cx][ix][rx]
                    LLb[rx,cx,ix,bx] = modelLL_bias(choices,target,bias)
                mLLib = np.argmax(LLb[rx,cx,ix,:])
                LLs[rx,cx,ix,modelnumber] = LLb[rx,cx,ix,mLLib]
                mp[rx,cx,ix,modelnumber,0] = biasrange[mLLib]

    # lapse model
    modelnumber = 4
    LLl = np.zeros((2,2,2,len(lapserange))) # (relevancy, context, congruency, parameterrange)
    for rx in [0,1]:             # relevancy
        for cx in [0,1]:         # context
            for ix in [0,1]:     # congruency
                for bx,lapse in enumerate(lapserange):
                    choices = action[cx][ix]
                    target = signals[cx][ix][rx]
                    LLl[rx,cx,ix,bx] = modelLL_lapse(choices,target,lapse)
                mLLil = np.argmax(LLl[rx,cx,ix,:])
                LLs[rx,cx,ix,modelnumber] = LLl[rx,cx,ix,mLLil]
                mp[rx,cx,ix,modelnumber,1] = lapserange[mLLil]

    # bias + lapse model
    modelnumber = 5
    LLbl = np.zeros((2,2,2,len(biasrange),len(lapserange))) # (relevancy, context, congruency, parameter1range, parameter2range)
    mLLibl = np.zeros((2,2,2,2),dtype=int) # (relevancy, context, congruency, parameter1argmax, parameter2argmax)
    for rx in [0,1]:             # relevancy
        for cx in [0,1]:         # context
            for ix in [0,1]:     # congruency
                for bx,bias in enumerate(biasrange):
                    for lx,lapse in enumerate(lapserange):
                        choices = action[cx][ix]
                        target = signals[cx][ix][rx]
                        LLbl[rx,cx,ix,bx,lx] = modelLL_biaslapse(choices,target,bias,lapse)
                mLLi = np.argmax(LLbl[rx,cx,ix,:,:])
                mLLibl[rx,cx,ix,:] = np.unravel_index(mLLi,(len(biasrange),len(lapserange)))
                LLs[rx,cx,ix,modelnumber] = LLbl[rx,cx,ix,mLLibl[rx,cx,ix,0],mLLibl[rx,cx,ix,1]]
                mp[rx,cx,ix,modelnumber,0] = biasrange[mLLibl[rx,cx,ix,0]]
                mp[rx,cx,ix,modelnumber,1] = lapserange[mLLibl[rx,cx,ix,1]]










    
    if globaldoplot or 0:

        fig,ax = plt.subplots(2,2,figsize=(2*8*1.2,2*8))
        # relevancy variable: rx = 0 contextual, or rx = 1 opposite
        for cx in [0,1]:        # context
            for ix in [0,1]:    # congruency

                label = ['visual context','audio context'][cx]+[' congruent',' incongruent'][ix]


                # context awaqre
                rx = 0 # relevancy

                axs = ax[0,0]
                mx = 3 # model number
                axs.plot(biasrange,LLb[rx,cx,ix,:],['-','--'][ix], lw=2, color=['navy','darkgreen'][cx], label=label)
                axs.scatter(mp[rx,cx,ix,mx,0],LLs[rx,cx,ix,mx],
                            s=100, color=['navy','darkgreen'][cx], marker=['o','x'][ix],label=None)
                axs.text(mp[rx,cx,ix,mx,0]+0.02,LLs[rx,cx,ix,4],'$\\beta$=%4.2f'%mp[rx,cx,ix,mx,0], color='grey', size='x-small', ha='left', va='center')

                axs.set_ylim(-1.5,0)
                axs.set_xlabel('$\\beta$ (bias model)')
                axs.set_title('context aware')


                axs = ax[0,1]
                mx = 4
                axs.plot(lapserange,LLl[rx,cx,ix,:],['-','--'][ix], lw=2, color=['navy','darkgreen'][cx], label=label)
                axs.scatter(mp[rx,cx,ix,mx,1],LLs[rx,cx,ix,mx],
                            s=100, color=['navy','darkgreen'][cx], marker=['o','x'][ix],label=None)
                axs.text(mp[rx,cx,ix,mx,1]+0.02,LLs[rx,cx,ix,mx],'$\\lambda$=%4.2f'%mp[rx,cx,ix,mx,1], color='grey', size='x-small', ha='left', va='center')

                axs.set_ylim(-1.5,0)
                axs.legend(fontsize=10)
                axs.set_xlabel('$\\lambda$ (lapse model)')
                axs.set_title('context aware')


                axs = ax[1,0]
                mx = 5
                axs.plot(biasrange,LLbl[rx,cx,ix,:,mLLibl[rx,cx,ix,1]],['-','--'][ix], lw=2, color=['navy','darkgreen'][cx], label=label)
                axs.scatter(mp[rx,cx,ix,mx,0],LLs[rx,cx,ix,mx],
                            s=100, color=['navy','darkgreen'][cx], marker=['o','x'][ix],label=None)
                axs.text(mp[rx,cx,ix,mx,0]+0.02,LLs[rx,cx,ix,mx],'$\\beta$=%4.2f'%mp[rx,cx,ix,mx,0], color='grey', size='x-small', ha='left', va='center')

                axs.set_ylim(-1.5,0)
                axs.legend(fontsize=10)
                axs.set_xlabel('$\\beta$ (bias+lapse model)')


                axs = ax[1,1]
                mx = 5
                axs.plot(lapserange,LLbl[rx,cx,ix,mLLibl[rx,cx,ix,0],:],['-','--'][ix], lw=2, color=['navy','darkgreen'][cx], label=label)
                axs.scatter(mp[rx,cx,ix,6,1],LLs[rx,cx,ix,mx],
                            s=100, color=['navy','darkgreen'][cx], marker=['o','x'][ix],label=None)
                axs.text(mp[rx,cx,ix,mx,1]+0.02,LLs[rx,cx,ix,mx],'$\\lambda$=%4.2f'%mp[rx,cx,ix,mx,1], color='grey', size='x-small', ha='left', va='center')

                axs.set_ylim(-1.5,0)
                axs.legend(fontsize=10)
                axs.set_xlabel('$\\lambda$ (bias+lapse model)')


        for av in [0,1]:
            for ah in [0,1,2,3,4,5]:
                axs = ax[av,ah]
                if av==0 and ah<4:
                    axs.axis('off')
                    continue

                axs.scatter(0,1, s=100, color='black', marker='o', label='contextual+%s congruent'%(['bias','lapse','bias+lapse','bias+lapse'][av*2+ah%2]))
                axs.scatter(0,1, s=100, color='black', marker='X', label='contextual+%s incongruent'%(['bias','lapse','bias+lapse','bias+lapse'][av*2+ah%2]))
                axs.legend(fontsize=10)
                axs.set_ylabel('LL')




        fig.suptitle(dn)

        fig.tight_layout()

        if globalsave or 0:
            fig.savefig(resultpath+'choice-p,LL-%s'%dn+ext)




    return LLs, meanresponses, mp











def behaviour_symmetry(dn,block,displaystatsonly=False):


    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot


    # load trials data

    g = preprocess.loadexperimentdata(dn)
    g['block']+=1
    g['success'] = g['punish']==False
    blv,bla = preprocess.getorderattended(dn)

    # apply behaviour masks
    action = ['go','nogo']
    congruency = ['congruent','incongruent']
    mask_clevers_list = [[],[]] # holds context dependent list
    mask_clevers = []           # holds indexed by original trial order
    for c in congruency:
        for a in action:
            mask_clever_contexts = get_mask_cleverness(dn, ma_threshold=0.5, visualfirst=True, go=a, congruency=c)
            for cx,mask_clever_context in enumerate(mask_clever_contexts):
                mask_clevers_list[cx].append(mask_clever_context)
            mask_clever = np.hstack( get_mask_cleverness(dn, ma_threshold=0.5, visualfirst=False, go=a, congruency=c) )
            mask_clevers.append(mask_clever)
    
    mask_clevers_list = [ np.vstack(mask_clevers_list[cx]).T for cx in [0,1]]
    mask_clevers = np.vstack(mask_clevers).T


    mask_contextuals = [ np.prod(mask_clevers_list[cx],axis=1) for cx in [0,1] ]
    # display the number of trials that has above threshold movingaverage for all 4 combinations of congruenct and action
    print(dn, 'V:%d/%d, A:%d/%d'%(np.sum(mask_contextuals[0]), len(mask_contextuals[0]), np.sum(mask_contextuals[1]), len(mask_contextuals[1])))
    if displaystatsonly: return
    mask = np.bool8(np.prod(mask_clevers, axis=1))








    # decode from full and from symmetric parts, and compare on the same plot
    symmetrylabels = ['entire session','successful periods','error periods']


    # setup stimulus:
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [   [ [ [], [45],    [] ],        [    [],[135],     [] ]  ], \
                            [ [ [], [],    [5000] ],        [    [],[],     [10000] ]  ], \
                            [ [ [blv[1]],  [],[] ], [ [bla[1]],   [],[] ]    ],     \
                            [ [True], [False] ], \
                        ]
    taskaspects = ['visual','audio','context','choice']
    columnnames = ['degree','freq','block','action']
    class1value = [45,5000,blv[1],True,0]
    
    
    # create the all trials, and just symmetriccal successful behaviour trial list
    triallist_all = preprocess.loadexperimentdata(dn, multimodalonly=True)
    triallist_all['block']+=1
    triallist_symmetric = triallist_all.iloc[mask,:]
    triallist_error = triallist_all.iloc[np.logical_not(mask),:]
    


    # collect neural data:
    wx = int((50*pq.ms/T['dt']).magnitude)  # width of the feature space in time for each neuron
    # wx = 0
    times = block.segments[0].analogsignals[0].times
    if wx>1: times = times[:-wx]
    n_timestamps = len(times)
    n_neurons = block.segments[0].analogsignals[0].shape[1]
    n_features = n_neurons
    n_targets = 1

    allcomplextrials = [ [[2,4], [], []] ]
    # (observations,times,features)
    neuralactivity_all = np.array(preprocess.collect_stimulusspecificresponses(block, allcomplextrials)[0])
    if wx>1:   # create wider feature space of several subsequent timepoints
        neuralactivity_all = neph.slidingwindow(neuralactivity_all,wx)
    neuralactivity_symmetric = neuralactivity_all[mask,:,:]
    neuralactivity_error = neuralactivity_all[np.logical_not(mask),:,:]


    # calculate decoders
    if recalculate:
        accuracies = np.zeros((n_timestamps,2,3,len(taskaspects),len(symmetrylabels)))
        coefs = np.zeros((n_timestamps,n_features,3,len(taskaspects),len(symmetrylabels)))
        for sx,(triallist,neuralactivity) in enumerate( zip( (triallist_all,triallist_symmetric, triallist_error),
                                                             (neuralactivity_all, neuralactivity_symmetric, neuralactivity_error) )  ):
            #

            for cx,task in enumerate(taskaspects):
                comparisonmask = []
                for comparison in comparisongroups[cx]:
                    if task != 'choice':
                        comparisonmask.append( \
                        ( (triallist['block'].isin(comparison[0])) | (len(comparison[0])==0) ) & \
                        ( (triallist['degree'].isin(comparison[1])) | (len(comparison[1])==0) ) & \
                        ( (triallist['freq'].isin(comparison[2])) | (len(comparison[2])==0) ) )
                    else:
                        comparisonmask.append( triallist['action'].isin(comparison) )

                targets = np.hstack([ triallist[columnnames[cx]][comparisonmask[i]]==class1value[cx] for i in [0,1]  ])
                predictors = np.vstack([neuralactivity[comparisonmask[i],:,:] for i in [0,1]])

                # fill out the targets axes, with dimensions = (observations,timepoints,features)
                targets = np.tile(targets, [1,len(times),1])
                targets = np.swapaxes(targets, 0, 2)

                accuracy,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'classification')
                
                # reshape to average over feature timewidth
                if wx>1:
                    coef = np.reshape(coef, [coef.shape[0], n_targets, wx, n_neurons, 3]).mean(axis=2)

                # (n_timestamps,train/test,stats,task,symmetry)
                accuracies[:,:,:,cx,sx] = accuracy.squeeze()
                # (n_timestamps,n_features,stats,task,symmetry)
                coefs[:,:,:,cx,sx] = coef.squeeze()

                

        pickle.dump((accuracies,coefs), open(cacheprefix+'symmetry/neural,VACC-all,symmetric-decoder,timecourse_%s.pck'%(dn), 'wb'))
    else:
        accuracies,coefs = pickle.load(open(cacheprefix+'symmetry/neural,VACC-all,symmetric-decoder,timecourse_%s.pck'%(dn), 'rb'))





    if doplot or globaldoplot:


        fig,ax = plt.subplots(len(symmetrylabels),len(taskaspects), figsize=(len(taskaspects)*8,len(symmetrylabels)*8) )
        taskcolors = ['navy','darkgreen','mediumvioletred','darkorange']


        for sx,symmetry in enumerate(symmetrylabels):
            for cx,task in enumerate(taskaspects):
                
                axs = ax[sx,cx]
                colors = ['grey',taskcolors[cx]]

                for k in [0,1]:
                    if k==0: continue

                    # (n_timestamps,train/test,stats,task,symmetry)
                    m = accuracies[:,k,0,cx,sx]
                    e = accuracies[:,k,2,cx,sx]
                    axs.plot(times, m, color=colors[k], lw=3, label=['train','test'][k])
                    axs.fill_between(times, m-e, m+e, color=colors[k], alpha=0.3)

                figs.setxt(axs)
                axs.set_ylim(0.45,1.05)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.5)


                if sx==0: axs.set_title(task)
                if cx==0: axs.set_ylabel('%s\naccuracy'%symmetry)
                # if cx==0 and sx==0: axs.legend(frameon=False)


        fig.suptitle(dn+' task variable decoders entire session vs. successfully symmetric behaviour trials'+\
                         '\n%d neurons, V:%d,%d, A:%d,%d trials'%(n_neurons,\
                          len(mask_contextuals[0]), np.sum(mask_contextuals[0]),\
                          len(mask_contextuals[1]), np.sum(mask_contextuals[1])  )  )

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'decode,VACC-all,symmetric,error_%s'%(dn)+ext)


















def behaviour_symmetry_context(dn,block,equalize=False,n_bootstrap=2,displaystatsonly=False):


    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot


    # load trials data
    g = preprocess.loadexperimentdata(dn, full=False, multimodalonly=True)
    g['block']+=1
    g['success'] = g['punish']==False
    blv,bla = preprocess.getorderattended(dn)


    # apply behaviour masks, to get a "clever"=True mask for all trials
    ma = 20
    ma_th = 0.5
    action = ['go','nogo']
    congruency = ['congruent','incongruent']
    mask_clevers_list = [[],[]] # holds context dependent list
    mask_clevers = []           # holds indexed by original trial order
    for c in congruency:
        for a in action:
            mask_clever_contexts = get_mask_cleverness(dn, ma=ma, ma_threshold=ma_th, visualfirst=False, go=a, congruency=c)
            for cx,mask_clever_context in enumerate(mask_clever_contexts):
                mask_clevers_list[cx].append(mask_clever_context)
            mask_clever = np.hstack( get_mask_cleverness(dn, ma=ma, ma_threshold=ma_th, visualfirst=False, go=a, congruency=c) )
            mask_clevers.append(mask_clever)
    
    mask_clevers_list = [ np.vstack(mask_clevers_list[cx]).T for cx in [0,1]]
    mask_clevers = np.vstack(mask_clevers).T


    mask_contextuals = [ np.prod(mask_clevers_list[cx],axis=1) for cx in [0,1] ]
    # display the number of trials that has above threshold movingaverage for all 4 combinations of congruenct and action
    trialslabel = '1st(%s):%d/%d, 2nd(%s):%d/%d'%(['V','A'][blv[1]==4],np.sum(mask_contextuals[0]), len(mask_contextuals[0]), ['A','V'][blv[1]==4], np.sum(mask_contextuals[1]), len(mask_contextuals[1]))
    print(dn, trialslabel)

    # mask variable for concatenated contexts, symmetric True for clever, antisymmetric True for clueless
    mask_symmetric_origin = np.bool8(np.prod(mask_clevers, axis=1))
    mask_antisymmetric_origin = np.logical_not(mask_symmetric_origin)
    print('trial counts concatenated sym: %d/%d, ant: %d/%d'%(np.sum(mask_symmetric_origin), len(mask_symmetric_origin), np.sum(mask_antisymmetric_origin), len(mask_antisymmetric_origin)))











    # decode 4 way same and cross from from symmetric and antisymmetric parts, and compare on randomized equal subsamples
    symmetrylabels = [['symmetric->symmetric','symmetric->antisymmetric'],['antisymmetric->antisymmetric','antisymmetric->symmetric']]

    
    
    # create the all trials, and just symmetrical symmetric behaviour trial list
    triallist = preprocess.loadexperimentdata(dn, multimodalonly=True)
    triallist['block']+=1
    
    # setup stimulus, and target labels for decoding:
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [blv[1], bla[1]]
    class1value = blv[1]
    classlabels = (triallist['block']==class1value)


    # collect neural data:
    wx = int((50*pq.ms/T['dt']).magnitude)  # width of the feature space in time for each neuron
    # wx = 0
    times = block.segments[0].analogsignals[0].times
    if wx>1: times = times[:-wx]
    n_timestamps = len(times)
    n_neurons = block.segments[0].analogsignals[0].shape[1]
    n_features = n_neurons
    n_targets = 1             # univariate targets

    allcomplextrials = [ [[2,4], [], []] ]
    # (observations,times,features)
    neuralactivity = np.array(preprocess.collect_stimulusspecificresponses(block, allcomplextrials)[0])
    if wx>1:   # create wider feature space of several subsequent timepoints
        neuralactivity = neph.slidingwindow(neuralactivity,wx)











    # calculate decoders
    if recalculate:

        np.random.seed(seed=0)

        accuracies = np.zeros((n_bootstrap,n_timestamps,2+2,3,2)) # (bootstraps,times,train-test-crosstrain-crosstest,stats,symmetrytrain)
        coefs = np.zeros((n_bootstrap,n_timestamps,n_features,3,2)) # (bootstraps,times,n_neurons,stats,symmetrytrain)

        for b in range(n_bootstrap):

            # randomly equalize train test trial set sizes and in both contexts
            if equalize:
                equalizelabel = 'equalized-'

                # collect first context symmetric, antisymmetric, then second context sy, as
                trials = [ np.copy(mask_contextuals[0]), 1-np.copy(mask_contextuals[0]), np.copy(mask_contextuals[1]), 1-np.copy(mask_contextuals[1]) ]
                trialnumbers = [ np.sum(trial) for trial in trials ]
                if b==0: print(trialnumbers)
                wh = np.argmin(trialnumbers)                
                minimumtrials = trialnumbers[wh]
                # print(dn,minimumtrials); return
                mask_equalizer = trials
                for sx in range(4):
                    if sx==wh: continue     # no need to resample the smallest fraction
                    indexlist = np.where(trials[sx]==1)[0]   # find active trials
                    discardindexlist = indexlist[np.random.permutation(len(indexlist))][:len(indexlist)-minimumtrials]   # subsample index of active trials
                    mask_equalizer[sx][discardindexlist] = 0        # do the subsampling by putting false into the mask for the not used trials
                # this will be (trials,symmetrytype)
                mask_equalizer_symmetric = np.hstack([mask_equalizer[0],mask_equalizer[2]]).astype(np.bool8)
                mask_equalizer_antisymmetric = np.hstack([mask_equalizer[1],mask_equalizer[3]]).astype(np.bool8)


                mask_symmetric = mask_symmetric_origin * mask_equalizer_symmetric
                mask_antisymmetric = mask_antisymmetric_origin * mask_equalizer_antisymmetric

            else:
                equalizelabel = ''

        





            
            
            # for sx,(triallist,neuralactivity) in enumerate( zip( triallist_comparisons, neuralactivity_comparisons)  ):
            for sx,(mask_train,mask_test) in enumerate(zip((mask_symmetric,mask_antisymmetric), (mask_antisymmetric, mask_symmetric))):
                    if b==0: print('crossdecoding',symmetrylabels[sx])
                    if sx==1: print('bootstrap %d/%d'%(b+1,n_bootstrap))

                    targets = classlabels
                    predictors = neuralactivity



                    # fill out the targets axes, with dimensions = (observations,timepoints,features)
                    targets = np.tile(targets, [1,len(times),1])
                    targets = np.swapaxes(targets, 0, 2)

                    accuracy,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'crossclassification', 'loo',\
                                                                        crossindicestrain=mask_train,crossindicestest=mask_test)
                    
                    # reshape to average over feature timewidth
                    if wx>1:
                        coef = np.reshape(coef, [coef.shape[0], n_targets, wx, n_neurons, 3]).mean(axis=2)


                    accuracies[b,:,:,:,sx] = accuracy.squeeze() # (bootstrap, times,train-test-crosstrain-crosstest,stats,symmetrytrain)
                    coefs[b,:,:,:,sx] = coef.squeeze() # (bootstrap, times,n_neurons,stats,symmetrytrain)




        # do for all trials as control    
        print('decoding all trials')
        targets = classlabels
        predictors = neuralactivity
        targets = np.tile(targets, [1,n_timestamps,1])
        targets = np.swapaxes(targets, 0, 2)
        accuracy,coef = nedi.get_linearregressionmultivariate(predictors,targets,times,'classification')
        if wx>1: coef = np.reshape(coef, [coef.shape[0], 1, wx, n_neurons, 3]).mean(axis=2)
        accuracyall = accuracy.squeeze()
        coefsall = coef.squeeze()




        pickle.dump((accuracies,coefs,accuracyall,coefsall), open(cacheprefix+'symmetry/neural,context-%ssymmetric,antisymmetric,cross-decoder-loo,boots,timecourse_%s.pck'%(equalizelabel,dn), 'wb'))
    else:
        if equalize: equalizelabel = 'equalized-'
        else: equalizelabel = ''

        accuracies,coefs,accuracyall,coefsall = pickle.load(open(cacheprefix+'symmetry/neural,context-%ssymmetric,antisymmetric,cross-decoder-loo,boots,timecourse_%s.pck'%(equalizelabel,dn), 'rb'))
        n_bootstrap = accuracies.shape[0]

    if n_bootstrap>1: loobootstraplabel = 'bootstrap-lo(p)o'
    else: loobootstraplabel = ''



    accuracies_m = accuracies.mean(axis=0)
    accuracies_e = accuracies.std(axis=0)/np.sqrt(n_bootstrap)



    if doplot or globaldoplot:

        # get baseline shuffle:
        n_resample = 10
        chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))
        chances_reduced = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d,reduced-%s.pck'%(n_resample,continuous_method),'rb'))

        fig,ax = plt.subplots(1,3, figsize=(3*10,1*8) )

        # plot all trials
        axs = ax[0]
        colors = ['grey','mediumvioletred']
        for k in [0,1]:
            if k==0: continue

            # (n_timestamps,train/test,stats,task,symmetry)
            m = accuracyall[:,k,0]
            e = accuracyall[:,k,2]
            axs.plot(times, m, color=colors[k], lw=3, label='entire session')
            axs.fill_between(times, m-e, m+e, color=colors[k], alpha=0.3)
        figs.setxt(axs)
        axs.set_ylim(0.45,1.05)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,0.5)
        figs.plottoaxis_chancelevel(axs,chances[dn])

        axs.legend(frameon=False)

        axs.set_ylabel('accuracy')
        
        
        
        # plot crossdecoding
        symmetrycolors = [['rebeccapurple','gold'],['darkorange','fuchsia']]
        for rx in range(2):             # train subset
            for sx in range(2):         # crosstest subset (same or cross)
                
                axs = ax[1+rx]
                colors = ['grey',symmetrycolors[rx][sx]]

                for k in [0,1]:
                    if k==0: continue             # exclude test

                    # (n_timestamps,train/test,stats,task,symmetry)
                    m = accuracies_m[:,k+sx*2,0,rx]
                    e =  np.sqrt(accuracies_m[:,k+sx*2,2,rx]**2 + accuracies_e[:,k+sx*2,2,rx]**2)
                    axs.plot(times, m, color=colors[k], lw=2, label=[None,symmetrylabels[rx][sx]][k])
                    axs.fill_between(times, m-e, m+e, color=colors[k], alpha=0.3)

            figs.plottoaxis_chancelevel(axs,chances_reduced[dn])
            figs.setxt(axs)
            axs.set_ylim(0.45,1.05)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)

            axs.legend(frameon=False)





        fig.suptitle(dn+' task variable decoders entire session vs. successfully symmetric behaviour trials'+\
                         '\n%d neurons, %s trials %s %s'%(n_neurons, trialslabel, equalizelabel[:-1], loobootstraplabel  )  )

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'decode-loo,boots,context-%ssymmetric,antisymmetric,cross_%s'%(equalizelabel,dn)+ext)







def moving_average_trials(success, masks, ma=20):
    # success is a n_trials array of correct = true, incorrect = false
    # mask is  (n_trials, {conggo,inconggo,congnogo,incongnogo}) array of booleans
    n_trials,n_masks = masks.shape
    slidingwindowperformance = np.zeros(masks.shape)
    for t in np.arange(n_trials):
        sl = slice( max(0,t-ma//2), min(t+ma//2,n_trials) )
        # get the success rate for the trials within the ma window with the index slice
        # but separately for each mask
        for imask in range(n_masks):
            # print(t,imask,sl,(((success==1) & masks[:,imask])[ sl ]))
            slidingwindowperformance[t,imask] = ((success==1) & masks[:,imask])[ sl ].sum()/sum(masks[:,imask][sl])
    return slidingwindowperformance


def consistent_masks(stimuli, success, relevant_idx=0):
    mask_congruent =  stimuli[:,0,:]==stimuli[:,1,:]
    mask_go = stimuli[:,relevant_idx,:]==1
    masks =  np.array([mask_congruent & mask_go, ~mask_congruent & mask_go, mask_congruent & ~mask_go, ~mask_congruent & ~mask_go])
    masks = masks.transpose(1,0,2)
    return masks


def segment_ma_consecutive(x):
    # find consecutive segments of x=1..
    # return the length of all segments
    consecutives = []
    if x[0]==1: flag = True
    else: flag = False
    c = 0
    for i in range(len(x)):
        if x[i]==1 and ~flag: flag = True
        if x[i]==1 and flag: c += 1
        if x[i]==0 and flag: consecutives.append(c); c = 0; flag = False
        if i==len(x)-1 and flag: consecutives.append(c)
    return consecutives


def  behaviour_generate_movingaveragestatistics():
    # generate trials: single context, visual and audio stimulus, choice with a bias term

    doplot = 1 or globaldoplot

    def model_param(stimulus, param, model):
        # bias model, get stimuli, make choice on beta parameter
        # stimulus: (n_trials, n_stimuli, n_repeats)
        p = np.zeros(stimulus.shape)
        if model=='mean':
            p[:] = param
        elif model=='bias':
            p[stimulus==1] = 1+param
            p[stimulus==0] = param
        elif model=='lapse':
            p[stimulus==1] = 1-param
            p[stimulus==0] = param
        np.clip(p,0,1,out=p)
        return p



    labels_models = ['mean','bias','lapse']
    n_params = 5
    param_ranges = [np.linspace(0,1,n_params), np.linspace(-1,1,n_params), np.linspace(0,0.5,n_params)]
    n_models = len(param_ranges)
    n_trials = 70
    n_modality = 2
    relevant_idx = 0
    n_repeats = 100000
    n_trialtypemasks = 4

    




    success_mas = []
    consecutive_mas = []
    stimuli = np.random.binomial(n=1,p=0.5,size=(n_trials,n_modality,n_repeats))
    for mx in range(n_models):
        # if mx>0: continue
        #  calculate successes and moving averages
        choices = np.zeros((n_trials,n_repeats,n_params))
        success = np.zeros((n_trials,n_repeats,n_params))
        success_ma = np.zeros((n_trials,n_trialtypemasks+1,n_repeats,n_params))
        for ux,u in enumerate(param_ranges[mx]):
            # if ux!=3:  continue
            ps = model_param(stimuli[:,relevant_idx,:], u, labels_models[mx])
            choices[:,:,ux] = np.random.binomial(n=1,p=ps,size=(n_trials,n_repeats))
            success[:,:,ux] = (stimuli[:,relevant_idx,:]==choices[:,:,ux])

            # get congruency from stimuli 1 and 2 equal or opposite
            masks = consistent_masks(stimuli, success)
            for rx in range(n_repeats):
                # success[masks[:,0,rx],rx,ux] = 1
                # success[masks[:,2,rx],rx,ux] = 1
                success_ma[:,:4,rx,ux] = moving_average_trials(success[:,rx,ux], masks[:,:,rx], ma=20)

        #  calculate consistent trials: where all four are logically true
        success_ma[:,4,:,:] = np.logical_and.reduce(success_ma[:,0:4,:,:]>0.5,axis=1)
        success_mas.append(success_ma)


        # count lengths of consecutive successes blocks
        consecutive_ma = []
        for ux,u in enumerate(param_ranges[mx]):
            c_ma = []
            for rx in range(n_repeats):
                # if ux!=3:  break
                c_ma.append( segment_ma_consecutive(success_mas[mx][:,4,rx,ux]) )
            # if mx==1 and ux==2: print(mx,ux,c_ma)
            consecutive_ma.append(c_ma)
        consecutive_mas.append(consecutive_ma)




    # count occurrences and lengths of consecutive successes

    num_ntrials_successes = np.zeros((n_models,n_params,n_trials))   # number of trials with t successes
    consecutive_counts = np.zeros((n_models,n_params,n_repeats,n_trials))

    for mx in range(n_models):
        # if mx>0: continue
        for ux in range(n_params):
            # if ux!=3: continue
            for rx in range(n_repeats):
                t = (success_mas[mx][:,4,rx,ux]).sum().astype(np.int32)
                if t>0: # otherwise the -1 index will be the last, and gives false count
                    num_ntrials_successes[mx,ux,t-1] += 1
                for c in consecutive_mas[mx][ux][rx]:         # use unique for single occurrence, and without unique for all occurrences
                    consecutive_counts[mx,ux,rx,c-1] += 1
                


    # probabilities
    cutpoint = 10       # criteria for mouse inclusion
    # probabilities of having n success_mas dependent on n
    prob_ntrials_successes = num_ntrials_successes/n_repeats

    # probability of having at least one consecutive successes of length n (from 1 to cutpoint through all trials)
    print(consecutive_counts.shape)
    print(consecutive_counts[1,2,:,:].sum(axis=0))
    # calculate the probability of having exactly 1, 2, etc. number of consecutive successes of length n
    prob_exactly_consecutive_length = np.zeros((n_models,n_params,n_trials,n_trials))
    for mx in range(n_models):
        # if mx>0: continue
        for ux in range(n_params):
            # if ux!=3: continue
            for t in range(n_trials):
                for rx in range(n_repeats):
                    c = consecutive_counts[mx,ux,rx,t].astype(np.int32)  # number of t length occurrences reported in repetition rx
                    if c>0:
                        prob_exactly_consecutive_length[mx,ux,t,c-1] += 1     # register, that at t length, there was c number of occurrences
    prob_atleastone_consecutive_length = (prob_exactly_consecutive_length.sum(axis=3))/n_repeats



    # mouse data
    dn = 'MT020_2'
    g = preprocess.loadexperimentdata(dn, full=False, multimodalonly=True)
    g['block']+=1
    g['success'] = g['punish']==False
    blv,bla = preprocess.getorderattended(dn)
    h = g[g['block']==blv[1]]
    stimuli_mouse = np.vstack([h['degree'].values==45.0, h['freq'].values==5000.0]).astype(np.int32).T
    choices_mouse = h['action'].values.astype(np.int32)
    success_mouse = h['success'].values.astype(np.int32)
    n_trials_mouse = len(choices_mouse)
    success_ma_mouse = np.zeros((n_trials_mouse,n_trialtypemasks+1))
    masks_mouse = consistent_masks(stimuli_mouse[:,:,np.newaxis], success_mouse[:,np.newaxis, np.newaxis])
    success_ma_mouse[:,:4] = moving_average_trials(success_mouse, masks_mouse[:,:,0], ma=20)
    success_ma_mouse[:,4] = np.logical_and.reduce(success_ma_mouse[:,0:4]>0.5,axis=1)
    consecutive_ma_mouse = segment_ma_consecutive(success_ma_mouse[:,4])

    print('mouse, consecutive blocks in visual context:',consecutive_ma_mouse)

    
    pickle.dump((n_trials,cutpoint,prob_ntrials_successes,prob_atleastone_consecutive_length,consecutive_ma_mouse),
                open(cacheprefix+'behaviour/probabilityconsistentperiods-n%d_%s'%(n_repeats,dn),'wb'))


    # plots
    if doplot:
    
        if 0:             # single sample ma-s 
            fig,axs = plt.subplots(n_models,n_params)
            for mx in range(n_models):
                for ux in np.arange(n_params):
                    ax = axs[mx,ux]
                    for kx in np.arange(n_trialtypemasks):
                        ax.plot(np.arange(n_trials)+1,success_mas[mx][:,kx,0,ux],
                                color=['lightseagreen','lightseagreen','red','red'][kx],ls=['-','--','-','--'][kx])
                    
                    
                    for expl in np.arange(n_trials)[success_mas[mx][:,4,0,ux]>0.5]:
                        ax.fill_between([expl-0.49999, expl+0.5],[-0.2,-0.2],[-0.1,-0.1],color='rebeccapurple',alpha=1)
                    

                    ax.set_title('%s=%4.2f'%(['<choice>','$\\beta$','$\\lambda$'][mx],param_ranges[mx][ux]),fontsize='small')
                    ax.set_ylim(-0.25,1.05)
                    ax.set_xticks([1,35,70])
                    if mx<n_models-1: ax.set_xticklabels([])
                    if mx==n_models-1: ax.set_xlabel('trial')
                    if ux>0: ax.set_yticklabels([])
                    if ux==0: ax.set_ylabel(labels_models[mx]+'\nfraction corr.')
                    figs.plottoaxis_chancelevel(ax,0.5)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
        

        if 0:            # probabilities of consistent consecutive block
            fig,axs = plt.subplots(n_models,n_params)
            for mx in range(n_models):
                for ux in np.arange(n_params):
                    ax = axs[mx,ux]

                    ax.plot(np.arange(n_trials)+1,prob_ntrials_successes[mx,ux,:], color='dodgerblue',lw=0.7,alpha=0.8)
                    ax.plot(np.arange(n_trials)+1,[prob_ntrials_successes[mx,ux,c:].sum() for c in np.arange(n_trials)], color='dodgerblue',lw=2,alpha=0.8)
                    ax.scatter(cutpoint,prob_ntrials_successes[mx,ux,cutpoint-1:].sum(), color='dodgerblue',alpha=0.8)

                    ax.plot(np.arange(n_trials)+1,prob_atleastone_consecutive_length[mx,ux,:], color='rebeccapurple',lw=0.7,alpha=0.8)
                    ax.plot(np.arange(n_trials)+1,[prob_atleastone_consecutive_length[mx,ux,c:].sum() for c in np.arange(n_trials)], color='rebeccapurple',lw=2,alpha=0.8)
                    ax.scatter(cutpoint,prob_atleastone_consecutive_length[mx,ux,cutpoint-1:].sum(), color='rebeccapurple',alpha=0.8)

                    ax.vlines(cutpoint,0,1,ls='--',lw=1,color='black',alpha=0.2)

                    ax.text(cutpoint+2,1.0,'$P(N\\geq 10)$ = %6.3f'%(prob_ntrials_successes[mx,ux,cutpoint:].sum()),
                            fontsize='xx-small',ha='left',va='top',color='dodgerblue')
                    ax.text(cutpoint+2,0.9,'$P(L\\geq 10)$ = %6.3f'%(prob_atleastone_consecutive_length[mx,ux,cutpoint:].sum()),
                            fontsize='xx-small',ha='left',va='top',color='rebeccapurple')

                    ax.set_ylim(-0.05,1.05)


                    ax.set_xticks([1,35,70])
                    if mx<n_models-1: ax.set_xticklabels([])
                    if mx==n_models-1: ax.set_xlabel('number of trials')
                    if ux>0: ax.set_yticklabels([])

                    ax.set_title('%s=%4.2f'%(['<choice>','$\\beta$','$\\lambda$'][mx],param_ranges[mx][ux]),fontsize='small')

                    if ux==0: ax.set_ylabel(labels_models[mx]+'\nprobability')

                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)




        if 1:
            fig,axs = plt.subplots(1,1)
            ax = axs  # production figure R1
            mx = 0    # mean bias model
            ux = 3    # mean bias = 0.75

            ax.plot(np.arange(n_trials)+1,[prob_ntrials_successes[mx,ux,c:].sum() for c in np.arange(n_trials)], color='dodgerblue',lw=2,alpha=0.8)
            ax.plot(np.arange(n_trials)+1,[prob_atleastone_consecutive_length[mx,ux,c:].sum() for c in np.arange(n_trials)], color='rebeccapurple',lw=2,alpha=0.8)

            num_success_mouse = sum(consecutive_ma_mouse)
            ax.scatter(cutpoint,prob_ntrials_successes[mx,ux,cutpoint-1:].sum(), color='dodgerblue',alpha=0.8)
            ax.scatter(num_success_mouse,prob_ntrials_successes[mx,ux,num_success_mouse-1:].sum(), color='dodgerblue',alpha=0.8)
            ax.scatter(cutpoint,prob_atleastone_consecutive_length[mx,ux,cutpoint-1:].sum(), color='rebeccapurple',alpha=0.8)
            ax.scatter(cutpoint,prob_atleastone_consecutive_length[mx,ux,num_success_mouse-1:].sum(), color='rebeccapurple',alpha=0.8)

            ax.vlines(cutpoint,0,1,ls='--',lw=1,color='black',alpha=0.2)

            ax.text(0.1,1.0,'$P(N\\geq 10)$ = %6.5f'%(prob_ntrials_successes[mx,ux,cutpoint:].sum()),
                    fontsize='xx-small',ha='left',va='top',color='dodgerblue',transform=ax.transAxes)
            ax.text(0.1,0.9,'$P(N\\geq 20)$ = %6.5f'%(prob_ntrials_successes[mx,ux,num_success_mouse:].sum()),
                    fontsize='xx-small',ha='left',va='top',color='dodgerblue',transform=ax.transAxes)
            ax.text(0.1,0.8,'$P(L\\geq 10)$ = %6.9f'%(prob_atleastone_consecutive_length[mx,ux,cutpoint:].sum()),
                    fontsize='xx-small',ha='left',va='top',color='rebeccapurple',transform=ax.transAxes)
            ax.text(0.1,0.7,'$P(L\\geq 20)$ = %6.9f'%(prob_atleastone_consecutive_length[mx,ux,num_success_mouse:].sum()),
                    fontsize='xx-small',ha='left',va='top',color='rebeccapurple',transform=ax.transAxes)

            # # combination of lengths we found
            # p = np.prod(np.array([prob_atleastone_consecutive_length[mx,ux,k-1:].sum() for k in consecutive_ma_mouse]))
            # ax.text(0.1,0.8,'$P(L\\geq 10)\cdot P(L\\geq 7)\cdot P(L\\geq 3)$ = %6.9f'%(p),
            #         fontsize='xx-small',ha='left',va='top',color='red',transform=ax.transAxes)
            # # probability of a single 17 length
            # p = prob_atleastone_consecutive_length[mx,ux,sum(consecutive_ma_mouse)-1:].sum()
            # ax.text(0.1,0.7,'$P(L\\geq 10+7+3)$ = %6.9f'%(p),
            #         fontsize='xx-small',ha='left',va='top',color='fuchsia',transform=ax.transAxes)


            ax.set_xticks([1,35,70])
            ax.set_xlabel('number of trials')

            ax.set_ylim(0,0.12)
            ax.set_yticks([0,0.05,0.1])
            ax.set_ylabel('probability')

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)



    return