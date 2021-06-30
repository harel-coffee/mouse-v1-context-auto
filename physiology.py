#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 10:16:45 2020

@author: mahajnal
"""










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




