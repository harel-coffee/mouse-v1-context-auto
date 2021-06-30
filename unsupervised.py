#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 10:18:07 2020

@author: mahajnal
"""


from config import *

import preprocess















# *************************************************************************
#                               dimensionreductions










def tcafactors(dn,block):
    
    method = 'tca'
    R = 10            # number of factors
    
    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot

    times = block.segments[0].analogsignals[0].times

    # collect neural data:
    data_alltrials  = [    [1,2,3,4],  [],    []      ],
    
    responses = preprocess.collect_stimulusspecificresponses(block,data_alltrials,'and')[0]    
    print(np.array(responses).shape, type(responses[0]))
    
    # collect runspeed data:
    runspeeds = preprocess.collect_stimulusspecificresponses(block,data_alltrials,'and',s=1)[0]
    runspeeds = np.array(runspeeds).squeeze()[:,T['stimstart_idx']:T['stimend_idx']].mean(axis=1)
    n_runbins = 100
    maxrunspeed = 60
    runspeed_idx = np.digitize(runspeeds, bins=np.linspace(0,maxrunspeed,n_runbins,endpoint=False))



    if method[:3] == 'tca':
        
        if recalculate:
            # tensor component analysis
            if method[-2:]=='nn': nn = True
            else: nn = False
            
            H = nedi.decompose_tca(np.array(responses),R,nn)
            pickle.dump(H,open(cacheprefix+'tca+/%s-factors%d_trials-H_%s'%(method,R,dn),'wb'))
            
        else:
            H = pickle.load(open(cacheprefix+'tca+/%s-factors%d_trials-H_%s'%(method,R,dn),'rb'))
        
        U = H.factors





    ixv45,ixv135,ixa45,ixa135,ixv5000,ixv10000,ixa5000,ixa10000 = preprocess.getstimulusidents(dn,block)
    ixblv, ixbla = preprocess.getcontextidents(dn,block)
    ixnarrow,ixbroad,ixsuperf,ixinput,ixdeep5,ixdeep6 = preprocess.getcelltypes(block)
    ixhit,ixmiss,ixcorrrej,ixfal = preprocess.assessbehaviouralperformance(dn)
    





    # assess the performance of each tensor component by decoders to specific experimental aspect

    decodercolors = ['navy','darkgreen','mediumvioletred','darkorange','darkred']
    decodernames = ['visual','audio','context','choice','withhold perf.']
    # find the performance of each factors by decoders
    classindices = [  [np.hstack( (ixa45,ixv45) ), np.hstack( (ixa135,ixv135) )], \
                      [np.hstack( (ixa5000,ixv5000) ), np.hstack( (ixa10000,ixv10000) )], \
                      [  ixblv, ixbla  ], \
                      [np.hstack( (ixhit,ixfal) ), np.hstack( (ixcorrrej,ixmiss) ) ], \
#                          [ ixhit, ixmiss ],\          # this cannot be done, as not enough misses for most mice
                      [ ixcorrrej, ixfal ]  ]     # visual, audio, context, choice, behavioural
    nct = len(classindices)
        
    powers = np.empty( (U.rank,nct) )                           # factors, taskaspects
    bestfactors = np.empty( (U.rank,nct), dtype='int16' )       # factors,taskaspects
    for ix,ixr in enumerate( range(U.rank) ):
        for k in range(nct):
            Z1 = U[0][classindices[k][0],ixr]      # 1st class neural latent
            Z2 = U[0][classindices[k][1],ixr]      # 2nd class neural latent
            Y1 = np.zeros(Z1.shape,dtype='int16')
            Y2 = np.ones(Z2.shape,dtype='int16')
            Z = np.hstack((Z1,Z2)).reshape(-1,1)
            Y = np.hstack((Y1,Y2))
            
            if len(Y1)>2 and len(Y2)>2:
                powers[ix,k] = nedi.classifierfortcalatents(Z,Y)
    for k in range(nct):
        bestfactors[:,k] = np.argsort(powers[:,k])[::-1]







    # plotting everything, each factor and tensor mode


    if doplot:
        runcolorlist = plt.cm.viridis( np.linspace(0, 1, n_runbins+1) )
        print(len(runcolorlist),max(runspeed_idx))

        subplotsize = 8
        fig,ax = plt.subplots(R,9,figsize=(9*subplotsize,R*subplotsize))
            
        for ix,ixr in enumerate( range(U.rank) ):
            
            axs = ax[ix,0]
            for k in range(nct):
                axs.bar((k-(nct-1)/2)/(nct+1),powers[ixr,k],width=1/(nct+1),color=decodercolors[k])
            if ix==0: axs.legend(decodernames)
            axs.xaxis.set_ticks([])
            axs.set_ylim(0.5,1.0)
            axs.set_ylabel('factor %d'%(ixr+1))#,factortitles[ix]))
            
    
            for ixtm, tensormode in enumerate(U):
    
                
                if ixtm==0:
                    axs = ax[ix,1]
                    axs.plot(ixv5000,tensormode[ixv5000, ixr],'o', color='sienna', label='a vis  5000')
                    axs.plot(ixv10000,tensormode[ixv10000, ixr],'o', color='darkorange', label='a vis 10000')
                    axs.plot(ixa5000,tensormode[ixa5000, ixr],'o', color='darkgreen', label='a aud  5000')
                    axs.plot(ixa10000,tensormode[ixa10000, ixr],'o', color='lime', label='a aud 10000')
                    axs.set_xlim([0,len(responses)+1])
                    axs.set_ylim([tensormode[:,ixr].min(),tensormode[:,ixr].max()])
                    figs.plottoaxis_chancelevel(axs)
    
                    axs = ax[ix,2]
                    axs.plot(ixv45,tensormode[ixv45, ixr],'bo', label='a vis  45')
                    axs.plot(ixv135,tensormode[ixv135, ixr],'co', label='a vis 135')
                    axs.plot(ixa45,tensormode[ixa45, ixr],'mo', label='a aud  45')
                    axs.plot(ixa135,tensormode[ixa135, ixr],'ro', label='a aud 135')
                    axs.set_xlim([0,len(responses)+1])
                    axs.set_ylim([tensormode[:,ixr].min(),tensormode[:,ixr].max()])
                    figs.plottoaxis_chancelevel(axs)
                    
                    axs = ax[ix,3]
                    axs.plot(ixhit,tensormode[ixhit, ixr],'o', color='rebeccapurple', label='hit')
                    axs.plot(ixmiss,tensormode[ixmiss, ixr],'o', color='plum', label='miss')
                    axs.plot(ixcorrrej,tensormode[ixcorrrej, ixr],'o', color='goldenrod', label='correct rejection')
                    axs.plot(ixfal,tensormode[ixfal, ixr],'o', color='yellow', label='false alarm')
                    axs.set_xlim([0,len(responses)+1])
                    axs.set_ylim([tensormode[:,ixr].min(),tensormode[:,ixr].max()])
                    figs.plottoaxis_chancelevel(axs)
                    
                    axs = ax[ix,4]
                    for n in range(tensormode.shape[0]):
                        axs.plot(n,tensormode[n, ixr],'o', color=runcolorlist[runspeed_idx[n]])
                    axs.set_xlim([0,len(responses)+1])
                    axs.set_ylim([tensormode[:,ixr].min(),tensormode[:,ixr].max()])
                    figs.plottoaxis_chancelevel(axs)
    
                if ixtm==1:
                    axs = ax[ix,5]
                    axs.plot(times,tensormode[:, ixr],'r',linewidth=2)
                    axs.set_xlim([T['starttime']+0*pq.ms,T['endtime']-0*pq.ms])
                    axs.set_ylim(-1.2,1.2)
                    figs.plottoaxis_chancelevel(axs)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    x = neph.removefrequencies(tensormode[:,ixr],[1.2],1000./T['dt'],filtertype='hp')
                    u_fx, u_Px = sp.signal.welch(x=x,fs=1000./T['dt'])
                    axs = ax[ix,6]
                    ax2 = axs.twinx()
                    axs.plot(u_fx[:20], u_Px[:20] ,color='mediumvioletred',linewidth=2)
                    axs.yaxis.set_ticks([])
                    #log power scale:
                    ax2.plot(u_fx[:20], u_Px[:20] ,color='darkgoldenrod',linewidth=2)
                    ax2.set_yscale('log')
                    ax2.yaxis.set_ticks([])
                    
                    
    
                elif ixtm==2:
                    continue
                    axs = ax[ix,7]
                    axs.bar(tensormode[:, ixr], color='lightcoral')
                    axs = ax[ix,8]
                    axs.bar(tensormode[:, ixr], color='mediumblue')
                    
                    continue
                    axs = ax[ix,7]
                    axs.bar(ixnarrow, tensormode[ixnarrow, ixr], color='lightcoral', label='narrow')
                    axs.bar(ixbroad, tensormode[ixbroad, ixr], color='lawngreen', label='broad')
                    axs.bar(tensormode.shape[0]+4, np.abs(tensormode[ixnarrow, ixr]).mean(axis=0), color='red', label='mean abs narrow')
                    axs.bar(tensormode.shape[0]+5, np.abs(tensormode[ixbroad, ixr]).mean(axis=0), color='darkgreen', label='mean abs broad')
                    
                    axs = ax[ix,8]
                    if len(ixdeep5)>0:
                        axs.bar(ixdeep5, tensormode[ixdeep5, ixr], color='mediumblue', label='deep5')
                        axs.bar(tensormode.shape[0]+4, np.abs(tensormode[ixdeep5, ixr]).mean(axis=0), color='royalblue', label='mean abs deep5')
                    if len(ixdeep6)>0:
                        axs.bar(ixdeep6, tensormode[ixdeep6, ixr], color='navy', label='deep6')
                        axs.bar(tensormode.shape[0]+5, np.abs(tensormode[ixdeep6, ixr]).mean(axis=0), color='blue', label='mean abs deep6')
                    if len(ixinput)>0:
                        axs.bar(ixinput, tensormode[ixinput, ixr], color='darkviolet', label='input')
                        axs.bar(tensormode.shape[0]+6, np.abs(tensormode[ixinput, ixr]).mean(axis=0), color='magenta', label='mean abs input')
                    if len(ixsuperf)>0:
                        axs.bar(ixsuperf, tensormode[ixsuperf, ixr], color='deepskyblue', label='superf')
                        axs.bar(tensormode.shape[0]+7, np.abs(tensormode[ixsuperf, ixr]).mean(axis=0), color='cyan', label='mean abs superf')
    
    
    
    
                if 1:
                    ax[R-1,1].set_xlabel('trial number')
                    ax[R-1,2].set_xlabel('trial number')
                    ax[R-1,3].set_xlabel('trial number')
                    ax[R-1,4].set_xlabel('trial number')
                    ax[R-1,5].set_xlabel('time along trial [s]')
                    ax[R-1,6].set_xlabel('signal frequency spectrum [Hz]')
                    ax[R-1,7].set_xlabel('neuron id (spike width)')
                    ax[R-1,8].set_xlabel('neuron id (tissue layer)')
                
                    ax[0,1].set_title('trials audio stimulus')
                    ax[0,2].set_title('trials visual stimulus')
                    ax[0,3].set_title('behavioural response')
                    ax[0,4].set_title('running speed')
                    ax[0,5].set_title('neural trajectory')
                    ax[0,6].set_title('n. traj. spectrum >1.2Hz')
                    ax[0,7].set_title('participation: neuron types')
                    ax[0,8].set_title('participation: layers')
                    
                    ax[0,1].legend();ax[0,2].legend();ax[0,3].legend()
                    
                    ax[0,7].legend();ax[0,8].legend()
            
            fig.suptitle(dn+' TCA, %s, %dms width, z-scored'%(continuous_method,T['bin']))
    
    
    
    
        save = 1
        if save:
            fig.savefig(resultpath+'%s-zsc,stim+bhv+r_l%d,dt%dms,window%dms_%s'%(method,R,T['dt'],T['bin'],dn) + '.png')
    
    
    if 0:
        fig,ax = plt.subplots(1,4,figsize=(36,4))
        if dn=='DT017':    S = [-1,-1, 1,-1]  # for DT017
        elif dn=='DT030':  S = [ 1,-1, 1, 1]  # for DT030

        for nct,taskaspects in enumerate(decodernames[:4]):     # concern only visual, audio, context and choice
            if nct<0:
                X = [  U[1][:,bestfactors[bfi,nct]] for bfi in [1,2]   ]
                X = np.mean(X ,axis=0)
            else:
                X = S[nct]*U[1][:,bestfactors[0,nct]] # trajectories only, and the best factor
                
            axs = ax[nct]
            axs.plot(times,X,linewidth=3,color=decodercolors[nct])
            axs.set_xlim([T['starttime']+0*pq.ms,T['endtime']-0*pq.ms])
            axs.set_ylim(-0.15,1.15)
            figs.plottoaxis_chancelevel(axs)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.setxt(axs)
            axs.set_yticks([0,1]);
#            axs.set_xlabel('[ms]')
            axs.set_title(decodernames[nct])
            axs.legend( ['factor %d, %s dec.acc. %4.2f'%(bestfactors[0,nct]+1, decodernames[nct], powers[bestfactors[0,nct],nct]) ] )
            if nct==0: axs.set_ylabel('trajectory tensor mode coeff.')
#        fig.suptitle(dn)
        
        save = 0
        if save:
            fig.savefig(resultpath+'3A-%s-tca,VACC-trajectories'%dn+ext)



































# ***************************
#            BAYESIAN





def latent_behaviour(dn):
    
    recalculate = 0 or globalrecalculate
    doplot = 1 or globaldoplot
    
    # data keys = ['start','duration','degree','block','freq','water','punish']
    data = preprocess.loadexperimentdata(dn)
    data = data.loc[data['block']==1,:]
    
    data['lick'] = ( (data['water']) &  (1-data['punish']) ) | \
                   ( (1-data['water']) &  (data['punish']) )
    
    n_trials = len(data)
    # print(data)
    
    
    # build model variables:
    # Y will contain the observables: stimuli and the reinforcer reward (reward is always last)
    # X will contain the task observables for prediction checking, not used during training
    
    # all variables, to learn correlations:
    # Y = np.array([ data['degree']==45, data['degree']==135, data['water']==True ], dtype=np.int16)
    # X = np.array([ data['lick']==True, data['lick']==False ], dtype=np.int16)
    
    # mutually exclusive 2-class setup just to learn the weight matrix
    Y = np.array([ data['degree']==45,data['water']==True, ], dtype=np.int16)
    X =  np.array([ data['lick']==True,], dtype=np.int16)
    
    if recalculate:
        Y_,X_,W_,b_ = nebay.latentcauses(X,Y)
        pickle.dump((Y_,X_,W_,b_),open(cacheprefix+'bayesian/latent,taskonly,visual_%s.pck'%(dn),'wb'))
    else:
        Y_,X_,W_,b_ = pickle.load(open(cacheprefix+'bayesian/latent,taskonly,visual_%s.pck'%(dn),'rb'))
        
    print(Y_.shape)
    print(X_.shape)
    print(W_)
    print(b_)
    
    
    if doplot:
        fig,ax = plt.subplots(2,3,figsize=(32,16))

        axs = ax[0,0]
        axs.plot(Y.T,'o')
        axs.set_ylabel('stimuli observation')

        axs = ax[1,0]
        axs.plot(Y_.T,'o')
        axs.set_ylabel('stimuli predictions')

        axs = ax[0,1]
        axs.plot(X.T,'o')
        axs.set_ylabel('choice observation')

        axs = ax[1,1]
        axs.plot(X_.T,'o')
        axs.set_ylabel('latent task choice')


        axs = ax[0,2]
        # axs.plot(X_[0,:]!=X_[1,:],'o',color='darkorange')
        axs.set_ylabel('latent 1 not = latent 2')

        axs = ax[1,2]
        axs.plot(np.mean((X_-X)**2,axis=0),'or')
        axs.set_ylabel('latent error X_-X')







    





def variationallatentgaussianprocess(dn,block,n_factors=5):

    counting = 'simulation'

    recalculate = 0 or globalrecalculate
    doplot = 1 or globaldoplot

    gp_char_timescale = ( 300*pq.ms ).rescale(block.segments[0].spiketrains[0][0].units).magnitude

    n_factors = 6
    
    # trials = [{'ID': i, 'y': y} for i, y in enumerate(sample['y'])]  # make trials
    # trials = [{'ID': tx, 'y': trial.spiketrains[0]} for tx,trial in enumerate(block.segments)]
    
    
    # vLGP requires dictionary format
    # exclude receptive field characterization 5th block
    
    # trials = [{'ID': tx, 'y': trial.analogsignals[0].magnitude, 'block': trial.annotations['block']}\
    #           for tx,trial in enumerate(block.segments) if trial.annotations['block']<=4]

    # there are two possibilities for input
    if counting=='ifr':
    #    use smoothed estimated instantenous firing rates:
        trials = [{'ID': tx, 'y': trial.analogsignals[0].magnitude, 'block': trial.annotations['block']}\
                  for tx,trial in enumerate(block.segments)\
                      if (trial.annotations['block'] in [2,4]) ]
                  
                  # and (trial.annotations['visual']==45) \
                  #     and (trial.annotations['audio']==5000)]
    
    elif counting=='spikes':
        #    or use the original spiketrains, in 1 ms blocks
        # y has to be (trajectory,neurons) for each trial
        trials = [ {'ID': tx, 'y': np.array([ neph.countspikes(trial.spiketrains[n],1*pq.ms,1*pq.ms,binsonly=True).magnitude \
                                    for n in range(len(trial.spiketrains)) ]).T, \
                    'block': trial.annotations['block'] }\
                  for tx,trial in enumerate(block.segments)   if (trial.annotations['block'] in [2,4]) ]
        
        n_factors = np.min(  (  n_factors, len(block.segments[0].spiketrains)  )  )


    elif counting=='test':
        dn = 'simulation'
        
        # # dynamics test:
        # a = np.array( [[-1,1],[2,2],[-2,5]] )      # 2 latent 3 neurons
        # x = np.vstack( [ np.arange(p_length)/p_length, np.sin(np.arange(p_length)) ] ).T[np.newaxis,:,:] + np.random.randn(n_trials,p_length,2)
        # y,h,rate = simulation.spike(x=x,a=a,b=np.array([1]))

        # stochastic test
        n_trials=20
        p_length = 100 *pq.ms
        mul_length = 1.25
        n_neuron = 3
        p_means = np.array([5,15,25])*pq.ms
        n_trajectory = (p_length/np.max(p_means)).magnitude.astype(np.int32)

        
        # trial: is a list of dictionaries, with ID and y in shape of (trajectory,neurons)
        trials = [{'ID': tx, 'y': np.array([ neph.countspikes_quantity(  np.random.poisson(p_means[n],size=(n_trajectory,1)).cumsum(axis=0)*pq.ms,\
                                                                         0*pq.ms, p_length*mul_length, 1*pq.ms  ) \
                                             for n in np.arange(n_neuron) ]).T }\
                  for tx in np.arange(n_trials) ]


    elif counting=='simulation':
        dn = 'simulation'
        np.random.seed(15)

        from vlgp import util,simulation

        ntrial = 50  # number of trials
        nbin = 6000  # number of time bins of each trial
        nneuron = 20 # number of neurons (spike trains)
        dim = 3 # latent dimension
        skip = 500
        lorenz = simulation.lorenz(skip +  ntrial * nbin, dt=5e-3, s=10, r=28, b=2.667, x0=np.random.random(dim))
        lorenz = sp.stats.zscore(lorenz[skip:, :]) 
        x = lorenz.reshape((ntrial, nbin, dim))  # latent dynamics in proper shape
        
        
        bias = np.log(15 / nbin)  # log base firing rate
        a = (np.random.rand(dim, nneuron) + 1) * np.sign(np.random.randn(dim, nneuron)) # loading matrix
        b = np.vstack((bias * np.ones(nneuron), -10 * np.ones(nneuron), -10 * np.ones(nneuron), -3 * np.ones(nneuron), 
                          -3 * np.ones(nneuron), -3 * np.ones(nneuron), -3 * np.ones(nneuron), -2 * np.ones(nneuron),
                          -2 * np.ones(nneuron), -1 * np.ones(nneuron), -1 * np.ones(nneuron)))  # regression weights
        
        y, _, rate = simulation.spike(x, a, b)
        sample = dict(y=y, rate=rate, x=x, alpha=a, beta=b)
        
        
        trials = [{'ID': i, 'y': y} for i, y in enumerate(sample['y'])]  # make trials
        
        
        np.savez('../../common/fluxturing/cache/data/lorentz.npz',y,x,a,b,bias)
        


    return
    

    A = np.array([trial['y'] for trial in trials])
    A = np.swapaxes(A,1,2)
    print('(trials,neurons, bins)', A.shape)



    if recalculate:
        vlgpfit = neba.vLGP(trials,n_factors=n_factors,gp_char_timescale=gp_char_timescale)
        pickle.dump(vlgpfit,open(cacheprefix+'vlgp/vlgp,f%d-%s_%s.pck'%(n_factors,dn,continuous_method),'wb'))
    else:
        vlgpfit = pickle.load(open(cacheprefix+'vlgp/vlgp,f%d-%s_%s.pck'%(n_factors,dn,continuous_method),'rb'))

    

    print(vlgpfit.keys())
    print(vlgpfit['config'].keys())
    print(vlgpfit['params'].keys())
    print(vlgpfit['trials'][0].keys())
    
    
    print('a:',vlgpfit['params']['a'].shape)
    print([(key,vlgpfit['params'][key]) for key in ['noise','sigma','omega','rank','gp_noise']])


    if doplot:
        
        # single trial plots:
        if 1:
            # triallist = [1,2,5,6]
            triallist = np.random.permutation(len(trials))[:4]

            # t = block.segments[0].analogsignals[0].times
            # t = np.arange(p_length.magnitude)*pq.ms*mul_length
            t = np.arange(len(vlgpfit['trials'][0]['y']))*pq.ms + T['starttime']


            U,S,V = np.linalg.svd(vlgpfit['params']['a'])
            print(U.shape,S.shape,V.shape)


    
            fig,ax = plt.subplots(1+n_factors+1,len(triallist)+1,figsize=((len(triallist)+1)*8,(1+n_factors+1)*8))
            
            for tx,tr in enumerate(triallist):
                
                # y = vlgpfit['trials'][tr]['y']
                mu = vlgpfit['trials'][tr]['mu']
                x = vlgpfit['trials'][tr]['x']
                v = vlgpfit['trials'][tr]['v']
                w = vlgpfit['trials'][tr]['w']
                
                # if tx==0: print('y mu x v w\n',y.shape,mu.shape,x.shape,v.shape,w.shape)
                
                # mu = vlgpfit['trials'][tr]['mu']  # extract posterior latent
                # W = np.linalg.lstsq(mu, x[tr,...], rcond=None)[0]
                mu = mu @ U.T
                
                
                
                
                # activity       
                axs = ax[0,tx]
                
                # axs.spy(y.T,aspect='auto')
                # axs.eventplot(block.segments[tx].spiketrains,color='red')
                
                axs.set_yticks([])
                figs.setxt(axs)
                figs.plottoaxis_stimulusoverlay(axs,T)
                if tx==0: axs.set_ylabel('%s\nspikes by neurons'%dn)
                axs.set_title('trial %d'%vlgpfit['trials'][tr]['ID'])
                
                
                
                # latents
                for lx in range(n_factors):
                    axs = ax[lx+1,tx]
                    
                    
                    # plt.plot(x[0, ...] + 2 * np.arange(3), color="b")
                    # plt.plot(mu + 2 * np.arange(3), color="r")
                    
                    # axs.plot(t, x[tr,:,lx], color='k' )
                    axs.plot(t, mu[:,lx], color='purple' )
                    # axs.set_ylim(-2000,2000)
                    if tx==0: axs.set_ylabel('latent %d'%(lx+1))
                    
                    figs.setxt(axs)
                    figs.plottoaxis_chancelevel(axs,0)
                    figs.plottoaxis_stimulusoverlay(axs,T)

        
                axs = ax[-1,tx].remove()
                axs = fig.add_subplot(1+n_factors+1,len(triallist)+1,(len(triallist)+1)*(1+n_factors)+1+tx, projection='3d')
                # axs.plot(x[tr,:,0],x[tr,:,1],x[tr,:,2], color='k' )
                axs.plot(mu[:,0],mu[:,1],mu[:,2],color='purple')
    
            
            # create trial averages
            mu_a = np.array([trial['mu'] for trial in vlgpfit['trials']])
            mu_m = mu_a.mean(axis=0)
            mu_e = 2*mu_a.std(axis=0)/np.sqrt(mu_a.shape[0])


            for lx in range(n_factors):
                axs = ax[lx+1,-1]
                axs.plot(t, mu_m[:,lx],lw=2,color='rebeccapurple')
                axs.fill_between(t,mu_m[:,lx]-mu_e[:,lx],mu_m[:,lx]+mu_e[:,lx], color='mediumvioletred',alpha=0.3)
                figs.setxt(axs)
                figs.plottoaxis_chancelevel(axs,0)
                figs.plottoaxis_stimulusoverlay(axs,T)

            save = 1 or globalsave
            if save:
                fig.savefig(resultpath+'vlgp,singletrials-f%d-%s'%(n_factors,dn)+ext)
            
            
        
        # trial average plots
        if 1:

            trial_ids = np.array([ trial['ID'] for trial in vlgpfit['trials'] ])
            block_ids = np.array([ trial['block'] for trial in vlgpfit['trials'] ],dtype=np.int16)
            ixv45,ixv135,ixa45,ixa135,ixv5000,ixv10000,ixa5000,ixa10000 = preprocess.getstimulusidents(dn,block,multimodalonly=True)
            hit,miss,correctrejection,falsealarm = preprocess.assessbehaviouralperformance(dn,modality='all',multimodalonly=True)
            
            taskaspects = ['visual','audio','context','choice']
            variablecolors = ['navy','darkgreen','mediumvioletred','orange']
            triallists = [  [np.concatenate((ixv45,ixa45)), np.concatenate((ixv135,ixa135))],\
                            [np.concatenate((ixa5000,ixv5000)), np.concatenate((ixa10000,ixv10000))],\
                            [np.concatenate((ixv45,ixv135)), np.concatenate((ixa45,ixa135))],\
                            [np.concatenate((hit,falsealarm)), np.concatenate((correctrejection,miss))] \
                          ]
            t = np.arange(len(vlgpfit['trials'][0]['y']))*pq.ms + T['starttime']
            mu_a = np.array([trial['mu'] for trial in vlgpfit['trials']])
            
            # rotate
            # U,S,V = np.linalg.svd(vlgpfit['params']['a'])
            # mu_a = mu_a @ U.T

            fig,ax = plt.subplots(n_factors,4,figsize=((4)*8,(n_factors)*8))
            
            for cx,(comparison,triallistpair) in enumerate(zip(taskaspects,triallists)):
                for lx in range(n_factors):
                    axs = ax[lx,cx]
                    for  cidx in [0,1]:          # response to given class
                        idx = np.array([ np.where(trial_ids==tr)[0][0] for tr in triallistpair[cidx] ])
                        mu_m = mu_a[idx,:,lx].mean(axis=0)
                        mu_e = 2*mu_a[idx,:,lx].std(axis=0)/np.sqrt(mu_a[idx,:,:].shape[0])
                        axs.plot(t, mu_m,linewidth=2,color=variablecolors[cx],alpha=1.-2./3.*cidx)
                        axs.fill_between(t,mu_m-mu_e,mu_m+mu_e,\
                                         color=variablecolors[cx],alpha=(1.-2./3.*cidx)/2.)
                    
                    figs.setxt(axs)
                    figs.plottoaxis_chancelevel(axs,0)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    if cx==0: axs.set_ylabel('latent %d'%(lx+1))
                    if lx==0: axs.set_title(comparison)
                
            fig.suptitle(dn+' vLGP, trial averaged latents')
                
            
            save = 1 or globalsave
            if save:
                fig.savefig(resultpath+'vlgp,variables,averaged-f%d-%s'%(n_factors,dn)+ext)






if __name__ == '__main__':
    print('and the little beasts are coming...')
