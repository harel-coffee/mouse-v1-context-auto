#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 10:13:00 2020

@author: mahajnal
"""

from config import *

import preprocess



# all routines here put all mice calculation results onto the same figure (aggregation)







def decoderstobehaviouralperformance(datanames):
    # correlate behavioural performance with task variable decoding accuracy for multiple mice

    examine = 'allexpcond'
    layeredfilename = 'all'
    stimulusspecificity = 'identical,go'
    
    n_mice = len(datanames)
    
    # offstim training interval starts:
    trainingoffsets = np.array([0,375,750,1225])*pq.ms
    trainingpoints = T['starttime'] + trainingoffsets
    
    # feature interval (along timepoints dimension)
    width = 100*pq.ms
    wx = int((width/T['dt']).magnitude)
    
    aucdeltat = 750*pq.ms
    
    
    miceperf = []    #   (hit %, miss %,    correct rejection %, false alarm %)
    miceauc = []     #   mice   x    modality    x     (on,off)    x   (diff,dec,offdec)
    for n,dn in enumerate(datanames):
        print(dn)


        # get the performance metrics for mice
        behaviour = []
        performance = []
        for m,modality in enumerate(['visual','audio','all']):
            behaviour.append(preprocess.assessbehaviouralperformance(dn,modality))
            
            performancecalc = np.asarray(  [ len(b) for b in behaviour[m] ]  )
            go = performancecalc[:2].sum()
            nogo = performancecalc[2:].sum()
            
            performance.append( performancecalc / np.array([  go,go,nogo,nogo  ]) )
            
        miceperf.append( performance )
#        print(dn,performance)
        
        
    



        # get the neural differentiability performance
        modalityauc = []
        for cx,comparison in enumerate(['visual','audio','context']):
            
            acrossdifference = pickle.load(open(cacheprefix+'continuous/responsedifference-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))

            dm,de = nedi.get_maxnormdifference(acrossdifference)
            dz = np.max((np.zeros(dm.shape),de),axis=0)       # clip below zero
            
            diffsignal = neo.AnalogSignal(dz,t_start=dm.t_start,sampling_period=dm.sampling_period,units=dm.units)
            
            decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # lower confidence bound


            for trx,trainingpoint in enumerate(trainingpoints):
                offstimdecoderfull = pickle.load(open(cacheprefix+'continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,comparison,trx),'rb'))
                t1 = int( ((-width+T['stimstarttime']-offstimdecoderfull[1].t_start)/T['dt']).magnitude)        # [1] is to use the test data
                t2 = int( ((-width+T['endtime']-offstimdecoderfull[1].t_start)/T['dt']).magnitude)
                
                print(trx,      trainingpoint, 'length', offstimdecoderfull[1].duration, 'length index:',len(offstimdecoderfull[1]), 'start', offstimdecoderfull[1].t_start,    'index:',  t1,t2)

                offstimdecoder = offstimdecoderfull[1][ t1:t2,:]

                # average out over the train starting points
                if trx==0:
                    offsignal = 1/len(trainingpoints) * (offstimdecoder[:,0])#-offstimdecoder[:,2])    # lower confidence bound
                else:
                    offsignal += 1/len(trainingpoints) * (offstimdecoder[:,0])#-offstimdecoder[:,2])

            # get AUC of the signals, which are the lower confidence bounds of differentiability
#            offsignal.t_start = offstimdecoderfull[1].t_start
            diffauc_on = neph.get_auctimenormalized(diffsignal,T['stimstarttime'],T['stimendtime'],'diff.on.'+comparison)
            decauc_on = neph.get_auctimenormalized(decsignal,T['stimstarttime'],T['stimendtime'])
            offauc_on = neph.get_auctimenormalized(offsignal,0.*pq.ms,3000.*pq.ms)

            diffauc_pre = neph.get_auctimenormalized(diffsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'],'diff.off.'+comparison)
            decauc_pre = neph.get_auctimenormalized(decsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])
            offauc_pre = neph.get_auctimenormalized(offsignal,-aucdeltat,0.*pq.ms)

            diffauc_post = neph.get_auctimenormalized(diffsignal,T['stimendtime'],T['stimendtime']+aucdeltat,'diff.off.'+comparison)
            decauc_post = neph.get_auctimenormalized(decsignal,T['stimendtime'],T['stimendtime']+aucdeltat)
            offauc_post = neph.get_auctimenormalized(offsignal,3000.*pq.ms,3000*pq.ms+aucdeltat)
            
            print(dn,comparison,offauc_on,offauc_pre,offauc_post)
            
            modalityauc.append( [[ diffauc_on,decauc_on,offauc_on ],\
                                 [ diffauc_pre,decauc_pre,offauc_pre ],\
                                 [ diffauc_post,decauc_post,offauc_post ]])

        miceauc.append(modalityauc)
        
#    print('shape',np.array(miceauc).shape)
#    print(np.array(miceperf).shape)
    
    
    for cx,comparison in enumerate(['visual','audio','context']):
        fig,ax = plt.subplots(3,6,figsize=(42,21))
        for n in range(n_mice):
            for ox,offlabel in enumerate(['on stim 3000 ms','pre stim %d ms'%aucdeltat, 'post stim %d ms'%aucdeltat]):
                for dx,diffmeasurelabel in enumerate(['maxn. rate diff.','pointwise. dec.','offstim tr. pred.dec.']):
                    for px,performancetype in enumerate(['GO hit','NOGO correct rejection']):
                        axs = ax[ox,dx*2+px]
#                        if datanames[n]>'ME100': color = 'k'; else: color = ''
                        axs.plot(  miceperf[n][cx][px*2],  miceauc[n][cx][ox][dx],'o', label=datanames[n]  )
                        
                        # draw behavioural accuracy: chance and max limits
                        axs.set_xlim([0.4,1.01]); axs.plot([0.5,0.5],[0,2],'r--',alpha=0.2); axs.plot([1,1],[0,2],'--',color='darkgreen',alpha=0.2)
                        
                        # draw neural differentiability limits, first is rate diff, second is decoder chance and max
                        if dx==0: axs.set_ylim([-0.05,2]); axs.plot([0,1],[0,0],'r--',alpha=0.2); 
                        else: axs.set_ylim([0.4,1.01]);  axs.plot([0,1],[0.5,0.5],'r--',alpha=0.2); axs.plot([0,1],[1,1],'--',color='darkgreen',alpha=0.2)
                        
                        if ox==0 and dx==0 and px==0: axs.set_ylabel('z.diff.95%% c.i.\n'+offlabel)
                        if ox==0 and dx>0 and px==0: axs.set_ylabel('dec. roc auc')
                        if ox>0 and dx==0 and px==0: axs.set_ylabel(offlabel)
                        if ox==0: axs.set_title(diffmeasurelabel+'\n'+performancetype)
                        if ox==2: axs.set_xlabel('behav. perf.')
                        
                        if px==0 and dx==0 and ox==0: axs.legend(loc=2)

        fig.suptitle(comparison)

    return






def decodersbehaviourwithtraining(datanames):
    # correlate all training and the recording session behavioural performance with task variable decoding accuracy

    datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    datanames.remove('DT008')
    n_mice = len(datanames)
    
    Xa = []
    Ya = []
    perf = np.zeros((n_mice,6))
    expertise = np.zeros(n_mice)

    for n,dn in enumerate(datanames):
        print(dn)
        a,b,b_s,_,_,_ = preprocess.loadtrainingbehaviouraldata(dn)       # a: all data, b: behav, b_s behav sem.
        perf[n,:] = [ b[1],b[3], b[5], 2*b_s[1], 2*b_s[3], 2*b_s[5] ]      # behaviour and its s.e.m.


        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        expertise[n] = acrossdecoder[1][ T['start_idx']:T['stimstart_idx'], 0 ].mean()

        Ya.append( [a[1],a[3], np.concatenate( [a[1],a[3]] ) ] )
        Xa.append( [ expertise[n]*np.ones(len(Ya[-1][0])), expertise[n]*np.ones(len(Ya[-1][1])), expertise[n]*np.ones(len(Ya[-1][2])) ] )


    Y = [ np.concatenate( [  Ya[n][d] for n in range(n_mice)  ] ) for d in range(3) ]
    X = [ np.concatenate( [  Xa[n][d] for n in range(n_mice)  ] ) for d in range(3) ]
    
#    behav = np.zeros((len(datanames),4))
#    
#    for n,dn in enumerate(datanames):
#        h,m,c,f = preprocess.assessbehaviouralperformance(dn)
#        behav[n,:] = np.array([ len(h),len(m),len(c),len(f) ])
#    
#    fraction_hitcorrrej = np.c_[behav[:,0]/(behav[:,:2].sum(axis=1)), behav[:,2]/(behav[:,2:].sum(axis=1))]
#    fraction_missfal = np.c_[behav[:,1]/(behav[:,:2].sum(axis=1)), behav[:,3]/(behav[:,2:].sum(axis=1))]
#    perf = fraction_missfal
    
    

        
    print(expertise.shape,perf.shape)
    
    
    # L, behaviour
#    norminv(fraction trials	with	correct	licking)	– norminv(fraction	trials	with	incorrect	licking
    colors = ['darkorange','firebrick','purple']
    ylabels = 'error fraction'

    fig, ax = plt.subplots(1,1,figsize=(10,10))
    for dx in [0,1,2]:     # for misses+fals do both separately
        axs = ax
        x = expertise
        y = perf[:,dx];
        
        x = np.array([ np.mean(Xa[n][dx]) for n in range(n_mice) ])
        y = np.array([ np.mean(Ya[n][dx]) for n in range(n_mice) ])
        
        axs.scatter(x,y,s=150,marker='o',color=colors[dx],alpha=0.8,label=['miss','false alarm','miss & false alarm'][dx])
        for n in range(n_mice): axs.errorbar(x[n],y[n],yerr=perf[n,dx+3],color=colors[dx],alpha=0.8)
#        if dx==1:
#            for n,dn in enumerate(datanames): axs.text(x[n]+0.01,y[n],dn,color='grey',fontsize='small',verticalalignment='center')

        axs.set_xlim(0.45,1.01)
#        axs.set_ylim(-0.04,1.04)
        axs.set_ylim(-0.04,0.54)
        if dx==1: figs.plottoaxis_crosshair(axs,0.5,0)
        axs.set_xlabel('context PRE', labelpad=0.1, verticalalignment='bottom')

        axs.set_xticks([0.5,1])
        axs.set_yticks([0,0.5])
        axs.set_yticklabels(['0','0.5'])
        

        l = sp.stats.linregress(X[dx],Y[dx])
        if l[3]<0.05:
            line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],color=colors[dx],linewidth=2)
        else:
            line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],'--',color=colors[dx],linewidth=2)
            
        xoffs = 0.05  #0.35
#        yoffs = y.mean()*1.4 - 0.06     #/2.5 + 0.35
        yoffs = [-0.032, 0.34, 0.07][dx]
        if l[3]<0.001:
            axs.text(line[0]+xoffs,yoffs,'p<%5.3f, c=%5.3f, R$^2$=%6.4f'%((0.001,l[2],l[2]**2)),color=colors[dx])
        else:
            axs.text(line[0]+xoffs,yoffs,'p=%5.3f, c=%5.3f, R$^2$=%6.4f'%((l[3],l[2],l[2]**2)),color=colors[dx])
        
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)

    axs.legend(frameon=False)
    axs.set_title('%d mice, %d trials'%(n_mice,np.sum(X[2])))
    axs.set_ylabel(ylabels, labelpad=0, verticalalignment='top')

    save = 0
    if save:
        fig.savefig(resultpath+'micetribe,training+recordedbehaviour-contextpre,multimodalonly'+ext)
    














def aggregatemice_relevantirrelevant_behaviour(datanames):
    # correlate behaviour with identity of contextually relevant vs. irrelevant modality stimulus

    doplot = 0 or globaldoplot

    taskaspects = ['gonogo-congruent,av','gonogo-conflict,av','gonogo-congruent,aa','gonogo-conflict,aa']


    print(datanames)

    if datanames[0][:2]=='AC': area = 'ACC'
    else: area = 'V1'


    for n,dn in enumerate(datanames):

        blv,bla = preprocess.getorderattended(dn)
        comparisongroups  = [\
                                [ [ [blv[1]], [45],  [5000] ], [ [blv[1]],   [135], [10000] ]  ], \
                                [ [ [blv[1]], [45], [10000] ], [ [blv[1]],   [135],  [5000] ]  ], \
                                [ [ [bla[1]], [45],  [5000] ], [ [bla[1]],   [135], [10000] ]  ], \
                                [ [ [bla[1]], [135], [5000] ], [ [bla[1]],    [45], [10000] ]  ], \
                            ]


        p = preprocess.get_conditioned_behaviour(dn,comparisongroups,taskaspects)

        p.reset_index(inplace=True)
        p.rename(columns={'index':'condition'},inplace=True)
        p.insert(0,'mouse',[dn,dn,dn,dn])
        p.set_index(['mouse','condition'],inplace=True)

        if n==0: P=p
        else: P = pd.concat([P,p],axis=0)

    print(P)
    P.xs(taskaspects[1],level='condition')



    if doplot:

        fig,ax = plt.subplots(1,1,figsize=(14*1,8*1))
        
        axs = ax

        for k in range(len(taskaspects)):
            print(taskaspects[k])
            D = P.xs(taskaspects[k],level='condition')
            print(D)
            axs.boxplot(x=D,positions=[k*4-0.5,k*4+0.5],whis=[5,95],widths=[1,1],notch=True,labels=['go','nogo'],bootstrap=1000)
        
        
        
        axs.set_title('\nattend visual                          attend audio\ncong              confl              cong              confl')

        axs.set_ylabel('%s, %d mice\ncongruent and conflicting\nbehavioural performance'%(area,len(datanames)))


        # fig.suptitle()


        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'behav-congruentconflicting_%s'%area+ext)
























def aggregatemice_singletrialneuralbehaviour_relevantirrelevantcongruentconflicting(datanames):
     # decode the relevant and the irrelevant stimuli, and congruent and conflicting stimuli trials separately, and for the two contexts

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot


    n_mice = len(datanames)


    taskaspects = ['visual,av','audio,av','visual,aa','audio,aa',\
                   'gonogo-congruent,av','gonogo-conflict,av','gonogo-congruent,aa','gonogo-conflict,aa']
    
    width = 50*pq.ms





    perftask_all = [[] for cx in range(len(taskaspects))]
    probasignals_tasks_all = [ [] for cx in range(len(taskaspects)) ]


    for n,dn in enumerate(datanames):

        # setup stimulus: 
        blv,bla = preprocess.getorderattended(dn)
        comparisongroups  = [\
                                [ [ [blv[1]], [45],      [] ], [ [blv[1]],   [135],      [] ]  ], \
                                [ [ [blv[1]],   [],  [5000] ], [ [blv[1]],      [], [10000] ]  ], \
                                [ [ [bla[1]], [45],      [] ], [ [bla[1]],   [135],      [] ]  ], \
                                [ [ [bla[1]],   [],  [5000] ], [ [bla[1]],      [], [10000] ]  ], \

                                [ [ [blv[1]], [45],  [5000] ], [ [blv[1]],   [135], [10000] ]  ], \
                                [ [ [blv[1]], [45], [10000] ], [ [blv[1]],   [135],  [5000] ]  ], \
                                [ [ [bla[1]], [45],  [5000] ], [ [bla[1]],   [135], [10000] ]  ], \
                                [ [ [bla[1]], [135], [5000] ], [ [bla[1]],    [45], [10000] ]  ], \
                            ]



        # assess performance

        perftask = preprocess.get_conditioned_behaviour_singletrial(dn,comparisongroups,taskaspects)

        for cx,comparison in enumerate(taskaspects):
            if n==0:
                perftask_all[cx] = perftask[cx]
            else:
                perftask_all[cx] = pd.concat((perftask_all[cx], perftask[cx]), axis=0)
        

        # assess neural
        probasignals_tasks = []
        for cx,comparison in enumerate(taskaspects):
            probasignals, acc = pickle.load(open(cacheprefix+'acc/relevantirrelevant+congruentconflict,singletrialneuralbehaviour_%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))
            probasignals_tasks.append(probasignals)
            if n==0:
                probasignals_tasks_all[cx] = probasignals
            else:
                probasignals_tasks_all[cx] = np.concatenate( (probasignals_tasks_all[cx],probasignals), axis=0 )



        # for i in range(len(taskaspects)):
        #     print(perftask[i].shape,probasignals_tasks[i].shape)



    
    for cx,comparison in enumerate(taskaspects):
        print(comparison, perftask_all[cx].shape,probasignals_tasks_all[cx].shape)



    


    if doplot:


        fig,ax = plt.subplots(4,2,figsize=(3*12,4*12))
        # we want to compare correct against error trials' class probability predictions
        # we have to do it for go and nogo trials separately (as they are the classes)

        for tx in [0,1]:                 # tasks
            for cx in [0,1]:             # context
                for mx in [0,1]:         # modality/congconfl
                    ix = tx*4+cx*2+mx
                    axs = ax[tx*2+mx,cx]
                    y_early = probasignals_tasks_all[ix].magnitude[:,150:175].mean(axis=1)
                    y_late = probasignals_tasks_all[ix].magnitude[:,300:450].mean(axis=1)
                    
                    for clx,signal in enumerate(['go','nogo']):    # choose go and nogo signals respectively
                        s = ['go','nogo'][clx]

                        mask_correct = perftask_all[ix][signal]==1     # correct trials
                        mask_error = perftask_all[ix][signal]==0       # error trials

                        axs.boxplot(x=[y_early[mask_correct],y_early[mask_error]],positions=[0+clx*2.5,6+clx*2.5],labels=['%s\nearly\ncorrect'%s,'%s\nearly\nerror'%s], notch=True)
                        axs.boxplot(x=[y_late[mask_correct],y_late[mask_error]], positions=[1+clx*2.5,7+clx*2.5],labels=['%s\nlate\ncorrect'%s,'%s\nlate\nerror'%s], notch=True)


                    axs.set_yticks([0,1])
                    axs.set_yticklabels(['P(go)=1','P(nogo)=1'])
                    figs.plottoaxis_chancelevel(axs,0.5)

                    if cx==0: axs.set_ylabel('predicted probability, neural decoder')


                    axs.set_title(taskaspects[tx*4+cx*2+mx],fontsize=20)

        fig.suptitle('%d ACC mice, probability of predicted stimuli representations in two contexts in correct and error trials'%n_mice)
        
        fig.tight_layout()

        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'relevantirrelevant+congruentconflict,singletrialneuralbehaviour_ACC_%s-%dms'%(continuous_method,T['dt'])+ext)

    return














def aggregatemice_singletrialneuralbehaviour_context(datanames):
    # decode the context, display probability for trials separately


    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot


    n_mice = len(datanames)


    taskaspects = ['context']
    classlabels = ['visual','audio']
    
    width = 50*pq.ms





    masks_all = [[] for cx in range(4)]
    perftask_all = [[] for cx in range(len(taskaspects))]
    probasignals_tasks_all = [ [] for cx in range(len(taskaspects)) ]


    for n,dn in enumerate(datanames):


        # setup stimulus: 


        blv,bla = preprocess.getorderattended(dn)
        comparisongroups  = [\
                                [ [ [blv[1]], [],      [] ], [ [bla[1]],   [],      [] ]  ], \
                            ]


        # create base mask for conditioned testing of decoder; always two element list for the two classes of context
        g = preprocess.loadexperimentdata(dn)
        g['block']+=1
        mask_base = [ g[ g['block'].isin([blv[1]]) ], g[ g['block'].isin([bla[1]]) ] ]
        # mask_base = g[ g['block'].isin([blv[1],bla[1]]) ]
        # create specific masks
        masks = [ [ ((mb['degree']==45) & (mb['freq']==5000))                  for mb in mask_base ],\
                [ ((mb['degree']==135) & (mb['freq']==10000))                for mb in mask_base ],\
                [ ((mask_base[0]['degree']==45) & (mask_base[0]['freq']==10000)),  ((mask_base[1]['degree']==135) & (mask_base[1]['freq']==5000)) ],\
                [ ((mask_base[0]['degree']==135) & (mask_base[0]['freq']==5000)),  ((mask_base[1]['degree']==45) & (mask_base[1]['freq']==10000)) ],\
                ]
        masks = np.array([  np.concatenate([mask[0],mask[1]]) for mask in masks        ])
        mask_labels = ['cong,go,correct','cong,go,error','cong,nogo,correct','cong,nogo,error',\
                    'confl,go,correct','confl,go,error','confl,nogo,correct','confl,nogo,error']

        for mx in range(4):
            masks_all[mx].extend(masks[mx])



        # assess performance
        perftask = preprocess.get_conditioned_behaviour_singletrial(dn,comparisongroups,taskaspects,classlabels=classlabels)


        for cx,comparison in enumerate(taskaspects):
            if n==0:
                perftask_all[cx] = perftask[cx]
            else:
                perftask_all[cx] = pd.concat((perftask_all[cx], perftask[cx]), axis=0)
        









        # assess neural
        probasignals_tasks = []
        for cx,comparison in enumerate(taskaspects):
            probasignals, acc = pickle.load(open(cacheprefix+'acc/singletrialneuralbehaviour_%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))
            if dn=='AC003': probasignals = 1 - probasignals
            probasignals_tasks.append(probasignals)
            if n==0:
                probasignals_tasks_all[cx] = probasignals
            else:
                probasignals_tasks_all[cx] = np.concatenate( (probasignals_tasks_all[cx],probasignals), axis=0 )



        # for i in range(len(taskaspects)):
        #     print(perftask[i].shape,probasignals_tasks[i].shape)


    masks_all = [  np.array(mask_all)  for mask_all in masks_all]

    n_trials = perftask_all[cx].shape[0]
    
    for cx,comparison in enumerate(taskaspects):
        print(comparison, perftask_all[cx].shape,probasignals_tasks_all[cx].shape,masks_all[0].shape,masks_all[1].shape)



    


    if doplot:

        if 0:

            fig,ax = plt.subplots(1,1,figsize=(1*17,1*14))
            # we want to compare correct against error trials' class probability predictions
            # we have to do it for both contexts separately (as they are the classes)

            ix = 0   # we only have context
            axs = ax


            y_pre = probasignals_tasks_all[ix].magnitude[:,0:150].mean(axis=1)
            y_early = probasignals_tasks_all[ix].magnitude[:,150:300].mean(axis=1)
            y_late = probasignals_tasks_all[ix].magnitude[:,300:450].mean(axis=1)
            
            for clx,signal in enumerate(classlabels):    # choose the two classes the decoder separated
                s = signal

                mask_correct = perftask_all[ix][signal]==1     # correct trials
                mask_error = perftask_all[ix][signal]==0       # error trials

                # axs.boxplot(x=[y_early[mask_correct],y_early[mask_error]],positions=[0+clx*2.5,6+clx*2.5],labels=['%s\nearly\ncorrect'%s,'%s\nearly\nerror'%s], notch=True)
                # axs.boxplot(x=[y_late[mask_correct],y_late[mask_error]], positions=[1+clx*2.5,7+clx*2.5],labels=['%s\nlate\ncorrect'%s,'%s\nlate\nerror'%s], notch=True)

                for yx,(y,lab) in enumerate(zip([y_pre,y_early,y_late],['pre','early','late'])):
                    axs.boxplot(x=[y[mask_correct],y[mask_error]],positions=[0+yx*7+clx*3,1.5+yx*7+clx*3],widths=0.8,labels=['%s\n%s\ncorrect'%(s,lab),'%s\n%s\nerror'%(s,lab)], notch=True)
                    y_logit = np.log(y/(1-y))
                    _,p = sp.stats.ttest_ind(y[mask_correct],y[mask_error])
                    _,p_logit = sp.stats.ttest_ind(y_logit[mask_correct],y_logit[mask_error])
                    axs.text(0.75+yx*7+clx*3,1.18,'prob t p=%5.4f\n odds t p=%5.4f'%(p,p_logit),horizontalalignment='center',verticalalignment='top',fontsize=15)
                    

            axs.set_yticks([0,1])
            axs.set_yticklabels([ 'P(%s)=1'%cs for cs in classlabels ])
            axs.set_ylim(-0.2,1.2)
            figs.plottoaxis_chancelevel(axs,0.5)

            if cx==0: axs.set_ylabel('predicted probability, neural decoder')


            axs.set_title(taskaspects[ix],fontsize=20)


            fig.suptitle('%d ACC mice, %d trials, probability of\npredicted context representations in correct and error trials'%(n_mice,n_trials))
            
            save = 0 or globalsave
            if save:
                fig.savefig(resultpath+'context,singletrialneuralbehaviour_ACC_%s-%dms'%(continuous_method,T['dt'])+ext)


        if 1:


            fig,ax = plt.subplots(1,1,figsize=(1*17,1*14))

            ix = 0   # we only have context
            # x = perftask[ix]
            y_pre = probasignals_tasks_all[ix].magnitude[:,0:150].mean(axis=1)
            y_early = probasignals_tasks_all[ix].magnitude[:,150:300].mean(axis=1)
            y_late = probasignals_tasks_all[ix].magnitude[:,300:450].mean(axis=1)


            # for mx,mask in enumerate(masks):
            axs = ax
            for clx,signal in enumerate(classlabels):    # choose the two classes the decoder separated
                s = signal


                mask_correct = perftask_all[ix][signal]==1     # correct trials
                mask_error = perftask_all[ix][signal]==0       # error trials

                print(mask_correct.shape,mask_error.shape, masks_all[0].shape)

                for yx,(y,lab) in enumerate(zip([y_pre,y_early,y_late],['pre','early','late'])):

                    ylist = []
                    for mask in masks_all:
                        ylist.extend( [ y[mask & mask_correct], y[mask & mask_error] ] )
                    parts = [0,1,3,4,6,7,9,10]
                    axs.boxplot(x=ylist,positions=[mx*0.6+yx*20+clx*10 for mx in parts], labels=None, widths=0.55, notch=True)
                    for mx,m in enumerate(parts):
                        axs.text(m*0.6+yx*20+clx*10,1.45,mask_labels[mx],rotation=90, horizontalalignment='center', verticalalignment='center',fontsize=10)
                        if mx%2==1:
                            y_logits = [np.log(ylist[mx-1]/(1-ylist[mx-1])), np.log(ylist[mx]/(1-ylist[mx])) ]
                            _,p = sp.stats.ttest_ind(ylist[mx-1],ylist[mx])
                            _,p_logit = sp.stats.ttest_ind(y_logits[0],y_logits[1])
                            axs.text(m*0.6-0.5+yx*20+clx*10,1.3,'prob t p=%5.4f\n odds t p=%5.4f'%(p,p_logit),rotation=90,horizontalalignment='center',verticalalignment='top',fontsize=7)
                    axs.text(4.5*0.6+yx*20+clx*10,-0.05,'attend %s\n%s'%(signal,lab),horizontalalignment='center', verticalalignment='center',fontsize=10)
                    
            axs.set_xticklabels([])
            axs.set_yticks([0,1])
            axs.set_yticklabels([ 'P(%s)=1'%cs for cs in classlabels ])
            axs.set_ylim(-0.1,1.6)
            figs.plottoaxis_chancelevel(axs,0.5)

            if cx==0: axs.set_ylabel('class prob.')


            axs.set_title(taskaspects[ix],fontsize=20)

            fig.suptitle('%d ACC mice, %d trials, probability of predicted\nrepresentations of context in congruent and conflicting trials'%(n_mice,n_trials))
            

            fig.tight_layout()

            save = 0 or globalsave
            if save:
                fig.savefig(resultpath+'context,singletrialneuralbehaviour-conditioned_ACC_%s-%dms'%(continuous_method,T['dt'])+ext)            

    return









def aggregatemice_subspaces_timeresolved_behaviour_ACC(datanames):
    # project activity to stimulus DV
    # show correct and error trials; are error trials closer to discrimination threshold than correct?
    # do this for all mice, and pool the trials together


    n_mice = len(datanames)
    n_neurons_sum = 0

    doplot = 0 or globaldoplot


    width = 50*pq.ms
    wx = int(width/T['dt'])



    Y_all = [  [  [] for j in range(2)        ] for i in range(8)  ]
    mask_correct_all = [  [  [] for j in range(2)        ] for i in range(8)  ]
    mask_error_all = [  [  [] for j in range(2)        ] for i in range(8)  ]


    for n,dn in enumerate(datanames):
        
        block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=True,recalculate=False)

        n_neuron = block.segments[0].analogsignals[0].shape[1]
        n_neurons_sum += n_neuron


        blv,bla = preprocess.getorderattended(dn)

        comparisongroups  = [\
                                [ [ [blv[1]], [45],      [] ], [ [blv[1]],   [135],      [] ]  ], \
                                [ [ [blv[1]],   [],  [5000] ], [ [blv[1]],      [], [10000] ]  ], \
                                [ [ [bla[1]], [45],      [] ], [ [bla[1]],   [135],      [] ]  ], \
                                [ [ [bla[1]],   [],  [5000] ], [ [bla[1]],      [], [10000] ]  ], \

                                [ [ [blv[1]], [45],  [5000] ], [ [blv[1]],   [135], [10000] ]  ], \
                                [ [ [blv[1]], [45], [10000] ], [ [blv[1]],   [135],  [5000] ]  ], \
                                [ [ [bla[1]], [45],  [5000] ], [ [bla[1]],   [135], [10000] ]  ], \
                                [ [ [bla[1]], [135], [5000] ], [ [bla[1]],    [45], [10000] ]  ], \
                            ]

        taskaspects = ['visual,av','audio,av','visual,aa','audio,aa',\
                    'gonogo-congruent,av','gonogo-conflict,av','gonogo-congruent,aa','gonogo-conflict,aa']






        activity = []   # this will be [taskaspects][class][trials][trajectory,neurons]   (list,list,list,analogsignal)


        # collect behaviour
        # perftasks = []
        # for cx,comparison in enumerate(taskaspects):
        perftask = preprocess.get_conditioned_behaviour_singletrial(dn,comparisongroups,taskaspects)
        # perftasks.append(perftask)


        # collect DBNVs, to project the activities collected above onto
        W = []
        for cx,comparison in enumerate(taskaspects):
            # collect activity
            activity.append(  preprocess.collect_stimulusspecificresponses(block, comparisongroups[cx]) )
            
            # collect decoders, already calculated, in this function we just display a different crossection:
            acrossdecoder = pickle.load(open(cacheprefix+'acc/stimulus,relevant-irrelevant_%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))
            w,w_mean = nedi.get_decoder_dv(acrossdecoder, n_neuron, n_preweightelements=7, timeaverage_start_dx=150, timeaverage_end_dx=150+150)
            W.append(w)



        W = np.array(W) # [taskaspects,neurons,trajectory]

        # print(len(acrossdecoder),len(acrossdecoder[0]))

        times = block.segments[0].analogsignals[0].times[:-wx]

        # find the orthogonal axes
        # Q,R = np.linalg.qr(W)     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix




        # extend each Y so it contains all trials for all mice
        for cx,comparison in enumerate(taskaspects):
                # project activity onto the DV

                action = ['go','nogo']

                for clx,signal in enumerate(action):    # choose go and nogo signals respectively

                    X = np.array(activity[cx][clx])                        # [trials,trajectory,neurons]
                    Y = np.array( [  X[:,t,:] @ W[cx,:,t]    for t in range(W.shape[2])     ] ).T         # [trials,trajectory]
                    # Y *= 2*(clx-0.5)        # symmetrise around the decision boundary
                    Y_all[cx][clx].extend(Y)
                    
                    P = perftask[cx][signal]
                    P = P[P.notna()]
                    mask_correct = P==1     # correct trials
                    mask_error = P==0       # error trials

                    mask_correct_all[cx][clx].extend(mask_correct)
                    mask_error_all[cx][clx].extend(mask_error)


    # print(len(Y_all),len(Y_all[0]),len(Y_all[0][0]),len(Y_all[0][0][0]),)




    if doplot:

    

        fig,ax = plt.subplots(2,len(taskaspects)//2, figsize=(len(taskaspects)//2*9,2*8))
        
        colors = [ ['darkgreen','red'],   ['dodgerblue','orange'] ]
        label_behav = ['correct','error']
        label_action = ['go','nogo']


        for cx,comparison in enumerate(taskaspects):

            axs = ax[cx//4,cx%4]
            
            for clx,signal in enumerate(label_action):    # choose go and nogo signals respectively
        
                # print(cx,clx,'    ',Y_aux.shape,Y.shape,mask_correct.shape,mask_error.shape)

                Y = np.array(Y_all[cx][clx])

                Y[:,:T['stimstart_idx']] = np.nan        # no valid DV projection possible before stimulus


                # print('X',X.shape, 'Y',Y.shape, sum(mask_correct), sum(mask_error))

                mask_correct = np.array(mask_correct_all[cx][clx])
                mask_error = np.array(mask_error_all[cx][clx])

                for ex,mask in enumerate([mask_correct,mask_error]):

                    m = Y[mask,:].mean(axis=0)
                    s = Y[mask,:].std(axis=0)
                    e = s/np.sqrt(Y[mask,:].shape[0])

                    # use this for sem and std
                    axs.plot(times,m,color=colors[clx][ex], label='%s %s %d'%(label_action[clx],label_behav[ex],sum(mask)))
                    axs.fill_between(times,m-e,m+e,color=colors[clx][ex], alpha=0.3)
                    axs.fill_between(times,m-s,m+s,color=colors[clx][ex], alpha=0.1)

                    # use this for 2 std lines
                    # axs.plot(times,m,color=colors[clx][ex], lw=2, label='%s %s %d'%(label_action[clx],label_behav[ex],sum(mask)))
                    # axs.plot(times,m+s,color=colors[clx][ex], alpha=0.5)
                    # axs.plot(times,m-s,color=colors[clx][ex], alpha=0.5)

            axs.set_ylim(-1.24,1.24)
            axs.set_xlim(times[0],times[-1])
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0)

            axs.legend(frameon=False)
            if cx==0 or cx==4:
                axs.set_ylabel(['relevant, irrelevant\nDV projection','congruent, conflicting\nDV projection'][cx//4])
            axs.set_title(comparison,fontsize=20)
            if cx>3: axs.set_xlabel('time from stimulus onset [ms]')
                

        fig.suptitle(' %d mice, %d neurons, activity projected onto DV at each timepoint, rel-irrel,cong-confl*contexts, shades: sem, std'%(n_mice,n_neurons_sum))


        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'ricc-dvprojections,timeresolved_allmice,decode_%s-%dms'%(continuous_method,T['dt'].magnitude)+ext)



    return

















def aggregatemice_representationmeasures_featuresvector(datanames):
    # collect features for each mice
    # show the correlation structure and pcs of features
    # plot a raw feature vs feature space

    
    plotcorr = 0                  # plot correlations of features
    plotpca = 0                 # plot pca directions
    plotfvf = 1                 # plot all mice in a feature vs. feature space
    
    save = 1                  # save a plot
    
    
    n_mice = len(datanames)
    
    layeredfilename = 'all'
    stimulusspecificity = 'all'
    
    aucdeltat=750*pq.ms
    windowwidth = 100*pq.ms
    timeaveragewidth = 100*pq.ms
    timeaverage_range = ((np.array([windowwidth+T['dt'], windowwidth+timeaveragewidth+T['dt']])/T['dt']).magnitude).astype('int16')



    powers,bestfactors = pickle.load(open(cacheprefix+'tca/latentfactorpowers,decoders-%dmice.pck'%(14),'rb')) # n_mice = 14
    # powers, bestfactors: mouse x factors x [visual,audio,context]


    taskaspects = ['visual','audio','context','choice']
    cueperiodcharacteristics = ['mean','std dev','slope','sq.r. of residual']
    
    
    # build the feature vector as from rate differences, decoders and dim. red. methods.
    features = []
    feature_names = []
    feature_color = []
    pca_colors = np.array([['navy','seagreen','chocolate','darkred'],\
                           ['mediumblue','darkseagreen','orange','crimson']])
    
    for n,dn in enumerate(datanames):
        print(dn)

        blv,bla = preprocess.getorderattended(dn)


        featurevector = []


        # first the number of neurons:
        block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=False)
        n_cell_narrow = np.sum(block.annotations['waveforms']==0)
        n_cell_broad = np.sum(block.annotations['waveforms']==1)
        n_cells = n_cell_broad + n_cell_narrow
        featurevector.extend( [n_cell_broad, n_cell_narrow, n_cells] )
        if n==0:
            feature_names.extend( ['# broad sp.','# narrow sp.','# neurons'] )
        
        examine = 'allexpcond'
        for cx,comparison in enumerate(taskaspects):
            acrossdifference = pickle.load(open(cacheprefix+'continuous/responsedifference-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
    

    
            dm,de = nedi.get_maxnormdifference(acrossdifference)
            dz = np.max((np.zeros(dm.shape),de),axis=0)       # clip below zero
            
            diffsignal = neo.AnalogSignal(dz,t_start=dm.t_start,sampling_period=dm.sampling_period,units=dm.units)
            
            decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # lower confidence bound

            diffauc_on = neph.get_auctimenormalized(diffsignal,T['stimstarttime'],T['stimendtime'],'diff.on.'+comparison)[0].magnitude
            decauc_on = neph.get_auctimenormalized(decsignal,T['stimstarttime'],T['stimendtime'])[0]

            diffauc_pre = neph.get_auctimenormalized(diffsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'],'diff.off.'+comparison)[0].magnitude
            decauc_pre = neph.get_auctimenormalized(decsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]

            diffauc_post = neph.get_auctimenormalized(diffsignal,T['stimendtime'],T['stimendtime']+aucdeltat,'diff.off.'+comparison)[0].magnitude
            decauc_post = neph.get_auctimenormalized(decsignal,T['stimendtime'],T['stimendtime']+aucdeltat)[0]

            featurevector.extend( [  diffauc_on,   decauc_on, \
                                   diffauc_pre,  decauc_pre, \
                                   diffauc_post, decauc_post   ])
            if n==0:
#                if not publish:
                if 0:
                    feature_names.extend(\
                            [ comparison+' diff on', comparison+' dec on', \
                            comparison+' diff pre', comparison+' dec pre', \
                            comparison+' diff post', comparison+' dec post' ]  )
                else:
                    feature_names.extend(\
                            [ comparison+' diff on', comparison+' decoder, ON', \
                            comparison+' diff pre', comparison+' decoder, PRE', \
                            comparison+' diff post', comparison+' decoder, POST' ]  )


        A = pickle.load(open(cacheprefix+'subspaces/angles,aspects,subspaces-%s_%s-%s-%dms%dms_%s.pck'%(examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn),'rb'))
            
        for cx1,comparison1 in enumerate(taskaspects):
            for cx2,comparison2 in enumerate(taskaspects):
                if cx2>cx1:
                    F_on   = neph.get_auctimenormalized( A[cx1,cx2], T['stimstarttime'],T['stimendtime'] )[0]/np.pi*180
                    F_pre  = neph.get_auctimenormalized( A[cx1,cx2], T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]/np.pi*180
                    F_post = neph.get_auctimenormalized( A[cx1,cx2], T['stimendtime'],T['stimendtime']+aucdeltat)[0]/np.pi*180
                    featurevector.extend( [ F_on, F_pre, F_post ] )
                    if n==0:
                        feature_names.extend( [ comparison1+'-'+comparison2+' on',\
                                                comparison1+'-'+comparison2+' pre',\
                                                comparison1+'-'+comparison2+' post'  ] )
            



        crossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,crosscontext_%s-%s.pck'%(dn,continuous_method),'rb'))
        ignoresignal = crossdecoder[1][:,0]
        attendsignal = crossdecoder[2][:,0]
        decignore = neph.get_auctimenormalized(ignoresignal,T['stimstarttime'],T['stimendtime'])[0]
        decattend = neph.get_auctimenormalized(attendsignal,T['stimstarttime'],T['stimendtime'])[0]
        featurevector.extend(  [  decattend - decignore  ]  )
        if n==0: feature_names.extend(  ['visual decoder, ignore-trained crosstest\nattend-ignore, ON']  )




        examine = 'attendignore'
        attentionstates = ['attend','ignore']
        for cx,comparison in enumerate([' visual',' audio']):
            diffauc_on = []; decauc_on = []
            for bx,bl in enumerate(attentionstates):
                acrossdifference = pickle.load(open(cacheprefix+'continuous/responsedifference-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,bl+comparison,stimulusspecificity),'rb'))
                acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,bl+comparison,stimulusspecificity),'rb'))
    
                dm,de = nedi.get_maxnormdifference(acrossdifference)
                dz = np.max((np.zeros(dm.shape),de),axis=0)       # clip below zero
                
                diffsignal = neo.AnalogSignal(dz,t_start=dm.t_start,sampling_period=dm.sampling_period,units=dm.units)
                
                decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # change to 2 for lower confidence bound or mean to the left
    
                diffauc_on.append(neph.get_auctimenormalized(diffsignal,T['stimstarttime'],T['stimendtime'],'diff.on.'+comparison)[0].magnitude)
                decauc_on.append(neph.get_auctimenormalized(decsignal,T['stimstarttime'],T['stimendtime'])[0])

            featurevector.extend( [  diffauc_on[0]-diffauc_on[1], decauc_on[0]-decauc_on[1]])
            if n==0:
#                if not publish:
                if 0:
                    feature_names.extend(\
                            [ attentionstates[0]+'-'+attentionstates[1]+comparison+' diff on', \
                              attentionstates[0]+'-'+attentionstates[1]+comparison+' dec on' ]  )
                else:
                    feature_names.extend(\
                            [ comparison[1:]+' difference\n'+attentionstates[0]+'-'+attentionstates[1]+', ON', \
                              comparison[1:]+' decoder\n'+attentionstates[0]+'-'+attentionstates[1]+', ON' ]  )




        for bl,blname in enumerate(['cue block 1', 'cue block 3']):
            blockorder = int(([blv[0],bla[0]][bl]-1)/2)
            for ptrx, ptname in enumerate([' pre ',' on ']):# ,trainingpoint in enumerate(trainingpoints):
                dtm,dtstd,dtsem,n_cuetrials = nedi.get_crossgradient(ptrx,bl,dn,timeaverage_range,continuous_method)
                m = dtm[:,0].mean()
                s = dtm[:,0].std()
#                print(m,len(dtm),len(np.arange(len(dtm))),'shape',dtm.shape)
                l = sp.stats.linregress(np.arange(len(dtm)),dtm[:,0])
                esr = np.sqrt(nedi.residuals(np.arange(len(dtm)),dtm[:,0],l[0],l[1]))
                featurevector.extend(   [m,s,l[0],esr]   )
                if n==0: feature_names.extend( [blname+ptname+'mean', blname+ptname+'std', blname+ptname+'slope',blname+ptname+'error'     ]   )



        # TCA powers, bestfactors: mouse x factors x [visual,audio,context]
        n_bestfactors = 3          # use only the first 3 factors
        for cx,comparison in enumerate(['visual','audio','context','choice']): 
#            print(len(featurevector))
            featurevector.extend( powers[n,bestfactors[n,:n_bestfactors,cx],cx])
#            print(len(featurevector))
            if n==0: feature_names.extend(  [ 'tca %s on #%d'%(comparison,nbf+1)  for nbf in range(n_bestfactors) ]  )



    
        # decode task variables from runspeed
        for cx,comparison in enumerate(['visual','audio','context','choice']): 
            acrossdecoder = pickle.load(open(cacheprefix+'phys/rundecodes-%s-%s-%s.pck'%(dn,continuous_method,taskaspects[cx]),'rb'))

            rundecsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # change to 2 for lower confidence bound or mean to the left
            rundecauc_pre = neph.get_auctimenormalized(rundecsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]
            rundecauc_early = neph.get_auctimenormalized(rundecsignal,T['stimstarttime'],T['stimendtime']/2+T['stimstarttime']/2)[0]
            rundecauc_late = neph.get_auctimenormalized(rundecsignal,T['stimendtime']/2+T['stimstarttime']/2,T['stimendtime'])[0]
            rundecauc_post = neph.get_auctimenormalized(rundecsignal,T['stimendtime'],T['stimendtime']+aucdeltat)[0]

            featurevector.extend( [  rundecauc_pre, rundecauc_early, rundecauc_late, rundecauc_post   ])
            if n==0:
                feature_names.extend(  [ comparison + ' run dec %s'%['pre','early','late','post'][i]\
                                                 for i in range(4)  ]  )

        # add each mice to the full features data matrix
        features.append(np.array(featurevector))


    
    
    print(np.array(features).shape, len(datanames), len(feature_names))
    print(feature_names)
    
    
    
    features = pd.DataFrame(features,index=datanames,columns=feature_names)
    feature_names = pd.Series(feature_names)
    
#    print('raw feature values')
#    print(features)
    
#    features = features.apply(sp.stats.zscore)
    
#    print('z-scored feature values')
#    print(features)
    


    features_pca,expl_var,pca_components = nedi.pca(features,n_components=min(4,n_mice))
    print(pca_components.shape)
    
    
    
    feature_corr = np.corrcoef(features.T)

    order = np.argsort( np.abs(feature_corr).mean(axis=0)  )[::-1]
    print(order)
    feature_corr_sorted = np.corrcoef(features.values[:,order].T)





    ticks = np.arange(len(feature_names),dtype='int16')
    
    if plotcorr:         # Correleation of features
        fig,ax = plt.subplots(1,2,figsize=(36,18))
        for fx in [0,1]:
    
            axs = ax[fx]
        
            display = [feature_corr,feature_corr_sorted][fx]
        
            axs.imshow(display,cmap=figs.getcorrcolormap())
        
        
            if fx==0:
                ticklabels = feature_names
                axs.set_title('feature correlations, sorted thematically')
            else:
                ticklabels = feature_names.values[order]
                axs.yaxis.tick_right()
                axs.set_title('feature correlations, sorted by <|corr.|>')
            axs.set_xticks(ticks)
            axs.set_xticklabels(ticklabels,rotation=90)
            axs.set_yticks(ticks)
            axs.set_yticklabels(ticklabels)
        fig.suptitle('Correlation of neural representation measures amongst %d mice'%n_mice)
        if save:
            fig.savefig(resultpath+'micetribe,reprfeat-space,corr_%s-%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
    
    
    
        
    if plotpca:      # principal components A
        fig,ax = plt.subplots(3,2,figsize=(32,48))
        pcacs = [0,1,2,3]    # to visualize from these list
        for fx in [0,2]:
        
            axs = ax[0,int(fx/2)]
            
            for n,mousename in enumerate(datanames):
                axs.scatter(features_pca[n,pcacs[fx]],features_pca[n,pcacs[fx+1]],s=150,marker='+')
                axs.text(features_pca[n,pcacs[fx]]+0.05,features_pca[n,pcacs[fx+1]]-0.02,mousename)
        #        axs.legend(datanames)
            axs.set_xlabel('features, principal component %d, expl. var. %4.2f'%(pcacs[fx]+1,expl_var[pcacs[fx]]))
            axs.set_ylabel('features, principal component %d, expl. var. %4.2f'%(pcacs[fx+1]+1,expl_var[pcacs[fx+1]]))
        fig.suptitle('Individual mice in the PCA-transformed feature space\ncontributions of features to the principal components')
    
    
        for fx in [0,1,2,3]:
            axs = ax[1+int(fx/2),fx%2]
            order_pca = np.argsort( np.abs(pca_components[pcacs[fx],:]) )[::-1]
            axs.bar( np.arange(len(feature_names)),np.abs(pca_components[pcacs[fx],order_pca]),\
                     color=pca_colors[(np.sign(pca_components[pcacs[fx],order_pca])==-1).astype('int16'),fx]  )
            axs.set_xticks(ticks)
            axs.set_xticklabels(feature_names[order_pca],rotation=90)
            axs.set_ylabel('coefficient for principal component %d     dark: +, light: -'%(fx+1))
    
        if save:
            fig.savefig(resultpath+'micetribe,reprfeat-space,pca_%s-%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
    
    
    if plotpca:      # principal components B
        pc_names=[]
        fig,ax = plt.subplots(1,1,figsize=(32,16))
        axs = ax
        for fx in [0,1,2,3]:
            axs.bar( np.arange(len(feature_names))+(fx-1.5)/5,pca_components[pcacs[fx]],\
                    width=0.19, color=pca_colors[(np.sign(pca_components[pcacs[fx],:])==-1).astype('int16'),fx] )
            pc_names.extend(['PC %2d'%(fx+1)])
        axs.legend(pc_names)
        axs.set_xticks(ticks)
        axs.set_xticklabels(feature_names,rotation=90)
        axs.set_ylim(-0.5,0.5)
        axs.set_ylabel('coefficient for principal components')
        fig.suptitle('Contributions of features to the principal components')
    
        if save:
            fig.savefig(resultpath+'micetribe,reprfeat-space,pca,coeffs_%s-%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
    



    if plotfvf:      # individual features vs. features
        sets = [6,7]#,3,4,5]
        sets = [0]
        sets = [8]
        sets = [1]
        
        for sx in sets:
            if sx==0:
                code = 'attend,ignore'
                x_names = ['context decoder, PRE']
                y_names = [ \
#                            ['visual dec on'],['attend-ignore visual dec on'],\
#                            ['audio dec on'],['attend-ignore audio dec on'] ]
                            ['context decoder, ON'],\
                            ['choice decoder, ON'],\
                            ['visual decoder, ON'],\
                            ['visual decoder\nattend-ignore, ON'],\
                            ['visual decoder, ignore-trained crosstest\nattend-ignore, ON']  ]
                colors = ['darkorange']
                figuregrouptitle = 'Pre-stimulus context decodability vs. attended context. Difference between attended and ignored decodability for each modality'
                if publish: colors = ['mediumvioletred','darkorange','navy','steelblue']; figuregrouptitle = None         #steelblue, mediumseagreen
                            
            elif sx==1:
                code = 'taskvariables'
                x_names = ['context dec pre']
                y_names = [ \
#                            ['visual dec on','tca visual on #1','tca visual on #2','tca visual on #3'],\
#                            ['audio dec on','tca audio on #1','tca audio on #2','tca audio on #3'],\
#                            ['context dec on','tca context on #1','tca context on #2','tca context on #3'],\
#                            ['choice dec on','tca choice on #1','tca choice on #2','tca choice on #3']   ]
                            ['visual dec on'],\
                            ['audio dec on'],\
                            ['context dec on'],\
                            ['choice dec on']   ]
                colors = ['darkorange','purple','mediumvioletred','hotpink']
                figuregrouptitle = 'Pre-stimulus context decodability vs. task conditions'
                if publish: colors = ['navy','darkgreen','mediumvioletred','darkorange','darkred']; figuregrouptitle=None
            
            elif sx==2:
                x_names = ['context dec pre']
                y_names = [ ['cue block 1 pre mean','cue block 3 pre mean','cue block 1 on mean','cue block 3 on mean'] ]#,\
#                            ['cue block 1 pre std','cue block 3 pre std','cue block 1 on std','cue block 3 on std'],\
#                            ['cue block 1 pre slope','cue block 3 pre slope','cue block 1 on slope','cue block 3 on slope'],\
#                            ['cue block 1 pre error','cue block 3 pre error','cue block 1 on error','cue block 3 on error'] ]
                colors = ['firebrick','darkslategray','mistyrose','paleturquoise']
                figuregrouptitle = 'Pre-stimulus context decodability vs. cue period single trial context change;\ntraining in multimodal blocks pre stimulus'
            
            elif sx==3:
                code = 'neurons,broad'
                x_names = ['# broad sp.']
                y_names = [ ['visual dec on','tca visual on #1','tca visual on #2','tca visual on #3'],\
                            ['audio dec on','tca audio on #1','tca audio on #2','tca audio on #3'],\
                            ['context dec on','tca context on #1','tca context on #2','tca context on #3'],\
                            ['choice dec on','tca choice on #1','tca choice on #2','tca choice on #3']   ]
                colors = ['darkorange','purple','mediumvioletred','hotpink']
                figuregrouptitle = 'Number of recorded broad spiking (excitatory) neurons vs. task conditions'

            elif sx==4:
                code = 'neurons,narrow'
                x_names = ['# narrow sp.']
                y_names = [ ['visual dec on','tca visual on #1','tca visual on #2','tca visual on #3'],\
                            ['audio dec on','tca audio on #1','tca audio on #2','tca audio on #3'],\
                            ['context dec on','tca context on #1','tca context on #2','tca context on #3'],\
                            ['choice dec on','tca choice on #1','tca choice on #2','tca choice on #3']   ]
                colors = ['darkorange','purple','mediumvioletred','hotpink']
                figuregrouptitle = 'Number of recorded broad spiking (excitatory) neurons vs. task conditions'

            elif sx==5:
                code = 'neurons'
                x_names = ['# neurons']
                y_names = [ ['visual dec on','tca visual on #1','tca visual on #2','tca visual on #3'],\
                            ['audio dec on','tca audio on #1','tca audio on #2','tca audio on #3'],\
                            ['context dec on','tca context on #1','tca context on #2','tca context on #3'],\
                            ['choice dec on','tca choice on #1','tca choice on #2','tca choice on #3']   ]
                colors = ['darkorange','purple','mediumvioletred','hotpink']
                figuregrouptitle = 'Number of recorded neurons vs. task conditions'

            elif sx==6:         # plot aspect decision normal angles in neural vectorspace
                x_names = ['context dec pre']
                y_names =[ ['visual-audio on'],\
                           ['visual-context on','audio-context on'],\
                           ['visual-choice on','audio-choice on'],\
                           ['context-choice on','context-choice pre', 'context-choice post']    ]
                colors = ['purple','seagreen','mediumvioletred','']
                figuregrouptitle = 'Decision normal angles from prestimulus context decodability'

            elif sx==7:         # plot aspect decision normal angles in neural vectorspace
                x_names = ['context-choice pre']
                y_names =[ ['visual-audio on'],\
                           ['visual-context on','audio-context on'],\
                           ['visual-choice on','audio-choice on'],\
                           ['context-choice on','context-choice post']    ]
                colors = ['purple','seagreen','mediumvioletred','']
                figuregrouptitle = 'Decision normal angles'

            elif sx==8:         # plot run speed decodability
                takechoice = True
                if takechoice:    # choose choice or context, here choice on
                    code = 'taskdecodingfromrunspeed,choicepost'
                    x_names = ['choice decoder, POST']
                    figuregrouptitle = 'Task variables decoded from runspeed vs. choice from neural on'
                else:                                           # here context pre
                    code = 'taskdecodingfromrunspeed,contextpre'
                    x_names = ['context decoder, PRE']
                    figuregrouptitle = 'Task variables decoded from runspeed vs. context from neural pre'
                trajectorygrouplabels = ['pre','early','late','post']
                y_names =[ ['%s run dec %s'%(comparison,trajectorygrouplabels[k])\
                                         for k in range(len(trajectorygrouplabels))]  \
                                     for comparison in taskaspects ]
                colors = ['navy','darkgreen','mediumvioletred','darkorange']


            xx = list(feature_names).index(x_names[0])
            if not publish: fig,ax = plt.subplots(1,len(y_names),figsize=(10*len(y_names),8))
            else: #fig,ax = plt.subplots(1,len(y_names),figsize=(9*len(y_names),7))
                fig,ax = plt.subplots(2,2,figsize=(12,11))
            if sx==8:
                fig,ax = plt.subplots(4,4,figsize=(32,32))
                
            for yfx, y_figgroup in enumerate(y_names):
                if len(y_names)>1:
                    if publish: axs = ax[int(yfx/2),yfx%2]
                    else: axs = ax[yfx]
                else: axs = ax
                for ynx, y_name in enumerate(y_figgroup):
                    alpha= 0.8
                    if not publish: color = colors[ynx]
                    else: color = colors[yfx]
                    if sx==8:
                        alpha = 1-yfx*0.15
                        color = colors[yfx]
                        axs = ax[yfx,ynx]
                    if yfx==1: datanamescorrected = datanames[:-3]
                    else: datanamescorrected = datanames
                    for n,mousename in enumerate(datanamescorrected):
                        yx = list(feature_names).index(y_name)
                        
                        
                        axs.scatter(features.iloc[n,xx],features.iloc[n,yx],s=150,marker='o',color=color,alpha=alpha)
#                        axs.text(features.iloc[n,xx]*1.05,features.iloc[n,yx]*0.98,mousename,color=colors[ynx],alpha=0.1)
                        if ynx==0 or sx==8: axs.text(features.iloc[n,xx]*1.02,features.iloc[n,yx]*0.995,mousename,color='k',alpha=0.3,fontsize=10)
#                        if ynx==0: axs.text(features.iloc[n,xx],tribe_stds[:,0,:,pox].min()*0.95,mousename,color='k',alpha=0.3,rotation=90,fontsize=10)
                        
                        
                        axs.set_xlabel(x_names[0])

                    l = sp.stats.linregress(features.iloc[:,xx],features.iloc[:,yx])
                    if l[3]<0.13:
                        line = np.linspace(start=features.iloc[:,xx].min()*0.8,stop=features.iloc[:,xx].max()*1.2,num=2)
                        axs.plot(line,l[0]*line+l[1],color=color,alpha=alpha,linewidth=2,label=y_name)
                    else:
                        line = np.linspace(start=features.iloc[:,xx].min()*0.8,stop=features.iloc[:,xx].max()*1.2,num=2)
                        axs.plot(line,l[0]*line+l[1],'--',color=color,alpha=alpha,linewidth=2,label=y_name)
                        
                    xoffs = 0.4; yoffs = 0
                    if sx in [3,4,5,7]: xoffs = 25
                    if sx==8: xoffs = 0.6; yoffs = 0.3
                    if l[3]<0.001:
                        axs.text(line[1]-xoffs,l[0]*line[1]+l[1]*1.02+yoffs,'p<%5.3f, $R^2$=%4.2f'%((0.001,l[2]**2)),color=color)
                    else:
                        axs.text(line[1]-xoffs,l[0]*line[1]+l[1]*1.02+yoffs,'p=%5.3f, $R^2$=%4.2f'%((l[3],l[2]**2)),color=color)
                    
                    if sx==8:
                        axs.legend()
                        axs.set_xlim(0.45,1.01); axs.set_ylim(0.45,1.01)
                        axs.set_xticks([0.5,1.0]); axs.set_yticks([0.5,1.0])
                        figs.plottoaxis_crosshair(axs,0.5,0.5)
                        if ynx==0: axs.set_ylabel(taskaspects[yfx])
                        if yfx==0: axs.set_title(trajectorygrouplabels[ynx])
                if not publish:
                    if not sx==8: axs.legend()
                    if sx==0: axs.set_title(['visual','attend-ignore visual','audio','attend-ignore audio','attend-ignore ignore crosstest'][yfx])
                    if sx==1 or sx in [3,4,5]: axs.set_title(taskaspects[yfx])
                    if sx==2: axs.set_ylabel('cue period prob. '+cueperiodcharacteristics[yfx])
                else:
                    figs.plottoaxis_notickmarks(axs)
                    if sx==0 and yfx==3: figs.plottoaxis_chancelevel(axs,0.0)# axs.set_yticks([0.]); axs.set_yticklabels('0')
                    axs.set_ylabel(y_name)
                
            fig.suptitle(figuregrouptitle)
            if save:
                if not publish:
                    fig.savefig(resultpath+'micetribe,reprfeat-space,feature-comparison,%s_%s-%dms%dms'%(code,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
                else:
                    if   sx==0: fig.savefig(resultpath+'7-contextdecpre,CCVVattendignore-tribe'+ext)
                    elif sx==1: fig.savefig(resultpath+'7A-contextdecpre,VACC-tribe'+ext)













def aggregatemice_subspaceangles(datanames,examine='allexpcond'):
    # this section draws the angles between aspects in "taskaspects" for multiple mice
    taskaspects = ['visual','audio','context','choice']

    lx1  = [0,   2, 2,   3, 3,  2]
    lx2  = [1,   0, 1,   0, 1,  3]
    xpos = [0,   1.5, 2.5,   4, 5,  6.5]
    colors = [ 'mediumturquoise',  'purple',      'deeppink',    'darkcyan',    'darkkhaki',    'maroon'  ]
    labels = [ 'visual\naudio', 'visual\ncontext','audio\ncontext',  'visual\nchoice','audio\nchoice', 'context\nchoice' ]

    # collect all the angles for all mice
    AT = np.zeros( (len(datanames),len(lx1)) )
    for n,dn in enumerate(datanames):
        A = pickle.load(open(cacheprefix+'subspaces/angles,aspects,subspaces-%s_%s-%s-%dms%dms_%s.pck'%(examine,'all',continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn),'rb'))
        for lx in range(len(lx1)):
            AT[n,lx] = A[lx1[lx],lx2[lx]][T['stimstart_idx']:T['stimend_idx']].mean()   # evoked only, on stimulus
    

    fig,ax = plt.subplots(1,1,figsize=(15,12))
    axs = ax
    box = axs.boxplot(x=AT,positions=xpos,whis=[5,95],labels=labels)
#    for e, c in zip(box['boxes'], colors):  e.set_color(c)

    
    axs.set_ylim([np.pi/4,np.pi/2])
#    figs.plottoaxis_chancelevel(axs,np.pi/2,'orthogonal')
    axs.set_yticks([np.pi/4,np.pi/2]); axs.set_yticklabels(['$45\degree$','$90\degree$'])
    axs.set_ylabel('decision normal angles')
                

#    fig.suptitle(dn+', decision normal angles in neural vectorspace, between aspects')
    
    save = 1
    if save:
        # for publication:
        fig.savefig(resultpath+'4A-aspects,decoders,subspaceangles'+ext)











def aggregatemice_spontaneousdecoder(datanames,examine='allexpcond'):
    # compare various prestimulus and on stimulus context, stimulus and choice representations

    n_mice = len(datanames)
    
    if examine=='allexpcond':
        arealabels = ['on stimulus','pre stimulus','across time']
        taskaspects = ['visual','audio','context','choice']
        crosslabels = ['context\nPRE->PRE','context\nPRE->ON']
        tasklabels=taskaspects
        notch=False
    elif examine=='attendignore':
        arealabels = ['on stimulus']
        taskaspects = ['attend visual','ignore visual','attend audio','ignore audio']
        tasklabels = ['attend\nvisual','ignore\nvisual','attend\naudio','ignore\naudio']
        notch=True
    
    stimulusspecificity = 'all'
    layeredfilename = 'all'

    trainingoffsets = np.array([0,375,750,1225])*pq.ms
    trainingpoints = T['starttime'] + trainingoffsets

    
    P = np.zeros( (n_mice,len(taskaspects),len(arealabels)) )
    S = np.zeros( (n_mice,len(taskaspects),len(arealabels)) )
    n_trials_av = np.zeros((n_mice)); n_trials_aa = np.zeros((n_mice))
    for n,dn in enumerate(datanames):
        offstimdecoderpre = pickle.load(open(cacheprefix+'continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,'context',1),'rb'))[1]
        offstimdecoderon = pickle.load(open(cacheprefix+'continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,'context',3),'rb'))[1]
        for cx,comparison in enumerate(taskaspects):
            decodertrajectory = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))[1] # use only test
            for ax,al in enumerate(arealabels):
                if al=='pre stimulus':
                    P[n,cx,ax] = decodertrajectory[T['start_idx']:T['stimstart_idx'],0].mean()
                elif al=='on stimulus':
                    P[n,cx,ax] = decodertrajectory[T['stimstart_idx']:T['stimend_idx'],0].mean()
                    S[n,cx,ax] = decodertrajectory[T['stimstart_idx']:T['stimend_idx'],2].mean() +\
                       decodertrajectory[T['stimstart_idx']:T['stimend_idx'],0].std()/np.sqrt(len(decodertrajectory[T['stimstart_idx']:T['stimend_idx'],0]))   # 1* s.e.m.
                elif al=='across time' and cx==0:
                    P[n,0,ax] = offstimdecoderpre[20:40,0].mean()
                    P[n,1,ax] = offstimdecoderon[20:40,0].mean()


    if examine=='allexpcond':
        fig = plt.figure(constrained_layout=False,figsize=(16,4))
        gs = fig.add_gridspec(1, 8)
        ax = []
        ax.append( fig.add_subplot(gs[0, 0:3]) )
        ax.append( fig.add_subplot(gs[0, 3:6]) )
        ax.append( fig.add_subplot(gs[0, 6:8]) )
#        fig,ax = plt.subplots(1,len(arealabels),figsize=(7.5*len(arealabels),5))
    
        for ex,al in enumerate(arealabels):
            if len(arealabels)>1: axs = ax[ex]
            else: axs = ax
            imax = 4
            if al=='across time': imax = 2
            
            ci = P[:,:imax,ex].std(axis=0)*2/np.sqrt(len(datanames))
            m = P[:,:imax,ex].mean(axis=0)
            ci = np.c_[m-ci,m+ci]
            
            if ex<2: labels=tasklabels; widths=0.5
            else: labels = crosslabels; widths=0.38
            axs.boxplot(x=P[:,:imax,ex],notch=notch,usermedians=m,conf_intervals=ci,whis=[5,95],labels=labels,widths=widths)

            if ex==0: axs.set_yticks([0.5,1.0])
            else: axs.set_yticks([])
            axs.set_ylim(0.475,1.)
            figs.plottoaxis_chancelevel(axs,0.5)
#            if ex==0: axs.set_ylabel('accuracy')

        fig.subplots_adjust(bottom=0.20)



    elif examine=='attendignore':
        if 0:
            # calculate confidence interval for accuracy, at 2 std
            block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=False)
            ixv45,ixv135,ixa45,ixa135,ixv5000,ixv10000,ixa5000,ixa10000 = preprocess.getstimulusidents(dn,block,multimodalonly=True)
            n_trials_av[n] = np.count_nonzero(ixv45)  +  np.count_nonzero(ixv135)
            n_trials_aa[n] = np.count_nonzero(ixa5000)  +  np.count_nonzero(ixa10000)
    
            S = 1.96 * np.sqrt(P * (1-P) /  0.25*np.c_[ n_trials_av, n_trials_aa, n_trials_aa, n_trials_av][:,:,None]   )
    
            n_colors = n_mice
            colorlist = plt.cm.viridis( np.linspace(0, 0.8, n_colors) )
            fig,ax = plt.subplots(2,2,figsize=(7.5*2,5*2))
            for lx in [0,1]:
                axs = ax[0,lx]
                for n in range(n_mice):
                    axs.plot(np.array([0,1])+n/80,P[n,lx:(lx+2),0],'o-',color=colorlist[n])
                    axs.errorbar(np.array([0,1])+n/80, P[n,lx:(lx+2),0], yerr=S[n,:2,0],color=colorlist[n], capsize=2 )
    #                if lx==0: axs.text(0, P[n,lx:(lx+2),0], datanames[n],color=colorlist[n])
    #                axs.fill_between([0,1], m-s*2, m+s*2, color=colorlist,alpha=1/3*alpha**2)#,label=label)
    #        if examine=='attendignore': axs.set_xticklabels(tasklabels,rotation=15)
    #        elif examine=='allexpcond': axs.set_xticklabels(tasklabels); axs.set_title(arealabels[ex])
                axs.set_yticks([0.5,1.0])
                axs.set_ylim(0.475,1.)
                figs.plottoaxis_chancelevel(axs,0.5)
                
                axs = ax[1,lx]
                AD = P[:,:lx,0]-P[:,lx:(lx+2),0]


        # create distribution for signed distance from hyperplane:
        # hyperplane = perpendicular to db;
        # get center of mass of all trials, M
        # run line through center of mass
        # add distances
        






    save = 0
    if save:
        # for publication:
        if examine=='allexpcond':
            fig.savefig(resultpath+'5-decoder,pre+onstimulus-tribe'+ext)
        elif examine=='attendignore':
            fig.savefig(resultpath+'2-decoder,attendignore-tribe'+ext)


#    for n,dn in enumerate(datanames):
#        n_neuron = block.segments[0].analogsignals[0].shape[1]
#        for cx,comparison in enumerate(taskaspects):
#            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
#            wx = int((len(acrossdecoder)-7)/n_neuron)
#            c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
#        
#        c_db = np.array(c_db)







def aggregatemice_crosscontextdecoder(datanames):



    correctonly = 0        # all trials



    taskaspects = ['attend visual','ignore visual']
    filename = ['crosscontext,attendattend','crosscontext,ignoreignore','crosscontext,attendignore','crosscontext,ignoreattend']

    times = np.arange(-1500,4510,10)[:596]*pq.ms

    labels = np.hstack((datanames,'all'))

    crossdecoders_all = []
    for n,dn in enumerate(datanames):

        crossdecoders = []  
        projections = []  #  (mice) x {trainedattend,trainedignore} x {prestim,onstim} x {testedignore,testedattend} x {go,nogo}
        variances = []  # (mice) x {train: attend,ignore} x {test: attend,ignore} x {class 1, class 2, class 1+2}

        for trainedwhere in [0,1]:
            trainidx = trainedwhere
            testidx = 1-trainedwhere
            filenameidx = np.array([[0,2],[3,1]])[trainidx,testidx]

            crossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,%s%s_%s-%s.pck'%(filename[filenameidx],['',',correctonly'][correctonly],dn,continuous_method),'rb'))
            crossdecoders.append(crossdecoder)

            projection,variance = pickle.load(open(cacheprefix+'continuous/%s%s,projections,variance_%s-%s.pck'%('cross',['',',correctonly'][correctonly],dn,continuous_method),'rb'))
            projections.append(projection)
            variances.append(variance)


            # print(dn,len(crossdecoder),len(crossdecoder[trainedwhere]),len(crossdecoder[trainedwhere][0]))
        
        crossdecoders_all.append(crossdecoders)



    fig,ax = plt.subplots(4,3,figsize=(3*1.41*8,4*8))


    # cross decoders
    for trainedwhere in [0,1]:
        M = []
        for n,dn in enumerate(datanames):
            axs = ax[0,trainedwhere]
            axs.plot(times,crossdecoders_all[n][trainedwhere][1][:,0], color='black')
            axs.plot(times,crossdecoders_all[n][trainedwhere][2][:,0], color='red')
            
            axs = ax[1,trainedwhere]
            m = crossdecoders_all[n][trainedwhere][1][:,0]-crossdecoders_all[n][trainedwhere][2][:,0]
            axs.plot(times,m, label=dn)
        
            M.append(m)

            axs = ax[2,trainedwhere]    
            axs.boxplot(positions=[n], x=m[150:450], notch=True, whis=[5,95], showfliers=False)
            
            axs = ax[3,trainedwhere]    
            axs.boxplot(positions=[n], x=m[150:200], notch=True, whis=[5,95], showfliers=False)


        M = np.array(M)

        axs = ax[2,trainedwhere]    
        axs.boxplot(positions=[n+2], x=np.vstack(M[:,150:450]), notch=True, whis=[5,95], showfliers=False)
        
        axs = ax[3,trainedwhere]    
        axs.boxplot(positions=[n+2], x=np.vstack(M[:,150:200]), notch=True, whis=[5,95], showfliers=False)


        ax[0,trainedwhere].set_ylim(0.45,1.05)
        figs.plottoaxis_stimulusoverlay(ax[0,trainedwhere],T)
        ax[0,trainedwhere].legend(['same','cross'],frameon=False)


        figs.plottoaxis_stimulusoverlay(ax[1,trainedwhere],T)


        figs.plottoaxis_chancelevel(ax[2,trainedwhere])
        ax[2,trainedwhere].set_ylim(-0.3,0.3)
        ax[2,trainedwhere].set_xticklabels([])


        figs.plottoaxis_chancelevel(ax[3,trainedwhere])
        ax[3,trainedwhere].set_ylim(-0.3,0.3)
        ax[3,trainedwhere].set_xticklabels(labels,rotation=60)






    ax[0,0].set_title('visual context')
    ax[0,1].set_title('audio context')
    ax[0,0].set_ylabel('accuracy')

    ax[1,0].set_ylabel('same$-$cross\naccuracy difference\nover time')

    ax[2,0].set_ylabel('same$-$cross\naccuracy difference\ndistribution, on stim')

    ax[3,0].set_ylabel('same$-$cross\naccuracy difference\ndistribution, on stim 0-0.5 sec')





    ax[1,2].set_ylabel('visual$-$audio\naccuracy difference\nover time')

    ax[2,2].set_ylabel('visual$-$audio\naccuracy difference\ndistribution, on stim')

    ax[3,2].set_ylabel('visual$-$audio\naccuracy difference\ndistribution, on stim 0-0.5 sec')










    
    # difference between same decoders
    M = []
    for n,dn in enumerate(datanames):
        axs = ax[1,2]
        m = crossdecoders_all[n][1][1][:,0]-crossdecoders_all[n][0][1][:,0]    # (animal)(trainblock)(tr,te,cte)(times,stats)
        axs.plot(times,m, label=dn)
        M.append(m)


        axs = ax[2,2]
        axs.boxplot(positions=[n], x=m[150:450], notch=True, whis=[5,95], showfliers=False)
        
        axs = ax[3,2]
        axs.boxplot(positions=[n], x=m[150:200], notch=True, whis=[5,95], showfliers=False)

        



    M = np.array(M)
    
    ax[1,2].set_title('visual$-$audio')

    axs = ax[2,2]
    axs.boxplot(positions=[n+2], x=np.vstack(M[:,150:450]), notch=True, whis=[5,95], showfliers=False)
    
    axs = ax[3,2]
    axs.boxplot(positions=[n+2], x=np.vstack(M[:,150:200]), notch=True, whis=[5,95], showfliers=False)



    figs.plottoaxis_stimulusoverlay(ax[1,2],T)



    figs.plottoaxis_chancelevel(ax[2,2])
    ax[2,2].set_ylim(-0.3,0.3)
    ax[2,2].set_xticklabels([])


    figs.plottoaxis_chancelevel(ax[3,2])
    ax[3,2].set_ylim(-0.3,0.3)
    ax[3,2].set_xticklabels(labels,rotation=60)


    S = np.vstack(M[:,150:200])
    t,p = sp.stats.ttest_1samp(S,0)
    print('timepoints: t',t,'p',p,'m+/-sem',S.mean(),'+/-', S.std()/np.sqrt(len(S)), 'std', S.std())
    R = M[:,150:200].mean(axis=1)
    t,p = sp.stats.ttest_1samp(R,0)
    print('mice:       t',t,'p',p,'m+/-sem',R.mean(),'+/-',R.std()/np.sqrt(len(R)),'std', R.std())


    ax[0,2].remove()


    fig.suptitle('visual decoding within and across contexts')
    fig.tight_layout()


    save = 0 or globaldoplot
    if save:
        fig.savefig(resultpath+'crosscontext,visual-tribe'+ext)

    return










def aggregatemice_acrosscontextcomparison(datanames):
    




    taskaspects = ['attend visual','ignore visual']
    filename = ['crosscontext,attendattend','crosscontext,ignoreignore','crosscontext,attendignore','crosscontext,ignoreattend']

    times = np.arange(-1500,4510,10)[:596]*pq.ms

    labels = np.hstack((datanames,'all'))

    acrossdecoders_all = []       # (mice)(taskaspects)(times,stats)
    for n,dn in enumerate(datanames):

        taskaspects = ['visual,av','visual,aa','audio,aa','audio,av']
        
        acrossdecoders = []
        for cx,comparison in enumerate(taskaspects):
            acrossdecoder = pickle.load(open(cacheprefix+'continuous/acrosscontextcomparison-%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))
            acrossdecoders.append(acrossdecoder)

            print(dn,cx,len(acrossdecoder))
        acrossdecoders_all.append(acrossdecoders)

    acrossdecoders_all = np.array(acrossdecoders_all)


    fig,ax = plt.subplots(4,2,figsize=(2*1.41*8,4*8))


    # cross decoders
    for mx in range(2):
        M = []
        for n,dn in enumerate(datanames):
            axs = ax[0,mx]
            axs.plot(times,acrossdecoders_all[n][mx*2][1][:,0], color='black')
            axs.plot(times,acrossdecoders_all[n][mx*2+1][1][:,0], color='red')
            
            axs = ax[1,mx]
            m = acrossdecoders_all[n][mx*2][1][:,0]-acrossdecoders_all[n][mx*2+1][1][:,0]
            axs.plot(times,m, label=dn)
        
            M.append(m)

            axs = ax[2,mx]
            axs.boxplot(positions=[n], x=m[150:450], notch=True, whis=[5,95], showfliers=False)
            
            axs = ax[3,mx]
            axs.boxplot(positions=[n], x=m[150:200], notch=True, whis=[5,95], showfliers=False)


        M = np.array(M)

        axs = ax[2,mx]
        axs.boxplot(positions=[n+2], x=np.vstack(M[:,150:450]), notch=True, whis=[5,95], showfliers=False)
        
        axs = ax[3,mx]
        axs.boxplot(positions=[n+2], x=np.vstack(M[:,150:200]), notch=True, whis=[5,95], showfliers=False)


        ax[0,mx].set_ylim(0.45,1.05)
        figs.plottoaxis_stimulusoverlay(ax[0,mx],T)
        ax[0,mx].legend(['attend','ignore'],frameon=False)



        figs.plottoaxis_stimulusoverlay(ax[1,mx],T)


        figs.plottoaxis_chancelevel(ax[2,mx])
        ax[2,mx].set_ylim(-0.3,0.3)
        ax[2,mx].set_xticklabels([])


        figs.plottoaxis_chancelevel(ax[3,mx])
        ax[3,mx].set_ylim(-0.3,0.3)
        ax[3,mx].set_xticklabels(labels,rotation=60)


        # stats
        print(['visual','audio'][mx], 'stats')
        S = np.vstack(M[:,150:200])
        t,p = sp.stats.ttest_1samp(S,0)
        print('timepoints: t',t,'p',p,'m+/-sem',S.mean(),'+/-', S.std()/np.sqrt(len(S)), 'std', S.std())
        R = M[:,150:200].mean(axis=1)
        t,p = sp.stats.ttest_1samp(R,0)
        print('mice:       t',t,'p',p,'m+/-sem',R.mean(),'+/-',R.std()/np.sqrt(len(R)),'std', R.std())






    ax[0,0].set_title('visual')
    ax[0,1].set_title('audio')
    ax[0,0].set_ylabel('accuracy')

    # ax[1,0].set_ylabel('same$-$cross\naccuracy difference\nover time')

    # ax[2,0].set_ylabel('same$-$cross\naccuracy difference\ndistribution, on stim')

    # ax[3,0].set_ylabel('same$-$cross\naccuracy difference\ndistribution, on stim 0-0.5 sec')





    # ax[1,2].set_ylabel('visual$-$audio\naccuracy difference\nover time')

    # ax[2,2].set_ylabel('visual$-$audio\naccuracy difference\ndistribution, on stim')

    # ax[3,2].set_ylabel('visual$-$audio\naccuracy difference\ndistribution, on stim 0-0.5 sec')










    fig.suptitle('visual decoding comparison across contexts')
    fig.tight_layout()


    save = 0 or globaldoplot
    if save:
        fig.savefig(resultpath+'acrosscontextcomparison,visual,audio-tribe'+ext)

    return
















def aggregatemice_crossgradientdecoder(datanames):
    # assess context in the initial and transition single modality cue blocks at the single trial level
    # create the cross gradient file: [ block (in original block order) ],
    #                                 ( {pre,on},trials,{dtm,dtmse,savgol} )
    
    datanames = np.array(datanames)
    n_mice = len(datanames)

    trialcompress = 30
    
    aucdeltat = 750*pq.ms     # for the order length

    crossgradients = 0.5*np.ones( (2,2,n_mice,100) )      # starting values for empty trials, this is raw trial numbers
    crossgradients2 = 0.5*np.ones( (2,2,n_mice,trialcompress) )      # starting values for empty trials, this is timecompressed
    blockorder = np.zeros((n_mice,2),dtype='int16')      # this will mean 0: attend visual first, 1: attend audio first
    decauc_pre = np.zeros((n_mice))
    

    
    
    ls = np.zeros((2,2,n_mice,2))
    
    for n,dn in enumerate(datanames):
    
        export_crossgradient = pickle.load(open(cacheprefix+'continuous/offstimgradient,cue,decodeprobs,savgol_%s-%s.pck'%(dn,continuous_method),'rb'))
#        export_crossgradients[n][bl][ptrx,:,0]
        
        blv,bla = preprocess.getorderattended(dn)
        
        for bl in [0,1]:          # block visual, audio
            n_cuetrials = export_crossgradient[bl].shape[1]
            blockorder[n,bl] = int(([blv[0],bla[0]][bl]-1)/2)
            for ptrx in [0,1]:    # pre, on
                crossgradients[blockorder[n,bl],ptrx,n,:n_cuetrials] = export_crossgradient[bl][ptrx,:,2]
                crossgradients2[blockorder[n,bl],ptrx,n,:] = np.interp( np.linspace(0,n_cuetrials-1,trialcompress), \
                              np.arange(n_cuetrials), export_crossgradient[bl][ptrx,:,2] )
                l = sp.stats.linregress(np.arange(trialcompress),crossgradients2[blockorder[n,bl],ptrx,n,:])
                ls[blockorder[n,bl],ptrx,n,:] = l[:2]


        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # lower confidence bound
        decauc_pre[n] = neph.get_auctimenormalized(decsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]        
        
                
    
    orders = np.argsort(ls[1,1,:,1])[::-1]     # cue 3, onstimulus,   all mice,    slope/intercept
    orders = np.argsort(ls[1,1,:,1]+ls[1,1,:,0]*15+crossgradients2[blockorder[n,1],1,n,-1]+\
                        crossgradients2[blockorder[n,0],1,n,-1]\
                        )[::-1]     # cue 3, prestimulus,   all mice,    slope/intercept
    orders = np.argsort(decauc_pre)[::-1]
    
    
    if not publish:
        blocklabels = ['cue 1','cue 3']
        timeperiodlabels = ['pre','on']
    else:
        blocklabels = ['initial','switched']
        timeperiodlabels = [' PRE',' ON']
        
    if not publish: fig, ax = plt.subplots(2,4,figsize=(32,24))
    else:
        fig = plt.figure(figsize=(14,5))
        gs = fig.add_gridspec(1, 13)
        ax = []
        ax.append( fig.add_subplot(gs[0, 0:3]) )
        ax.append( fig.add_subplot(gs[0, 3:6]) )
        ax.append( fig.add_subplot(gs[0, 6:9]) )
        ax.append( fig.add_subplot(gs[0, 9:12]) )
        ax.append( fig.add_subplot(gs[0, 12:13]) )
        
    for blo in [0,1]:          # block cue 1, cue 3
        for ptrx in [0,1]:    # pre, on
            
            G = crossgradients[blo,ptrx,orders,:]
            
            if not publish:
                axs = ax[blo,ptrx]
                axs.imshow(G,vmin=0.0,vmax=1.0,cmap='RdBu',aspect='auto')
                figs.plottoaxis_notickmarks(axs)
    #            axs.set_yticks(range(n_mice)); axs.set_yticklabels(datanames[orders]);
    #            axs.set_xticks([])
                
                if ptrx==0: axs.set_ylabel(timeperiodlabels[ptrx]);
                axs.set_title(blocklabels[blo])

                axs = ax[blo,ptrx+2]
            else:
                axs = ax[blo*2+ptrx]

            G = crossgradients2[blo,ptrx,orders,:]
            
            im = axs.imshow(G,vmin=0.0,vmax=1.0,cmap='RdBu',aspect='auto')
            figs.plottoaxis_notickmarks(axs)
#            axs.set_ylabel(timeperiodlabels[ptrx]);
#            axs.set_yticks(range(n_mice)); axs.set_yticklabels(datanames[orders]);
#            axs.set_xticks([])

            axs.set_title(blocklabels[blo]+timeperiodlabels[ptrx])
            
#    fig.subplots_adjust(right=0.95)
#    axs = fig.add_axes([0.97, 0.15, 0.02, 0.7])
#    fig.colorbar(im, cax=ax[4], label='context probability'  )
#    ax[4].set_yticks([0,0.5,1.0]); ax[4].set_yticklabels(['0','0.5','1'])


    save = 0
    if save:
        if not publish:
            fig.savefig(resultpath+'micetribe,cuetrials,crossgradient-precontextdecsort'+ext)
        else:
            # for publication:
            fig.savefig(resultpath+'6-tribecolors-cuetrials,crosstestgradient'+ext)










def aggregatemice_cueperiod_neuraltobehaviour(datanames):
    # relate context representation and behaviour in the cueing period

    # subselect only clever animals
#    datanames = ['DT008','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']



    datanames = np.array(datanames)
    n_mice = len(datanames)

    trialcompress = 30
    n_max_trials = 30
    
    aucdeltat = 750*pq.ms     # for the order length

    crossgradients = 0.5*np.ones( (2,2,n_mice,n_max_trials) )      # starting values for empty trials, this is raw trial numbers
    crossgradients2 = 0.5*np.ones( (2,2,n_mice,trialcompress) )      # starting values for empty trials, this is timecompressed
    blockorder = np.zeros((n_mice,2),dtype='int16')      # this will mean 0: attend visual first, 1: attend audio first
    cx_decauc_pre = np.zeros((n_mice))       # this is context PRE: timeaveraged decoder accuracy
    cuetrainings = []
    cue_ests = []

    
    
    ls = np.zeros((2,2,n_mice,2))
    
    for n,dn in enumerate(datanames):
    
        export_crossgradient = pickle.load(open(cacheprefix+'continuous/offstimgradient,cue,decodeprobs,savgol_%s-%s.pck'%(dn,continuous_method),'rb'))
#        export_crossgradients[n][bl][ptrx,:,0]
        
        blv,bla = preprocess.getorderattended(dn)
        
        for bl in [0,1]:          # block visual, audio
            n_cuetrials = export_crossgradient[bl].shape[1]
            blockorder[n,bl] = int(([blv[0],bla[0]][bl]-1)/2)
            for ptrx in [0,1]:    # pre, on
                crossgradients[blockorder[n,bl],ptrx,n,:n_cuetrials] = export_crossgradient[bl][ptrx,:n_max_trials,2]
                crossgradients2[blockorder[n,bl],ptrx,n,:] = np.interp( np.linspace(0,n_cuetrials-1,trialcompress), \
                              np.arange(n_cuetrials), export_crossgradient[bl][ptrx,:,2] )
                l = sp.stats.linregress(np.arange(trialcompress),crossgradients2[blockorder[n,bl],ptrx,n,:])
                ls[blockorder[n,bl],ptrx,n,:] = l[:2]


        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # lower confidence bound
        cx_decauc_pre[n] = neph.get_auctimenormalized(decsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]        


        _,_,_,_,cuetraining,cue_est = preprocess.loadtrainingbehaviouraldata(dn)
#        cuetraining is (init/trans)(sessions,h+c/m+f,trials)=(2)(8,2,30)
        
        
        # create the regression predictor lists:
        # we need at each mouse: a list of performances from sessions separately for each trial
        # grouped by initial and transitional periods
        # (we only reshape over first and last few trials later)
        # [mice][init/trans,sessions,h+c/m+f,trials]           =(14)(8,2,30)

        cuetrainings.append(cuetraining)
        cue_ests.append(cue_est)

    cue_ests = np.swapaxes(np.array(cue_ests),0,1)

    print(crossgradients.shape)
    print(cx_decauc_pre.shape)
    print(cue_ests.shape)
    


    # find out the optimal correlates of transition
    # averages of first N and last N trials
    window = 5         # averaging at the begining and end of the cue periods (30 trials restricted)
    
    neural =    np.zeros((2,2,n_mice,2,2))     # block,preon,mouse,firstN-lastN,stats
    behaviour = np.zeros((2,n_mice,2,2))       # block,mouse,firstN-lastN,stats
    for n in range(n_mice):
        neural[:,:,n,0,0] = crossgradients[blockorder[n,:],:,n,:window].mean(axis=2)
        neural[:,:,n,1,0] = crossgradients[blockorder[n,:],:,n,-window:].mean(axis=2)
        neural[:,:,n,0,1] = crossgradients[blockorder[n,:],:,n,:window].std(axis=2)/np.sqrt(window)
        neural[:,:,n,1,1] = crossgradients[blockorder[n,:],:,n,-window:].std(axis=2)/np.sqrt(window)
    behaviour[:,:,0,0] = cue_ests[:,:,0,:window].mean(axis=2)
    behaviour[:,:,1,0] = cue_ests[:,:,0,-window:].mean(axis=2)
    behaviour[:,:,0,1] = cue_ests[:,:,0,:window].std(axis=2)/np.sqrt(window)
    behaviour[:,:,1,1] = cue_ests[:,:,0,-window:].std(axis=2)/np.sqrt(window)
    

    
    # collect regression behavioural data as long lists for each animal (i.e. at each context PRE value)
    initials = []
    transitions = []
    for n,dn in enumerate(datanames):
        initials.append( [ np.reshape(cuetrainings[n][0][:,0,:window],-1),\
                           np.reshape(cuetrainings[n][0][:,0,-window:],-1) ]   )
        transitions.append( [ np.reshape(cuetrainings[n][1][:,0,:window],-1),\
                              np.reshape(cuetrainings[n][1][:,0,-window:],-1) ] )



    print(neural.shape)
    print(behaviour.shape)





    if 1:  # plot various measures of neural-behavaioural correlates of transition
        fig, ax = plt.subplots(4,4,figsize=(32,32))
        colors = [['red','blue'],['darkred','navy'],['purple'],['forestgreen']]
        blocklabel = ['initial','transition']
        trialgrouplabel = ['first %d'%window,'last %d'%window]
        measurelabel = ['neural','behaviour']
        for dff in [0,1]:
            for ptrx in [0,1]:
                for bl in [0,1]:          # block initial and transition
                    period = [initials,transitions][bl]         # shortcut call
                    for nbix in [0,1]:    # neural and behavioural
                        axs = ax[2*dff+nbix, ptrx*2+bl]
                        for i in [0,1]:          # first and last, on the same figure
                            if nbix==1 and ptrx==1: x = neural[bl,1,:,i,0]
                            if dff==1 and i==1: continue
                            else: x = cx_decauc_pre
                            
                            # plot the overall binomial p values as points:

                            X = np.concatenate(  [ [  x[n] for d in range( len(period[n][i]) ) ] \
                                               for n in range(n_mice)  ]  )

                            if dff==0:
                                y = [neural[bl,ptrx,:,i,0], behaviour[bl,:,i,0]][nbix]
                                Y = np.concatenate(  [ [  period[n][i][d] for d in range( len(period[n][i]) ) ] \
                                                   for n in range(n_mice)  ]  )
                            elif dff==1:

                                y = [  neural[bl,ptrx,:,1,0] - neural[bl,ptrx,:,0,0],\
                                       behaviour[bl,:,1,0]   - behaviour[bl,:,0,0]     ][nbix]
                                Y = np.concatenate(  [ [  period[n][1][d]-period[n][0][d] for d in range( len(period[n][0]) ) ] \
                                                   for n in range(n_mice)  ]  )


#                    dY = (  [ (period[n][1] - period[n][0]).mean()  \
#                                       for n in range(n_mice)  ]  )

                                
                            axs.scatter(x,y,s=150,color=colors[2*dff+nbix][i],label=trialgrouplabel[i])
        
                            
    
                            l = sp.stats.linregress(X,Y)
                            if l[3]<0.05:
                                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                                axs.plot(line,l[0]*line+l[1],color=colors[2*dff+nbix][i],linewidth=2)
                            else:
                                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                                axs.plot(line,l[0]*line+l[1],'--',color=colors[2*dff+nbix][i],linewidth=2)
                                
                            xoffs = 0.0  #0.35
                            yoffs = [0.5,0.2][nbix]+0.05*i
                            if l[3]<1e-5:
                                axs.text(line[0]+xoffs,yoffs,'p<%s, slope=%5.3f, R$^2$=%6.4f'%(('10e$^{-5}$',l[2],l[2]**2)),color=colors[2*dff+nbix][i])
                            elif l[3]<0.001:
                                axs.text(line[0]+xoffs,yoffs,'p<%5.3f, slope=%5.3f, R$^2$=%6.4f'%((0.001,l[2],l[2]**2)),color=colors[2*dff+nbix][i])
                            else:
                                axs.text(line[0]+xoffs,yoffs,'p=%5.3f, slope=%5.3f, R$^2$=%6.4f'%((l[3],l[2],l[2]**2)),color=colors[2*dff+nbix][i])
        
                        axs.set_xlim(0.45,1.05)
                        axs.set_ylim([0.,0][nbix],1.05)
                        figs.plottoaxis_crosshair(axs,0.5,None)
    
                        axs.legend(frameon=False)
                        axs.set_title('%s %s %s'%(blocklabel[bl], measurelabel[nbix], [['PRE','ON'][ptrx],''][nbix]))
                        axs.set_xlabel('context PRE')
                        axs.set_ylabel(measurelabel[nbix])
                        if nbix==1 and ptrx==1: 
                            axs.set_xlabel('neural ON')
                            axs.set_title('%s'%blocklabel[bl])
#                            axs.set_xlim(0.0,1.05)
                        if dff==1:
                            axs.set_ylim(-0.5,0.5)

        save = 0
        if save:
            fig.savefig(resultpath+'micetribe,cuetrials,crossgradient-behaviour,regression'+ext)




    if 1: # plot group differences
        fig, ax = plt.subplots(2,4,figsize=(32,16))
        colors = [['magenta','lightgreen'],['purple','darkgreen']]
        blocklabel = ['initial','transition']
        trialgrouplabel = ['first %d'%window,'last %d'%window]
        measurelabel = ['neural','behaviour']
        for ptrx in [0,1]:
            for bl in [0,1]:          # block initial and transition
                period = [initials,transitions][bl]         # shortcut call
                for nbix in [0,1]:    # neural and behavioural
                    if nbix==1 and ptrx==1: continue
                    axs = ax[nbix, ptrx*2+bl]
                    dy = [neural[bl,ptrx,:,:,0], behaviour[bl,:,:,0]][nbix] @ np.array([-1,1])
                    dys = [neural[bl,ptrx,:,:,1], behaviour[bl,:,:,1]][nbix] @ np.array([1,1])
                    
                    print(len(period[n][1]))
                    dY = (  [ (period[n][1] - period[n][0]).mean()  \
                                       for n in range(n_mice)  ]  )
                    dYs = (  [ (period[n][1] - period[n][0]).std()/np.sqrt( len(period[n][1]) )  \
                                       for n in range(n_mice)  ]  )
                    
                    
                    if nbix==0: dY = dy; dYs = dys
                    
                    
                    axs.bar(x=np.arange(n_mice),height=dY,color=colors[nbix][bl],label='%s - %s'%(trialgrouplabel[1],trialgrouplabel[0]))
                    for n in range(n_mice): axs.errorbar(x=n,y=dY[n],yerr=dYs[n],color='darkgrey')
                    
                    axs.legend(frameon=False)
                    axs.set_title('%s %s %s'%(blocklabel[bl], measurelabel[nbix], ['PRE','ON'][ptrx]))
                    axs.set_ylabel(measurelabel[nbix])


#                    axs.set_xlim(0.45,1.05)
#                    axs.set_ylim([0.,0][nbix],1.05)
#                    figs.plottoaxis_crosshair(axs,0.5,None)
    
                    axs.set_xticks(np.arange(n_mice))
                    axs.set_xticklabels(datanames,rotation=45,fontsize='x-small')
                    axs.legend(frameon=False)
                    axs.set_title('%s %s %s'%(blocklabel[bl], measurelabel[nbix], [['PRE','ON'][ptrx],''][nbix]))
#                    axs.set_xlabel('mice')
                    axs.set_ylabel(measurelabel[nbix])

        save = 0
        if save:
            fig.savefig(resultpath+'micetribe,cuetrials,crossgradient-behaviour,differences'+ext)




# test mice for orientation priors
def aggregatemice_crossorentation(datanames, blocks, T):
    
    n_mice = len(datanames)
    
    labels = ['0°-90°','30°-120°','60°-150°','45°-135°']
    

    
    projections_all = []
    variances_all = []
    pvals_all = []
    CI = []
    for n,dn in enumerate(datanames):
        print('Collecting projections ',dn,'...')
        block = blocks[n]
        proj_oris = []
        var_oris = []
        chi_oris = []
        for trainidx in [2,3,4,5]:
            testidx = 6
            v_aux,p_aux = decoder_crosscontexttest(dn,block,retvar=True,trainidx=trainidx,testidx=testidx)
            # dimensions: projections[0][preon][trlx][class]
            proj_oris.append(   np.r_[p_aux[0][0][1][0],p_aux[0][0][1][1]]   )
#            var_oris.append(   np.r_[v_aux[1][2],np.std(proj_oris[-1])]    )
            var_oris.append(    v_aux[1][2]    )
            chi_oris.append(   np.r_[ v_aux[1][2] * np.sqrt((len(proj_oris[-1])-1)/sp.stats.chi2.isf(0.025,df=(len(proj_oris[-1])-1) )),  \
                                      v_aux[1][2] * np.sqrt((len(proj_oris[-1])-1)/sp.stats.chi2.isf(0.975,df=(len(proj_oris[-1])-1) )) ]  )
        projections_all.append(proj_oris)
        variances_all.append(var_oris)
        CI.append(chi_oris)
        aux = []
        for trx in [0,1,2]:
            s,p = sp.stats.bartlett(proj_oris[trx],proj_oris[3])
            aux.append(p)
        pvals_all.append(aux)
    P = np.array(projections_all)  
    V = np.array(variances_all)
    CI = np.array(CI)
    p = np.array(pvals_all)
    
    print(P.shape,V.shape,CI.shape,p.shape)
    print(V)
    print(CI)
    

    orders = [0,1,2,3]
    orders = [1,3,2,0]
    colors = ['darkred','darkcyan','purple','darkgreen']
    fig,ax = plt.subplots(1,1,figsize=(32,15))
    axs = ax
    for trxidx,trx in enumerate(orders):
        axs.bar(x=np.arange(n_mice)*1.25+(trxidx)/6, height=V[:,trx], width=1/6,\
                color=colors[trx], alpha=0.5, tick_label=datanames )
    axs.legend([ labels[i] for i in orders] )
    for trxidx,trx in enumerate(orders):
        axs.errorbar(x=np.arange(n_mice)*1.25+(trxidx)/6, y=V[:,trx],\
                    yerr=[ V[:,trx]-CI[:,trx,0], -V[:,trx]+CI[:,trx,1] ],\
                    capsize=2,capthick=2,ls='',color=colors[trx],lw=2)

    axs.set_xlabel('mice')
    axs.set_yticks([0,0.1,0.2,0.3])
    axs.set_ylabel('std with 95% c.i. [SU]')
    axs.set_title('st.dev.-s of multimodal block activity\nprojected onto 5th block decoder vectors (colors) ')

    save = 1
    if save:
        fig.savefig(resultpath+'micetribe,crossorientation-projections,variance'+ext)



        





def aggregatemice_layercelltypecontributions(datanames,examine='attendignore'):
    # compare layer and celltype contributions to various task variable representations with decoders

    #'allexpcond' or 'attendignore'
    n_mice = len(datanames)
    
    stimulusspecificity = 'all'               # marginal on all other conditions
    if examine=='allexpcond':
        taskaspects = ['visual','audio','context','choice']
        taskfilenames = taskaspects
    elif examine=='attendignore':
        taskaspects = ['attend visual','ignore visual','attend audio','ignore audio']
        taskfilenames = ['avis','ivis','aaud','iaud']
    
    trajectorylength = 596
    block = preprocess.loaddatamouse(datanames[0],T,continuous_method,recalculate=False)
    times_diff = block.segments[0].analogsignals[0].times
    times = times_diff[:trajectorylength]

    layerlist = [0,1,2,3,4]       # choose which layer to focus on: 0 all, 1 superf, 2 input, 3 deep5, 4 deep6
#    layerlist = [0]
    celltypelist = [0,1]        # choose which spike width to use: 0 narrow, 1 broad
    leaveorusegroup = ['leave','only','random']

    # last dimension is mean,s.e.m.
    dec_all = np.zeros( ( len(taskaspects),1,1,n_mice,trajectorylength,3 ) )
    dec_leave = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,n_mice,trajectorylength,3 ) )
    dec_leave_random = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,n_mice,trajectorylength,3 ) )
    dec_only = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,n_mice,trajectorylength,3 ) )

    diff_all = np.zeros( ( len(taskaspects),1,1,n_mice,trajectorylength+5,3 ) )
    diff_leave = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,n_mice,trajectorylength+5,3 ) )
    diff_leave_random = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,n_mice,trajectorylength+5,3 ) )
    diff_only = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,n_mice,trajectorylength+5,3 ) )

    
    a_dec_leave = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,trajectorylength,3 ) )
    a_dec_only = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,trajectorylength,3 ) )
    a_diff_leave = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,trajectorylength+5,3 ) )
    a_diff_only = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,trajectorylength+5,3 ) )

    a_dec_leave_random = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,trajectorylength,3 ) )
    a_diff_leave_random = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1,trajectorylength+5,3 ) )

    
#    dec_leave.fill(np.nan); dec_only.fill(np.nan); diff_leave.fill(np.nan); diff_only.fill(np.nan);
    
    dec_N = np.zeros( ( len(taskaspects),len(celltypelist),len(layerlist)-1) )


    colors_layers = np.array([['lightcoral','deeppink','red','darkred'],\
                              ['deepskyblue','darkviolet','blue','navy']])   # superf, input, deep5, deep6
#    colors_layers_ignore = np.array([['maroon','rebeccapurple','red','darkred'],\
#                              ['deepskyblue','darkviolet','blue','navy']])   # superf, input, deep5, deep6


    for cx,comparison in enumerate(taskaspects):
        for cidx in celltypelist:
            for layered in layerlist:
                if layered==0 and cidx>0: continue    # don't run over celltypes, if it's all layers and celltypes, it's by definition all celltypes in the first run
                for n,dn in enumerate(datanames):
                    # check if there is a such neuron group at all present in the given mouse
                    
                    # sort out which layers and celltypes to leave out or only use
                    if layered==0:
                        layeredfilename = 'all'
                        acrossdifference = pickle.load(open(cacheprefix+'continuous/responsedifference-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
                        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
                        dec_all[cx,0,0,n,:,0] = acrossdecoder[1][:trajectorylength,0].squeeze()           # 0->[0,2]
                        dm,de = nedi.get_maxnormdifference(acrossdifference)
                        diff_all[cx,0,0,n,:,0] = dm.squeeze()
                    else:
                        lidx = layered-1
                        if not os.path.isfile(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,   ',%s,%d-%s,%s'%(leaveorusegroup[0],lidx+1,layernames[lidx],spikewidthnames[cidx]),\
                                     dn,continuous_method,comparison,stimulusspecificity)):
                            continue
                        dec_N[cx,cidx,lidx] += 1
#                        print(n,dn,['narrow','broad'][cidx],layernames[lidx],dec_N[cx,cidx,lidx] )
                        for gx in range(len(leaveorusegroup)):
                            if gx==1: continue   # forget "only"
                            layeredfilename = ',%s,%d-%s,%s'%(leaveorusegroup[gx],layered,layernames[lidx],spikewidthnames[cidx])
                            acrossdifference = pickle.load(open(cacheprefix+'continuous/responsedifference-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
                            acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
                            if gx==0:     # leave
                                dec_leave[cx,cidx,lidx,n,:,:] = acrossdecoder[1][:,:].squeeze()
                                dm,de = nedi.get_maxnormdifference(acrossdifference)
                                diff_leave[cx,cidx,lidx,n,:,0] = dm.squeeze()
                                diff_leave[cx,cidx,lidx,n,:,2] = de.squeeze()
                                a_dec_leave[cx,cidx,lidx,:,:] += dec_leave[cx,cidx,lidx,n,:,:]
                                a_diff_leave[cx,cidx,lidx,:,:] += diff_leave[cx,cidx,lidx,n,:,:]
                            elif gx==1:    # only
                                dec_only[cx,cidx,lidx,n,:,:] = acrossdecoder[1][:,:].squeeze()
                                dm,de = nedi.get_maxnormdifference(acrossdifference)
                                diff_only[cx,cidx,lidx,n,:,0] = dm.squeeze()
                                diff_only[cx,cidx,lidx,n,:,2] = de.squeeze()
                                a_dec_only[cx,cidx,lidx,:,:] += dec_only[cx,cidx,lidx,n,:,:]
                                a_diff_only[cx,cidx,lidx,:,:] += diff_only[cx,cidx,lidx,n,:,:]
                            elif gx==2:    # randomized leave-group-out, to compare randomized to
                                dec_leave_random[cx,cidx,lidx,n,:,:] = acrossdecoder[1][:,:].squeeze()
                                dm,de = nedi.get_maxnormdifference(acrossdifference)
                                diff_leave_random[cx,cidx,lidx,n,:,0] = dm.squeeze()
                                diff_leave_random[cx,cidx,lidx,n,:,2] = de.squeeze()
                                a_dec_leave_random[cx,cidx,lidx,:,:] += dec_leave_random[cx,cidx,lidx,n,:,:]
                                a_diff_leave_random[cx,cidx,lidx,:,:] += diff_leave_random[cx,cidx,lidx,n,:,:]
    
    a_dec_leave /= dec_N[:,:,:,None,None]
    a_diff_leave /= dec_N[:,:,:,None,None]
    a_dec_only /= dec_N[:,:,:,None,None]
    a_diff_only /= dec_N[:,:,:,None,None]
    a_dec_leave_random /= dec_N[:,:,:,None,None]
    a_diff_leave_random /= dec_N[:,:,:,None,None]

    print(dec_all.shape,dec_leave.shape,dec_only.shape)
#    print(dec_leave.max())
#    print(dec_only.max())


    if 0:
        for cx,comparison in enumerate(taskaspects):
            fig,ax = plt.subplots(4,4,figsize=(36,30))
            for cidx in celltypelist:
                for lidx in [0,1,2,3]:
                    
                    
    #                D_leave = (dec_leave[cx,cidx,lidx,:,:,0] - dec_all[cx,0,0,:,:,0]).sum(axis=0,keepdims=True)/dec_N[cx,cidx,lidx]
    #                D_only = (dec_only[cx,cidx,lidx,:,:,0] - dec_all[cx,0,0,:,:,0]).sum(axis=0,keepdims=True)/dec_N[cx,cidx,lidx]
    #                
    #
    #                R_leave = (diff_leave[cx,cidx,lidx,:,:,0] - diff_all[cx,0,0,:,:,0]).sum(axis=0,keepdims=True)/dec_N[cx,cidx,lidx]
    #                R_only = (diff_only[cx,cidx,lidx,:,:,0] - diff_all[cx,0,0,:,:,0]).sum(axis=0,keepdims=True)/dec_N[cx,cidx,lidx]
    #
    #                S_leave = R_leave + diff_all[cx,0,0,:,:,0].mean(axis=0)
    #                S_only = R_only + diff_all[cx,0,0,:,:,0].mean(axis=0)
    
                    R_all = diff_all[cx,0,0,:,:,:].mean(axis=0)
                    R_all_se = diff_all[cx,0,0,:,:,:].std(axis=0)/np.sqrt(n_mice)*2
                    R_leave = a_diff_leave[cx,cidx,lidx,:,:]
                    R_only = a_diff_only[cx,cidx,lidx,:,:]
                    R_leave_random = a_diff_leave_random[cx,cidx,lidx,:,:]
                    
                    D_all = dec_all[cx,0,0,:,:,:].mean(axis=0)
                    D_leave = a_dec_leave[cx,cidx,lidx,:,:]
                    D_only = a_dec_only[cx,cidx,lidx,:,:]
                    D_leave_random = a_dec_leave_random[cx,cidx,lidx,:,:]
                    
    
    
                    axs = ax[lidx,1-cidx]      # just to display broad spiking first and narrow spiking second column
                    
                    
    #                axs.plot(times,dec_all[cx,0,0,:,:,0].T,color='darkorange',label='all')
    #                axs.plot(times,dec_leave[cx,cidx,lidx,:,:,0].T,color='firebrick',label='leave')
    #                axs.plot(times,dec_only[cx,cidx,lidx,:,:,0].T,color='g',label='only')
                    
    #                axs.plot(times_diff,diff_all[cx,0,0,:,:,0].mean(axis=0).T,color='darkorange',label='all')
                    
    #                axs.plot(times_diff,R_all[:,0],color='k',linewidth=1,alpha=0.5,label='all cells')
                    axs.plot(times_diff,R_leave[:,0],'--',color=colors_layers[cidx,lidx],linewidth=3,alpha=0.67,label='%s %s left out'%(layernames[lidx],spikewidthnames[cidx]))
                    axs.plot(times_diff,R_leave_random[:,0],color=colors_layers[cidx,lidx],linewidth=3,alpha=1.0,label='5x randomized same # cells left out')
    
    #                axs.fill_between(times_diff, (R_leave[:,0]-R_leave[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
    #                                    (R_leave[:,0]+R_leave[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
    #                                    alpha=0.04,color=colors_layers[cidx,lidx])
    #                axs.fill_between(times_diff, (R_leave_random[:,0]-R_leave_random[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
    #                                    (R_leave_random[:,0]+R_leave_random[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
    #                                    alpha=0.04,color=colors_layers[cidx,lidx])
    
                    axs.fill_between(times_diff, (R_leave[:,0]-R_all_se[:,0]),\
                                        (R_leave[:,0]+R_all_se[:,0]),\
                                        alpha=0.04,color=colors_layers[cidx,lidx])
                    axs.fill_between(times_diff, (R_leave_random[:,0]-R_all_se[:,0]),\
                                        (R_leave_random[:,0]+R_all_se[:,0]),\
                                        alpha=0.04,color=colors_layers[cidx,lidx])
    
    #                    axs.plot(times_diff,R_only.T)#,color=colors_layers[lidx],linewidth=2)
                    
                    axs.set_ylim(-0.5,4)
                    axs.set_xlim(-1500,4500)
                    figs.plottoaxis_chancelevel(axs,0)
                    figs.plottoaxis_stimulusoverlay(axs,T)
    
                    if lidx==0: axs.set_title('difference\n'+spikewidthnames[cidx])#+' N %d '%dec_N[cx,cidx,lidx])
                    if cidx==0: axs.set_ylabel(layernames[lidx])
    #                axs.legend(['all','superf','input','deep5','deep6'])
                    axs.legend()
    
    
    
                    axs = ax[lidx,1-cidx+2]
    
    #                axs.plot(times,D_all[:,0],color='k',linewidth=1,alpha=0.5,label='all cells')
                    axs.plot(times,D_leave[:,0],'--',color=colors_layers[cidx,lidx],linewidth=1,alpha=0.67,label='%s %s left out'%(layernames[lidx],spikewidthnames[cidx]))
                    axs.plot(times,D_leave_random[:,0],color=colors_layers[cidx,lidx],linewidth=1,alpha=1.,label='5x randomized same # cells left out')
    #                    axs.plot(times,D_only.T)#,color='darkseagreen')
                    axs.fill_between(times, (D_leave[:,0]-D_leave[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
                                        (D_leave[:,0]+D_leave[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
                                        alpha=0.1,color=colors_layers[cidx,lidx])
                    axs.fill_between(times, (D_leave_random[:,0]-D_leave_random[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
                                        (D_leave_random[:,0]+D_leave_random[:,2]/np.sqrt(dec_N[cx,cidx,lidx])),\
                                        alpha=0.3,color=colors_layers[cidx,lidx])
    
                    axs.set_ylim(0.395,1.05)
                    axs.set_xlim(-1500,4500)
                    figs.plottoaxis_chancelevel(axs,0.5)
                    figs.plottoaxis_stimulusoverlay(axs,T)
    #                axs.legend(['all','superf','input','deep5','deep6'])
                    axs.legend()
                    
                    if lidx==0: axs.set_title('decoder\n'+spikewidthnames[cidx])#+' N %d '%dec_N[cx,cidx,lidx])
    
                    
            fig.suptitle('%s aspect representation by celltypes in layers, %d mice'%(comparison,n_mice))
            
            save = 1
            if save:
                fig.savefig(resultpath+'micetribe,celltypeslayers-%s_%d-%s_%s-%dms%dms'%(examine,cx+1,comparison,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
    
    




    # show contributions by the difference between random leave and specific leave 
    if 1:
        layerlist = [0,1,2,3]
        R = np.zeros(   ( len(taskaspects),len(celltypelist),len(layerlist), trajectorylength+5)   )
        D = np.zeros(   ( len(taskaspects),len(celltypelist),len(layerlist), trajectorylength)   )
        savvywidth = 201 # 57,137,201
        S = np.zeros(   ( len(taskaspects),len(celltypelist),len(layerlist), trajectorylength)   )
        
        for cx,comparison in enumerate(taskfilenames):
            for cidx in celltypelist:
                for lidx in layerlist:
                    R[cx,cidx,lidx,:] = ( a_diff_leave_random[cx,cidx,lidx,:,0] - a_diff_leave[cx,cidx,lidx,:,0] )
                    D[cx,cidx,lidx,:] = ( a_dec_leave_random[cx,cidx,lidx,:,0] - a_dec_leave[cx,cidx,lidx,:,0] )
                    S[cx,cidx,lidx,:] = sp.signal.savgol_filter(D[cx,cidx,lidx,:], savvywidth, 5)
                    
        fig,ax = plt.subplots(4,2,figsize=(18,36))
        for cidx in celltypelist:
            for lidx in layerlist:
                axs = ax[lidx,1-cidx]           
                #colors_layers[cidx,lidx]
                axs.plot(times,S[0,cidx,lidx,:],color='darkred',linewidth=2,alpha=0.8)
                axs.plot(times,S[1,cidx,lidx,:],color='navy',linewidth=2,alpha=0.8)
                axs.plot(times,D[0,cidx,lidx,:],color='darkred',linewidth=1,alpha=0.2)
                axs.plot(times,D[1,cidx,lidx,:],color='navy',linewidth=1,alpha=0.2)
                axs.legend(taskaspects[:2])
                axs.set_ylim(-0.06,0.06)
                axs.set_xlim(-1500,4500)
                figs.plottoaxis_chancelevel(axs,0.0)
                figs.plottoaxis_stimulusoverlay(axs,T)
                axs.set_title(spikewidthnames[cidx])#+' N %d '%dec_N[cx,cidx,lidx])
                if cidx==1: axs.set_ylabel(layernames[lidx])

                continue
                axs = ax[lidx,1-cidx+2]
                axs.plot(times,S[0,cidx,lidx,:]-S[1,cidx,lidx,:],color=colors_layers[cidx,lidx],linewidth=2,alpha=1.0)
                axs.set_ylim(-0.06,0.06)
                axs.set_xlim(-1500,4500)
                figs.plottoaxis_chancelevel(axs,0.0)
                figs.plottoaxis_stimulusoverlay(axs,T)
                axs.set_title(spikewidthnames[cidx])#+' N %d '%dec_N[cx,cidx,lidx])
#                if cidx==1: axs.set_ylabel(layernames[lidx])




        fig.suptitle('Differences between attended and ignored context in contributions of celltypes in layers, %d mice\ntimepointwise difference between leave-out specific group and leave-out 5x randomized same # cells\npositive values mean that celltypes in layer contribute more than average to decoding'%(n_mice))
                     #Two right columns: difference of contributions to visual decodability between attended and ignored context


        save = 1
        if save:
            fig.savefig(resultpath+'micetribe,celltypeslayers,contributiontrajectory-%s_%dsv_%s-%dms%dms'%(examine,savvywidth,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)






    
    

    if 0:
        layerlist = [0,1,2,3]
        R = np.zeros(   ( len(taskaspects),len(celltypelist),len(layerlist) )   )
        D = np.zeros(   ( len(taskaspects),len(celltypelist),len(layerlist) )   )
        for cx,comparison in enumerate(taskfilenames):
            for cidx in celltypelist:
                for lidx in layerlist:
                    
                    R[cx,cidx,lidx] = -(a_diff_leave[cx,cidx,lidx,:,0] - a_diff_leave_random[cx,cidx,lidx,:,0]).mean()
                    D[cx,cidx,lidx] = -(a_dec_leave[cx,cidx,lidx,:,0] - a_dec_leave_random[cx,cidx,lidx,:,0]).mean()
                    
                    
                        
        fig,ax = plt.subplots(4,4,figsize=(36,36))
        for cidx in celltypelist:
            for lidx in layerlist:
                axs = ax[lidx,1-cidx]
                axs.bar([0,1,2,3],R[:,cidx,lidx],color=colors_layers[cidx,lidx],tick_label=taskaspects)
                axs.set_xticklabels(taskaspects,rotation=15)
                figs.plottoaxis_chancelevel(axs,0.0)
                axs.set_ylim(-0.15,0.15)
                if lidx==0: axs.set_title('rate\n'+spikewidthnames[cidx])
                if cidx==0: axs.set_ylabel(layernames[lidx])
                
                axs = ax[lidx,1-cidx+2]
                axs.bar([0,1,2,3],D[:,cidx,lidx],color=colors_layers[cidx,lidx],tick_label=taskaspects)
                axs.set_xticklabels(taskaspects,rotation=15)
                figs.plottoaxis_chancelevel(axs,0.0)
                axs.set_ylim(-0.02,0.02)
                if lidx==0: axs.set_title('decoder\n'+spikewidthnames[cidx])

    
        fig.suptitle('Aspects representation by celltypes in layers, %d mice\nauc differences between leave-out specific group and leave-out 5x randomized same # cells\npositive values mean that celltypes in layer contribute to decoding'%(n_mice))

        save = 1
        if save:
            fig.savefig(resultpath+'micetribe,celltypeslayers-%s_aucall_%s-%dms%dms'%(examine,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)













def aggregatemice_layertypecontributions_firingratestatistics(datanames,examine='allattendallignore'):
    # compare layer and celltype contributions to various task variable representations with firing rate comparisons
    
    n_mice = len(datanames)
    layerlist = [1,2,3,4]       # choose which layer to focus on: 0 all, 1 superf, 2 input, 3 deep5, 4 deep6
    celltypelist = [0,1]        # choose which spike width to use: 0 narrow, 1 broad
    
    showallcells = n_mice==1

    if showallcells: rates = []; ratese = []
    
    # collect across mice
    for n,dn in enumerate(datanames):
        block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=False,recalculate=False)
#        print(dn,block.annotations['idents'])
        
        blv,bla = preprocess.getorderattended(dn)
        if examine=='allexpcond':
            comparisongroups  = [ \
                                    [ [ [2,4], [45],        [] ],    [ [2,4],[135],        [] ] ],\
                                    [ [ [2,4], [],        [5000]  ], [ [2,4],[],        [10000] ] ],\
                                    [ [ [blv[1]],   [],     [] ],    [ [bla[1]],[],        [] ] ],\
                                    [ [], [] ]\
                                ]
            taskaspects = ['visual','audio','context','choice']
            logical = 'and'
            colors = ['navy','darkgreen','mediumvioletred','orange']
            classlabels = [['  45°','135°'],[' 5 kHz','10 kHz'],['attend visual','attend audio'],['lick','withhold lick']]
        elif examine=='attendignore':
            comparisongroups  = [ \
                                    [ [ [blv[1]], [45],        [] ],  [ [blv[1]],[135],        [] ] ],\
                                    [ [ [bla[1]], [45],        [] ],  [ [bla[1]],[135],        [] ] ],\
                                    [ [ [bla[1]],   [],    [5000] ],  [ [bla[1]],[],      [10000] ] ],\
                                    [ [ [blv[1]],   [],    [5000] ],  [ [blv[1]],[],      [10000] ] ]\
                                ]
            taskaspects = ['attend visual','ignore visual','attend audio','ignore audio']
            logical = 'and'
            colors = ['blue','navy','green','darkgreen']
        elif examine=='allattendallignore':
            comparisongroups  = [ \
                                    [ [ [blv[1]], [],        [] ] ],\
                                    [ [ [bla[1]], [],        [] ] ]\
                                ]
            taskaspects = ['attend','ignore']
            logical = 'and'
            colors = ['darkred','navy']
        elif examine=='context':
            comparisongroups  = [ \
                                    [ [ [blv[1]], [],        [] ],  [ [bla[1]], [],        [] ] ],\
                                ]
            taskaspects = ['context']
            logical = 'and'
            




        


        if n==0: # initialize; for {m,e} "_t" pooling trials, "_" aggregates across mice
            if examine=='context' or examine=='allexpcond':
                wave_ids = []
                layer_ids = []
                features_attend = []
                features_ignore = []

            if showallcells:             # single mouse and channels:
                rss = np.empty(  (n_mice, len(taskaspects),len(layerlist),len(celltypelist) ), dtype=object )
                rates = np.empty(  (n_mice, len(taskaspects),len(layerlist),len(celltypelist) ), dtype=object )
                ratese = np.empty(  (n_mice, len(taskaspects),len(layerlist),len(celltypelist) ), dtype=object )
            # all mouse and all channels
            rs = np.empty(  ( len(taskaspects),len(layerlist),len(celltypelist) ), dtype=object )
            rscl =  np.empty(  ( len(taskaspects),2,len(layerlist),len(celltypelist) ), dtype=object )
            m = pq.Quantity( np.zeros( ( len(taskaspects),len(layerlist),len(celltypelist),block.segments[0].analogsignals[0].shape[0]  ) ), units=block.segments[0].analogsignals[0].units)
            e = m.copy()
            m_cl = pq.Quantity( np.zeros( ( len(taskaspects),2,len(layerlist),len(celltypelist),block.segments[0].analogsignals[0].shape[0]  ) ), units=block.segments[0].analogsignals[0].units)
            e_cl = m_cl.copy()
            n_cells = np.zeros(  ( len(layerlist),len(celltypelist) )  )# keep count of each celltype and layer
            

        
        
        
        
        # create the features vectors for each cell across mice for clustering
        if examine=='context':
            layer_ids.extend(block.annotations['idents'])
            wave_ids.extend(block.annotations['waveforms'])
            responses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[0],logical=logical)
            for clx in range(len(responses)):
                response = np.array(responses[clx])
                
                bl = response[:,T['start_idx']:T['stimstart_idx'],:].mean(axis=(0,1))
                ev = response[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=(0,1))
                ev_diff = ev-bl
                ev_rel = ev/bl
                tr_e = (response[:,T['stimstart_idx']:T['stimstart_idx']+30,:].mean(axis=(0,1)) ) 
                tr_p = (response[:,T['stimend_idx']:T['stimend_idx']+30,:].mean(axis=(0,1)))
                tr_e_rel = (response[:,T['stimstart_idx']:T['stimstart_idx']+30,:].mean(axis=(0,1)) - ev) / (ev-bl)
                tr_p_rel = (response[:,T['stimend_idx']:T['stimend_idx']+30,:].mean(axis=(0,1)) - ev) / (ev-bl)
                F = [bl, ev, ev_diff, ev_rel, tr_e, tr_p, tr_e_rel, tr_p_rel]
                if clx==0: features_attend.extend(  np.array(F).T  )
                else:      features_ignore.extend(  np.array(F).T  )

        
        
        
        
        # assess the trial average firing rates:    
        onlyresponsive = True
        if onlyresponsive: # select only responsive neurons:
            m_b,s_b,m_e,s_e = neph.getfiringratestatistics(block,T)   # get mean and std of baseline and evoked
            responsivecells = np.logical_and( np.abs(m_b + s_b*0.2 - m_e) > 0, (m_e+m_b)/2 > 0.5)
        else:   # or select all neurons:
            responsivecells = np.ones((block.segments[0].analogsignals[0].shape[1]),dtype='int16')
        
                        

        for cidx in celltypelist:
            for layered in layerlist:
                lidx=layered-1
                lmask = np.logical_and(block.annotations['idents']==layered, block.annotations['waveforms']==cidx) 
                mask = np.where(  np.logical_and( lmask, responsivecells ) )[0]  # combine histology presence and responsiveness
                if len(mask)<1: continue  # skip if no cells with these histological type
                for cx,comparison in enumerate(taskaspects):
                    
                    # collect neural responses
                    if not comparison=='choice': # visual, audio, context:
                        response = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],logical=logical)
                    else:  # choice:
                        response = preprocess.collect_stimulusspecificresponses_choice(block,dn)
                    
                    if examine=='allexpcond':
                        cxclset = [0,1]
                    else:
                        cxclset = [0]
                    
                    for cxcl in cxclset:
                    
                        if showallcells:
                            if rss[n,cx,lidx,cidx] == None: rss[n,cx,lidx,cidx] = []
                            aux = np.array(response[cxcl])[:,:,lmask]
                            rss[n,cx,lidx,cidx].extend(aux)    # add to the individual mouse, per-channel trial averages
    
                        if rs[cx,lidx,cidx]==None: rs[cx,lidx,cidx] = []
                        if rscl[cx,cxcl,lidx,cidx]==None: rscl[cx,cxcl,lidx,cidx] = []
    
                        aux = np.array(response[cxcl])[:,:,lmask].swapaxes(2,1).reshape((-1,block.segments[0].analogsignals[0].shape[0]))
                        rs[cx,lidx,cidx].extend(aux)        # add to the across mice responses for all channels all trial averages into m+/-e
                        
                        rscl[cx,cxcl,lidx,cidx].extend(aux)
                        

                # and increase the counter for number of cells available in this celltype and layer
                
                n_cells[lidx,cidx] += len(mask)
                


    # create the trial averages:
    for cx,comparison in enumerate(taskaspects):
        for cidx in celltypelist:
            for layered in layerlist:
                lidx = layered-1
                if showallcells:
                    for n in range(n_mice):
                        if rss[n,0,lidx,cidx]!=None:
                            aux = pq.Quantity(rss[n,cx,lidx,cidx], units=block.segments[0].analogsignals[0].units)
                            rates[n,cx,lidx,cidx] = aux.mean(axis=0)
                            ratese[n,cx,lidx,cidx] = aux.std(axis=0)/np.sqrt(aux.shape[0])*2
                
                if rs[0,lidx,cidx]!=None:
                    aux = pq.Quantity(rs[cx,lidx,cidx], units=block.segments[0].analogsignals[0].units)
                    m[cx,lidx,cidx,:] = aux.mean(axis=0)
                    e[cx,lidx,cidx,:] = aux.std(axis=0)/np.sqrt(aux.shape[0])*2
                    
#                print(cx,lidx,cidx,'set:   ',cxcl,cxclset,cxcl in cxclset,rscl[0,0,lidx,cidx]!=None)
                if rscl[0,0,lidx,cidx]!=None:
                    for cxcl in cxclset:
                        aux = pq.Quantity(rscl[cx,cxcl,lidx,cidx], units=block.segments[0].analogsignals[0].units)
                        m_cl[cx,cxcl,lidx,cidx,:] = aux.mean(axis=0)
                        e_cl[cx,cxcl,lidx,cidx,:] = aux.std(axis=0)/np.sqrt(aux.shape[0])*2
                    


    if 0 and showallcells:
        pickle.dump(rss,open(cacheprefix+'phys/rate,trajectories,celltypes,allmice-%s.pck'%(continuous_method),'wb'))

        
        
        

    # create the clustering
#    if examine=='context':
#        import sklearn.cluster as skcl
#        cla = skcl.KMeans(n_clusters=4)
    wave_ids = np.array(wave_ids)
    layer_ids = np.array(layer_ids)
    features = np.array([features_attend,features_ignore])

    


    # plots
    
    if 1:       # individual rate displays
        layernamesroman=['Layer I+II/III','Layer IV','Layer V','Layer VI']
        if showallcells: fig,ax = plt.subplots(4,6,figsize=(48,24))
        elif examine!='allexpcond': fig,ax = plt.subplots(4,2,figsize=(16,24))
        elif examine=='allexpcond': fig,ax = plt.subplots(4,8,figsize=(16*4,24))
        for cx,comparison in enumerate(taskaspects):
            for cidx in celltypelist:
                for layered in layerlist:
                    lidx=layered-1
                    
        
                    if examine=='allexpcond':
                        for clcx in cxclset:
                            axs = ax[lidx,1-cidx+cx*2]
                            figs.plottoaxis_plottrajectory(axs,response[0][0].times,m_cl[cx,clcx,lidx,cidx,:],e_cl[cx,clcx,lidx,cidx,:],\
                                                           colors[cx],linewidth=2,alpha=(1-0.4*clcx)/(4-cidx)*3)
                        if cidx==1 and lidx==0:
                            axs.legend(classlabels[cx])
                            if cx==0: axs.text(-1000,4.5,'2 s.e.m.')
                            if publish:
                                figs.labelaxis(axs,y=1.12,panel=['A    visual grouped trials','B    audio grouped trials','C    context grouped trials','D    choice grouped trials'][cx])

                        axs.set_xlim(T['starttime']+200*pq.ms,T['endtime']-200*pq.ms) # -1500,4500
                        if cidx==0: axs.set_ylim(4,18); axs.set_yticks([5,10,15])
                        else: axs.set_ylim(2,6); axs.set_yticks([2,4,6])
#                        figs.plottoaxis_chancelevel(axs)
                        figs.plottoaxis_stimulusoverlay(axs,T)
                        if lidx==0: axs.set_title(spikewidthnames[cidx])
                        figs.setxt(axs)
                        if cidx==1 and cx==0: axs.set_ylabel(layernamesroman[lidx]+'\n[Hz]')
                        axs.spines['right'].set_visible(False)
                        axs.spines['top'].set_visible(False)
                        
                    else:
                        axs = ax[lidx,1-cidx]
                        figs.plottoaxis_plottrajectory(axs,response[0][0].times,m[cx,lidx,cidx,:],e[cx,lidx,cidx,:],\
                                                       colors[cx],linewidth=2,label=taskaspects[cx])
                        axs.legend(taskaspects)
                        figs.plottoaxis_chancelevel(axs)
                        
                    if cx==1 and examine!='allexpcond':
                        axs.set_xlim(T['starttime'],T['endtime']) # -1500,4500
                        if continuous_method=='count':
                            if cidx==0: axs.set_ylim(0,20)
                            else: axs.set_ylim(0,10)
    #                    axs.set_ylim(-0.5*pq.Hz,np.max(m+e)+0.5*pq.Hz)
                        figs.plottoaxis_stimulusoverlay(axs,T)
    #                    figs.plottoaxis_stimulusoverlay(axs,T,offphase=[-6000*pq.ms,-3000*pq.ms])
    #                    figs.plottoaxis_stimulusoverlay(axs,T,offphase=[6000*pq.ms,9000*pq.ms])
                        if cidx==1: axs.set_ylabel(layernames[lidx])
                        axs.set_title(spikewidthnames[cidx]+' $n_{cells}$=%d'%(n_cells[lidx,cidx]))
                        figs.setxt(axs)
    
    
    
                    if showallcells:
                        axs = ax[lidx,(1-cidx)*2+2+cx]
                        for n,dn in enumerate(datanames):
                            if type(rates[n,cx,lidx,cidx])==pq.Quantity:
                                figs.plottoaxis_plottrajectories(axs,response[0][0].times,\
                                        rates[n,cx,lidx,cidx], ratese[n,cx,lidx,cidx],['Reds','Blues'][cx],linewidth=2,alpha=0.33)
                        axs.legend([taskaspects[cx]])
                        axs.set_xlim(T['starttime'],T['endtime']) # -1500,4500
                        if continuous_method=='count':
                            if cidx==0: axs.set_ylim(0,45)
                            else: axs.set_ylim(0,25)
                        figs.plottoaxis_stimulusoverlay(axs,T)
                        figs.plottoaxis_chancelevel(axs)
                        axs.set_title(spikewidthnames[cidx])
                        figs.setxt(axs)
    
    
    
        if showallcells:    
            fig.suptitle(datanames[0]+'  '+taskaspects[0]+', '+taskaspects[1]+'  trial averaged rates +/- 2 s.e.m. [Hz], %s, %d ms\naveraged over each cell groups (columns 1-2)\nall cells individually (columns 3-6), shades correspond to same cells'%(continuous_method,T['bin']))
        elif examine!='allexpcond':
            fig.suptitle(taskaspects[0]+', '+taskaspects[1]+'  trial averaged rates  +/- 2 s.e.m. [Hz], %s, %d ms, $n_{mouse}$=%d'%(continuous_method,T['bin'],n_mice))
    
        
    
        save = 0
        if save:
            if showallcells:
                fig.savefig(resultpath+'%s,celltypeslayers-%s_firingrate_%s-%dms%dms'%(datanames[0],examine,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
            else:
                fig.savefig(resultpath+'micetribe,celltypeslayers-%s_firingrate,responsive_%s-%dms%dms'%(examine,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
            





    if 0 and examine=='context':
        colors_layers = np.array([['lightcoral','deeppink','red','darkred'],\
                                  ['deepskyblue','darkviolet','blue','navy']])   # # narrow,broad / superf, input, deep5, deep6
        labels = ['narrow superf', 'narrow input','narrow deep5','narrow deep6',\
                  'broad superf', 'broad input','broad deep5','broad deep6']
        axislabels = ['baseline','evoked','evoked-baseline difference','evoked/baseline relative',\
                      'onstim transient','poststim transient','onstim transient relative','poststim transient relative']
        
        fig,ax = plt.subplots(int(features.shape[2]/2),2,figsize=(20,40))
        for px in range(2):
            for fx in range(int(features.shape[2]/2)):
                axs = ax[fx,px]
                
                for cidx in range(2):
                    for lidx in range(4):
                        axs.plot( features[px, (layer_ids==lidx) & (wave_ids==cidx), fx*2 ],  \
                                  features[px, (layer_ids==lidx) & (wave_ids==cidx), fx*2+1 ],  \
                                  'o', color=colors_layers[cidx,lidx]   )
                axs.legend(labels)
                axs.set_title(['attend visual','ignore visual'][px])
                axs.set_xlabel(axislabels[fx*2])
                axs.set_ylabel(axislabels[fx*2+1])



















def aggregatemice_decoderaccuracyaverages(datanames,examine='allexpcond'):
    # average out accuracies between mice
    
    n_mice = len(datanames)
    colorlist = plt.cm.viridis( np.linspace(0, 0.8, n_mice) )
    
    stimulusspecificity = 'all'               # marginal on all other conditions
    taskaspects = ['visual','audio','context','choice']
    trajectorylabels = ['pre','stim 0-1.5s','stim 1.5s-3s','post']
    trajectoryendidx = np.array([  0, T['stimstart_idx'], (T['stimstart_idx']+T['stimend_idx'])/2, T['stimend_idx'], T['end_idx']   ], dtype='int16')
    trajectorypointcenters = trajectoryendidx[:-1] + (trajectoryendidx[1])/2
    trajectorypointcentertimes = trajectorypointcenters*T['dt']+T['starttime']
    trajectorylength = 550#596
    trajectoryendidx[3] = trajectorylength

    block = preprocess.loaddatamouse(datanames[0],T,continuous_method,recalculate=False)
    times_diff = block.segments[0].analogsignals[0].times
    times = times_diff[:596]
    
    
    stimulusspecificity = 'all'
    layeredfilename = 'all'
    
    
    
    avaccs = np.zeros((len(taskaspects),len(trajectorylabels),n_mice))
    avaccs_se = np.zeros((len(taskaspects),len(trajectorylabels),n_mice))
    
    fig,ax = plt.subplots(3,4,figsize=(36,24))
    for cx,comparison in enumerate(taskaspects):
        acrossdecoders = []
        for n,dn in enumerate(datanames):
            acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            acrossdecoders.append(acrossdecoder)
            for tx,t in enumerate(trajectoryendidx[1:]):
                avaccs[cx,tx,n] = acrossdecoder[1][trajectoryendidx[tx]:trajectoryendidx[tx+1],0].squeeze().mean(axis=0)
                avaccs_se[cx,tx,n] = acrossdecoder[1][trajectoryendidx[tx]:trajectoryendidx[tx+1],2].squeeze().mean(axis=0)*2
                


        axs = ax[2,cx]
        axs.boxplot(x=avaccs[cx,:,:].squeeze().T, positions=trajectorypointcentertimes.magnitude,\
                    widths=trajectorypointcenters[1]*T['dt'].magnitude*0.3,notch=False, whis=[5,95]  )
        if cx==0: axs.set_ylabel('(C)\ndistribution of mean accuracies')
        
        ip_n = [2,6,7,10]#np.random.permutation(n)
        for n,dn in enumerate(datanames):
            for dx in range(3):
                axs = ax[dx,cx]
                if dx==0: axs.set_title(comparison)
                if cx==0 and dx==0: axs.set_ylabel('(A)\ndecoder trajectories, example single mice')
                
                if dx==0:
                    if n in ip_n[:5]:
                        figs.plottoaxis_plottrajectory(axs,times,\
                                                   acrossdecoders[n][1][:596,0].squeeze(),\
                                                   acrossdecoders[n][1][:596,2].squeeze(),\
                                                   colorlist=colorlist[n], label=dn)
    
                if dx==1: axs.errorbar(trajectorypointcentertimes+30*(n-n/2)*pq.ms,avaccs[cx,:,n],yerr=avaccs_se[cx,:,n],\
                             capsize=3, alpha=0.5, linewidth=1, linestyle='',marker='o', label=dn, color=colorlist[n])
    
                if cx==0 and dx==1: axs.set_ylabel('(B)\ndecoder acc. time-averages')
    
                if cx==0 and dx<2: axs.legend(fontsize='x-small',frameon=False)
                axs.set_ylim(0.399,1.001)
                axs.set_xlim(-1500,6500)
                figs.setxt(axs,[0,3000])
                if n==0: figs.plottoaxis_stimulusoverlay(axs,T)
                if n==0: figs.plottoaxis_chancelevel(axs,0.5)
    
    
    fig.suptitle('Visual and audio identity, context and choice decoder performance across animals and across time\n'+\
                 '(A) decoder dynamics of example mice trajectories\n(B) time-averaged at important timepoints from trial onset\n'+\
                 '(C) distribution of the mean values from (B): orange: median, box: 25-75 perc., whiskers: 5-95 perc.')
    
    save = 0
    if save:
        fig.savefig(resultpath+'micetribe,decoderperformances_%s-%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
        








def aggregate_spontaneousprojections(datanames):
    # show activitiy projections to various decision vector subspaces for all mice during spontaneous activity (pre stimulus)
    examine = 'allexpcond'
    layeredfilename = 'all'
    print(datanames)
    n_mice = len(datanames)
    
    taskaspects = ['visual','audio','context','choice',\
                       'char,45-135','char,90-180','char,270-180','char,60-150','char,120-210']
    angleoffset = []
    dbprojections = [0,4,5,6,7,8]      # visual and visual in characterization block

    tribe_stds = np.zeros((n_mice,len(dbprojections),3,2))     #  [n_mice,comparisons,classes,pre-onstim]
    for n,dn in enumerate(datanames):
        projections = pickle.load(open(cacheprefix+'subspaces/projections,character-%s_%s-%s-%dms%dms_%s.pck'%(examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn),'rb'))
        for cx in range(len(dbprojections)):
            for clx in range(2):
                for pox in range(2):
                    tribe_stds[n,cx,clx,pox] = projections[cx][clx][pox].std()
                    tribe_stds[n,cx,2,pox] = np.hstack( (projections[cx][0][pox],projections[cx][1][pox] )).std()
        
       
    
    print(tribe_stds.shape)
    
    
    
    if 0:
        classcolors = ['darkcyan','mediumvioletred','navy']
        fig,ax = plt.subplots(4,len(dbprojections)-2,figsize=((len(dbprojections)-2)*9,4*9))
        for cx,tx in enumerate(dbprojections[2:]):
            for pox in range(2):
    
                axs = ax[pox,cx]
                axs.plot(tribe_stds[:,0,2,pox],tribe_stds[:,cx+2,2,pox],'o',alpha=1-clx/0.86,color=classcolors[2])
                for n,dn in enumerate(datanames):
                    axs.text(tribe_stds[n,0,clx,pox]+0.02,tribe_stds[n,cx+2,clx,pox]-0.01,dn,color='k',alpha=0.3,fontsize=10)
                axs.set_ylabel(taskaspects[tx])
                axs.set_xlabel('task,45-135')
                if cx==0: axs.legend(['both classes'])
                axs.set_xlim(0,6)
                axs.set_ylim(0,1.65)
    
                x = tribe_stds[:,0,2,pox]; y = tribe_stds[:,cx+1,2,pox]
                l = sp.stats.linregress(x,y)
    #            if l[3]<0.05:
                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                axs.plot(line,l[0]*line+l[1],color='darkorange',linewidth=2)
                axs.text(3.3,0.2,'p=%4.2f, $R^2$=%4.2f'%((l[3],l[2]**2)),color='darkorange')
    
    
    
    
                
                
                axs = ax[pox+2,cx]
                for clx in range(2):
                    axs.plot(tribe_stds[:,0,clx,pox],tribe_stds[:,cx+2,clx,pox],'o',alpha=1-clx/0.86,color=classcolors[clx])
                    for n,dn in enumerate(datanames):
                        axs.text(tribe_stds[n,0,clx,pox]+0.02,tribe_stds[n,cx+2,clx,pox]-0.01,dn,color='k',alpha=0.3,fontsize=10)
                axs.set_ylabel(taskaspects[tx])
    #            if cx==0:
    #                axs.set_ylabel('%s\n%s'%(['spontaneous','evoked'][pox],taskaspects[0]))
                axs.set_xlabel('task,45-135')
                if cx==0: axs.legend(['class 1','class 2'])
                axs.set_xlim(0,6)
                axs.set_ylim(0,1.65)
                
                
                x = tribe_stds[:,0,:,pox].ravel(); y = tribe_stds[:,cx+1,:,pox].ravel()
                l = sp.stats.linregress(x,y)
    #            if l[3]<0.05:
                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                axs.plot(line,l[0]*line+l[1],color='darkorange',linewidth=2)
                axs.text(3.3,0.2,'p=%4.2f, $R^2$=%4.2f'%((l[3],l[2]**2)),color='darkorange')
                # line[1]-2,l[0]*line[1]+l[1]-0.02,   <- original x,y points



    if 1:
        fig,ax = plt.subplots(2,1,figsize=(12,18))
        for pox in range(2):
            axs = ax[pox]
#            for cx,tx in enumerate(dbprojections[2:]):
            axs.boxplot(x=tribe_stds[:,:,2,pox],positions=np.arange(len(dbprojections)))
            axs.set_xticklabels(np.array(taskaspects)[np.array(dbprojections)],rotation=15)
            axs.set_ylim(0.3,1.5)









def aggregatemice_subspaces_behaviour(datanames):
    # select attend condition
    # find the visual dbnv
    # project activities onto that vector
    # label them based on animal choice
    # see if they are better or worse discriminable
    
    doplot = 0 or globaldoplot
    dump = 1
    
    
    depth_idx = int((1500*pq.ms / T['dt'] ).magnitude)      # pre and on stimulus averaging depth
    
    # taskcontexts = ['attendvisual','ignorevisual']
    # taskcontexts = ['visualattend','visualignore']
    # taskcontexts = ['00']
    # taskcontexts = ['visualattend','audioattend']
    taskcontexts = ['visualattend']
    
    for tcx,taskcontext in enumerate(taskcontexts):
        # taskaspects = [taskcontext,'visual']
        taskaspects = [taskcontext,'context']
        # taskaspects = ['context']
    
        # triallistindex = [1]           # use this for visual: which task aspect to use the classes for as data points (2,3 = context visual go and nogo)
        # dbprojections = [0,2,4,5,6,9,10]            # visual-context-choice for 3D, +visual attend, visual ignore
        dbprojections = [0]                          # attend/ignore, project to visual dbnv 0
        # dbprojections = [1]                          # attend/ignore, project to context dbnv 1
    
    
    #        colors = ['navy','darkturquoise','darkred','fuchsia']
        colors = ['b','c','m','r']
        # basiscolors = ['navy','blue','dodgerblue','fuchsia','sienna','gold','darkorange']
        # basisvectors = [0,1]      # indexes of dbprojections, dbprojections was:
        labels = [\
                   ['attended visual dec. normal','attended visual choice',\
                    'attend visual,  45, hit','attend visual,  45, miss','attend visual, 135, corr.rej.','attend visual, 135, false alarm'], \
                   ['','',\
                    ' hit',' miss',' correct rejection',' false alarm'] \
                 ]
        
        
        # go over mice
        projections_all = []
    
        n_mice = len(datanames)
        behavfrac = np.zeros((n_mice,2))   # (mice, 45 135 successes)
        for n,dn in enumerate(datanames):
            block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=True,recalculate=False)
    
            n_trajectory,n_neuron = block.segments[0].analogsignals[0].shape
            
            blv,bla = preprocess.getorderattended(dn)
        
            comparisongroups =  [\
                                   [ [ [blv[1]], [45],    [] ], [ [blv[1]], [135],     [] ] ], \
                                   [ [ [bla[1]], [],    [5000] ], [ [bla[1]], [],     [10000] ] ], \
                                   [ [ [blv[1]], [],    [] ],   [ [bla[1]], [],     [] ] ], \
                                ]
            # remove the not focused attend/ignore, so that comparisongroups will be proper
            del comparisongroups[1-tcx]        # delete the other context not used in this tcx iteration
            # del comparisongroups[:2]
            print(comparisongroups)
            
        
        
            # collect responses
            prestimactivity = []       # this will have for each comparisongroup a [taskaspect][classes][trials,neurons] list of matrices
            onstimactivity = []
            for cx,comparison in enumerate(taskaspects):
    
                # print('comparison',cx,comparisongroups[cx])
    
                stimulusIDgroups = comparisongroups[cx]
                # collect neural responses; first two     45 hit 135 correct rejection, next two 45 miss, 135 false alarm; then reorder for hmcf

                acrossresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups,correctonly=1)
                acrossresponses.extend( preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups,erroronly=1) )
                acrossresponses = [acrossresponses[i] for i in [0,2,1,3]]  # reorder for h,m,cr,fa
        
                auxpre = []
                auxon = []
                for clx in range(4):   # h,m,c,f
                    if dn=='DT008' and cx==11: # this mouse doesn't have characterization
                        auxpre.append(np.zeros((1,1)))
                        auxon.append(np.zeros((1,1)))
                    else:
                        if len(acrossresponses[clx])>0:
                            auxpre.append(np.array( acrossresponses[clx] )[:,:T['stimstart_idx'],:].mean(axis=1).squeeze())
                            # auxon.append(np.array( acrossresponses[clx] )[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1).squeeze())
                            auxon.append(np.array( acrossresponses[clx] )[:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,:].mean(axis=1).squeeze())
                        else:
                            auxpre.append(np.zeros((0,n_neuron)))
                            auxon.append(np.zeros((0,n_neuron)))
                prestimactivity.append(auxpre)
                onstimactivity.append(auxon)
                # print(tcx,taskcontext,cx,comparison,onstimactivity[cx][0].shape,onstimactivity[cx][1].shape)
                    
            # collect DBNVs, to project the activities collected above onto
            c_db = []
            for cx,comparison in enumerate(taskaspects):
                acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
                wx = int((len(acrossdecoder)-7)/n_neuron)
                c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
                
            c_db = np.array(c_db) # [comparisongroup,neurons,trajectory,stats]
        
        
            # mean context decision normal vector with averaging over EA:
            c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,:].mean(axis=2)
            c_db_means_matrix = np.vstack( [c_db_means[dbprojections[i]][:,0] for i in range(len(dbprojections))] ).T
        
        
            # find the orthogonal axes
            Q,R = np.linalg.qr(c_db_means_matrix)     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
        
        
            dbnbasiscomps = np.dot(Q.T,c_db_means_matrix) / np.linalg.norm(c_db_means_matrix,axis=0,keepdims=True)         # orthonormbasis x dbbasistypes
        
        
            # transform the activity vectors to the new orthonormal basis
            projections = []  #{visual,context,choice} x {pre,onstim} x {go hit, go miss, nogo cr, nogo falseal}
            for cx,_ in enumerate(taskaspects):  # go over activities to project them
                for coordx,_ in enumerate(dbprojections):        # go over the desired projection coordinates: visual/audio, choice, context; cx is the taskaspects index, coords is the Q coordinate in order of orthonormali
                    # print('projection index, dbnv number',coordx,cx,'from dbprojections',dbprojections[:3])
                    projected_prestim = []
                    projected_onstim = []
                    for clx in range(4):  # h,m,c,f # go over the classes and success
                        # print(clx,len(onstimactivity[cx][clx]),Q[:,coordx].shape)
                        # for trlx in triallistindex: # go over the activity groups: here context attend, and context attend other
                        # project pre- and onstimulus activities separately
                        # print(clx,trlx)
                        # print('pre',len(prestimactivity), 'on',len(onstimactivity),)
                        # print('pre',len(prestimactivity[trlx]),len(prestimactivity[trlx][clx]),'on',len(onstimactivity[trlx]),len(onstimactivity[trlx][clx]))
                        projected_prestim.append(  np.atleast_1d(  np.dot(   prestimactivity[cx][clx],Q[:,coordx]  ) ) )
                        projected_onstim.append(  np.atleast_1d(  np.dot(   onstimactivity[cx][clx],Q[:,coordx] ) ) )
                    projections.append( [projected_prestim,projected_onstim] )
            # print(len(projections),len(projections[0]),len(projections[0][1][3]))
            # projections = np.array(projections)
            projections_all.append(projections)
    
        
        activities_all = []     # this will be [n_mice][comparison][performance] = 15 x 4 x 4
        activities_aggregate = [   [  [],[],[],[]   ] for i in range(len(taskaspects))  ]
        for n,dn in enumerate(datanames):
            print('collecting results', dn)
            # ixv45,ixv135,ixa45,ixa135,ixv5000,ixv10000,ixa5000,ixa10000 = preprocess.getstimulusidents(dn,block,multimodalonly=True)
            # ixhit,ixmiss,ixcorrrej,ixfal = preprocess.assessbehaviouralperformance(dn)
            
            # print(len(projections_all),len(projections_all[n]),len(projections_all[n][0]),len(projections_all[n][0][1]))
            # print('full',projections_all[n][0][1])
            # print([ len(projections_all[n][0][1][i]) for i in [0,1,2,3] ])

            activities_comparison = []
            for cx,comparison in enumerate(taskaspects):
            # for cx,_ in enumerate(dbprojections):        # go over the desired projection coordinates: visual/audio, choice, context; cx is the taskaspects index, coords is the Q coordinate in order of orthonormali
                # print(cx,len(projections_all[n][cx][1][0]),len(projections_all[n][cx][1][1]))

                
                preon = 1 # onstim
                # preon = 0 # pre stim, for context only!




                # if tcx in [0]:   # attend
                #     _,ixgohit,_ = np.intersect1d(ixv45,ixhit,return_indices=True)
                #     _,ixnogocorrrej,_ = np.intersect1d(ixv135,ixcorrrej,return_indices=True)
                #     _,ixgomiss,_ = np.intersect1d(ixv45,ixmiss,return_indices=True)
                #     _,ixnogofal,_ = np.intersect1d(ixv135,ixfal,return_indices=True)
                # elif tcx in [1]:  # ignore
                #     _,ixgohit,_ = np.intersect1d(ixa45,ixhit,return_indices=True)
                #     _,ixnogocorrrej,_ = np.intersect1d(ixa135,ixcorrrej,return_indices=True)
                #     _,ixgomiss,_ = np.intersect1d(ixa45,ixmiss,return_indices=True)
                #     _,ixnogofal,_ = np.intersect1d(ixa135,ixfal,return_indices=True)
    
                # behavfrac[n,0] = len(ixgohit)/len(ixv45)
                # behavfrac[n,1] = len(ixnogocorrrej)/len(ixv135)
                
                # h/h+m,   c/c+f:
                behavfrac[n,0] = len(projections_all[n][cx][preon][0])/(len(projections_all[n][cx][preon][0])+len(projections_all[n][cx][preon][1]))
                behavfrac[n,1] = len(projections_all[n][cx][preon][2])/(len(projections_all[n][cx][preon][2])+len(projections_all[n][cx][preon][3]))



                # collect activities for display
        
        
                # flip
                # print(projections_all[n][cx][preon][0])
                # print(ixgohit)
                # print(projections_all[n][cx][preon][0][ixgohit])
                ix = np.argmax(np.abs(projections_all[n][cx][preon][0]))
                sign = np.sign(projections_all[n][cx][preon][0][ix])
                
                # single mice
                activity = [ sign*projections_all[n][cx][preon][i] for i in [0,1,2,3] ]
                activities_comparison.append(activity)
                
                # concatenate all trials for all mice
                for k in [0,1,2,3]: # h,m,c,f
                    activities_aggregate[cx][k].extend(activity[k])
        
            activities_all.append(activities_comparison)
        
        
        
        if dump:
            pickle.dump((activities_all,activities_aggregate),open(cacheprefix+'subspaces/dbnvprojections-visual,behaviourdifference-%s-%dms.pck'%(continuous_method,T['dt']),'wb'))
            
        
        
        # plot each mouse separately
        if doplot:
            n_panels = int(np.ceil(np.sqrt(n_mice)))
            fig,ax = plt.subplots(n_panels,n_panels,figsize=(12*n_panels,8*n_panels))
            
            cx = 0    # choose the context conditioned block (0), or context (1) for displaying projected activities
            lx = 1     # labels for the conditions  0 visual        1 context
            dprimes = np.zeros((n_mice,2))   #    (n_mice x {45,135})
            for n,dn in enumerate(datanames):
                print(dn,'activity len',[len(activities_all[n][cx][i]) for i in [0,1,2,3]])
        
                # print(dn,'activity vector@ hmcf',activities_all[n][cx])
    
                # find significantly different mice with behaviour
                s45,p45 = sp.stats.ttest_ind(activities_all[n][cx][0], activities_all[n][cx][1], equal_var=False)
                s135,p135 = sp.stats.ttest_ind(activities_all[n][cx][2], activities_all[n][cx][3], equal_var=False)
        
                dprimes[n,:] = [ nedi.dprime_normal(activities_all[n][cx][0], activities_all[n][cx][1]),\
                                 nedi.dprime_normal(activities_all[n][cx][2], activities_all[n][cx][3]) ]
            
                axs = ax[int(np.floor(n/n_panels)),n%n_panels]
                
            
                for k in range(len(activity)):
                    axs.hist(activities_all[n][cx][k],bins=10,color=colors[k],label=taskcontext+labels[lx][2+k],alpha=0.5)
                # axs.set_xlim(-2,2)
                axs.plot([0,0],[0,axs.get_ylim()[1]],'--',color='grey')
                axs.legend(frameon=False)
                # axs.set_title('%s; 45: p=%5.4f, 135: p=%5.4f'%(dn,p45,p135,))
                axs.set_title('%s; p=%5.4f, p=%5.4f'%(dn,p45,p135,))
        
        
            fig.suptitle(taskaspects[dbprojections[0]])
            save = 0
            if save:
                fig.savefig(resultpath+'visualdiscrimination,behaviourconditioned,%s-singlemice_%s_%dms%dms'%(taskcontext,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)







        # plot all mice aggregate    
        if doplot:
        
            # find out statistically if performance makes a difference in discrimination peaks
            cx = 0    # only one projection
            s45,p45 = sp.stats.ttest_ind(activities_aggregate[cx][0], activities_aggregate[cx][1], equal_var=False)
            s135,p135 = sp.stats.ttest_ind(activities_aggregate[cx][2], activities_aggregate[cx][3], equal_var=False)
            # s45,p45,s135,p135 = 0,0,0,0
        
        
            # selection = np.abs(dprimes)<1
            fig,ax = plt.subplots(2,2,figsize=(24,24))
            for d in [0,1]:       # separate 45 and 135 degrees
                axs = ax[0,d]
                for k in d*2+np.array([0,1]):
                    # print(k,'hmcf'[k]+':',len(activities_aggregate[cx][k]))
                    axs.hist(activities_aggregate[cx][k],bins=20,color=colors[k],label=taskcontext+labels[lx][2+k],alpha=0.5)
                    m = np.mean(activities_aggregate[cx][k])
                    e = 2*np.std(activities_aggregate[cx][k])/np.sqrt(len(activities_aggregate[cx][k]))
                    axs.plot([m,m],[0,axs.get_ylim()[1]],lw=3,color=colors[k])
                    axs.plot([m-e,m+e],[10+k,10+k],lw=7,color=colors[k])
                axs.set_xlim(-3,3)
                axs.plot([0,0],[0,axs.get_ylim()[1]],'--',color='grey')
                axs.legend(frameon=False)
                # axs.set_title('%d mice; %d: p=%5.4f'%(n_mice,[45,135][d],(p45,p135)[d]))
                axs.set_title('%d mice; p=%5.4f'%(n_mice,(p45,p135)[d]))
            
                axs = ax[1,0]
                axs.hist( dprimes[:,d], color=colors[2*d], bins=5, label=['45 hit-miss','135 corr.rej.-fal.al.'][d], alpha=0.5 )
            axs.legend(frameon=False)
            axs.set_title('%d mice; signed d-prime'%(n_mice))
            axs.set_xlabel('signed d-prime')
            axs.set_ylabel('# of mouse')
            
            axs = ax[1,1]
            for k in [0,1]:
                x = dprimes[:,k]
                y = behavfrac[:,k]
    
                axs.scatter(x,y,marker='o',s=150,color=colors[2*k],alpha=0.8,label=['hit/45','corr.rej./135'][k])
                # for n,dn in enumerate(datanames):
                #     axs.text(x,y,dn,fontsize='x-small',color=colors[2*k],alpha=0.3)
        
        
                l = sp.stats.linregress(x,y)
                print(k,l)
                if l[3]<0.05:
                    line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                    axs.plot(line,l[0]*line+l[1],color=colors[2*k],linewidth=2)
                else:
                    line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                    axs.plot(line,l[0]*line+l[1],'--',color=colors[2*k],linewidth=2)
                    
                xoffs = 1
            #        if sx in [3,4,5,7]: xoffs = 25
                if l[3]<0.05:
                    axs.text(line[1]-xoffs,0.47+k*0.1,'p<%5.3f, $\\rho$=%4.2f'%((0.001,l[2])),color=colors[2*k])
                else:
                    axs.text(line[1]-xoffs,0.47+k*0.1,'p=%5.3f, $\\rho$=%4.2f'%((l[3],l[2])),color=colors[2*k])
            
            axs.legend(frameon=False)
            axs.set_xlabel('signed d-prime')
            axs.set_ylabel('behaviour success')
            
            fig.suptitle(taskaspects[dbprojections[0]])
            
            save = 0
            if save:
                fig.savefig(resultpath+'visualdiscrimination,behaviourconditioned-projected,%s-%s-aggregate_%s_%dms%dms'%(taskaspects[1],taskcontext,continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)







    return













def aggregatemice_subspace_projectiontodbnvs_behaviour(datanames):
    # project multiple mouse activity onto standardized task variable coordinates found by decoders

    n_mice = len(datanames)    
    projected_dynamics_all = []      #  (mice)(taskaspects)(classes)(trials,trajectory,dbnvs)
    differences = [] # (mice)(class)(mean/s.e.m)(trajectory)
    for n,dn in enumerate(datanames):
        projected_dynamics,V,B = pickle.load( open(cacheprefix+'subspaces/projections,dbnv-visual,context_%s.pck'%(dn),'rb'))
        projected_dynamics_all.append(projected_dynamics)


        dx = 0     # we are intersted only in visual dbnv projections in the two contexts and behavioural choices
        cx_av = 0  # attend visual
        cx_iv = 1  # ignore visual
        aix_h = 0  # used indices for correct trials only: hit, correct rejectoin
        aix_c = 2

        differences.append([])
        for aix in [aix_h,aix_c]:
    
            # mean
            trial_av = np.array(projected_dynamics[cx_av][aix]).mean(axis=0)[:,dx]
            trial_iv = np.array(projected_dynamics[cx_iv][aix]).mean(axis=0)[:,dx]
            
            nf = np.max(np.hstack((np.abs(trial_av),np.abs(trial_iv)))) # normalization factor
            print(dn,nf)
            trial_av /= nf
            trial_iv /= nf
        
            # s.e.m.
            trial_av_s = np.array(projected_dynamics[cx_av][aix]).std(axis=0)[:,dx]/nf
            trial_iv_s = np.array(projected_dynamics[cx_iv][aix]).std(axis=0)[:,dx]/nf
            trial_av_e = trial_av_s / np.sqrt(len(trial_av_s))
            trial_iv_e = trial_iv_s / np.sqrt(len(trial_iv_s))
            
            
            
            differences[-1].append(    [trial_av - trial_iv, (trial_av_s + trial_iv_s)*2, (trial_av_e + trial_iv_e)*2 ] )

    n_trajectory = len(differences[n][0])
    times = np.arange(T['offsettime'],T['endtime']+T['dt'],T['dt'])



    fig, ax = plt.subplots(1,2,figsize=(32,12))
    
    for aix in [0,1]:
    
        axs = ax[aix]
        
        for n,dn in enumerate(datanames):
            
            axs.plot(times,differences[n][aix][0],color=[0,0,(n+1.)/n_mice],lw=2,alpha=0.5,label=dn)
            
            axs.fill_between(times,differences[n][aix][0]-differences[n][aix][2],\
                                   differences[n][aix][0]+differences[n][aix][2],\
                                   color=np.array([0,0,(n+1.)/n_mice]),alpha=0.05)
    
        
        axs.set_title(['hits','correct rejections'][aix])
        axs.set_ylabel('contextual difference index')
        axs.set_xlabel('time from stimulus onset')
    
        axs.set_ylim(-1,1)
        axs.set_xlim(times[0],times[-1])
        axs.set_yticks(np.arange(-1,1.1,0.5))
        axs.legend(frameon=False,ncol=2)
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs)
        
        

    fig.suptitle('difference between two contexts using only hit and correct rejection trials\n(trial average and 2 s.e.m. diff. / max || ) of activities projected onto visual dbnv')

    save = 0
    if save or globalsave:
        fig.savefig(resultpath+'tribe,contextdifference,projectonto,dbnv,visual_%s_%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)



    return








def aggregatemice_subspace_anglesbetweendbnvs(datanames):
    # aggregate for all mice angles between decision vector directions
    
    taskaspects = ['visual','audio','context','choice','choice,av','choice,aa']

    # pairlists = [[[0,2]],[[0,3],[0,4],[0,5]], [[2,3],[2,4],[2,5]]   ]
    pairlists = [[[0,2]], [[1,2]],  [[2,4],[2,5]],    [[0,4],[0,5]],     [[1,4],[1,5]]    ]
    titles = ['visual-context','audio-context','choice-context','visual-choice','audio-choice']
    
    singlecolors = [['dodgerblue'],['darkgreen'],['purple','seagreen'],['navy','darkgreen'],['navy','darkgreen']]
    labels = [['vi-cx'],['au-cx'],['cx-ch(av)','cx-ch(aa)'],['vi-ch(av)','vi-ch(aa)'],['au-ch(av)','au-ch(aa)'],]

    # this is how angles were calculated
    # angles = np.zeros((len(taskaspects),len(taskaspects),n_trajectory))
    # angles_highres[cxr,t,t_c] = (taskvariables,n_trajectory,n_trajectory) as angles

    
    angles_all = []
    angles_highres_all = []
    for n, dn in enumerate(datanames):
        # angles_aspects = []
        # angles_highres_aspects = []
        angles = pickle.load(open(cacheprefix+'subspaces/angles,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_highres = pickle.load(open(cacheprefix+'subspaces/angles,highres,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        # angles_aspects.append(angles)
        # angles_highres_aspects.append(angles_highres)
        angles_all.append(angles)
        angles_highres_all.append(angles_highres)

    angles_all = np.array(angles_all)

    times = np.arange(0,angles_all.shape[3])*T['dt'] + T['offsettime']



    fig,ax = plt.subplots(1,len(pairlists),figsize=(len(pairlists)*8,1*8))
    
    for plx,pairlist in enumerate(pairlists):
        axs = ax[plx]
        
        for px,pair in enumerate(pairlist):
            for n,dn in enumerate(datanames):
                x = neph.smooth(angles_all[n,pair[0],pair[1]],kernelwidth=6)
                axs.plot(times,x,lw=0.8,color=singlecolors[plx][px],alpha=0.2)

            m = angles_all[:,pair[0],pair[1]].mean(axis=0)
            e = angles_all[:,pair[0],pair[1]].std(axis=0)/np.sqrt(len(datanames))
            
            axs.plot(times,m,color=singlecolors[plx][px],lw=2,label=labels[plx][px])
            axs.fill_between(times,m-2*e,m+2*e,color=singlecolors[plx][px],alpha=0.2)


        
        axs.legend(frameon=False)
        axs.set_xlim(T['offsettime']+200*pq.ms,T['endtime']-200*pq.ms)
        figs.setxt(axs)
        axs.set_yticks([0,45,90,135,180])
        axs.set_ylim(0,180)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,90)

        axs.set_title(titles[plx])

    
    fig.suptitle('decision vector angles between stimuli, context, choice in the two contexts, mean and 2 s.e.m. of %d mice'%len(datanames))
    

    save = 0
    if save or globalsave:
        fig.savefig(resultpath+'tribe,anglesDBNV,viaucxchavchaa_%s_%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)


    

    
    
    return








def aggregatemice_subspace_contextrotationdecay(datanames):
    # show context representation shift speed with a linear fit to tested accuracies across timepoints


    trainingpoints = np.array([-1500, -100, 0, 1500, ])*pq.ms
    taskaspects = ['visual','audio','context','choice']



    # offstimdecoders_all  as   (mice,taskvariables,trainingpoint,n_trajectory)    
    # angles_highres_all   as (nice,taskvariables,n_trajectory,n_trajectory) as angles
    offstimdecoders_all = []
    angles_highres_all = []
    for n,dn in enumerate(datanames):

        offstimdecoders = []
        for cx in [0,2]:         # context only
           comparison = taskaspects[cx]
           offstimdecoder = []
           for trx,trainingpoint in enumerate(trainingpoints):
               offstimdecoder.append(pickle.load(open(cacheprefix+'continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,comparison,trx),'rb')))
           offstimdecoders.append(offstimdecoder)
        offstimdecoders_all.append(offstimdecoders)

        angles_highres = pickle.load(open(cacheprefix+'subspaces/angles,highres,alongDBNVs-VACC_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_highres_all.append(angles_highres)





    fig,ax = plt.subplots(2,4,figsize=(4*12,2*8) )

    colors = ['navy','mediumvioletred']
    skip = 20
    for trxax,trx in enumerate([0,1,2,3]):
        for cx in [0,1]:
        
            axs = ax[cx,trxax]
            accs = []
            for n,dn in enumerate(datanames):
                times = offstimdecoders_all[n][cx][trx][1].times[skip:-skip]       # T['stimstart_idx']:T['stimend_idx']
                accuracy = offstimdecoders_all[n][cx][trx][1][ : ,0 ]
                accuracy = neph.smooth(accuracy,kernelwidth=10)[skip:-skip]
                accs.append(accuracy)
                
                axs.plot(times,accuracy,color=colors[cx],lw=0.8,alpha=0.5,label=dn)
            accs = np.array(accs)
            
            m = accs.mean(axis=0)
            e = accs.std(axis=0)/np.sqrt(len(datanames))
            
            axs.plot(times,accs.mean(axis=0),color=colors[cx],lw=4)
            axs.fill_between(times,m-2*e,m+2*e,color=colors[cx],alpha=0.2)
            
            
            figs.setxt(axs)
            axs.set_xlim(-1500*pq.ms+skip*T['dt'],4500*pq.ms-skip*T['dt'])
            axs.set_ylim(0.3,1)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)
    
            
            # axs.legend(frameon=False);
            
            axs.set_title('trained at %d ms'%trainingpoints[trx])
            axs.set_ylabel(['visual','context'][cx]+'\nroc auc at each timepoint')
            axs.set_xlabel('time from trial onset [ms]')
                
    
    fig.suptitle('decoding test across timecourse, trained at 4 different timepoints, %d mice'%len(datanames))
    

    save = 0
    if save or globalsave:
        fig.savefig(resultpath+'tribe,crossdecoding,timepoints,vi,cx_%s_%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)




    
    
    return








def aggregatemice_calculatedecayconstant(datanames):
    # fit decoder accuracy cross-tests with a decaying exponential:
    # c * np.exp(m * x) + b
    # -m is decay constant: proportion of accuracy decaying to its 1/e ~ 1/3 
    # average lifetime tau = 1/m
    # halflife, T1/2 = ln (2)/m   = 
    
    doplot = 0
    recalculate = 1
    
    if recalculate:

        # this collects crossdecoder tests:
        crossdecoder_matrix_allaspects_all = []   #   (mice)(taskaspects)(trainingpoints)(testtimecourse,stats)
        print('Decay constant\nloading accuracies...')
        for n,dn in enumerate(datanames):
            crossdecoder_matrix_allaspects = pickle.load(open(cacheprefix+'continuous/crossdecodematrix-%s-%dms_%s.pck'%(continuous_method,T['dt'],dn),'rb'))
            crossdecoder_matrix_allaspects_all.append(crossdecoder_matrix_allaspects)
        crossdecoder_matrix_allaspects_all = np.array(crossdecoder_matrix_allaspects_all)
    
    
    
        
        print('fitting exponentials...')
        # forward decay rate;  modeled as     exp(-t/tau - c) + B    =   A * exp(-t/tau) + B,      A = exp(-c)
        # log(a[index:index+window])
        sh = crossdecoder_matrix_allaspects_all.shape
        # fig,axs = plt.subplots(sh[1],sh[0],figsize=(n_mice*8,sh[1]*8))
        fdr_window = 50
        fdr_skip_window=0
        fdr = np.zeros( (sh[0],sh[1],sh[2],3) )     # (mice,taskaspects,trainingpoints,{tau,multiplier,additive})
        bdr = np.zeros( (sh[0],sh[1],sh[2],3) )     # (mice,taskaspects,trainingpoints,{tau,multiplier,additive})
        for n in range(sh[0]):
            print(datanames[n])
            for cx in range(sh[1]):
                for ix in range(sh[2]-fdr_window-fdr_skip_window):
                    #fdr[n,cx,ix] = neph.decay(crossdecoder_matrix_allaspects_all[n,cx,ix,ix:ix+fdr_window,0])
                    m,c,b = neph.decay(crossdecoder_matrix_allaspects_all[n,cx,ix,ix+fdr_skip_window:ix+fdr_skip_window+fdr_window,0])
                    tau = -m
                    # A = np.exp(-c)
                    A = c
                    fdr[n,cx,ix,:] = (tau,A,b)
                    # print(ix,ix+fdr_window,\
                    #  crossdecoder_matrix_allaspects_all.shape,\
                    #  crossdecoder_matrix_allaspects_all.shape[2]-ix+fdr_window,\
                    #  -ix,-ix-fdr_window,\
                    #  crossdecoder_matrix_allaspects_all[n,cx,-ix,-1-ix-fdr_skip_window:-1-ix-fdr_skip_window-fdr_window:-1,0].shape\
                    #      )
                    m,c,b = neph.decay(crossdecoder_matrix_allaspects_all[n,cx,-ix,-1-ix-fdr_skip_window:-1-ix-fdr_skip_window-fdr_window:-1,0])
                    tau = -m
                    # A = np.exp(-c)
                    A = c           # do backward adjusting the time
                    bdr[n,cx,ix,:] = (tau,A,b)
    
    #            if np.any(fdr[n,cx,:,0]<0) or (fdr[n,cx,:,0]>2000):  fdr[n,cx,:,0] = np.nan
    
        # for n in range(sh[0]):
        #     for cx in range(sh[1]):
        #         axs[cx,n].plot(fdr[n,cx,:,0])
        #         if cx==0: axs[cx,n].set_title(datanames[n])
        #         if n==0: axs[cx,n].set_ylabel(['visual','audio','context','choice'][cx])
        
        # fig.suptitle('charateristic persistence time, tau:  exp(-t/tau)')
        # fig.suptitle('charateristic wavelength, beta:  exp(-$ \\beta \\cdot $ t)')
    
    
        pickle.dump((fdr,bdr),open(cacheprefix+'continuous/crossdecodematrix,decayconstant-%s-%dms.pck'%(continuous_method,T['dt']),'wb'))

    else:
        
        fdr,bdr = pickle.load(open(cacheprefix+'continuous/crossdecodematrix,decayconstant-%s-%dms.pck'%(continuous_method,T['dt']),'rb'))
        

    if doplot:
        mice = ['DT014','DT017','DT019','ME103']
        
        fig,ax = plt.subplots(3,len(mice),figsize=(len(mice)*8,3*8))
        for nx,dn in enumerate(mice):
            n = datanames.index(dn)
            times = np.arange(fdr.shape[2])*10-1500
            for coeffidx in [0,1,2]:
        
                axs = ax[coeffidx,nx]
                
                axs.plot(times,fdr[n,0,:,coeffidx],lw=2,color='navy')
                axs.plot(times,fdr[n,2,:,coeffidx],lw=2,color='mediumvioletred')
                # axs.set_xticks([0,3000,600])
                # axs.set_xticklabels(np.array([stx,stx+300,stx+600])-150)
                if coeffidx==0:
                    # axs.set_ylim(0,0.05)
                    axs.set_title(dn)
                else:
                    axs.set_ylim(-3,3)
                    figs.plottoaxis_chancelevel(axs,0)
                    
                figs.plottoaxis_stimulusoverlay(axs,T)
                


        
    return







def aggregatemice_lowdimensionnullspacecontext(datanames):


    doplot = 1 or globaldoplot

    datanames = ['ME110','ME113','DT009','DT014','DT017','DT018','DT019','DT021','DT022','MT020_2']
    datanames_lowdimcontext = ['DT014','DT022','MT020_2']


    # indices_lowdimcontext = [  datanames.index(lowdim) for lowdim in datanames_lowdimcontext   ]



    # load full- and reduced-space crosstime decoder accuracies

    crossdecoders_all = []

    for n,dn in enumerate(datanames):     # change to datanames!!!
        
        crossdecoders = []
        
        exptyp = ',allcontexts'
        crossdecoder_matrix_allaspects = pickle.load(open(cacheprefix+'continuous/crossdecodematrix%s,contextcond,cv-%s-%dms_%s.pck'%(exptyp, continuous_method,T['dt'],dn),'rb'))
        crossdecoders.append(crossdecoder_matrix_allaspects[2])           # 2 is context


        if dn in datanames_lowdimcontext:
            for n_nsrd in [1,3]:  # nullspace recursion depth       1, 3, 10
                exptyp = ',recurrentnullspacecontext%d'%n_nsrd
                crossdecoder_matrix_allaspects = pickle.load(open(cacheprefix+'continuous/crossdecodematrix%s,contextcond,cv-%s-%dms_%s.pck'%(exptyp, continuous_method,T['dt'],dn),'rb'))
                crossdecoders.append(crossdecoder_matrix_allaspects[0])        # this question has only a single reply - unlike ',allcontexts' above

        
        crossdecoders_all.append(crossdecoders)


    print(len(crossdecoders_all), crossdecoders_all[0][0].shape)

    # baseline random decoder labels        
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    if doplot:








        if 0: # show explicitly side bby side the non-blocked animals with original, subtracting 1, and 3 dims of context DVs.
            fig,ax = plt.subplots(3,len(datanames_lowdimcontext),figsize=(len(datanames_lowdimcontext)*8,3*8))


            for nx,dn in enumerate(datanames_lowdimcontext):
                for sx,s in enumerate([0,1,3]):
                    axs = ax[sx,nx]

                    cmap = figs.getcorrcolormap('correlation')
                    cf = axs.pcolormesh(crossdecoders_all[nx][sx][:,:,0],vmin=0,vmax=1,cmap=cmap)
                    
                    axs.set_aspect('equal')
                    
                    ticks=[150,450]
                    ticklabels=['0 ms','3000 ms']
                    axs.set_xticks(ticks)
                    axs.set_xticklabels(ticklabels)                
                    axs.set_yticks(ticks)
                    axs.set_yticklabels(ticklabels,rotation=90)                
                    
                    if sx==0: axs.set_title(dn)
                    
                    axs.set_xlabel('test timecourse')

                    axs.set_ylabel(['','%d dims subtracted\n'%s][nx==0]+'train timecourse')

                    # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                    fig.colorbar(cf,ax=axs,ticks=np.arange(0,1.1,0.25))









        if 1:     # publication story
            # A example animal 2 non-blocked, B context nullspace recurrent timecourse, example animal 2, C example animal 2, with the first context DV subtracted, 
            # D context nullspace recurrent deletion for all animals, mean acc remaining, E 1st subspace drop in context vs. drop between on and pre original

            fig,ax = plt.subplots(2,3,figsize=(3*8,2*8))


            panels = ['I','K']
            dn = 'MT020_2'
            n = datanames.index(dn)
            ssids = [0,1]        # subspace code
            

            for i,si in enumerate(ssids):

                mapres = np.arange(0.5,1.,0.01)
                cmap = plt.cm.PuRd(mapres)
                cmap[mapres<chances[dn],:] = np.array([0.94,0.94,0.94,1.0])
                cmap = clrs.ListedColormap(cmap, name='PuRdE', N=cmap.shape[0])

                axs = ax[0,2*i]
                
                x = crossdecoders_all[n][si][:,:,0]
                cf = axs.pcolormesh(x,vmin=0.5,vmax=1,cmap=cmap)
                
                axs.set_aspect('equal')
                
                ticks=[150,450]
                ticklabels=['0 ms','3000 ms']
                axs.set_xticks(ticks)
                axs.set_xticklabels(ticklabels)                
                axs.set_yticks(ticks)
                axs.set_yticklabels(ticklabels,rotation=90)                
                
                axs.set_title('%s accuracy'%['full space','nullspace of first context'][i],fontsize='small')
                
                axs.set_xlabel('test timecourse')
                axs.set_ylabel('train timecourse')
            
                # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                fig.colorbar(cf,ax=axs,ticks=np.arange(0.5,1.1,0.25))
            
            
                figs.labelaxis(axs,panels[i])


            



            panel = 'J'
            axs = ax[0,1]
            comparison = 'context'
            colors = ['darkorange','darkred']

            acrossdecoder_nullspaces,ranks = pickle.load(open(cacheprefix+'subspaces/nullspace,recurrent,decodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
            # times = acrossdecoder_nullspace[1].analogsignals[0].times

            n_neurons = len(ranks)+1
            for px in range(n_neurons):
                if px==0:
                    acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
                else:
                    acrossdecoder = acrossdecoder_nullspaces[px-1]
                    colors[1]= np.array([1,0,0.5])*(px/n_neurons/2)+np.array([0.25,0,0.12])
                
                a = neph.smooth(acrossdecoder[1][:,0],mode='same')
                axs.plot(acrossdecoder[1].times,a,color=colors[px>0],lw=[2,1][px>0])

            figs.setxt(axs)
            axs.set_xlim(-1300,4200)
            axs.set_ylim(0.45,1.05)
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,chances[dn])

            axs.set_title('recurrent context nullspaces',fontsize='small')
            axs.set_ylabel('context accuracy')

            axs.spines['right'].set_visible(False)
            axs.spines['top'].set_visible(False)
            figs.labelaxis(axs,panel)







            panel = 'L'
            axs = ax[1,0]
            comparison = 'context'
            colors = ['darkorange','darkred']

            for n,dn in enumerate(datanames):

                acrossdecoder_nullspaces,ranks = pickle.load(open(cacheprefix+'subspaces/nullspace,recurrent,decodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
                # times = acrossdecoder_nullspace[1].analogsignals[0].times

                n_neurons = len(ranks)+1
                x = []
                for px in range(n_neurons):
                    if px==0:
                        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
                    else:
                        acrossdecoder = acrossdecoder_nullspaces[px-1]
                        colors[1]= np.array([1,0,0.5])*(px/n_neurons/2)+np.array([0.25,0,0.12])
                    
                    x.append( acrossdecoder[1][:,0].mean(axis=0) )

                    axs.plot(px,x[-1],'o',markerfacecolor=colors[px>0], markeredgecolor=colors[px>0],markersize=3)
                
                axs.plot(np.arange(n_neurons),x,color='mediumvioletred',alpha=0.3,lw=1)
            

            axs.set_ylim(0.45,0.8)
            axs.set_yticks([0.5,0.8])
            m_c = np.array(list(chances.values())).mean()
            e_c = np.array(list(chances.values())).std()/np.sqrt(len(datanames))
            figs.plottoaxis_chancelevel(axs,m_c+e_c)

            axs.set_xlabel('subtracted dimensions')
            axs.set_ylabel('context accuracy')


            axs.spines['right'].set_visible(False)
            axs.spines['top'].set_visible(False)
            figs.labelaxis(axs,panel)










            # drop in mean first context decoder accuracy   vs.   drop between mean pre and on accuracy
            panel = 'M'
            axs = ax[1,1]
            comparison = 'context'
            color = 'darkmagenta'



            x = []
            y = []
            neurons = []
            for n,dn in enumerate(datanames):

                # get the mean first context accuracy drop:
                
                acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
                acrossdecoder_nullspaces,ranks = pickle.load(open(cacheprefix+'subspaces/nullspace,recurrent,decodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
                # times = acrossdecoder_nullspace[1].analogsignals[0].times

                dx = acrossdecoder[1].magnitude[:,0].mean(axis=0) - acrossdecoder_nullspaces[0][1].magnitude[:,0].mean(axis=0)

                x.append(dx)


                # get average pre vs on accuracy drop
                skipoffdiag = 10
                C = crossdecoders_all[n][0][:,:,0]
                c_pre = C[:T['stimstart_idx'],:T['stimstart_idx']]
                c_on = C[T['stimstart_idx']:T['stimend_idx'],T['stimstart_idx']:T['stimend_idx']]
                for c in (c_pre,c_on):
                    sz = c.shape[0]
                    for imax in range(skipoffdiag):
                        for i in range(sz-imax):
                            c[i+imax,i] = np.nan
                            c[i,i+imax] = np.nan
                c_cross = np.hstack( [ C[T['stimstart_idx']:T['stimend_idx'],:T['stimstart_idx']].ravel(), \
                                       C[:T['stimstart_idx'], T['stimstart_idx']:T['stimend_idx'] ].ravel() ] )

                
                dy = np.nanmean(c_pre)/2 + np.nanmean(c_on)/2  - np.mean(c_cross)
                
                y.append( dy )
            
                n_neurons = len(ranks)+1
                neurons.append(n_neurons)

                if n==0: print('name, neurons, nullspace 1, offdiag   > dx, dy ')
                print(dn, n_neurons, acrossdecoder_nullspaces[0][1].magnitude[:,0].mean(axis=0), np.nanmean(c_pre)/2 + np.nanmean(c_on)/2, chances[dn], '  >', dx, dy)
            

            axs.plot(x,y,'o',markerfacecolor=color, markeredgecolor=color)
            
            for nx,n in enumerate(neurons):
                axs.text(x[nx]-0.001,y[nx]+0.001,'%s %d'%(datanames[nx],n),fontsize=9)


            removes = [7,2]
            for k in removes:
                x.pop(k)
                y.pop(k)
                datanames.pop(k)
            print(datanames)

            x = np.array(x)
            y = np.array(y)
            l = sp.stats.linregress(x,y)
            print(k,l)
            line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],color='grey',linewidth=0.5)
            




            axs.set_xlabel('accuracy loss in first context nullspace')
            axs.set_ylabel('accuracy difference\nacross preon    -   within preon')


            axs.spines['right'].set_visible(False)
            axs.spines['top'].set_visible(False)
            figs.labelaxis(axs,panel)





        ax[1,2].remove()


        fig.tight_layout()




        save = 1 or globalsave
        if save:
            fig.savefig('../../../publish/journals/journal2020spring/figures/'+'Fig4_addition_lowdimcontextpersistent'+ext)




    return

















def aggregatemice_attendignore(datanames, examine='attendignore'):
    # compare attended and ignored condition stimulus representations

    timedepth = 1000*pq.ms      # how deep we want to average over
    timedepth_idx = int(timedepth/T['dt'])
    
    n_mice = len(datanames)    
        
    taskaspects = ['attend visual','ignore visual','attend audio','ignore audio']
    taskaspectsshort = ['att.vis.','ign.vis.','att.au.','ign.au.']
    logical = 'and'
    stimulusspecificity = 'all'
    layeredfilename = 'all'



    
    m = np.zeros((n_mice,len(taskaspects),2))
    s = np.zeros((n_mice,len(taskaspects),2))
    tr = np.zeros((n_mice,len(taskaspects),2))
    se = np.zeros((n_mice,len(taskaspects),2))
    r = np.zeros((n_mice,len(taskaspects)))
    
    roc_m = np.zeros((n_mice,len(taskaspects),T['stimend_idx']-T['stimstart_idx']))
    roc_s = np.zeros((n_mice,len(taskaspects),T['stimend_idx']-T['stimstart_idx']))
    acc_d_a = np.zeros((n_mice,int(len(taskaspects)/2)))
    acc_d_i = np.zeros((n_mice,int(len(taskaspects)/2)))
    
    for n,dn in enumerate(datanames):
        blv,bla = preprocess.getorderattended(dn)
        comparisongroups  = [ \
                                [ [ [blv[1]], [45],        [] ],  [ [blv[1]],[135],        [] ] ],\
                                [ [ [bla[1]], [45],        [] ],  [ [bla[1]],[135],        [] ] ],\
                                [ [ [bla[1]],   [],    [5000] ],  [ [bla[1]],[],      [10000] ] ],\
                                [ [ [blv[1]],   [],    [5000] ],  [ [blv[1]],[],      [10000] ] ]\
                            ]


        block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=True,recalculate=False)        
        
        
        n_timecourse, n_neurons = block.segments[0].analogsignals[0].shape
        # collect decision boundary normal vectors
    
        c_db = []          # vector components (coefficients) of the decision normal in the activity space
        onstimactivity = []     # neural activity
        roc = []             # roc auc accuracies
        trajectories = []
        for cx,comparison in enumerate(taskaspects):

            stimulusIDgroups = comparisongroups[cx]
            # collect neural responses
            if not comparison=='choice': # visual, audio, context:
                acrossresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups,logical)
            # acrossresponses is [class][trials][trajectory,neurons]
            aux = []
            for clx in [0,1]:    # classes (go and nogo signal)
                aux.append(np.array( acrossresponses[clx] )[:,T['stimstart_idx']:T['stimstart_idx']+timedepth_idx,:].mean(axis=1).squeeze())
                # aux is (trials,neurons)        # print(clx,aux[clx].shape)
            onstimactivity.append(aux)
            # onstimactivity is [comparisongroup][class][trials,neurons])
    

            # collect decoder coefficients
            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            wx = int((len(acrossdecoder)-7)/n_neurons)
            c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )

            # and collect roc auc accuracies
            roc.append(     neph.get_auctimenormalized(acrossdecoder[1][:,0],T['stimstarttime'],T['stimstarttime']+timedepth)            )
            
            trajectories.append( acrossdecoder[1][T['stimstart_idx']:T['stimend_idx'],:]  )
#            print(len(acrossdecoder),len(acrossdecoder[1]),len(acrossdecoder[1][0]),type(acrossdecoder[0]))
            

            
        c_db = np.array(c_db) # [comparisongroup,neurons,trajectory,stats]

        c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimstart_idx']+timedepth_idx,0].mean(axis=2)            # this is a comparison group by neuron by   stats matrix, we take only the mean value ,0]
        norms_c_db = np.linalg.norm(c_db_means,axis=1,keepdims=True)

#        print(c_db_means.shape)
#        print(norms_c_db.shape)
#        print(len(onstimactivity), len(onstimactivity[0]), onstimactivity[0][0].shape)

        

        for cx,comparison in enumerate(taskaspects):
            for clx in [0,1]:   # go over the classes
                 y = np.dot(   onstimactivity[cx][clx], c_db_means[cx] / norms_c_db[cx] )
                 m[n,cx,clx] = y.mean()
                 s[n,cx,clx] = y.std()
                 tr[n,cx,clx] = len(y)
            r[n,cx] = roc[cx][0]
            roc_m[n,cx,:] = trajectories[cx][:,0].squeeze() # mean accuracy dynamics
            roc_s[n,cx,:] = trajectories[cx][:,2].squeeze() # s.e.m.
#        print(m)
#        print(s)
#        print(tr)
    se = 2*s/np.sqrt(tr)
#        print(se)
    rse = 2*  np.sqrt( r*(1-r) / tr.sum(axis=2) )
    
    for cxcat in [0,1]:
        for n in range(n_mice):
            acc_d_a[n,cxcat], acc_d_i[n,cxcat] = neph.get_aucpositivenegative(\
                           roc_m[n,cxcat*2,:],roc_m[n,cxcat*2+1,:],\
                           roc_s[n,cxcat*2,:],roc_s[n,cxcat*2+1,:])



    # plots

    colorlist = plt.cm.viridis( np.linspace(0, 0.8, n_mice) )

    if not publish: fig,ax = plt.subplots(1,5,figsize=(45,9)); z=4
    else: fig,ax = plt.subplots(1,1,figsize=(14,7)); z=0
    
    b = 80
    x0 = np.array([0,1,3,4])
    
    if not publish:
        axs = ax[0]
        for n in range(n_mice):
            for ccatx in [0,2]:
                x = x0[ccatx:ccatx+2]
                y = np.abs(m[n,ccatx:ccatx+2,0]-m[n,ccatx:ccatx+2,1])
                s = se[n,ccatx:ccatx+2,:].sum(axis=1)
                axs.plot(x+n/b, y, 'o-', color=colorlist[n] )
                axs.errorbar(x+n/b, y=y, yerr=s, color=colorlist[n] )
            axs.set_xticks(x0)
            axs.set_xticklabels(taskaspectsshort)
        axs.legend(datanames)
        axs.set_title('class center of mass difference')
    
    
        axs = ax[1]
        for n in range(n_mice):
            for ccatx in [0,2]:
                x = x0[ccatx:ccatx+2]
                y = np.sqrt(2)*sp.stats.norm.cdf(r[n,ccatx:ccatx+2])
                s = np.sqrt(2)*sp.stats.norm.cdf(r[n,ccatx:ccatx+2]+rse[n,ccatx:ccatx+2]) - y
                axs.plot(x+n/b, y, 'o-', color=colorlist[n] )
                axs.errorbar(x+n/b, y=y, yerr=s, color=colorlist[n] )
            axs.set_xticks(x0)
            axs.set_xticklabels(taskaspectsshort)
    #    axs.legend(datanames)
        axs.set_title('d`')

    
        axs = ax[2]
        for n in range(n_mice):
            for ccatx in [0,2]:
                x = x0[ccatx:ccatx+2]
                y = r[n,ccatx:ccatx+2]
                s = rse[n,ccatx:ccatx+2]
                axs.plot(x+n/b, y, 'o-', color=colorlist[n] )
                axs.errorbar(x+n/b, y=y, yerr=s, color=colorlist[n] )
            axs.set_xticks(x0)
            axs.set_xticklabels(taskaspectsshort)
        axs.set_ylim(0.435,1)
        figs.plottoaxis_chancelevel(axs,0.5)
    #    axs.legend(datanames)
        axs.set_title('accuracy roc auc')
    
    
        axs = ax[3]
        for ccatx in [0,2]:
            x = x0[ccatx]
            da = (r[:,ccatx] - rse[:,ccatx]) - (r[:,ccatx+1] + rse[:,ccatx+1])
            di = (r[:,ccatx+1] - rse[:,ccatx+1]) - (r[:,ccatx] + rse[:,ccatx])
            ca = np.sum(da>0) / n_mice
            ci = np.sum(di>0) / n_mice
            c0 = 1 - ca - ci
            print(ca,c0,ci)
            axs.bar([x-0.5,x,x+0.5], [ca,c0,ci], width=0.3, color=['navy','darkgreen'][int(ccatx/2)])
        axs.set_ylim(0,1)
        axs.set_xticks([-0.5,0,0.5,2.5,3,3.5])
        axs.set_xticklabels(['att.','=\nvisual','ign.','att.','=\naudio','ign.'] )
        axs.set_title('change categories')


    
    if not publish: axs = ax[z]
    else: axs = ax
    
    for ccatx in [0,1]:
        x = ccatx*2
        
        for n in range(n_mice):
            axs.plot(x+n/b-0.25,acc_d_a[n,ccatx]+acc_d_i[n,ccatx],'o',color=colorlist[n])      #['navy','darkgreen'][ccatx])
        
        ml = (acc_d_a[:,ccatx]+acc_d_i[:,ccatx]).mean()
        cl = (acc_d_a[:,ccatx]+acc_d_i[:,ccatx]).std()/np.sqrt(n_mice)*2
        axs.boxplot(x=acc_d_a[:,ccatx]+acc_d_i[:,ccatx],positions=[x+0.25],notch=True,whis=[0,95],\
                    usermedians=[ml],conf_intervals=[[ml-cl,ml+cl]])
        
    axs.set_xlim(-1,3)
    figs.plottoaxis_chancelevel(axs,0.0)
    axs.set_xticks([0,2])
    axs.set_xticklabels(['visual','audio'])
    axs.set_yticks([0,0.1])
    axs.set_yticklabels(['0','0.1'])
    axs.set_title('accuracy timecourse significant difference: attended to  -  ignored',fontsize=20)
    axs.set_ylabel('accuracy [RoC AUC]')
    


        
    if not publish: fig.suptitle('attend - ignore')
    save = 0
    if save:
        if not publish:
            fig.savefig(resultpath+'micetribe,attendignore,cm+dp+roc+p_%s-%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
        else:
            fig.savefig(resultpath+'2-decoder,attendignore'+ext)













def aggregatemice_subspaces_PCA(datanames):
    # calculate PCA with the covariance along trials with averaged activity in spontaneous and evoked activity

    #    n_neurons = block.segments[0].analogsignals[0].shape[1]


    pcadimlist = [1,2,3,4,5,7,10,15,30,60]
    # pcadimlist = [1,2,3,4,5,7,10,15,20]
    
    aggregateuntil = 5   # for comparisons, that exists for all mice

    pcadimlist = pcadimlist[:aggregateuntil]

    taskaspects =     taskaspects = ['visual','audio','context','choice']
    colors = ['navy','darkgreen','mediumvioletred','orange']
    

    tribe_ea_all = []
    tribe_sa_all = []
    tribe_sa_stat = []
    tribe_ea_stat = []
    

    # pcadimlist.append(n_neurons)
    #         for r in pcadimlist: # go for all pca dims and the last one as the actual number of neurons
    #             acrossdecoders.append(acrossdecoder)       # collect for plot
    #             acrossrandomdecoders.append(acrossrandomdecoder)    # collect for plot


    
    # fig, ax = plt.subplots(4,4,figsize=(24,24))
    
    for frx,pcamean in enumerate(['sa','ea']):
    # where to run the PCA on

        cux_sa = []
        cux_ea = []
        for cx,comparison in enumerate(taskaspects):
            aux_sa = []
            aux_ea = []
            for rx,r in enumerate(pcadimlist):
                tribe_sa = []
                tribe_ea = []
                # tribe_sa.append([])
                for n,dn in enumerate(datanames):
                    print(n,dn)
                    if dn=='DT030': continue
        
                    if not rx+1==len(pcadimlist):
                        acrossdecoder,acrossrandomdecoder = pickle.load(open(cacheprefix+'pca/responsedecodes+O,pcadim%d,%s-%s-%s-%s.pck'%(r,pcamean,dn,continuous_method,comparison),'rb'))


                        fn = cacheprefix+'pca/responsedecodes,pcadim%d,%s-%s-%s-%s.pck'%(r,pcamean,dn,continuous_method,taskaspects[cx])
                        if os.path.isfile(fn):
                            acrossdecoder = pickle.load( open(fn,'rb') )
        
                            tribe_sa.append( acrossdecoder[1][T['start_idx']:T['stimstart_idx'],0].mean() )
                            tribe_ea.append( acrossdecoder[1][T['stimstart_idx']:T['stimend_idx'],0].mean() )
                        mk = 'o'
                        al = 0.2
                    else:    # use the best without PCA
                        fn = cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all')
                        acrossdecoder = pickle.load( open(fn,'rb') )
    
                        tribe_sa.append( acrossdecoder[1][T['start_idx']:T['stimstart_idx'],0].mean() )
                        tribe_ea.append( acrossdecoder[1][T['stimstart_idx']:T['stimend_idx'],0].mean() )
                        mk = '+'
                        al = 0.8
                        
                            
                
                tribe_sa = np.array(tribe_sa)
                tribe_ea = np.array(tribe_ea)
                
                tribe_sa_stat.append(  [tribe_sa.mean(), tribe_sa.std(), len(tribe_sa) ] )
                tribe_ea_stat.append(  [tribe_ea.mean(), tribe_ea.std(), len(tribe_ea) ] )
                
                # axs = ax[cx,frx*2]
                # axs.plot( r, tribe_sa[None,:], mk, color=colors[cx], alpha=al )
                # axs.set_xlim(0,pcadimlist[-1]+1)
                # axs.set_ylim(0.45,1.01)
                # figs.plottoaxis_chancelevel(axs,0.5)
                # if frx==0: axs.set_ylabel(taskaspects[cx])
                # if cx==0: axs.set_title('pca %s -> sa'%pcamean)
                
                # axs = ax[cx,frx*2+1]
                # axs.plot( r, tribe_ea[None,:], mk, color=colors[cx], alpha=al  )
                # axs.set_xlim(0,pcadimlist[-1]+1)
                # axs.set_ylim(0.45,1.01)
                # figs.plottoaxis_chancelevel(axs,0.5)
                # if cx==0: axs.set_title('pca %s -> ea'%pcamean)
    
                if rx<aggregateuntil:
                    aux_sa.append(tribe_sa)
                    aux_ea.append(tribe_ea)
        
            cux_sa.append(aux_sa)
            cux_ea.append(aux_ea)
        tribe_sa_all.append(cux_sa) # taskaspects x pcadimlist x mice x ->sa,->ea
        tribe_ea_all.append(cux_ea)

                
      
        
    if 1:
        fn = cacheprefix+'pca/responsedecodes,pca-%s.pck'%(continuous_method)
        pickle.dump((tribe_sa_all,tribe_ea_all), open(fn,'wb') )
        
        



    if 0:
        print(np.array(tribe_sa_stat).shape)
        tribe_sa_stat = np.array(tribe_sa_stat).reshape((2,4,len(pcadimlist),3))
        tribe_ea_stat = np.array(tribe_ea_stat).reshape((2,4,len(pcadimlist),3))
        tribe_sa_stat_o = tribe_sa_stat[:,:,-1,:]
        tribe_ea_stat_o = tribe_ea_stat[:,:,-1,:]
        tribe_sa_stat = tribe_sa_stat[:,:,:-1,:]
        tribe_ea_stat = tribe_ea_stat[:,:,:-1,:]
        pcadimlist = pcadimlist[:-1]

        fig, ax = plt.subplots(4,4,figsize=(24,24))
        
        for frx,pcamean in enumerate(['sa','ea']):
        # where to run the PCA on:
    
            for cx,comparison in enumerate(taskaspects):
                    
                axs = ax[cx,frx*2]
                axs.plot( pcadimlist, tribe_sa_stat[frx,cx,:,0], linewidth=2,color=colors[cx], alpha=1.0 )
                axs.fill_between( pcadimlist, tribe_sa_stat[frx,cx,:,0]-tribe_sa_stat[frx,cx,:,1]*2/np.sqrt(tribe_sa_stat[frx,cx,:,2]),\
                                              tribe_sa_stat[frx,cx,:,0]+tribe_sa_stat[frx,cx,:,1]*2/np.sqrt(tribe_sa_stat[frx,cx,:,2]), color=colors[cx], alpha=0.6 )
                axs.fill_between( pcadimlist, tribe_sa_stat[frx,cx,:,0]-tribe_sa_stat[frx,cx,:,1]*2,\
                                              tribe_sa_stat[frx,cx,:,0]+tribe_sa_stat[frx,cx,:,1]*2, color=colors[cx], alpha=0.2 )
                axs.set_xlim(0,pcadimlist[-1]+1)
                axs.set_ylim(0.45,1.01)
                figs.plottoaxis_chancelevel(axs,0.5)
                if frx==0: axs.set_ylabel(taskaspects[cx])
                if cx==0: axs.set_title('pca %s -> sa'%pcamean)
                
                axs = ax[cx,frx*2+1]
                axs.plot( pcadimlist, tribe_ea_stat[frx,cx,:,0], linewidth=2,color=colors[cx], alpha=1.0 )
                axs.fill_between( pcadimlist, tribe_ea_stat[frx,cx,:,0]-tribe_ea_stat[frx,cx,:,1]*2/np.sqrt(tribe_ea_stat[frx,cx,:,2]),\
                                              tribe_ea_stat[frx,cx,:,0]+tribe_ea_stat[frx,cx,:,1]*2/np.sqrt(tribe_ea_stat[frx,cx,:,2]), color=colors[cx], alpha=0.6 )
                axs.fill_between( pcadimlist, tribe_ea_stat[frx,cx,:,0]-tribe_ea_stat[frx,cx,:,1]*2,\
                                              tribe_ea_stat[frx,cx,:,0]+tribe_ea_stat[frx,cx,:,1]*2, color=colors[cx], alpha=0.2 )
                axs.set_xlim(0,pcadimlist[-1]+1)
                axs.set_ylim(0.45,1.01)
                figs.plottoaxis_chancelevel(axs,0.5)
                if cx==0: axs.set_title('pca %s -> ea'%pcamean)
        
            
            if cx==3: axs.set_xlabel('# pca dims')
    
    
    
        fig.suptitle('Accuracies of decoders using only the best principal components\npca performed on the mean FR in spontaneous and evoked activity')
        save = 0
        if save:
            fig.savefig(resultpath+'pcareduceddecoders_all,distrib-%s-%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
    


    if 0:
        
        tribe_sa = np.array(tribe_sa_all)
        tribe_ea = np.array(tribe_ea_all)
        print(tribe_sa.shape)
        max_pc = 4

        fig, ax = plt.subplots(max_pc,4,figsize=(8*5,8*max_pc))
        
        for pcx in range(max_pc):
            for sx,place in enumerate(['sa','ea']):
        
                axs = ax[pcx,2*sx]
                axs.boxplot(x=tribe_sa[sx,:,pcx,:].T,notch=True)
                if pcx==0: axs.set_title('pca %s -> sa'%place)
                axs.set_xticklabels(taskaspects)
                if sx==0:
                    ax[pcx,sx].set_ylabel('# pc %d'%(pcx+1))
                axs.set_xlim(0.5,4.5)
                axs.set_ylim(0.475,1.01)
                figs.plottoaxis_chancelevel(axs,ch=0.5)
                
            
                axs = ax[pcx,2*sx+1]
                axs.boxplot(x=tribe_ea[sx,:,pcx,:].T,notch=True)
                if pcx==0: axs.set_title('pca %s -> ea'%place)
                axs.set_xticklabels(taskaspects)
                if sx==0:
                    ax[pcx,sx].set_ylabel('# pc %d'%(pcx+1))
                axs.set_xlim(0.5,4.5)
                axs.set_ylim(0.475,1.01)
                figs.plottoaxis_chancelevel(axs,ch=0.5)

        fig.suptitle('Accuracies of decoders using only the best principal components\npca performed on the mean FR in spontaneous and evoked activity\n%d mice'%(len(datanames)-1))
        save = 0
        if save:
            fig.savefig(resultpath+'pcareduceddecoders_1-4pc-%s-%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)







def aggregatemice_subspaces_PCA_saveonly(datanames):
    
    tribedata_sa = []
    tribedata_ea = []
    for n,dn in enumerate(datanames):
        block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=False)
        

        # load these for all contexts:       pcadimlist,p_sa,p_sa_r,p_ea,p_ea_r 
        
        data = subspaces_PCA(dn,block,returnall=True,preon=0)
        tribedata_sa.append(data)
        data = subspaces_PCA(dn,block,returnall=True,preon=1)
        tribedata_ea.append(data)

    
    pickle.dump(tribedata_sa, open(cacheprefix+'pca/tribedata-pca,sa,rand+O','wb') )
    pickle.dump(tribedata_ea, open(cacheprefix+'pca/tribedata-pca,ea,rand+O','wb') )










def aggregatemice_subspaces_PCA_explainedvariance(datanames):
    # quantify explained variance of pcs
    # plot the principal components' explained variance and their ordered cumulative ratio
    # establish a reference for total number of neurons measured
    
    
    n_subs = 10
    pcagroups = [      [[2,4],[],[]]   ]         # use the multimodal blocks only as activities
    n_mice = len(datanames)
    
    titles = ['total var.','expl. var.','cumulative expl. var.','expl. var. ratio','cumulative expl. var. ratio']
    n_pc = 4
    
    n_cov = 5
    singularvalues = np.zeros((n_mice*n_subs,n_pc,n_cov+1))

    for n,dn in enumerate(datanames):
        
        block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=True,recalculate=False)
        responses_full = preprocess.collect_stimulusspecificresponses(block,pcagroups,'and')
        
        for s in range(n_subs):
            
            # repeatedly take smaller than all number of neurons randomly
            n_neurons = responses_full[0][0].shape[1]
            print(n_neurons)
            n_subsample_neurons = int(np.floor(np.random.rand()*(n_neurons-n_pc)))+n_pc
            mask = np.random.permutation(n_neurons)[:n_subsample_neurons]
            responses =  np.array( responses_full[0] )[:,:,mask]
            
            
            # n_trials = len(responses)
            n_trials,n_trajectory,n_neurons = responses.shape
            X = responses[:,T['start_idx']:T['end_idx'],:].reshape(-1,n_neurons)
            Y,E,V = nedi.pca(X,n_components=n_neurons)
    
            print('PCA eigenvalues',E)
            print('neurons %d, PCs %d,'%(n_neurons,2),'    PCA transform shape',V.shape,'    transformed shape',Y.shape)
    
            
            
            totalvariance = np.sum(E)
            explainedvariance = E[:n_pc]
            cumulativeexplainedvariance = np.cumsum(E[:n_pc])
            explainedvarianceratio = explainedvariance / totalvariance
            cumulativeexplainedvarianceratio = np.cumsum(explainedvarianceratio)
            
            print(dn,  explainedvariance, explainedvarianceratio)
            print('     ',cumulativeexplainedvariance, cumulativeexplainedvarianceratio)
            
            
            singularvalues[n*n_subs+s,:,0] = np.repeat(totalvariance,n_pc)
            singularvalues[n*n_subs+s,:,1] = explainedvariance[:n_pc]
            singularvalues[n*n_subs+s,:,2] = cumulativeexplainedvariance[:n_pc]
            singularvalues[n*n_subs+s,:,3] = explainedvarianceratio[:n_pc]
            singularvalues[n*n_subs+s,:,4] = cumulativeexplainedvarianceratio[:n_pc]
            singularvalues[n*n_subs+s,:2,5] = n_trials,n_neurons

    





    # plot figure
    fig,ax = plt.subplots(1+n_pc,n_cov,figsize=(8*n_cov,7*(1+n_pc)))
    
    
    for k in range(n_cov):
        axs = ax[0,k]
        for n in range(n_mice*n_subs):
            cs = n/n_mice/n_subs
            
            if n_subs==1: label = datanames[n%n_subs]+' N %d'%singularvalues[n,1,-1]
            else: label=None
        
            axs.plot(np.arange(1,n_pc+1),singularvalues[n,:,k].T,color=[0.5*cs,1*cs,0.5+0.5*cs],\
                     label=label)
        
        if k==0 and n_subs==1: axs.legend(frameon=False)

        axs.set_xlabel('# PC')
        axs.set_title(titles[k])
        axs.set_ylabel('variance [SU]')




        
        # plot
    
    
    for pc_idx,pc in enumerate(np.arange(n_pc)+1):
        for k in range(n_cov):       # number of neurons vs. variances:
            axs = ax[1+pc_idx,k]

        
            axs.plot(singularvalues[:,1,-1],singularvalues[:,pc-1,k],'o',color='darkgreen')
            if n_subs==1:
                for n in range(n_mice*n_subs): 
                    axs.text(singularvalues[n,1,-1]-2,singularvalues[n,pc-1,k]-0.04,datanames[n],color='grey',alpha=0.5)

    

            # try to fit

            def exp_param(x, t, y):
                return        x[2] + x[1] * t**x[0]    -    y
            # def exp_param(x, t, y):
            #     return        t**x[0]    -    y

            x0 = np.array([-0.5,1,0])      # initialize the parameters
            # perform fit, with how neuron number would determine variances (total, exp, cumexp)
            res_robust = sp.optimize.least_squares(exp_param, x0, loss='soft_l1', f_scale=0.01,\
                                                   args=(singularvalues[:,1,-1],singularvalues[:,pc-1,k]) )
            m,c,b = res_robust.x
            # m = res_robust.x
        

        
            # plot fitted power functions
            x = np.arange(1,50+1)
            # axs.plot(x,c*np.sqrt(x)+b,'k:',alpha=0.5,label='$\sqrt{N}$')
            # axs.plot(x,c*x+b,'k--',alpha=0.5,label='$N$')
            # axs.plot(x,c*x**2+b,'k-',alpha=0.5,label='$N^2$')
            axs.plot(x,c*x**m+b,color='firebrick',alpha=0.5,label='%4.2f*N$^{%4.2f}$+%4.2f'%(c,m,b))
            
            
            # if k==0 and pc_idx==0:
            axs.legend(frameon=False,loc='upper left')
            
            
            if pc_idx==n_pc-1: axs.set_xlabel('N neurons')
            if k==0: axs.set_title('%s'%(titles[k]))
            else:    axs.set_title('%d PC %s'%(pc,titles[k]))
            if k<=2:
                axs.set_ylabel('variance [SU]')
                axs.set_ylim(0,np.max(singularvalues[:,0,0])*1.1)
            else:
                axs.set_ylabel('ratio')
                axs.set_ylim(-0.05,1.05)
    


    # fig.tight_layout()



    save = 0 or globalsave
    if save:
        if n_subs==1: fig.savefig(resultpath+'pca_total,explained-variance-asymptotic'+ext)
        else:         fig.savefig(resultpath+'pca_total,explained-variance-asymptotic,neuronsubsampling'+ext)










def aggregatemice_subspaces_PCA_signalvariance(datanames):
    # show contributions of signal and noise variance for various task variables
    
    # total variance = signal variance + noise variance
    # COV yi,yj = Ex[ Covy yi yj | x ] + Varx( Ey[yi|x] * Ey[yj|x] )

    n_subs = 10
    pcagroups = [      [[2,4],[],[]]   ]         # use the multimodal blocks as activities
    n_mice = len(datanames)
    
    taskaspects = ['visual','audio','context','choice']
    n_tasks = len(taskaspects)
    
    titles = ['total var.','expl. var.','cumulative expl. var.','expl. var. ratio','cumulative expl. var. ratio']
    n_pc = 4
    

    # C = [ np.zeros((n_tasks,n_pc,5)) for n in range(n_mice) ]  # last dim:  av.noi.cov., var.signal av.,  :n_pc, :all, concat(mean)
    C = []   # (mice)(n_tasks,60,5) last dim:  av.noi.cov., var.signal av.,  :n_pc, :all, concat(mean)
    n_neurons = np.zeros((n_mice),dtype=np.int16)
    
    
    for n,dn in enumerate(datanames):
        
        # C will be = np.zeros(( n_tasks,n_pc,5 ))  # avg. noise covariance. and signal-average variance, and total variances
        correlations,n_neurons[n] = pickle.load(open(cacheprefix+'subspaces/noise+signal-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'rb'))
        C.append(correlations)


    order = np.argsort(n_neurons)
    labels = [ dn+'\n  %d'%n_neurons[n] for n,dn in enumerate(datanames) ]

    if 0:
        fig,ax = plt.subplots(4,1,figsize=(24*1,6*4))
        pos = np.arange(n_mice)[order]
        width = 0.2
        for cx,comparison in enumerate(taskaspects):
            axs = ax[cx]
            bottoms = np.zeros((n_mice,5)) # will hold previous bar tops for stacked bars
            
            for n,dn in enumerate(datanames):
                n_pc_display = n_neurons[n]
                
                axs.bar(pos[n]-2*width,C[n][cx,0,4],width=width,color='black',align='edge')      # total variance
                axs.bar(pos[n]-1*width,C[n][cx,0,2],width=width,color='grey',align='edge')       # total variance for this 4 pc
                
                for px in range(n_pc_display):
        
                    axs.bar(pos[n],C[n][cx,px,0],width=width, bottom=bottoms[n,0],color=np.array([1,0,0])*(px+1)/n_pc_display,align='edge',edgecolor='white')  # noise
                    axs.bar(pos[n]+1*width,C[n][cx,px,1],width=width, bottom=bottoms[n,1],color=np.array([0,1,0])*(px+1)/n_pc_display,align='edge',edgecolor='white')   # signal
                    
                    bottoms[n,:]+=C[n][cx,px,:]   # raise the stack
            
            
            axs.set_xticks(np.arange(n_mice))
            axs.set_xticklabels([ labels[order[n]] for n in range(n_mice) ],rotation=30)
            axs.set_ylabel(comparison+'\nVar [SU$^2$]')
            
            
            if cx==0:
                axs.text(-1,47,'total variance for all # pc = # n, concat',color='black')
                axs.text(-1,44,'total variance for all # pc = # n',color='grey')
                # axs.text(-1,44,'total variance for 4 pc',color='grey')
                axs.text(-1,41,'average noise variance, pcs stacked',color='darkred')
                axs.text(-1,38,'variance of average signals, pcs stacked',color='darkgreen')
        
        
    
        # fig.tight_layout()
        fig.suptitle('Total variance vs. average noise variance and condition-variance of average signal\nfor 4 task variables, mice are listed with increasing recorded number of neurons')
    
    
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'covariance-noise,signal-timeaveraged,allpcs'+ext)


    if 1:
        # in this plot, we show how the first few PCs contain all signal variance of all available signal variance
        fig,ax = plt.subplots(4,1,figsize=(4*6,4*6))
        
        pos = np.arange(n_mice)#[order]
        width = 0.4

        for cx,comparison in enumerate(taskaspects):
            axs = ax[cx]

            bottoms = np.zeros((n_mice)) # will hold previous bar tops for stacked bars

            for n,dn in enumerate(datanames):
                # n_pc_display = n_neurons[n]
                n_pc_display = n_pc

                C_signal = C[n][cx,:n_pc_display,1]/(C[n][cx,0,2])
                
                # C_all_n_pc = C[n][cx,:n_pc_display,1].sum()
                # C_cumulative = 0
                for px in range(n_pc_display):
                    # C_cumulative += C[n][cx,px,1]
                    # if n==0: print(px,n_pc_display)
                    axs.bar(pos[n],C_signal[px],width=width, bottom=bottoms[n],\
                            color=np.array([0,1,0])*(px+1)/n_pc_display,edgecolor='white')
                    # axs.plot(n_neurons[n],C[n][cx,px,1]/C_all_n_pc,'o',color=np.array([0,1,0])*(px+1)/n_pc_display)
                    bottoms[n] += C_signal[px]
                    
            axs.set_xticks(pos)
            # axs.set_xticklabels([ labels[order[n]] for n in range(n_mice) ],rotation=30)
            axs.set_xticklabels(labels)
            axs.set_ylabel(comparison+'\nVar [SU$^2$]')


        fig.suptitle('Condition-variance of average signal: #PCs / 4 PCs total variance \nfor 4 task variables, mice are listed with increasing recorded number of neurons')


        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'covariance-signal,relative,totalvariance,4pc-timeaveraged'+ext)

        
    return








def aggregatemice_layertypecontributions_predictiveglm(datanames):
    # predict mean neural activity from task variables, and compare with decoders
    
    taskaspects = ['visual','audio','context','choice']


    n_mice = len(datanames)
    ctgroup_all = np.zeros(8)
    G = np.zeros((8,len(taskaspects)))
    B = np.zeros((8,len(taskaspects)))
    M = np.zeros(n_mice)
    cxpre = np.zeros(n_mice)
    
    modeldisplay = np.array([4,5,6,7,12,13,14,15,16,17,18],dtype='int16')
    modeldisplaylabels = ['* -vi -r','* -au -r','* -cx -r','* -ch -r',\
                          ' *-vi +r','* -au +r','* -cx +r','* -ch +r',\
                          '* -r','* +r', 'r']
    contextcomp_idx = [ [2,8], [6,9] ]       # first two is without run, last two is with run
    
    for n,dn in enumerate(datanames):

        block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=False)

#        n_neurons = block.segments[0].analogsignals[0].shape[1]

        W,mm = pickle.load(open(cacheprefix+'glm/glm_predictivecoeffs_%s'%(dn),'rb'))
        M[n] = mm

        c_db_means = pickle.load(open(cacheprefix+'subspaces/subspacecontrib,coeffs,cxvich-%s'%dn,'rb'))


        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        cxpre[n] = neph.get_auctimenormalized(acrossdecoder[1][:,0],T['starttime'],T['stimstarttime'])[0]



        ctgroups,ctcolors,ctlabels = preprocess.getcelltypesinlist(block,withnewline=False)

        for i in range(8):
            for cx in range(len(taskaspects)):
                if len(ctgroups[i])>0:
                    G[i,cx] += np.mean(np.abs(W[ctgroups[i],cx]))
                    B[i,cx] += np.mean(np.abs(c_db_means.T[ctgroups[i],cx]))

        
        mc,r,aic,dev = pickle.load(open(cacheprefix+'glm/glm_predictive,modelselection_%s'%(dn),'rb'))
        if n==0:
            R = np.zeros((n_mice,len(modeldisplay)))
            A = np.zeros((n_mice,len(modeldisplay)))
            D = np.zeros((n_mice,len(modeldisplay)))
        R[n,:] = r[modeldisplay]
        A[n,:] = aic[modeldisplay]
        D[n,:] = dev[modeldisplay]


    
    
    if 1:
        fig, ax = plt.subplots(4,2,figsize=(28,24))
    
        for cx in range(len(taskaspects)):

            axs = ax[cx,0]
            for i in range(8):
                axs.bar(i,B[i,cx],color=ctcolors[i])
#                axs.errorbar(i,B[i,cx],ctgroups[i])          # get the 1/sqrt N
            axs.set_ylabel(taskaspects[cx]+'\n$\sum$| | mapping coeffs')

            axs = ax[cx,1]
            for i in range(8):
                axs.bar(i,G[i,cx],color=ctcolors[i])
            axs.set_ylabel('$\sum$| | predictor coeffs')
    
    
        ax[0,0].legend(ctlabels)
        ax[0,0].set_title('Decoding')
        ax[0,1].set_title('Predicting')
    
        for j in [0,1]:
            ax[3,j].set_xlabel('cell types and layers')
            ax[3,j].set_xticks([]); ax[3,j].set_xticklabels([])
    
        fig.suptitle('Logistic regression decoding task variables and GLM predicting neural mean firing rate, aggregate over %d mice\n'%n_mice+\
                     '$\sum$| | of test coefficients for task variables by cell type and layer, cell # normalized\n')
            
        save = 0
        if save:
            fig.savefig(resultpath+'micetribe_glm-predictive+decoder_coefficients,neurons-visual,audio,context,choice'+ext)

    
    if 1:
        fig,ax = plt.subplots(1,2,figsize=(36,12))
        
        axs = ax[0]
        colorlist = plt.cm.viridis( np.linspace(0, 0.8, n_mice) )
        axs.bar(np.arange(n_mice),M,color=colorlist)
        axs.set_xticks(np.arange(n_mice)); axs.set_xticklabels(datanames,rotation=45)
        axs.set_ylabel('explained variance')
        
        axs = ax[1]
        axs.hist(M)
        axs.plot([M.mean(),M.mean()],[0,5])
        axs.set_xlabel('explained variance, mean=%4.3f'%M.mean())
        
        fig.suptitle('GLM predicting neural mean firing rate\nexplained total variance across mice, using visual, audio, context and choice predictors')
        save = 0
        if save:
            fig.savefig(resultpath+'micetribe_glm-predictive_coefficients,neurons-visual,audio,context,choice'+ext)



    if 1:
        colorlist = plt.cm.viridis( np.linspace(0, 0.8, n_mice) )

        fig,ax = plt.subplots(1,3,figsize=(36,24))

        axs = ax[0]
#        for i,ml in enumerate(modeldisplaylabels):
        for n in range(n_mice):
            axs.plot(R[n,:],'-o',color=colorlist[n],alpha=0.7)
        axs.set_ylim(0,0.6)
        axs.grid(True)
        axs.legend(datanames)
        axs.set_xticks(np.arange(len(modeldisplay))); axs.set_xticklabels(modeldisplaylabels,rotation=45)
        axs.set_title('explained variance')
        
        axs = ax[1]
#        for i,ml in enumerate(modeldisplaylabels):
        for n in range(n_mice):
            axs.plot(A[n,:],'-o',color=colorlist[n],alpha=0.7)
        axs.grid(True)
        axs.set_xticks(np.arange(len(modeldisplay))); axs.set_xticklabels(modeldisplaylabels, rotation=45)
        axs.set_title('akaike information criterion')

        axs = ax[2]
        for n in range(n_mice):
            axs.plot(D[n,:],'-o',color=colorlist[n],alpha=0.7)
        axs.grid(True)
        axs.set_xticks(np.arange(len(modeldisplay))); axs.set_xticklabels(modeldisplaylabels, rotation=45)
        axs.set_title('marginal log likelihood')


        
        fig.suptitle('Model Selection for GLM: *: vi,au,cx,ch;  +/-r include omit run, -vi omit visual ')
        
        save = 0
        if save:
            fig.savefig(resultpath+'micetribe_glm-predictive_modelselection-visual,audio,context,choice,run'+ext)



    if 1:    # expert test with model likelihoods
        fig,ax = plt.subplots(1,2,figsize=(24,12))
        
        x = cxpre
        
        l_titles = ['aic difference','marginal log likelihood difference']
        
        for lx,P in enumerate([A,D]):
            
            
                y = np.c_[ P[ :, contextcomp_idx[0][1] ] - P[ :, contextcomp_idx[0][0] ],\
                           P[ :, contextcomp_idx[1][1] ] - P[ :, contextcomp_idx[1][0] ] ]
        
                labels = ['LL cx+ - cx-   r- ','LL cx+ - cx-   r+ '] 
                colors = ['darkred','orange']
        
                axs = ax[lx]
    
                for m in range(2):
                    for n,dn in enumerate(datanames):
                        
                        axs.scatter(x[n],y[n,m],s=150,marker='o',color=colors[m],alpha=0.8)
                        axs.text(x[n]*1.02,y[n,m]*0.995,datanames[n],color='k',alpha=0.3,fontsize=10)
            
                    l = sp.stats.linregress(x,y[:,m])
                    if l[3]<0.13:
                        line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                        axs.plot(line,l[0]*line+l[1],color=colors[m],linewidth=2,label=labels[m])
                    else:
                        line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                        axs.plot(line,l[0]*line+l[1],'--',color=colors[m],linewidth=2,label=labels[m])
                        
                    xoffs = 0.4
        #            if sx in [3,4,5,7]: xoffs = 25
                    if l[3]<0.001:
                        axs.text(line[1]-xoffs,l[0]*line[1]+l[1]*1.02,'p<%5.3f, $R^2$=%4.2f'%((0.001,l[2]**2)),color=colors[m])
                    else:
                        axs.text(line[1]-xoffs,l[0]*line[1]+l[1]*1.02,'p=%5.3f, $R^2$=%4.2f'%((l[3],l[2]**2)),color=colors[m])
                axs.set_xlabel('context pre dec acc.')
                axs.set_ylabel(l_titles[lx])
                
                axs.legend()


        fig.suptitle('%d mice, expertism/engagement (context pre) vs. FR predictive model strength of context predictor\n'%len(datanames)+\
                     'include/omit context: cx+/cx-  and run: r+/r-;  AIC and MLL')
        
        save = 0
        if save:
            fig.savefig(resultpath+'micetribe_glm-predictive_expert,loglike-visual,audio,context,choice,run'+ext)








def aggregatemice_contextpreonsubpop(datanames):
    # test if without putative inhibitory cells the excitatory only cells provide similar context representation


    n_mice = len(datanames)

    dispersions = []
    datanames_dispersion = []
    n_th = 1         # min number of neurons per mice for aggregate displays
    
    fs = np.ceil(np.sqrt(n_mice)).astype(np.int16)
    cmap = figs.getcorrcolormap('correlation')
    fig,ax = plt.subplots(fs,fs,figsize=(fs*12,fs*8))
    
    comparison = 'context'
    for n,dn in enumerate(datanames):
        print(dn)
        blv,bla = preprocess.getorderattended(dn)
        block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=False)
        n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
        celltypes = block.annotations['celltypes']
        
        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))

        c_db_matrix = np.zeros((n_neurons,n_trajectory))       # prepare the dbnv vector holder

        wx = int((len(acrossdecoder)-7)/n_neurons)
        coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
        # mean of on stimulus
        c_db_matrix = coeff[:,:,0].T      # neurons times trajectory
        times = block.segments[0].analogsignals[0].times[:-wx]


        # order the neurons by average activity during stimulus presentation
        if 0:    # change to true if want to separate inhibitory and excitatory neurons
            order1 = np.argsort(c_db_matrix[T['stimstart_idx']:T['stimend_idx'],celltypes==1].mean(axis=0))
            order2 = np.argsort(c_db_matrix[T['stimstart_idx']:T['stimend_idx'],celltypes==0].mean(axis=0))
            cellmatrix1 = c_db_matrix[:,celltypes==1][:,order1]
            cellmatrix2 = c_db_matrix[:,celltypes==0][:,order2]
            print(cellmatrix1.shape,cellmatrix2.shape)
            cellmatrix = np.hstack( (cellmatrix1,cellmatrix2) )
            changepoint = np.sum(celltypes==1).astype(np.int16)
        else:
            order = np.argsort(c_db_matrix[T['stimstart_idx']:T['stimend_idx'],:].mean(axis=0))
            cellmatrix = c_db_matrix[:,order]

        axs = ax[n//fs,n%fs]
        # axs.plot(times,c_db_matrix)
        # cf = axs.pcolormesh(times,np.arange(n_neurons+1),c_db_matrix[:,order].T,vmin=-0.5,vmax=+0.5,cmap=cmap)
        cf = axs.pcolormesh(times,np.arange(n_neurons+1),cellmatrix.T,vmin=-0.5,vmax=+0.5,cmap=cmap)
        fig.colorbar(cf,ax=axs,ticks=np.arange(-1,1.1,0.5))
        # axs.hlines(changepoint,T['starttime'],T['endtime'],color='white',lw=2)
        axs.vlines( (T['stimstarttime'],T['stimendtime']),0,n_neurons,color='white',lw=2)
        # axs.set_ylim(0,n_neurons)
        axs.set_yticks([])
        figs.setxt(axs)
        axs.set_title(dn)
        
        # add the dispersions for all neurons:
        m_o = c_db_matrix[T['stimstart_idx']:T['stimend_idx'],:].mean(axis=0)
        m_p = np.r_[ c_db_matrix[0:T['stimstart_idx'],:],  c_db_matrix[T['stimend_idx']:,:]  ]
        m_p = m_p.mean(axis=0)
        
        if n_neurons>n_th:
            dispersions.append(np.c_[m_p,m_o])
            datanames_dispersion.append(dn)
        

        
    fig.suptitle('context decoder weights (color) timecourse (horizontal axis),\nneurons (vertical axis) separated as excitatory and inhibitory cells\nordered by mean "ON" weight')


    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'tribe_contextdecodercoeffs-celltypesplit-colormesh'+ext)










    if 1:

        H = np.concatenate([ np.abs(dispersions[n].reshape((-1,1))) for n in range(len(datanames_dispersion)) ])
        x_res = np.arange(0,0.501,0.01)
        hg,_ = np.histogram(H, bins=x_res, weights=H)
        hg /= np.sum(hg)
        chg = np.cumsum(hg)
        w = np.where(chg>0.1)
        print(w[0][0],x_res[w[0][0]])
        
    
        # th = 0.075         # minimum value of decoder coefficient for aggregate displays
        th = x_res[w[0][0]]
    
        fig,ax = plt.subplots(1,6,figsize=(6*10,10))
        axs = ax[0]
        for n,dn in enumerate(datanames_dispersion):
            axs.plot(dispersions[n][:,0],dispersions[n][:,1],'o',alpha=0.5,label=dn+' %d'%(dispersions[n].shape[1]))
            axs.set_xticks(np.arange(-1,1.1,0.25))
            axs.set_xlim(-0.5,0.5)
            axs.set_yticks(np.arange(-1,1.1,0.25))
            axs.set_ylim(-0.5,0.5)
            
            figs.plottoaxis_crosshair(axs)
            
            axs.set_xlabel('pre')
            axs.set_ylabel('on')
        axs.legend(frameon=False)
    
    
        rectangle = plt.Rectangle((-th, -th), 2*th, 2*th, color='grey',lw=1.5,fill=False)
        axs.add_artist(rectangle)
    
            
        axs = ax[1]
        ad = np.concatenate(dispersions,axis=0)
        axs.hist2d(ad[:,0],ad[:,1],bins=np.linspace(-1,1.001,50))
        axs.set_xlabel('pre')
        axs.set_ylabel('on')
        axs.set_xticks(np.arange(-1,1.1,0.25))
        axs.set_xlim(-0.5,0.5)
        axs.set_yticks(np.arange(-1,1.1,0.25))
        axs.set_ylim(-0.5,0.5)
        figs.plottoaxis_crosshair(axs,color='white')
        
        # rectangle = plt.Rectangle((-th, -th), 2*th, 2*th, lw=2, color='red',fill=False)
        # axs.add_artist(rectangle)
        axs.vlines([-th,th],-0.5,0.5,color='white')
        axs.hlines([-th,th],-0.5,0.5,color='white')
    
        axs.text(0.3,0.04,'th=%4.2f, 10%% weigthed perc.'%th)
    
        s_label = ['constant','sign change','on-off']
        axs.text(-0.4,0.04,s_label[0],color='white')
        axs.text(-0.4,0.0,s_label[2],color='white')
        axs.text(-0.4,-0.04,s_label[1],color='white')
    
    
        axs = ax[2]
        counts = np.array(\
                 [ np.sum( (ad[:,0]>=th) & (ad[:,1]>=th) ),\
                   np.sum( (ad[:,0]<=-th) & (ad[:,1]<=-th) ),\
                   np.sum( (ad[:,0]>=th) & (ad[:,1]<=-th) ),\
                   np.sum( (ad[:,0]<=-th) & (ad[:,1]>=th) ),\
                   np.sum( (  np.abs(ad[:,0])<th) & (ad[:,1]>=th) ),\
                   np.sum( (  np.abs(ad[:,0])<th) & (ad[:,1]<=-th) ),\
                   np.sum( (  ad[:,0]>=th) & (np.abs(ad[:,1])<th) ),\
                   np.sum( (  ad[:,0]<=-th) & (np.abs(ad[:,1])<th) ),\
                       ] )
    
        pos = np.arange(len(counts))
        axs.bar( pos,counts )
        axs.set_xticks(pos)
        axs.set_xticklabels(['pre+ on+','pre- on-','pre+ on-','pre- on+','pre0 on+','pre0 on-','pre+ on0','pre- on0'],rotation=90)
        axs.set_ylabel('number of cells weight>%4.2f'%th)
        

        axs = ax[3]
        counts = np.array(\
                [  np.sum( (np.abs(ad[:,0])>th ) &  (np.abs(ad[:,1])>th )     ),\
                   np.sum( np.logical_or( (ad[:,0]>=th) & (ad[:,1]<=-th), (ad[:,0]<=-th) & (ad[:,1]>=th) )    ),\
                   np.sum( np.logical_or(  (np.abs(ad[:,0])<th) & (np.abs(ad[:,1])>=th),  (np.abs(ad[:,0])>=th) & (np.abs(ad[:,1])<th) )   ),\
                       ] )
        counts = counts/np.sum(counts)*100
        pos = np.arange(len(counts))
        axs.bar( pos,counts )
        axs.set_xticks(pos)
        axs.set_xticklabels(['constant','sign change','on-off'])
        axs.set_ylabel('percent of units')
        

        
        axs = ax[4]      # histogram to establish noise ceiling
        for i in [0,1]:
            H = []
            H = np.concatenate([ dispersions[n][:,i] for n in range(len(datanames_dispersion)) ])
            axs.hist(H, bins=np.arange(-0.5,0.501,0.05), weights=H, color=['red','blue'][i],alpha=0.5,label=['pre','on'][i])
        axs.vlines([th,-th],-3,3,color='k')
        axs.legend(frameon=False)
        
        axs = ax[5]      # histogram to establish noise ceiling all pre and on weights
        H = np.concatenate([ np.abs(dispersions[n].reshape((-1,1))) for n in range(len(datanames_dispersion)) ])
        axs.hist(H, bins=np.arange(0,0.501,0.025), weights=H, color='darkseagreen')
        axs.vlines([th],-3,3,color='k')
        
        # print('q percentiles 0.1:', w  )
        # lH = np.log(H)
        # lm = np.mean(lH)
        # ls = np.std(lH)
        # axs.plot(np.arange(0,0.501,0.025),sp.stats.lognorm.pdf(np.arange(0,0.501,0.025),np.exp(ls),np.exp(lm)))
        
        
        fig.suptitle('All neurons, mean pre vs. on context decoder weights; mice neurons>%d'%n_th)
        
        # fig.tight_layout()
        
        save = 0 or globalsave
        if save:
            if n_th==1:
                fig.savefig(resultpath+'tribe_contextdecodercoeffs-meanperneuron'+ext,bbox_inches='tight')
            else:
                fig.savefig(resultpath+'tribe_contextdecodercoeffs-meanperneuron,onlyNgt%d'%n_th+ext,bbox_inches='tight')









def aggregatemice_behavioursymmetrycontext(datanames):
    # show symmetrically and assymetrically behaving mouse
    # cross decode between parts

    n_mice = len(datanames)
    accuracylist = []
    accuracyalllist = []
    for dn in datanames:
        # accuracy (times,train-test-crosstrain-crosstest,stats,symmetrytrain)
        accuracies,coefs,accuracyall,coefsall = pickle.load(open(cacheprefix+'symmetry/neural,context-equalized-symmetric,antisymmetric,cross-decoder-loo,boots,timecourse_%s.pck'%(dn), 'rb'))
        accuracylist.append(accuracies)
        accuracyalllist.append(accuracyall)
    accuracylist = np.array(accuracylist) # (mice,bootstrap,times,train-test-crosstrain-crosstest,stats,symmetrytrain)
    accuracyalllist = np.array(accuracyalllist) # (mice,bootstrap,times,train-test,stats)
    n_bootstrap = accuracylist.shape[1]


    n_times = accuracyalllist.shape[1]

    mask_on = np.zeros(n_times, dtype=np.bool8)
    mask_on[T['stimstart_idx']:T['stimend_idx']] = True
    mask_off = np.ones(n_times, dtype=np.bool8)
    mask_off[T['stimstart_idx']:T['stimend_idx']] = False





    # load baseline chance levels for subsampled number of trials ('reduced'), matching to minimum symmetry number of trials
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))
    chances_reduced = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d,reduced-%s.pck'%(n_resample,continuous_method),'rb'))
    






    # figure
    fig,ax = plt.subplots(3,2, figsize=(2*10,3*8) )




    for mx,mask in enumerate((mask_off, mask_on)):
        n_subtrialtimes = np.sum(mask)


        # plot all trials
        axs = ax[0,mx]
        

        m = accuracyalllist[:,mask,1,0].mean(axis=1)
        e = np.sqrt( (accuracyalllist[:,mask,1,0].std(axis=1)/np.sqrt(n_subtrialtimes))**2 + accuracyalllist[:,mask,1,2].mean(axis=1)**2 )

        axs.bar(x=np.arange(n_mice)*6-2, height=m, yerr=e, color='grey', label='all', alpha=0.8)


        symmetrycolors = [['rebeccapurple','gold'],['darkorange','fuchsia']]

        xranges = np.arange(n_mice)*6-2
        for n in range(n_mice):
            chs = chances[datanames[n]]
            figs.plottoaxis_chancelevel( axs, chs, xs=[ xranges[n]-0.5, xranges[n]+0.5] )
        
        
        
        # plot crossdecoding
        symmetrylabels = [['symmetric->symmetric','symmetric->asymmetric'],['antisymmetric->asymmetric','asymmetric->symmetric']]
        symmetrycolors = [['rebeccapurple','gold'],['darkorange','fuchsia']]
        for rx in range(2):             # train subset
            for sx in range(2):         # crosstest subset (same or cross)
                
                color = symmetrycolors[rx][sx]


                # (n_timestamps,train/test,stats,task,symmetry)
                m = accuracylist[:,:,mask,1+sx*2,0,rx].mean(axis=(1,2))

                e = np.sqrt( (accuracylist[:,:,mask,1+sx*2,0,rx].std(axis=(1,2))/np.sqrt(n_bootstrap*n_subtrialtimes))**2 + \
                              accuracylist[:,:,mask,1+sx*2,2,rx].mean(axis=(1,2))**2 )
                print(m.shape,e.shape)
                axs.bar(x=np.arange(n_mice)*6-2+1+rx*2+sx, height=m, yerr=e, color=color, label=symmetrylabels[rx][sx])


            xranges = np.arange(n_mice)*6
            for n in range(n_mice):
                chs = chances_reduced[datanames[n]]
                figs.plottoaxis_chancelevel( axs, chs, xs=[ xranges[n]-1, xranges[n]+2.5] )



        axs.set_xticks(np.arange(n_mice)*6)
        axs.set_xticklabels(datanames, rotation=60)

        axs.set_ylim(0.45,1.05)
        figs.plottoaxis_chancelevel(axs,0.5)
        axs.legend(frameon=False)

        axs.set_ylabel('time averaged accuracy')
        axs.set_title(['off stimulus','on stimulus'][mx])





        # compare same-test between sym and asym
        axs = ax[1,mx]

        sx = 0

        m_sym = accuracylist[:,:,mask,1,0,0].mean(axis=1).T
        m_asym = accuracylist[:,:,mask,1,0,1].mean(axis=1).T
        e_sym = accuracylist[:,:,mask,1,2,0].mean(axis=1).T
        e_asym = accuracylist[:,:,mask,1,2,1].mean(axis=1).T

        d = m_sym-m_asym#-(e_sym+e_asym)

        P = axs.violinplot(d,  quantiles=np.tile([0.025, 0.15866, 0.84134, 0.975 ],(n_mice,1)).T,\
                           showextrema=False, showmedians=True)         # color='mediumturquoise',
        # P = axs.boxplot(d,  whis=[ 2.5, 97.5 ], notch=True, bootstrap=1000)
        axs.set_xticks(np.arange(n_mice)+1)
        axs.set_xticklabels(datanames, rotation=60)
        figs.plottoaxis_chancelevel(axs,0)
        
        # for n in range(n_mice):
        #     axs.plot([n+0.5,n+1.5],)

        axs.set_ylabel('accuracy difference\nsymmetric-asymmetric')







        # random chance controlled time comparison (3rd and 4th rows)


        axs = ax[2,mx]
        sx = 0
        
        
        # random criteria sym only:
        timedmasks = [ np.logical_and(mask, accuracylist[n,:,:,1,0,0].mean(axis=0)>chances_reduced[datanames[n]] ) for n in range(n_mice) ]
        # random criteria sym and asym as well:
        # timedmasks = [ np.logical_and( mask, np.logical_and(accuracylist[n,:,:,1,0,0].mean(axis=0)>chances_reduced[datanames[n]],\
        #                                                    accuracylist[n,:,:,1,0,1].mean(axis=0)>chances_reduced[datanames[n]]) ) \
        #                  for n in range(n_mice) ]

        print([ sum(timedmask) for timedmask in timedmasks])


        m_sym = [ accuracylist[n,:,timedmasks[n],1,0,0].mean(axis=1) for n in range(n_mice)]
        m_asym = [ accuracylist[n,:,timedmasks[n],1,0,1].mean(axis=1) for n in range(n_mice)]
        e_sym = [ accuracylist[n,:,timedmasks[n],1,2,0].mean(axis=1)  for n in range(n_mice)]
        e_asym = [ accuracylist[n,:,timedmasks[n],1,2,1].mean(axis=1) for n in range(n_mice)]

        d = [ m_sym[n]-m_asym[n]    for n in range(n_mice) ]        #-(e_sym+e_asym)

        # P = axs.violinplot(d,  quantiles=np.tile([0.025, 0.15866, 0.84134, 0.975 ],(n_mice,1)).T,\
        #                    showextrema=False, showmedians=True)         # color='mediumturquoise',
        P = axs.boxplot(d,  whis=[ 2.5, 97.5 ], notch=True, bootstrap=1000)
        axs.set_xticks(np.arange(n_mice)+1)
        axs.set_xticklabels(datanames, rotation=60)
        figs.plottoaxis_chancelevel(axs,0)



        axs.set_ylabel('above random accuracy difference\nsymmetric-asymmetric')




    fig.tight_layout()
    
    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'behaviour,symmetry,context,invariant'+ext)






















def aggregatemice_chancelevel(datanames, reducedtrials=False):
    # calculate randomized chance level significance boundary around 0.5 accuracies
    
    n_mice = len(datanames)
    n_resample = 10
    taskaspects = ['visual','audio','context','choice']
    if reducedtrials: reducedlabel = ',reduced'
    else: reducedlabel = ''

    chances = dict()
    
    for n,dn in enumerate(datanames):

        acrossdecoders = [] # get together all random runs
        for cx,comparison in enumerate(taskaspects):
            acrossdecoders_taskaspect = pickle.load(open(cacheprefix+'subspaces/shuffledecoders,resampled-%s-full,r%d%s-%s_%s.pck'%(comparison,n_resample,reducedlabel,continuous_method,dn),'rb'))
            acrossdecoders.extend(acrossdecoders_taskaspect)


        
        
        D = np.array([acrossdecoders[rx][1] for rx in range(n_resample*len(taskaspects))])
        m = D.mean(axis=0)                        #  mean across resample randomized trial labels
        e = 2*D.std(axis=0)/np.sqrt(n_resample*len(taskaspects))   # 2 sem across resample randomized trial labels
        # mean + sem of cv-means over resamples
        # plus 1 sem around mean of cv-means over resamples
        
        c = m[:,0]+m[:,2]/2+e[:,0]/2+e[:,2]/2
        chances[dn] = c.mean()
        chances[dn] += c.std()/np.sqrt(len(c))

    print(chances)

    fn = 'subspaces/chances,allmice,resampled-full,r%d%s-%s.pck'%(n_resample,reducedlabel,continuous_method)
    print('saving to ', fn)
    pickle.dump(chances, open(cacheprefix+fn,'wb'))









def aggregatemice_numberofhighlowperformancetrials(datanames):
    from physiology import get_mask_cleverness
    n_mice = len(datanames)
    print('Calculating number of high and low performance trials.')
    trialnumbers = np.zeros((n_mice,4))
    fractioncorrect = np.zeros((n_mice,4))        # all
    fractioncorrect_in = np.zeros((n_mice,4))     # incongruent nogos
    for n,dn in enumerate(datanames):
        print(n,dn)
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
        # mask_clevers = np.vstack(mask_clevers).T


        mask_contextuals = [ np.prod(mask_clevers_list[cx],axis=1).astype(bool) for cx in [0,1] ]
        # display the number of trials that has above threshold movingaverage for all 4 combinations of congruenct and action
        trialslabel = '1st(%s):%d/%d, 2nd(%s):%d/%d'%(['V','A'][blv[1]==4],np.sum(mask_contextuals[0]), len(mask_contextuals[0]), ['A','V'][blv[1]==4], np.sum(mask_contextuals[1]), len(mask_contextuals[1]))
        trialnumbers[n,:] = np.array([np.sum(mask_contextuals[0]), len(mask_contextuals[0])-np.sum(mask_contextuals[0]), np.sum(mask_contextuals[1]), len(mask_contextuals[1])-np.sum(mask_contextuals[1]) ])
        fractioncorrect[n,:] = np.array([np.mean(g['success'][g['block']==2][mask_contextuals[0]]), 1-np.mean(g['success'][g['block']==2][mask_contextuals[0]]),\
                                         np.mean(g['success'][g['block']==4][mask_contextuals[1]]),1-np.mean(g['success'][g['block']==4][mask_contextuals[1]]) ])
        fractioncorrect_in[n,:] = np.array([np.mean(g['success'][g['block']==2][mask_clevers_list[0][:,3]]), 1-np.mean(g['success'][g['block']==2][mask_clevers_list[0][:,3]]),\
                                         np.mean(g['success'][g['block']==4][mask_clevers_list[1][:,3]]),1-np.mean(g['success'][g['block']==4][mask_clevers_list[1][:,3]]) ])


        # for i in range(2):
        #     for k in range(4):
        #         m = np.mean(g['success'][g['block']==(i+1)*2][mask_clevers_list[i][:,k]])
        #         print(i,k,m)


    pickle.dump(trialnumbers, open(cacheprefix+'behaviour/numtrials-n%d-highlowperformance'%n_mice,'wb'))
    pickle.dump(fractioncorrect, open(cacheprefix+'behaviour/fraccorrect-n%d-highlowperformance'%n_mice,'wb'))
    pickle.dump(fractioncorrect_in, open(cacheprefix+'behaviour/fraccorrect,in-n%d-highlowperformance'%n_mice,'wb'))






















def aggregatemice_modelLLs(datanames, onlyboth=True):
    from physiology import behaviour_symmetry_highperformance
    n_mice = len(datanames)

    LL = np.zeros((n_mice,2,2,3,6))  # n_mice, context, congruency, gonogoboth, model
    meanresponses = np.zeros((n_mice,2,2,3))
    mp = np.zeros((n_mice,2,2)) # n_mice, context, congruency
    for dx,dn in enumerate(datanames):
        LL[dx,:,:,:,:], meanresponses[dx,:,:,:], mp[dx,:,:] = behaviour_symmetry_highperformance(dn)

    
    pickle.dump(LL, open(cacheprefix+'behaviour/LL-strategies-highformance.pck','wb'))
    

    # plot LLs (log likelihoods of models)
    # make a subplot for each context and congruency
    # and plot a barplot for each model
    # as a mean and errorbars between mice



    # plot constants
    labels_contexts = ['visual','auditory']
    labels_congruency = ['congruent','incongruent']
    labels_models = [ 'contextual','opposite','lick bias','contextual w/lick bias','opposite w/lick bias','squeezed contextual']
    colors_models = ['purple','red','darkorange','seagreen','gold','fuchsia']


    # plot mean and all individual mouse

    LLm = np.mean(LL,axis=0)
    LLe = np.std(LL,axis=0)/np.sqrt(n_mice)

    n_models = LL.shape[4]


    if globaldoplot or 0:
            

        fig,ax = plt.subplots(1+n_mice,2+4,figsize=((4+2)*6,(n_mice+2)*6),sharex=True,sharey=True)


        # overal marginal over mice and congruency context etc.
        axs = ax[0,0]
        for i in range(n_models):
            axs.bar(x=[i-0.33,i,i+0.33], height=LLm.mean(axis=(0,1))[:,i], yerr=LLe.sum(axis=(0,1))[:,i], color=colors_models[i], alpha=1,
                    width=0.8/3, label=labels_models[i], hatch=[None,None,'.'] )
        axs.legend(frameon=False)
        axs.set_title('mean total')
        axs.set_ylabel('LL')
        axs.set_ylim(-1.5,0)

        # incongruent only overal marginal over mice and context etc.
        axs = ax[0,1]
        for i in range(n_models):
            axs.bar(x=[i-0.33,i,i+0.33], height=LLm[:,1::2,:,:].mean(axis=(0,1))[:,i], yerr=LLe[:,1::2,:,:].sum(axis=(0,1))[:,i],
                    width=0.8/3, color=colors_models[i], hatch=[None,None,'.'])
        axs.set_title('mean incongruent total')
        


        # marginal over congruency context etc. for each mouse, and then only for incongruents
        for n in range(n_mice):
            axs = ax[n+1,0]
            for i in range(n_models):
                axs.bar(x=[i-0.33,i,i+0.33],height=LL[n,:,:,:,:].mean(axis=(0,1))[:,i],
                        width=0.8/3, color=colors_models[i], alpha=0.8, hatch=[None,None,'.'])
                axs.set_title(datanames[n])

            if n==n_mice: axs.set_xlabel('models')

            axs = ax[n+1,1]
            for i in range(n_models):
                axs.bar(x=[i-0.33,i,i+0.33],height=LL[n,:,1::2,:,:].mean(axis=(0,1))[:,i],
                        width=0.8/3, color=colors_models[i], alpha=0.8, hatch=[None,None,'.'])

            if n==n_mice: axs.set_xlabel('models')




        # plot each context and congruency individually for all mice averaged, and all mice separately
        for cx in [0,1]:
            for ix in [0,1]:

                axs = ax[0,cx*2+ix+2]
                for i in range(n_models):
                    axs.bar(x=[i-0.33,i,i+0.33],height=LLm[cx,ix,:,i], yerr=LLe[cx,ix,:,i],
                            width=0.8/3, color=colors_models[i], hatch=[None,None,'.'])
                axs.set_title('%s %s'%(labels_contexts[cx],labels_congruency[ix]))


                for n in range(n_mice):
                    axs = ax[n+1,cx*2+ix+2]
                    for i in range(n_models):
                        axs.bar(x=[i-0.33,i,i+0.33],height=LL[n,cx,ix,:,i],
                                width=0.8/3, color=colors_models[i], alpha=0.5, hatch=[None,None,'.'])

                        # annotate lick bias:
                        for gx in [0,1,2]: axs.text([1.67,2,2.33][gx],LL[n,cx,ix,gx,2], '%4.2f'%meanresponses[n,cx,ix,gx],
                                                    verticalalignment='top',horizontalalignment='center', fontsize=8)
                        # annotate squeezed contextual:
                        axs.text(5.33,LL[n,cx,ix,2,5], '%4.2f'%mp[n,cx,ix], verticalalignment='top',horizontalalignment='center', fontsize=8)
                    if cx==0 and ix==0: axs.set_title(datanames[n])

                    if n==n_mice: axs.set_xlabel('models')


        # fig.plot_title('model LLs over mice n=%d'%n_mice)

        if globalsave or 0:
            # fig.savefig(resultpath+'tribe-modelLLs-highperformance'+ext)
            fig.savefig('tribe-modelLLs-highperformance'+ext)
















