#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 10:27:02 2020

@author: mahajnal
"""

from matplotlib import animation
from config import *
import preprocess


def compareacrossconditions_maxnormdiffdecoder(dn,block,examine='allexpcond'):
    # group trials into two classes along different task variables
    # compare groups with a firing rate difference metric
    # compare groups with decoders: crossvalidated predictive accuracy
    # save decoder learnings
    # make figures

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot
    
#    layerlist = [1,2,3,4]       # choose which layer to focus on: 0 all, 1 superf, 2 input, 3 deep5, 4 deep6
    layerlist = [0]
    celltypelist = [0,1]        # choose which spike width to use: 0 narrow, 1 broad
    useorleavegroup = 2
    leaveorusegroup = ['leave','only','random']       # leave out groups; random is baseline left-out average
    randomizedleftoutaverage_nruns = 5

    width = 50*pq.ms               # larger input sizes do not help
    if continuous_method=='count': width = 1000*pq.ms

    blv,bla = preprocess.getorderattended(dn)
    if examine=='allexpcond':
        comparisongroups  = [ \
                              [  [[ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                                 [[ [blv[1]], [45],    [5000] ], [ [blv[1]],[135],     [5000] ]  ] ],\
                              [  [[ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                                 [[ [bla[1]],  [45],[5000] ], [ [bla[1]],   [45],[10000] ]    ] ],\
                              [  [[  blv,  [],[] ],   [ bla,   [], [] ] ],\
                                 [[ [blv[1]],  [45],[5000] ], [ [bla[1]],   [45],[5000] ]   ] ], \
                              [[],[]] \
                            ]  # note that choice is determined by performance, not by block setup, see inline later
        horizontals = ['visual','audio','context','choice']
        logical = 'and'
        trialconditions = ['all']#,'identical,go']
        
    elif examine=='gonogo':
        comparisongroups  = [ \
                              [  [ [ [blv[1]], [45], [], [bla[1]],[],[5000] ] ,    [ [blv[1]], [135], [], [bla[1]], [], [10000] ] ]  ] \
                            ]
        horizontals = ['go-nogo']
        logical = 'or' # got this 6 components needed in the comparisongroup smallest container
        trialconditions = ['']

    elif examine=='attendignore':
        comparisongroups  = [ \
                              [  [ [ [blv[1]], [45],        [] ],  [ [blv[1]],[135],        [] ] ],\
                                 [ [ [blv[1]], [45],    [5000] ],  [ [blv[1]],[135],    [5000] ]  ] ],\
                              [  [ [ [bla[1]], [45],        [] ],  [ [bla[1]],[135],        [] ] ],\
                                 [ [ [bla[1]], [45],    [5000] ],  [ [bla[1]],[135],    [5000] ] ]  ],\
                              [  [ [ [bla[1]],   [],    [5000] ],  [ [bla[1]],[],      [10000] ] ],\
                                 [ [ [bla[1]], [45],    [5000] ],  [ [bla[1]],[45],    [10000] ] ]  ],\
                              [  [ [ [blv[1]],   [],    [5000] ],  [ [blv[1]],[],      [10000] ] ],\
                                 [ [ [blv[1]], [45],    [5000] ],  [ [blv[1]],[45],    [10000] ] ]  ]\
                            ]
        horizontals = ['attend visual','ignore visual','attend audio','ignore audio']
        logical = 'and'
        trialconditions = ['all']#,'identical,go']
        
        
    elif examine=='character':
        comparisongroups  = [ \
                              [  [ [ [blv[1]], [45],        [] ],  [     [blv[1]],[135],    [] ] ]  ],\
                              [  [ [ [5],    [45],        [] ],  [        [5],[135],    [] ] ]  ],\
                              [  [ [ [5],    [90],        [] ],  [        [5],[180],    [] ] ]  ],\
                              [  [ [ [5],   [270],        [] ],  [        [5],[180],    [] ] ]  ]\
                              ]
                                 
        
        horizontals = ['task,45-135','char,45-135','char,90-180','char,270-180']
        logical = 'and'
        trialconditions = ['all']#,'identical,go']
        

    else: print('examination type not found');   return



    print('layers',layerlist,'logical',logical)
    print(comparisongroups)
    print(horizontals)

    for cidx in celltypelist:
        if layerlist[0]==0 and cidx>0: continue    # don't run over cells, if it's all layers and celltypes
        for layered in layerlist:
    
            # sort out which layers and celltypes to leave out or only use
            
            # print('layers:', block.annotations['celltypes'], 'celltypes:', block.annotations['waveforms'])
    
            if layered>0:
                lidx = layered
                layeredfilename = ',%s,%d-%s,%s'%(leaveorusegroup[useorleavegroup],lidx,layernames[lidx-1],spikewidthnames[cidx])
                layeredtitle = leaveorusegroup[useorleavegroup]+': '+layernames[lidx-1]+', '+spikewidthnames[cidx]
                print(layeredtitle)
                if useorleavegroup==0:    # leave this group
                    lmask = np.where(    np.logical_not(np.logical_and(block.annotations['celltypes']==lidx,block.annotations['waveforms']==cidx)) )[0]
                elif useorleavegroup==1:   # use only thist group
                    lmask = np.where(    np.logical_or(block.annotations['celltypes']==lidx, block.annotations['waveforms']==cidx) )[0]
                elif useorleavegroup==2:
                    max_neurons = np.sum(np.logical_not(np.logical_and(block.annotations['celltypes']==lidx,block.annotations['waveforms']==cidx)))
                    lmask = np.random.permutation(len(block.annotations['celltypes']))[:max_neurons]
                if (useorleavegroup==1 and len(lmask)==0) or (useorleavegroup==0 and len(lmask)==len(block.annotations['celltypes'])):
                    print(dn,'layer+celltype empty: ',layeredtitle)
                    continue
                
            else:
                lmask = np.arange(block.segments[0].analogsignals[0].shape[1])
                layeredfilename = 'all'
                layeredtitle = 'all'
                


    
            # setup stimulus:
            maxcolumns = len(trialconditions)*len(horizontals)

            if doplot:
                if not publish:
                    fig,ax = plt.subplots(3,maxcolumns, figsize=(11*maxcolumns,30))
                else:
                    fig,ax = plt.subplots(2,2,figsize=(15,10))
    
            for cx,comparison in enumerate(horizontals):
                for sx,stimulusspecificity in enumerate(trialconditions):
                    print(cx,sx,comparison,stimulusspecificity,comparisongroups[cx][sx])
                    
                    stimulusIDgroups = comparisongroups[cx][sx]
        
            
                    # collect neural responses
                    if not comparison=='choice': # visual, audio, context:
                        acrossresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups,logical)
                    else:  # choice:
                        acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
                    
    
                
                    print(len(acrossresponses[0][0]),type(acrossresponses[0][0]))
                    if layered:
                        for j,response in enumerate(acrossresponses):
                            for k,asig in enumerate(response):
                                acrossresponses[j][k] = acrossresponses[j][k][:,lmask]
                
            
                    # calculate difference metrics
                    # usually only use after else; special treatment for randomized leave-group-out runs useorleavegroup=2
                    if recalculate:
                        if len(layerlist)>1 and useorleavegroup==2:
                            for runs in range(randomizedleftoutaverage_nruns):
                                aux = nedi.get_responsedifference(acrossresponses)
                                if runs==0:
                                    acrossdifference = [ aux[a]/randomizedleftoutaverage_nruns for a in range(3) ]
                                else:
                                    for a in np.arange(3):
                                        acrossdifference[a] += aux[a]/randomizedleftoutaverage_nruns
                        else: # generally you need this:
                            acrossdifference = nedi.get_responsedifference(acrossresponses)
                        # now save
                        pickle.dump(acrossdifference,open(cacheprefix+'continuous/responsedifference-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'wb'))
                    else:
                        acrossdifference = pickle.load(open(cacheprefix+'continuous/responsedifference-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            
                    if recalculate:
                        if len(layerlist)>1 and useorleavegroup==2:
                            for runs in range(randomizedleftoutaverage_nruns):
                                aux = nedi.get_responsedecoding(acrossresponses,width=width)
                                if runs==0:
                                    acrossdecoder = [ aux[a]/randomizedleftoutaverage_nruns for a in range(3) ]
                                else:
                                    for a in np.arange(3):
                                        acrossdecoder[a] += aux[a]/randomizedleftoutaverage_nruns
                        else:  # generally you need this:
                            acrossdecoder = nedi.get_responsedecoding(acrossresponses,width=width)
                        # now save
                        pickle.dump(acrossdecoder,open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'wb'))
                    else:
                        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
                    

        
                    
                # plot metrics
                    if doplot:
                        if not publish:
                            if maxcolumns>1: axs = ax[0,cx*len(trialconditions)+sx]
                            else: axs = ax[0]
                            figs.plottoaxis_difference(acrossdifference,axs)
                            figs.plottoaxis_stimulusoverlay(axs,T)
                            figs.plottoaxis_chancelevel(axs)
                            if cx==0 and sx==0: axs.legend();axs.set_ylabel('z-score difference')
                            axs.set_title('%s bw/%s cond.'%(comparison,stimulusspecificity))
                            
                            
                            if maxcolumns>1: axs = ax[1,cx*len(trialconditions)+sx]
                            else: axs = ax[1]
                            figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs)       # plot the performance
                            figs.plottoaxis_stimulusoverlay(axs,T)
                            figs.plottoaxis_chancelevel(axs,0.5,'chance level')
                            if cx==0 and sx==0: axs.legend(); axs.set_ylabel('roc auc at each timepoint')
                            
                            if maxcolumns>1: axs = ax[2,cx*len(trialconditions)+sx]
                            else: axs = ax[2]
                            figs.plottoaxis_decoderangles(acrossdecoder[2:7],axs)       # plot 5 of the subspace directions
                            figs.plottoaxis_chancelevel(axs,np.pi/2,'orthogonal')
                            figs.plottoaxis_chancelevel(axs,0.0,'parallel')
                            axs.set_ylim([0-0.1,np.pi/2+0.1])
                            figs.plottoaxis_stimulusoverlay(axs,T)
                            axs.set_yticks([0,np.pi/2]); axs.set_yticklabels(['$0$','$\pi/2$'])
                            if cx==0 and sx==0: axs.legend(); axs.set_ylabel('angles between\ndecision boundary normal and principal components')
                            axs.set_xlabel('time from trial onset [ms]')
                        
                        
                        else:     # for publication to CCN2019
                            captionsigns = ['A','B','C','D']
                            if examine=='allexpcond':
                                decodercolors = ['navy','darkgreen','mediumvioletred','darkorange']
                            elif examine=='attendignore':
                                decodercolors = ['darkred','navy','orangered','darkcyan']
                            axs = ax[int(cx/2), cx%2]
                            figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs,colorlist=['',decodercolors[cx]],plottrain=False)       # plot the test performance
                            figs.plottoaxis_stimulusoverlay(axs,T)
                            figs.plottoaxis_chancelevel(axs,0.5,'chance level')
                            figs.setxt(axs)
                            axs.set_yticks([0.5,1.0])# axs.set_yticklabels([0.5,1.0])
                            axs.text(-2440,1.04,captionsigns[cx],fontsize=24,fontweight='bold')
                            if cx>1: axs.set_xlabel('[ms]')
                            axs.set_ylabel(comparison+' accuracy')
                            
                            if examine=='attendignore':
                                axs.set_title(horizontals[cx])
                        
    
        

            if doplot:
                if not publish:
                    fig.suptitle(dn+', %s, %s, %d/%d neurons, timecourse of differentiability metrics:\n'%(examine,layeredtitle,len(lmask),len(block.annotations['celltypes']))+\
                                 'channelwise maxnorm rate diff. (green),\n '+\
                                 'logreg. short %dms window decoder (blue+orange),\n'%width.magnitude+\
                                 'decision boundary direction to principal components (purple to green)')
                    
                    save = 0 or globalsave
                    if save:
                        fig.savefig(resultpath+'difference+decoder+angles_all-VACC-%s_%s_%s-w%dms-%dms%dms_%s'%(examine,layeredfilename,continuous_method,width.magnitude,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                
                else: # publication to CCN2019
    #                fig.suptitle(dn)
                    
                    save = 0 or globalsave
                    if save:
                        if examine=='allexpcond':
                            fig.savefig(resultpath+'3-%s-dec,VACC-timecourse'%dn+ext)
                        elif examine=='attendignore':
                            fig.savefig(resultpath+'2A-%s-dec,attendignore-timecourse'%dn+ext)










def decoder_celltypes(dn,block):

    recalculate = 0 or globalrecalculate
    doplot = 1

    width = 50*pq.ms

    blv,bla = preprocess.getorderattended(dn)
    comparisongroup  = [    [  blv,  [],[] ],   [ bla,   [], [] ]   ]
    comparison = 'context'
    celltypes = ['narrow','broad']

    # create list of neurons:
    print('neuron list: ', block.annotations['celltypes'])     # 0 narrow spiking, 1 broad spiking
    celllists = [ block.annotations['celltypes']==k for k in [0,1] ]

    # get responses of all neurons, select from them later
    acrossresponses_full = preprocess.collect_stimulusspecificresponses(block,comparisongroup)


    if doplot:
        fig,ax = plt.subplots(1,2,figsize=(2*10,1*10))
    
    
    for idx,celllistmask in enumerate(celllists):
        if True:#sum(celllistmask)>=8:
            print('%s: calculating for %d %s cells.'%(dn,sum(celllistmask),celltypes[idx]))
            if sum(celllistmask)==0: continue         # if there is not of this type, leave
            acrossresponses = []
            for aix in [0,1]:
                acrossresponses.append([  trial[:,celllistmask]  for trial in acrossresponses_full[aix]  ])


            if recalculate:
                acrossdecoder = nedi.get_responsedecoding(acrossresponses,width=width)
                pickle.dump(acrossdecoder,open(cacheprefix+'continuous/celltyperestricted-%s,%s-%s-%s-.pck'%(comparison,celltypes[idx],dn,continuous_method),'wb'))
            else:
                acrossdecoder = pickle.load(open(cacheprefix+'continuous/celltyperestricted-%s,%s-%s-%s-.pck'%(comparison,celltypes[idx],dn,continuous_method),'rb'))

            if doplot:
                axs = ax[idx]
                figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs)       # plot the performance
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.5)
                figs.setxt(axs)
                axs.set_ylabel('context accuracy')
                axs.set_title('%s %d/%d'%(celltypes[idx],sum(celllistmask),len(celllistmask)))
    if doplot:
        fig.suptitle(dn+'    context decoder only from selected celltypes')
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'celltyperestricted,%s-%s'%(comparison,dn)+ext)










def comparisonacrossconditions_choice(dn,block):
    
    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot
    

    width = 50*pq.ms               # larger input sizes do not help


    blv,bla = preprocess.getorderattended(dn)
    
    contextual_blocks = [blv[1],bla[1]]

    taskaspects = ['choice,av','choice,aa']

    print('contextual choice')
    for cx,comparison in enumerate(taskaspects):
        print(cx,comparison)
        acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn,\
                               onlyinblock=[contextual_blocks[cx]])
        if recalculate:
            acrossdecoder = nedi.get_responsedecoding(acrossresponses,width=width)
            pickle.dump(acrossdecoder,open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'wb'))
        else:
            acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
    

    
    
    
    return


    


def angles_between_decoders(dn,block):
    
    # display the trajectory of the angles between the decoder DBNVs


    doplot = 0 or globaldoplot
    dump = 1


    n_neuron = block.segments[0].analogsignals[0].shape[1]


    taskaspects = ['visual','audio','context','choice']
    # taskaspects = ['context','choice','choice,av','choice,aa']
    # taskaspects = ['visual','audio','context','choice','choice,av','choice,aa']
    
    # find the decision normal vectors
    c_db = []          # vector components (coefficients) of the decision normals in the activity space

    for cx,comparison in enumerate(taskaspects):

        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))

        numsignal = 7       # normal non-cross decoder
        wx = int((len(acrossdecoder)-numsignal)/n_neuron)     
        c_db.append(  np.reshape(np.array(acrossdecoder[numsignal:]), (wx,n_neuron,acrossdecoder[numsignal].shape[0],acrossdecoder[numsignal].shape[1]) ).mean(axis=0)    )
            
    c_db = np.array(c_db) # [{v,a,cx,ch},neurons,trajectory,stats]
    
    n_trajectory = c_db.shape[2]
    # for the fixed comparisons, how much to look for before and after stim onset
    depth_idx = int((100*pq.ms/T['dt']).magnitude)
    shift_idx = int((300*pq.ms/T['dt']).magnitude)

    # print(c_db.shape, '[{v,a,cx,ch},neurons,trajectory,stats]')

    # find angles for each timepoint
    #              to other dbnv-s (angles)
    #              to pre, and post same dbnv (tofixed)
    #              each t to t angle matrix same dbnv (highres)
    angles = np.zeros((len(taskaspects),len(taskaspects),n_trajectory))
    angles_tofixed = np.zeros((3,len(taskaspects),n_trajectory))
    angles_highres = np.zeros((len(taskaspects),n_trajectory,n_trajectory))
    for cxr,cr in enumerate(taskaspects):

        # first aspect to its own fixed angle
        fixed_pre = np.mean(c_db[cxr,:,T['stimstart_idx']-depth_idx:T['stimstart_idx'],0],axis=1)
        fixed_on = np.mean(c_db[cxr,:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,0],axis=1)
        fixed_lateon = np.mean(c_db[cxr,:,T['stimstart_idx']+shift_idx:T['stimstart_idx']+shift_idx+depth_idx,0],axis=1)
        for t in range(n_trajectory):
            for px,fixed in enumerate((fixed_pre,fixed_on,fixed_lateon)):
                aux = c_db[cxr,:,t,0].T @ fixed
                aux = aux / np.linalg.norm(c_db[cxr,:,t,0]) / np.linalg.norm(fixed)
                angles_tofixed[px,cxr,t] = np.arccos(aux)
            
            # with high resolution at each timepoint to each timepoint
            for t_c in range(n_trajectory):
                aux = c_db[cxr,:,t_c,0].T @ c_db[cxr,:,t,0]
                aux = aux / np.linalg.norm(c_db[cxr,:,t_c,0]) / np.linalg.norm(c_db[cxr,:,t,0])
                angles_highres[cxr,t,t_c] = np.arccos(aux)
                
        
        # now to the cross angles between aspects
        for cxc,cc in enumerate(taskaspects):
            for t in range(n_trajectory):
                aux = c_db[cxr,:,t,0].T @ c_db[cxc,:,t,0]
                aux = aux / np.linalg.norm(c_db[cxr,:,t,0]) / np.linalg.norm(c_db[cxc,:,t,0])
                angles[cxr,cxc,t] = np.arccos(aux)
    angles = angles * 180/np.pi
    angles_tofixed = angles_tofixed * 180/np.pi
    angles_highres = angles_highres * 180/np.pi


    if dump:
        # saving angles_highres   as (taskvariables,n_trajectory,n_trajectory) as angles
        pickle.dump(angles_highres,open(cacheprefix+'subspaces/angles,highres,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'wb'))
        pickle.dump(angles,open(cacheprefix+'subspaces/angles,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'wb'))


    # now we don't need to average over the trajectory, but rather show the angles at all timepoints:


    times = block.segments[0].analogsignals[0].times[:n_trajectory]
    fixed_colors = ['dodgerblue','purple','darkred']

    

    if doplot and 1:
        fig,ax = plt.subplots(len(taskaspects),len(taskaspects),figsize=(8*len(taskaspects),8*len(taskaspects)))
        
        
        for cxr,cr in enumerate(taskaspects):
            for cxc,cc in enumerate(taskaspects):
                axs = ax[cxr,cxc]
    
                if cxr==cxc:
                    # axs.axis('off')
                    # continue
                    for fx in range(angles_tofixed.shape[0]):
                        axs.plot(times,angles_tofixed[fx,cxr,:],color=fixed_colors[fx])
                    axs.legend(['pre -100ms..0ms','on 0..100ms','on 300..400ms'],frameon=False)
                
                else:
                    axs.plot(times,angles[cxr,cxc,:])
                
                axs.set_xlim([times[0],times[-1]])
                
                axs.set_ylim(0,180)
                axs.set_yticks(np.arange(0,181,30))
                
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.setxt(axs)
                figs.plottoaxis_chancelevel(axs,ch=90)
                
                # axs.set_title(cr+' $\cdot$ '+cc)
                axs.set_title(cc)
                axs.set_ylabel(cr)
    
        fig.suptitle('%s, angles between DBNVs'%dn)
        
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'angles,DBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn)+ext)
    



    
    if doplot and 1:
        fig,ax = plt.subplots(1,len(taskaspects),figsize=(11*len(taskaspects),8*1))
        
        
        for cx,aspect in enumerate(taskaspects):
                axs = ax[cx]
                
                cmap = figs.getcorrcolormap('correlation')
                cf = axs.pcolormesh(angles_highres[cx,:,:],vmin=0,vmax=180,cmap=cmap)
                
                axs.set_aspect('equal')
                
                ticks=[150,450]
                ticklabels=['0 ms','3000 ms']
                axs.set_xticks(ticks)
                axs.set_xticklabels(ticklabels)                
                axs.set_yticks(ticks)
                axs.set_yticklabels(ticklabels,rotation=90)                
                
                axs.set_title(aspect)

                # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                fig.colorbar(cf,ax=axs,ticks=np.arange(0,181,30))

    
        # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        # fig.colorbar(cf,ax=plt.axes([0.85, 0.1, 0.075, 0.8]))

        fig.suptitle('%s, angles within DBNVs between different timepoints'%dn)
        
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'angles,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn)+ext)
    
    










def decoder_driftsearch(dn,block):
    # looks for potential drift artefacts _after_ spike sorted data; can guide curation iteratively
    # https://github.com/MouseLand/Kilosort2/wiki/3.-More-on-drift-correction
    width = 500*pq.ms             # this is the width for the deocders
    filtercutoff = 90*pq.min
#    wx = int((width/T['dt']).magnitude)        # this is the width of time averaging
    wx = int((1500*pq.ms/T['dt']).magnitude)
    ma_l = 30
    blv,bla = preprocess.getorderattended(dn)
    # in the order they appear in the session, then in visual audio order
    comparisongroups  = [ \
                          [   [ [1,2],   [],    [] ], [ [3,4],  [],     [] ] ],\
                          [   [ [2],   [],    [] ], [ [4],  [],     [] ] ],\
                          [   [ [blv[1]],   [],    [] ], [ [bla[1]],  [],     [] ] ],\
                        ]
    taskaspects = ['context']
    

    design = False         # use the time-averaged pre for filter design; otherwise use the preprocess.removedrift() output
    uselogodds = False
    standardize = True
    remove_cells_driftedtonoise = True
    
    
    

    for cx,comparison in enumerate(taskaspects):

        print(dn,' drift search in %s...'%comparison)
    
        # train and test with CV and get single trial probabilities
        across_responses = preprocess.collect_stimulusspecificresponses(block, comparisongroups[cx])
        across_responses_multionly = preprocess.collect_stimulusspecificresponses(block, comparisongroups[1])
        classnumbers = [len(across_responses[0]),len(across_responses[1])]
        blocknumbers = [len(across_responses_multionly[0]),len(across_responses_multionly[1])]
        blocknumbers = [ classnumbers[0]-blocknumbers[0], blocknumbers[0], classnumbers[1]-blocknumbers[1], blocknumbers[1] ]
        bn_c = np.r_[0,np.cumsum(blocknumbers)]

        if 1:
            # get the cross validated probabilities
            acrossdecoders = []
            for preon in ['pre','on']:
                acrossdecoders.append(nedi.get_decoderprobability(across_responses,T=T,width=width,preon=preon))
            
            pickle.dump(acrossdecoders,open(cacheprefix+'phys/driftsearchdecodes,%s_%s-%s.pck'%(comparison,dn,continuous_method),'wb'))
        else:
            acrossdecoders = pickle.load(open(cacheprefix+'phys/driftsearchdecodes,%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))



#    print(responses.shape,len(responses),type(responses[15]))
#    fr_pre = []; fr_on = []
    responses = np.concatenate(across_responses)
    fr_pre = responses[:,T['stimstart_idx']-wx:T['stimstart_idx'],:].mean(axis=1)
    fr_on  = responses[:,T['stimstart_idx']:T['stimstart_idx']+wx,:].mean(axis=1)






    
    # collect original raw Hz firing rates (mainly for design purposes):
    if standardize==False:
        block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=False,recalculate=False)
    
    
    
    
    
    # filter the raw data and find the noise-drifter cells held in 'mask' as 0, and mask_lsq for the least squares estimate
    block_filt,mask,mask_lsq = preprocess.removedrift(dn,T,filtercutoff=filtercutoff)
    
    
    
    
#    across_responses_filt = preprocess.collect_stimulusspecificresponses(block_filt, comparisongroups[cx])    
#    across_responses_filt = across_responses.copy()      # don't filter, just remove the unwanted cells below
    print('Mask',len(mask),np.sum(mask),mask)
        
    if remove_cells_driftedtonoise:
        ch_remove = np.where(1-mask)[0]+1          # channels here are numbered according to the silicon probe
        # remove them from both original and filtered data
        block_rm = neph.removechannelsfromsession(block,ch_remove)
        block_filt_rm = neph.removechannelsfromsession(block_filt,ch_remove)
        
    # give which version of filtering and removing to use:  block, block_filt, block_rm, block_filt_rm
    block_final = block_filt_rm
        


    across_responses_filt = preprocess.collect_stimulusspecificresponses(block_final, comparisongroups[0])    


    
    
    if 1:
        # get the cross validated probabilities
        acrossdecoders_filt = []
        for preon in ['pre','on']:
            acrossdecoders_filt.append(nedi.get_decoderprobability(across_responses_filt,T=T,width=width,preon=preon))
        
        pickle.dump(acrossdecoders_filt,open(cacheprefix+'phys/driftremoveddecodes,%s_%s-%s.pck'%(comparison,dn,continuous_method),'wb'))
    else:
        acrossdecoders_filt = pickle.load(open(cacheprefix+'phys/driftremoveddecodes,%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))
        
        
        
        

    responses_filt = np.concatenate(across_responses_filt)

    # use the filtered data file to estimate the filtered average prestimulus firing rates
    fr_pre_filt = responses_filt[:,T['stimstart_idx']-wx:T['stimstart_idx'],:].mean(axis=1)
    

    n_trials, n_neurons = fr_pre.shape
    n_neurons_filt = fr_pre_filt.shape[1]
    neuronsinplay = np.zeros((n_neurons,5))


#    u_Px_pre = np.zeros((int((n_trials/2)+1),n_neurons))
#    u_Px_pre_filt = np.zeros((int((n_trials/2)+1),n_neurons))
    
    # old cell order with just end - start
#    n_compare_trials = 10
#    cellorder = np.argsort( np.abs( ( fr_pre[:n_compare_trials,:] - fr_pre[-n_compare_trials:] ).mean(axis=0))/fr_pre.std(axis=0) )



    # frequency spectrum
    cut_trials = 20  #   high pass filter how many trials above which slow waves to remove
    sampling_trials = 1 #       one trial sampling
    for n in range(n_neurons_filt):

        if design:
        # this part is for design; alternatively you should use the actual filtered raw data based calculations for fr_pre_filt
    
            if n==0: fr_pre_filt = np.zeros((n_trials,n_neurons))    # initialize
            fr_pre_filt[:,n] = neph.removefrequencies(fr_pre[:,n],[1/cut_trials],1/sampling_trials,filtertype='hp')
            if n==0:
                _,f = sp.signal.welch(x=fr_pre[:,n],fs=1/sampling_trials)
                u_Px_pre = np.zeros((len(f),n_neurons))
                u_Px_pre_filt = np.zeros((len(f),n_neurons))
                
            u_fx_pre, u_Px_pre[:,n] = sp.signal.welch(x=fr_pre[:,n],fs=1/sampling_trials)
            u_fx_pre_filt, u_Px_pre_filt[:,n] = sp.signal.welch(x=fr_pre_filt[:,n],fs=1/sampling_trials)
        
        

        # establish criteria for excluding cells
#        v = np.abs(np.std(fr_pre_filt[:bn_c[2],n]) - np.std(fr_pre_filt[bn_c[2]:,n]))        # compare variance around the detrended mean to detect zero drift
#        v = v / ( np.std(fr_pre_filt[:bn_c[2],n]) + np.std(fr_pre_filt[bn_c[2]:,n]) ) * 2
            v = np.abs(fr_pre_filt[:,n])
            l = sp.stats.linregress(np.arange(n_trials)+1,v)
            neuronsinplay[n,:] = l   # p value for regression of standard deviations around the moving mean to detect off-on cells
    if not design: neuronsinplay = mask_lsq


    



    # find trends with least squares    
    lsq = np.zeros((n_neurons,5,5))      # neurons,{1,2,3,4,all-blocks},lsqsparams
    lsq_filt = np.zeros((n_neurons,5,5))    # lsq post filter


    n_filt = -1
    for n in range(n_neurons):
        if mask[n]: n_filt += 1
        for c in [0,1,2,3,4]:
            if c<4:
                l = sp.stats.linregress( np.arange(bn_c[c],bn_c[c+1])+1, fr_pre[ bn_c[c]:bn_c[c+1] , n ]  )
                if mask[n]:
                    l_filt = sp.stats.linregress( np.arange(bn_c[c],bn_c[c+1])+1, fr_pre_filt[ bn_c[c]:bn_c[c+1] , n_filt ]  )
            else:
                l = sp.stats.linregress( np.arange(n_trials)+1, fr_pre[ : , n ]  )
                if mask[n]:
                    l_filt = sp.stats.linregress( np.arange(n_trials)+1, fr_pre_filt[ : , n_filt ]  )
            lsq[n,c,:] = l
            if mask[n]:
                lsq_filt[n_filt,c,:] = l_filt



    cellorder = np.argsort(lsq[:,4,0] * (lsq[:,4,3]<0.05) )        # times if significant at all slope 1th, p 4th
    cellorder = np.arange(n_neurons)
    
    nsq = int(np.ceil(np.sqrt(n_neurons)))    # number of displays in one direction






                        






    


    
    plotlist = [0,1]
    if 1:
        isvisualfirst = int((blv[0]-1)/2)    # check if visual is the first two blocks
#        i = 1-isvisualfirst       

        if 0 in plotlist: fig,ax = plt.subplots(1,2,figsize=(32,8))            # decoders
        if 1 in plotlist: fig_rate,ax_rate = plt.subplots(nsq,nsq,figsize=(nsq*10,nsq*6))      # firing rates
        if 2 in plotlist: fig_freq,ax_freq = plt.subplots(nsq,nsq,figsize=(nsq*10,nsq*6))      # requency spectra, before and after
        
        colors = [ ['dodgerblue','blue','forestgreen','mediumseagreen'], ['navy','mediumblue','darkgreen','darkcyan']     ]
        labels = [ '%s'%['V','Va','A','Av'][i-1] for i in [ blv[0],blv[1],bla[0],bla[1] ] ]
        for preon in [0,1]:
        
            
            # first draw the log odds space decoder probabilities
            if 0 in plotlist: 
                axs = ax[preon]
                
                for dx,acrossdecoder in enumerate([acrossdecoders, acrossdecoders_filt]):
    
                    x = acrossdecoder[preon][:,1-isvisualfirst]         # display class probabilities as attend visual
                    if uselogodds:
                        x = np.log(x/(1-x))           # uncomment to show in log odds space
                    
                    x_ma = sp.signal.savgol_filter(x, 29, 5)
                    cli = [[0,1,2,3],[2,3,0,1]][isvisualfirst]
                    for c in [0,1,2,3]:
                        color = [colors[preon][cli[c]],['fuchsia','mediumvioletred','gold','darkorange'][cli[c]]][dx]
                        axs.plot(np.arange(blocknumbers[c])+bn_c[c],x[bn_c[c]:bn_c[c+1]],'-o',color=color,alpha=0.4, label=labels[c]+[' raw',' filt'][dx])
                        axs.plot(np.arange(blocknumbers[c])+bn_c[c],x_ma[bn_c[c]:bn_c[c+1]],'-',lw=3,color=color, alpha=0.6)
            
                
                axs.set_xticks(np.r_[1,bn_c[1:]])
                axs.set_xlim(0,len(x)+1)
                if not uselogodds: axs.set_ylim(-0.01,1.01)
    #            axs.set_ylim(-2,2)
                figs.plottoaxis_chancelevel(axs,0.5*(1-uselogodds))
                axs.set_title(['pre','on'][preon])
                axs.legend(frameon=False)
        






            
            # draw firing rates on another fig:
            if 1 in plotlist: 
                yoffs = -2*preon + 1
                n_filt = -1
                for n in range(n_neurons):
                    if mask[n]: n_filt += 1
                    axs = ax_rate[int(np.floor(cellorder[n]/nsq)),cellorder[n]%nsq]
                    
    #                x = [fr_pre,fr_on][preon][:,n]                        # draw rates with pre and on on the same subplot
                    x = [fr_pre,fr_pre_filt][preon][:,[n,n_filt][preon]]                   # draw pre and filtered pre on the same subplot
                    x_ma = sp.signal.savgol_filter(x, 29, 5)
    #                if n==0: label=; else: label=None
                    for c in [0,1,2,3]:
                        axs.plot(np.arange(blocknumbers[c])+bn_c[c],yoffs+x[bn_c[c]:bn_c[c+1]],'-',color=colors[preon][cli[c]],alpha=0.5)
                        axs.plot(np.arange(blocknumbers[c])+bn_c[c],yoffs+x_ma[bn_c[c]:bn_c[c+1]],'-',lw=3,color=colors[preon][cli[c]], label=['pre ','pre filt '][preon]+labels[c])
    
    
                    
                    if preon==0:
                        if mask[n]:
                            if mask_lsq[n,3]<1e-1: color = 'k'
                            else: color = 'darkgreen'
                        else: color = 'red'
                        axs.set_title('#%d/%d; p=%.1e'%(cellorder[n]+1,n+1,neuronsinplay[n,3]),color=color)
                        if cellorder[n]==0: axs.legend(frameon=False,ncol=2)
#                    else:
#                        if not mask[n]: continue
                    lsqn = [lsq[n],lsq_filt[n_filt]][preon]
                    for c in [0,1,2,3,4]:
                        if c<4:
                            triallist = np.array([bn_c[c],bn_c[c+1]])
                        else:
                            triallist = np.array([1,n_trials])
                        lt = ['-','--'][lsqn[c,3]>0.05]
                        axs.plot(triallist,yoffs+triallist*lsqn[c,0]+lsqn[c,1],lt,lw=3,\
                                 color=[['fuchsia','gold'],['mediumvioletred','darkorange']][preon][c==4])
                        
    
                    if preon==1:
                        axs.set_xticks(np.r_[1,bn_c[1:]])
                        axs.set_xlim(0,len(x)+1)
#                        axs.set_ylim(-3*fr_pre[:,n].std()*1.2-1,3*fr_pre[:,n].std()*1.2+1)
                        axs.set_ylim(axs.get_ylim())
                        figs.plottoaxis_chancelevel(axs,1)
                        figs.plottoaxis_chancelevel(axs,-1)
    
                for n in np.arange(n_neurons,nsq**2):
                    if preon==1:
                        axs = ax_rate[int(np.floor(n/nsq)),n%nsq]
                        fig_rate.delaxes(axs)
    

            
            
            
            
            
            
            
            # draw frequency spectrum detrending on yet another figure:
            if 2 in plotlist: 
                if preon==1: continue
                for n in range(n_neurons):
                    axs = ax_freq[int(np.floor(cellorder[n]/nsq)),cellorder[n]%nsq]
                    
                    axs.plot(u_fx_pre,u_Px_pre[:,n],lw=2,color='seagreen',alpha=0.8)
                    axs.plot(u_fx_pre,u_Px_pre_filt[:,n],lw=2,color='darkgreen',alpha=0.8)
    
                    axs.set_xticks(1/np.array([n_trials,n_trials/2,n_trials/4]))
                    axs.set_xticklabels(['%d'%n_trials,'%d'%int(n_trials/2),'%d'%int(n_trials/4)])
    #                axs.set_xlim(0,len(x)+1)
                    
                    axs.set_title('#%d, neuron id: %d'%(cellorder[n]+1,n+1))
                
    
                for n in np.arange(n_neurons,nsq**2):
                    if preon==0:
                        axs = ax_freq[int(np.floor(n/nsq)),n%nsq]
                        fig_freq.delaxes(axs)




        
        if 0 in plotlist: 
            if not uselogodds: fig.suptitle(dn+' stationarity of decoder class probabilities')
            else: fig.suptitle(dn+' stationarity of decoder class log odds')
        
        if 1 in plotlist: fig_rate.suptitle(dn+' stationarity of firing rates'+['',', standardized'][standardize])
        if 2 in plotlist: fig_freq.suptitle(dn+' stabilizing firing rates: high pass detrending')
        
        
        
        save = 0
        if save:
            if 0 in plotlist:
                if not uselogodds: fig.savefig(resultpath+'driftsearch,%dms_%s-%dms%dms_%s'%(width.magnitude,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                else: fig.savefig(resultpath+'driftsearch,logodds,%dms_%s-%dms%dms_%s'%(width.magnitude,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
            if standardize:
                if 1 in plotlist: fig_rate.savefig(resultpath+'driftsearch,firingrates,highpassdetrended,standardized,%dms_%s-%dms%dms_%s'%(width.magnitude,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                if 2 in plotlist: fig_freq.savefig(resultpath+'driftsearch,frequencies,highpassdetrending,standardized,%dms_%s-%dms%dms_%s'%(width.magnitude,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
            else:
                if 1 in plotlist: fig_rate.savefig(resultpath+'driftsearch,firingrates,highpassdetrended%dms_%s-%dms%dms_%s'%(width.magnitude,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                if 2 in plotlist: fig_freq.savefig(resultpath+'driftsearch,frequencies,highpassdetrending,%dms_%s-%dms%dms_%s'%(width.magnitude,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)








        






















    




def decoder_crosstest(dn,block):
    # compare decoder performance across spontaneous=stimulus off   versus   evoked=stimulus on    periods during trials
    # trains in stim on, test in stim off and vice versa
    # context is the main question, others are control or other questions depending on pre or post stimulus comparison

    # figs.rasterplot(block)
    
    recalculate = 0 or globalrecalculate
    doplot = 1 or globaldoplot

    # setup stimulus:
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [   [ [ [blv[1]], [45],    [5000] ], [ [blv[1]],[135],     [5000] ]  ], \
                            [ [ [bla[1]],  [45],[5000] ], [ [bla[1]],   [45],[10000] ]    ],    \
                            [ [ [blv[1]],  [45],[5000] ], [ [bla[1]],   [45],[5000] ]    ],     \
                            [[],[]] \
                        ]
    taskaspects = ['visual','audio','context','choice']
    
    windowwidth = 100*pq.ms

    # trainingoffsets = np.array([0,375,750,1225])*pq.ms
    # trainingoffsets = np.array([0, 500, 1000, 1400])*pq.ms
    trainingoffsets = np.array([-1500, -100, 0, 1500, ])*pq.ms
    # trainingpoints = T['starttime'] + trainingoffsets
    trainingpoints = trainingoffsets
    
    if doplot:
        fig,ax=plt.subplots(4,len(trainingpoints),figsize=(32,24))

    for cx,comparison in enumerate(taskaspects):
        # if cx not in [2]: continue
        for trx,trainingpoint in enumerate(trainingpoints):
       
            print(cx,trx,comparison,trainingpoint,comparisongroups[cx])
            
            stimulusIDgroups = comparisongroups[cx]

            # collect neural responses
            if not comparison=='choice': # visual, audio, context:
                offresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups)
            else:  # choice:
                offresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)

            # offresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups)

            
            
            
    
    
            # calculate crosstest over stimulus off and on

            if recalculate:
                offstimdecoder = nedi.get_offdecoding(offresponses,T,width=windowwidth,trainingpoint=trainingpoint)
                pickle.dump(offstimdecoder,open(cacheprefix+'continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,comparison,trx),'wb'))
            else:
                offstimdecoder = pickle.load(open(cacheprefix+'continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,comparison,trx),'rb'))
    
#
            print(offstimdecoder[0].times[0],offstimdecoder[0].times[-1])
            
        # plot metrics
            if doplot:
                axs = ax[cx,trx]
                figs.plottoaxis_offdecoderrocauc(offstimdecoder,axs,T)
                figs.plottoaxis_stimulusoverlay(axs,T)
                axs.fill_between([ offstimdecoder[0].times[0],offstimdecoder[0].times[-1] ],[0,1],[0,1],color='dodgerblue',alpha=0.9)
                figs.plottoaxis_chancelevel(axs,0.5,'chance level')
                if cx==0: axs.set_title('train %d-%d ms'%(trainingpoint,trainingpoint+windowwidth))
                if trx==0: axs.legend(); axs.set_ylabel(comparison+'\nroc auc at each timepoint')
                if cx==len(taskaspects)-1: axs.set_xlabel('time from trial onset [ms]')
            

    if doplot:
        fig.suptitle(dn+', timecourse of time-crosstests:\n'+\
                     'logistic regression, %d ms window decoder trained at %d points'%(windowwidth.magnitude,len(trainingpoints)))
        
        save = 0
        if save:
            fig.savefig(resultpath+'offon,crosstest_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)











def decoder_crosstest_highres(dn,block):
    # like crosstest, but training points are used for all timepoints
    # returns a matrix of trainpoints times testpoints
    # here using all possible multimodal trials (non-conditioned on other variables)

    # figs.rasterplot(block)
    
    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot

    # setup stimulus:
    blv,bla = preprocess.getorderattended(dn)
    # comparisongroups  = [   [ [ [blv[1]], [45],    [5000] ], [ [blv[1]],[135],     [5000] ]  ], \
    #                         [ [ [bla[1]],  [45],[5000] ], [ [bla[1]],   [45],[10000] ]    ],    \
    #                         [ [ [blv[1]],  [45],[5000] ], [ [bla[1]],   [45],[5000] ]    ],     \
    #                         [[],[]] \
    #                     ]
    # taskaspects = ['visual','audio','context','choice']

    # for both contexts
    # comparisongroups  = [   [ [ [2,4], [45],    [] ], [ [2,4],[135],     [] ]  ], \
    #                         [ [ [2,4], [],    [5000] ], [ [2,4],[],     [10000] ]  ], \
    #                         [ [ [blv[1]],  [],[] ], [ [bla[1]],   [],[] ]    ],     \
    #                         [[],[]] \
    #                     ]
    # taskaspects = ['visual','audio','context','choice']
    taskorder = [0,1,2,3]
    exptyp = ''
    


    # for separately conditioning on the two contexts
    comparisongroups  = [ \
                            [ [ [blv[1]], [45],    [] ],   [ [blv[1]],[135],     [] ]  ], \
                            [ [ [blv[1]], [],    [5000] ], [ [blv[1]],[],     [10000] ]  ], \
                            [[blv[1]],[]], \
                            [ [ [bla[1]], [45],    [] ],   [ [bla[1]],[135],     [] ]  ], \
                            [ [ [bla[1]], [],    [5000] ], [ [bla[1]],[],     [10000] ]  ], \
                            [[bla[1]],[]] \
                        ]
    taskaspects = ['visual,av','audio,av','choice,av',\
                   'visual,aa','audio,aa','choice,aa']
    taskorder = [0,3,1,4,2,5]
    exptyp = ',avaa'



    windowwidth = 10*pq.ms # sliding for train test indexing
    width = 50*pq.ms       # feature width of neural activity

    # specific training points
    # trainingoffsets = np.array([0,375,750,1225])*pq.ms
    # trainingoffsets = np.array([0, 500, 1000, 1400])*pq.ms
    # trainingoffsets = np.array([-1500, -100, 0, 1500, ])*pq.ms

    # sliding window for all points
    trainingoffsets = np.arange(-1500*pq.ms,4500*pq.ms-width/1.1,windowwidth)*pq.ms
    # trainingoffsets = np.arange(4400*pq.ms,4500*pq.ms-width/1.1,windowwidth)*pq.ms
    trainingpoints = trainingoffsets
    

    if recalculate:
        crossdecoder_matrix_allaspects = [] #  (taskaspects)(trainingpoints)(testtimecourse,stats)
        for cx,comparison in enumerate(taskaspects):
            print(dn,cx,comparison,comparisongroups[cx],trainingpoints[0],':',windowwidth,':',trainingpoints[-1],'n tr:',len(trainingpoints))
            
            
            crossdecoder_matrix = []
            for trx,trainingpoint in enumerate(trainingpoints):
           
                
                stimulusIDgroups = comparisongroups[cx]
    
                # collect neural responses
                if not comparison[:6]=='choice': # visual, audio, context:
                    offresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups)
                else:  # choice:
                    offresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn,onlyinblock=stimulusIDgroups[0])
    
                # offresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups)
        
                # calculate crosstest over timecourse
                crossdecoder = nedi.get_offdecoding(offresponses,T,width=width,windowwidth=windowwidth,trainingpoint=trainingpoint)
                crossdecoder_matrix.append(crossdecoder)   # only the test is returned, (timepoints,{mean,2std,2sem});   to check the training score, get it from other deocders

            crossdecoder_matrix_allaspects.append(crossdecoder_matrix)

        crossdecoder_matrix_allaspects = np.array( crossdecoder_matrix_allaspects )
        print('crossdecoder matrix shape:   (tasks)(traintimes)(testtimepoints)(stats)   ',crossdecoder_matrix_allaspects.shape)
        # pickle.dump(crossdecoder_matrix_allaspects,open(cacheprefix+'continuous/crossdecodematrix,cv-%s-%dms_%s.pck'%(continuous_method,T['dt'],dn),'wb'))
        pickle.dump(crossdecoder_matrix_allaspects,open(cacheprefix+'continuous/crossdecodematrix,contextcond,cv-%s-%dms_%s.pck'%(continuous_method,T['dt'],dn),'wb'))
    else:
        # crossdecoder_matrix_allaspects = pickle.load(open(cacheprefix+'continuous/crossdecodematrix,cv-%s-%dms_%s.pck'%(continuous_method,T['dt'],dn),'rb'))
        crossdecoder_matrix_allaspects = pickle.load(open(cacheprefix+'continuous/crossdecodematrix,contextcond,cv-%s-%dms_%s.pck'%(continuous_method,T['dt'],dn),'rb'))
        print('crossdecoder matrix shape:   (tasks)(traintimepoints)(testtimepoints)(stats)   ',crossdecoder_matrix_allaspects.shape)
    # print(len(crossdecoder_matrix_allaspects),len(crossdecoder_matrix_allaspects[0]),crossdecoder_matrix_allaspects[0][0].shape)

    # crossdecoder_matrix_allaspects = crossdecoder_matrix_allaspects.swapaxes(1,2)
    
    n_neurons = block.segments[0].analogsignals[0].shape[1]

    if doplot:
        # fig,ax = plt.subplots(2,4,figsize=(12*4,20*1))
        fig,ax = plt.subplots(2,len(taskaspects),figsize=(8*len(taskaspects),8*2))
        
        
        for cx,aspect in enumerate(taskaspects):
                axs = ax[0,cx]
                
                cmap = figs.getcorrcolormap('correlation')
                cf = axs.pcolormesh(crossdecoder_matrix_allaspects[taskorder[cx],:,:,0],vmin=0,vmax=1,cmap=cmap)
                
                axs.set_aspect('equal')
                
                ticks=[150,450]
                ticklabels=['0 ms','3000 ms']
                axs.set_xticks(ticks)
                axs.set_xticklabels(ticklabels)                
                axs.set_yticks(ticks)
                axs.set_yticklabels(ticklabels,rotation=90)                
                
                axs.set_title(taskaspects[taskorder[cx]])
                
                axs.set_xlabel('test timecourse')
                axs.set_ylabel('train timecourse')

                # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                fig.colorbar(cf,ax=axs,ticks=np.arange(0,1.1,0.25))
                
                
                
                
                # now the same with surface plot:
                ax[1,cx].remove()
                ax[1,cx] = fig.add_subplot(2,len(taskaspects),len(taskaspects)+1+cx,projection='3d')
                axs = ax[1,cx]

                skip = 50
                smoother = np.ones((skip,skip)); smoother /= np.sum(smoother)
                smoothz = sp.signal.convolve2d(crossdecoder_matrix_allaspects[taskorder[cx],:,:,0],smoother,mode='valid')
                trainlength,testlength = smoothz.shape

                # trainlength = len(crossdecoder_matrix_allaspects[0,:,0,0])
                # testlength = len(crossdecoder_matrix_allaspects[0,0,:,0])
                xmesh,ymesh = np.meshgrid(np.arange(0,testlength),np.arange(0,trainlength))

                axs.plot_surface(xmesh,ymesh,smoothz, vmin=0,vmax=1,cmap=cmap, edgecolor='none')
                # axs.plot_wireframe(xmesh,ymesh,smoothz)

                axs.view_init(elev = 30, azim=60+180)

                # axs.invert_xaxis()
                # axs.invert_yaxis()
                ticks=[150-skip,450-skip]
                axs.set_xticks(ticks); axs.set_xticklabels(ticklabels)
                axs.set_yticks(ticks); axs.set_yticklabels(ticklabels)
                axs.set_zticks(np.arange(0,1.1,0.25))

            
                axs.set_zlim(0.4,1)

                axs.set_xlabel('test'); axs.set_ylabel('train');axs.set_zlabel('accuracy')


        fig.suptitle(dn+', %d neurons, timecourse of time-crosstests:\n'%n_neurons+\
                      'logistic regression 5 x cv, %d ms window decoder trained and tested at each timepoint'%(width.magnitude)+\
                       '\ncv accuracy: colors and elevation, spike sorting: %s'%continuous_method   )
        
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'crossdecoder,cv_matrix%s_%s-%dms%dms_%s'%(exptyp,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)



    return











def decoder_crosstestgradient(dn,block):
    # gradient over single modality block to show switching from one context to the other:
    # initial (first single modality) and transition (second single modality) blocks
    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot
    
    # load training sessions behaviour data, and success/error estimates by session-averaging
    print(dn)
    _,_,_,_,cuetrainig,cue_est = preprocess.loadtrainingbehaviouraldata(dn,recalculate=True)

    

#        figs.rasterplot(block)

    # setup stimulus:
    blv,bla = preprocess.getorderattended(dn)
#    comparisongroups_train  = [ [ [blv[1]],  [45],[5000] ], [ [bla[1]],   [45],[5000] ]    ] # go signals only
    comparisongroups_train  = [ [ [blv[1]],  [],[] ], [ [bla[1]],   [],[] ]    ] # all
    windowwidth = 100*pq.ms

    trainingoffsets = np.array([0, 500, 1000, 1400])*pq.ms
    trainingpoints = T['starttime'] + trainingoffsets

    timeaveragewidth = 100*pq.ms
#    timeaverage_range_off = ((np.array([windowwidth+T['dt'], 2*windowwidth+T['dt']])/T['dt']).magnitude).astype('int16')
#    timeaverage_range_on = ((np.array([T['stimstarttime']-T['starttime'], T['stimstarttime']+windowwidth-T['starttime']])/T['dt']).magnitude).astype('int16')
    timeaverage_range = ((np.array([windowwidth+T['dt'], windowwidth+timeaveragewidth+T['dt']])/T['dt']).magnitude).astype('int16')
#    print('ta-on',timeaverage_range)

    # collect neural responses
    # train in multimodal blocks
    stimulusIDgroups = comparisongroups_train
    offresponses_train = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups)
    
    # test in single modality cueing blocks
    offresponses_test = [ [  trial.analogsignals[0] for trial in block.segments
                           if trial.annotations['block']==blv[0]   ],\
                          [  trial.analogsignals[0] for trial in block.segments
                           if trial.annotations['block']==bla[0]   ]   ]
#    classes_test = np.array(\
#                        [ [    trial.annotations['visual']==135  for trial in block.segments 
#                        if trial.annotations['block']==blv[0]     ] ,\
#                    [    trial.annotations['audio']==10000  for trial in block.segments 
#                        if trial.annotations['block']==bla[0]     ]         ]       )
    classes_test = [np.zeros(   (len(offresponses_test[0])),dtype='int16'   ), np.ones(   (len(offresponses_test[1])),dtype='int16'   )]

    print(np.array(offresponses_test).shape, np.array(classes_test).shape)
#    fig,ax=plt.subplots(len(offresponses_test[0]),len(trainingpoints),figsize=(36,4*len(offresponses_test[0])))
#    print((classes_test))

    if 1:
        for bl in [0,1]:
            n_cuetrials = len(offresponses_test[bl])
    
            side_n_cuetrials = int(np.ceil(np.sqrt(n_cuetrials)))
#            fig,ax=plt.subplots(side_n_cuetrials,side_n_cuetrials,figsize=(6*side_n_cuetrials,6*side_n_cuetrials))
        
                
            for trx,trainingpoint in enumerate(trainingpoints):
                # if trx not in [3]: continue
                print(trx,trainingpoint,comparisongroups_train)
        
                # calculate crosstest over stimulus off and on
        
                if recalculate:
                    offstimdecoderlist = nedi.get_offgradientdecoding(offresponses_train,offresponses_test[bl],classes_test[bl],T,width=windowwidth,trainingpoint=trainingpoint)
                    pickle.dump(offstimdecoderlist,open(cacheprefix+'continuous/offstimgradientdecodes,allstim,b%d-%s-%s-t%d.pck'%(bl,dn,continuous_method,trx),'wb'))
                else:
                    offstimdecoderlist = pickle.load(open(cacheprefix+'continuous/offstimgradientdecodes,allstim,b%d-%s-%s-t%d.pck'%(bl,dn,continuous_method,trx),'rb'))
        
                
            # plot metrics
#                for sx,trialdecoder in enumerate(offstimdecoderlist):      #[offstimdecoderlist[0], offstimdecoderlist[1],offstimdecoderlist[10],offstimdecoderlist[11], offstimdecoderlist[20],offstimdecoderlist[21] ]):
#                    j=int(np.floor(sx/side_n_cuetrials))
#                    k=sx%side_n_cuetrials
#                    axs = ax[j,k]
#                    figs.plottoaxis_offdecoderrocauc(trialdecoder,axs,T)
#                    figs.plottoaxis_stimulusoverlay(axs,T)
#                    figs.plottoaxis_chancelevel(axs,0.5,'chance level')
#                    if k==0 and j==0: axs.legend(); axs.set_ylabel('probability at each timepoint')
#                    if j==side_n_cuetrials-1: axs.set_xlabel('time from trial onset [ms]')
#                    axs.set_title('trial %d'%(sx+1))
#                    axs.set_ylim([0,1])
#                    
#        
#            fig.suptitle(dn+', timecourse of differentiability metrics during single modality cue block:\n'+\
#                         'logreg. long window decoder trained only off stimulus')
#            
#            save=0
#            if save:
#                fig.savefig(resultpath+'offon,crosstestgradient,allcuetrials_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
        
    

    if doplot:

        bp = preprocess.assessbehaviouralperformancefullarray(dn)
        bp_vals = np.array([1,-1,1,-1])

        if not publish:
            fig,ax=plt.subplots(4,2,figsize=(32,28))
        else:
            fig,ax=plt.subplots(2,2,figsize=(24,16))           # this is for publication without performance
        
        titles = ['pre stimulus','on stimulus']
        colorlists = [ ['darkorange','chocolate'], ['tomato','firebrick'] ]
        perfcolorlist = [ 'seagreen','orange','mediumblue','orangered'   ]   # hit miss correct rejection false alarm
        

        # create the cross gradient file: [ block (in original block order) ],
        #                                 ( {pre,on},trials,{dtm,dtmse,savgol} )
        export_crossgradient = []
        for bl in [0,1]:
            blockorder = int(([blv[0],bla[0]][bl]-1)/2)
            for ptrx in [0,1]:# ,trainingpoint in enumerate(trainingpoints):
                dtm,dtstd,dtsem,n_cuetrials = nedi.get_crossgradient(ptrx,bl,dn,timeaverage_range,continuous_method)
                dtm_c_s = sp.signal.savgol_filter(dtm[:,0], 29, 5)    # window width, and polynomial order
                if ptrx==0:
                    export_crossgradient.append(       np.zeros((2,n_cuetrials,3))    )
                export_crossgradient[bl][ptrx,:,0] = dtm[:,0]
                export_crossgradient[bl][ptrx,:,1] = dtsem[:,0]
                export_crossgradient[bl][ptrx,:,2] = dtm_c_s

                axs = ax[ ptrx+3-1, blockorder]     # display by block order, not visual and then audio
                figs.plottoaxis_offgradientdecoder(n_cuetrials,dtm,dtsem,dtm_c_s,axs,colorlist=colorlists[ptrx])
                figs.plottoaxis_chancelevel(axs,ch=0.5,label='chance')    # +', %d ms'%trainingpoints[ptrx].magnitude
                axs.set_yticks([0,0.5,1])
                axs.set_xlim([0.1,n_cuetrials+0.9])
                if blockorder==0: axs.set_ylabel(titles[ptrx]+'\nlogistic regr. context probability')
                if ptrx==0:
                    if bl==0: axs.set_title('attend to visual, cue block %d'%(blv[0]))
                    elif bl==1: axs.set_title('attend to audio, cue block %d'%(bla[0]))


            if publish:  continue
            # display recording session performance:    hit miss correct rejection false alarm
            axs = ax[0, blockorder]
            bpx = bp[blockorder*2][0]
            figs.plottoaxis_performancebarplot(n_cuetrials,bp_vals,bpx,axs,perfcolorlist)
            axs.set_xlim([0.1,n_cuetrials+0.9])
            if blockorder==0: axs.set_ylabel('performance\nhit (green) - miss (orange)\ncorr.rej.(blue) - false al.(red) ')
            if blockorder==1: axs.set_xlabel('trial number')
            if bl==0: axs.set_title('attend to visual, cue block %d'%(blv[0]))
            elif bl==1: axs.set_title('attend to audio, cue block %d'%(bla[0]))
#            axs.legend([ ])
            
            
            # display training session performance
            axs = ax[1, blockorder]
            axs.plot(np.arange(30)+1,cue_est[bl,0,:],color='darkblue')
            axs.fill_between(np.arange(30)+1,cue_est[bl,0,:]-cue_est[bl,2,:],cue_est[bl,0,:]+cue_est[bl,2,:],\
                             color='darkblue',alpha=0.2)
            axs.set_ylim(-0.04,1.04)
            axs.set_yticks([0,0.5,1])
            axs.set_title('cueing during training session, cue block %d: %s'%(bl*2+1,['initial','transition'][bl]))
            axs.set_ylabel('session average p')

            
    #        ax.legend('raw','smoothed')
        
        
        fig.suptitle(dn+', log.lin. classifier, offstimulus-trained on multimodality blocks (2,4)\ntimecourse-averaged singletrial probabilities in preceding single modality cue blocks (1,3)')


        if 1:
            pickle.dump(export_crossgradient,open(cacheprefix+'continuous/offstimgradient,cue,decodeprobs,savgol_%s-%s.pck'%(dn,continuous_method),'wb'))

        save = 1
        if save:
            if not publish:
                fig.savefig(resultpath+'offon,crosstestgradient,trainingbehaviour,%dms_%s-%dms%dms_%s'%(timeaveragewidth.magnitude,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
            else:
                fig.savefig(resultpath+'6-%s-cuetrials,crosstest'%dn+ext)





















def decoder_crosscontexttest(dn,block,retvar=False):                     #trainidx=None,testidx=None

    # test visual discrimination ability in the two contexts
    # test relevant and test irrelevant stimulus modality in terms of the context

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot
    
    width = 50*pq.ms
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [ \
                          [   [ [blv[1]],   [45],    [] ], [ [blv[1]],  [135],     [] ] ],\
                          [   [ [bla[1]],   [45],    [] ], [ [bla[1]],  [135],     [] ] ], \
                          [   [ [5],    [90,270],    [] ], [ [5],     [0,180],     [] ] ],\
                          [   [ [5],    [30,210],    [] ], [ [5],     [120,300],     [] ] ],\
                          [   [ [5],    [60,240],    [] ], [ [5],     [150,330],     [] ] ],\
                          [   [ [5],    [45,225],    [] ], [ [5],     [135,315],     [] ] ],\
                          [   [ [2,4],      [45],    [] ], [ [2,4],     [135],     [] ] ], \
                        ]
    taskaspects = ['attend visual','ignore visual',\
                   'OT block 0°-90°','OT block 30°-120°','OT block 60°-150°','OT block 45°-135°','all multimodal blocks 45°-135°']
    filename = ['crosscontext,attendattend','crosscontext,ignoreignore','crosscontext,attendignore','crosscontext,ignoreattend',\
                'crossorientation90','crossorientation30','crossorientation60','crossorientation45']
    labelprojections = [ ['ignore','attend'], ['ignore','attend'],\
                         ['  0°- 90°',' 45°-135°'], [' 30°-120°',' 45°-135°'], [' 60°-150°',' 45°-135°'],\
                         [' 45°-135°',' 45°-135°'] ]        # this last is the test: multimodal
    
    crosslabels = [['attend','attend'],['ignore','ignore'],['attend','ignore'],['ignore','attend']]
    preonlabels = ['PRE','ON']
    classlabels = [' 45°','135°']
    
    correctonly = 0      # use only correct behaviour trials, if True
    
    
    # print(trainidx,testidx)
    # if trainidx==None or testidx==None:
        # 2 symmetric, then 2 cross
        # trainidx = 0; testidx = 0; # train attend visual test attend visual
        # trainidx = 1; testidx = 1; # train ignore visual test ignore visual
        # trainidx = 0; testidx = 1; # train attend visual test ignore visual
        # trainidx = 1; testidx = 0; # train ignore visual test attend visual
        # print('None!!')
    
        # receptive field characterization block
        # trainidx = 2; testidx = 6       # train block 5 90,270/0,180 test all visual
        # trainidx = 3; testidx = 6       # train block 5 30,210/120,300 test all visual
        # trainidx = 4; testidx = 6       # train block 5 60,240/150,330 test all visual
        # trainidx = 5; testidx = 6       # train block 5 45,225/135,315 test all visual

    # if trainidx!=1 and testidx!=1: return
    
    
    # run through both crosstrain situation: trained on attend, and trained on ignore
    crossdecoders = []
    projections = []  #  {trainedattend,trainedignore} x {prestim,onstim} x {testedignore,testedattend} x {go,nogo}

    for trainedwhere in [0,1]:
        trainidx = trainedwhere
        testidx = 1-trainedwhere
    
    
        filenameidx = np.array([[0,2],[3,1]])[trainidx,testidx]
        numsignal = 7+int(trainidx!=testidx)
    
        print('Doing %s...'%filename[filenameidx])
    
        # train on ignore 75%, test on ignore 25% and attend 25%
        responses_trained = preprocess.collect_stimulusspecificresponses(block, comparisongroups[trainidx], correctonly=correctonly)
        responses_predicted = preprocess.collect_stimulusspecificresponses(block, comparisongroups[testidx], correctonly=correctonly)
        
    
            
        if recalculate:
            if trainidx!=testidx:    # cross
                crossdecoder = nedi.get_crosscontextdecoding(responses_trained,responses_predicted,width=width)
            else:                    # ident (same as in the initial sliding decoders)
                crossdecoder = nedi.get_responsedecoding(responses_trained,width=width)
            pickle.dump(crossdecoder,open(cacheprefix+'continuous/responsedecodes,%s%s_%s-%s.pck'%(filename[filenameidx],['',',correctonly'][correctonly],dn,continuous_method),'wb'))
        else:
            crossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,%s%s_%s-%s.pck'%(filename[filenameidx],['',',correctonly'][correctonly],dn,continuous_method),'rb'))
        
        crossdecoders.append(crossdecoder)
        print(len(crossdecoder),numsignal,np.array(crossdecoder[numsignal:]).shape,type(crossdecoder[numsignal]),crossdecoder[numsignal].shape)
    
    
        # calculate projections onto ignore (or 5th block, so the trained) decision boundary normal vectors (coming from training):
    
        # collect responses
        auxisa = []
        auxasa = []
        auxiea = []
        auxaea = []
        for clx in range(2):  # go over classes in both evoked and spontaneous
            auxisa.append(np.array( responses_trained[clx] )[:,T['start_idx']:T['stimstart_idx'],:].mean(axis=1).squeeze())
            auxasa.append(np.array( responses_predicted[clx] )[:,T['start_idx']:T['stimstart_idx'],:].mean(axis=1).squeeze())
            auxiea.append(np.array( responses_trained[clx] )[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1).squeeze())
            auxaea.append(np.array( responses_predicted[clx] )[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1).squeeze())
        onstimactivity = [[auxiea,auxaea]]         # {trainedignore} x {testignore, testattend} x {trials} x {cells}
        prestimactivity = [[auxisa,auxasa]]         # {trainedignore} x {testignore, testattend} x {trials} x {cells}
    
    
        n_neuron = block.segments[0].analogsignals[0].shape[1]
        
        # find the decision normal vectors
        c_db = []          # vector components (coefficients) of the decision normal in the activity space
        
        wx = int((len(crossdecoder)-numsignal)/n_neuron)     # here due to cross tests, we have +1 signal, so not 7, but 8 is the start position of the coefficients
        c_db.append(  np.reshape(np.array(crossdecoder[numsignal:]), (wx,n_neuron,crossdecoder[numsignal].shape[0],crossdecoder[numsignal].shape[1]) ).mean(axis=0)    )
            
        c_db = np.array(c_db) # [1,neurons,trajectory,stats]        (here only 1 in the first dimension, only one taskaspect: train on ignore)
    
        # mean decision normal vector with averaging over SA and EA:
        c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=2)  # use entire stimulus time          # this is a comparison group by neuron by   stats matrix
        norms_c_db = np.linalg.norm(c_db_means,axis=0)
    
    #    findneo
    
    
        dbprojections = [0]         # only use the first two, which is the visual in ignore
    
    
        # ingredients:
        # activity projected onto the basis:
        # projections[cx,sx,gx]      # {visual,context,choice} x {pre,onstim} x {attend go, attend nogo, ignore go, ignore nogo}
        # decision normal vectors as basis:
        # c_db_means[dbprojections[0]][:,0]      # visual db normal
    
        
        # create orthonormal basis from the list of db normals above
        c_db_means_matrix = np.vstack( [c_db_means[dbprojections[i]][:,0] for i in range(len(dbprojections))] ).T
        #  use only the mean-> [:,0] from the distribution
    
        Q = c_db_means_matrix
    
    
    
        # transform the activity vectors to the new orthonormal basis
        # projections:   {trainedattend,trainedignore} x {prestim,onstim} x {testedignore,testedattend} x {go,nogo}
        for cx,comparisonid in enumerate(dbprojections):        # go over the desired projection coordinates: visual/audio, choice, context
            projected_prestim = []
            projected_onstim = []
            for clx in [0,1]:   # go over the classes (45 and 135 degrees in this instance), and times evoked spontaneous
                auxpre = []
                auxon = []
                for trlx in [0,1]: # go over the activity groups: tested on ignore, tested on attend
                    auxpre.append(  np.dot(   np.array(prestimactivity[cx][clx][trlx]),Q[:,cx] ) )
                    auxon.append(  np.dot(   np.array(onstimactivity[cx][clx][trlx]),Q[:,cx] ) )
                projected_prestim.append(auxpre)
                projected_onstim.append(auxon)
            projections.append( [projected_prestim,projected_onstim] )
    #    projections = np.array(projections)
        





    # collect variance of the activities projected onto decision boundary normal vector
    variances = []
    preon = 0 # spontaneous only
    for twidx in [0,1]: #  trained where: attend or ignore
        variance = []
        for trlx in [0,1]: # test on attend or ignore
            var = ( [ np.std( projections[twidx][preon][trlx][0]), np.std( projections[twidx][preon][trlx][1] ),\
                      np.std( np.r_[projections[twidx][preon][trlx][0],projections[twidx][preon][trlx][1]]  )  ] )
            variance.append(var)
        variances.append(variance)
    variances = np.array(variances)     # {train: attend,ignore} x {test: attend,ignore} x {class 1, class 2, class 1+2}


    # save projections and variance
    if 0:
        pickle.dump((projections,variances),open(cacheprefix+'continuous/%s%s,projections,variance_%s-%s.pck'%('cross',['',',correctonly'][correctonly],dn,continuous_method),'wb'))






    # FIGURES
    
    # plot metrics
    if 0:
        fig,ax = plt.subplots(1,2,figsize=(32,12))
        for twidx in [0,1]:                    # left train on attend test on attend and crossignore,  right vv.
            axs = ax[twidx]
            figs.plottoaxis_decoderrocauc(crossdecoders[twidx][:3],axs,plotcrosstest=True)       # plot the performance
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5,'chance level')
            axs.legend()
            axs.set_title(dn+', trained %s, testing %s '%(crosslabels[twidx+2][0],crosslabels[twidx+2][1]))
            
            save = 0
            if save:
                fig.savefig(resultpath+'%stest_%s-%dms%dms_%s'%(filename[filenameidx],continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
        





    if 0:
        fig,ax = plt.subplots(2,2,figsize=(24,12))
        preon = 1 # plot only on stimulus,
        # rows -> trained on attend or ignore
        # columns -> tested on attend or ignore
        for twidx in [0,1]:
            for trlx in [0,1]:
                axs = ax[twidx,trlx]          # use projections tested on: {trainedwhere} x {onstim,prestim} x {ignore,attend} x {45,135}
                axs.hist(projections[twidx][preon][trlx][0],bins=15,color='darkgreen',alpha=0.8-0.3*twidx,label='%s'%(classlabels[0]))
                axs.hist(projections[twidx][preon][trlx][1],bins=15,color='darkred',alpha=0.8-0.3*twidx,  label='%s'%(classlabels[1]))
                axs.set_title('%s     tested (projected) %s '%(['direct','cross'][trlx],crosslabels[twidx+trlx*2][1]))
                if trlx==0:
                    axs.set_ylabel('trained '+crosslabels[twidx+trlx*2][0])
                
                axs.legend()
                axs.set_xlim(-1,1)
        fig.suptitle(dn+' projections to %s DBNV'%taskaspects[trainidx]) 
        save = 0
        if save:
            fig.savefig(resultpath+'%s,project_%s-%dms%dms_%s'%(filename[filenameidx],continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)



    # plot projections # old
    # if 0:
    #     fig,ax = plt.subplots(2,2,figsize=(24,12))
    #     for preon in [0,1]:
    #         for trlx in [0,1]:
    #             axs = ax[preon,trlx]          # use projections tested on: {trainedwhere} x {onstim,prestim} x {ignore,attend} x {45,135}
    #             axs.hist(projections[0][preon][trlx][0],bins=15,color='darkgreen',alpha=0.8-0.3*preon,label='%s'%(crosslabels[filenameidx][trlx]))
    #             axs.hist(projections[0][preon][trlx][1],bins=15,color='darkred',alpha=0.8-0.3*preon,  label='%s'%(crosslabels[filenameidx][trlx]))
    #             if preon==0:
    #                 axs.set_title('activites from '+crosslabels[filenameidx][trlx])
    #             if trlx==0:
    #                 axs.set_ylabel(preonlabels[preon])
                
    #             axs.legend()
    #             axs.set_xlim(-1,1)
    #     fig.suptitle(dn+' projections to %s DBNV'%taskaspects[trainidx]) 
    #     save = 0
    #     if save:
    #         fig.savefig(resultpath+'%s,project_%s-%dms%dms_%s'%(filename[filenameidx],continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)


    
    if retvar: return variances, projections











def decoder_withincontext_saveprojections(dn,block):

    width = 50*pq.ms
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [ \
                          [   [ [blv[1]],[45],    [] ], [ [blv[1]],[135],     [] ] ],\
                          [   [ [bla[1]],[45],    [] ], [ [bla[1]],[135],     [] ] ] \
                        ]
    taskaspects = ['attend visual','ignore visual']



    # go over attend and then ignore
    responses_attend = preprocess.collect_stimulusspecificresponses(block, comparisongroups[0])
    responses_ignore = preprocess.collect_stimulusspecificresponses(block, comparisongroups[1])
            
    acrossdecoderattend = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('attendignore','all',dn,continuous_method,taskaspects[0],'all'),'rb'))
    acrossdecoderignore = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('attendignore','all',dn,continuous_method,taskaspects[1],'all'),'rb'))



    # calculate projections onto decision boundary normal vectors (coming from training) for both conditions

    # collect responses
    auxa = []
    auxi = []
    for clx in range(2):  # go over classes
        auxa.append(np.array( responses_attend[clx] )[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1).squeeze())
        auxi.append(np.array( responses_ignore[clx] )[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1).squeeze())
    onstimactivity = [[auxa,auxi]]         # {trainedwithin} x {testattend,testignore} x {trials} x {cells}

#    print(len(onstimactivity),len(onstimactivity[0]),len(onstimactivity[0][0]),len(onstimactivity[0][0][0]))

    n_neuron = block.segments[0].analogsignals[0].shape[1]
    
    # find the ignore condition decision normal vectors
    c_db = []          # vector components (coefficients) of the decision normal in the activity space
    
    for cx,comparison in enumerate(taskaspects):
        acrossdecoder = [acrossdecoderattend,acrossdecoderignore][cx]
        wx = int((len(acrossdecoder)-7)/n_neuron)    # pca starts only at idx 7 = 8th position: train test 0-1 and angles 2-6
        c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
            
    c_db = np.array(c_db) # [attend-ignore,neurons,trajectory,stats]        (here only 1 in the first dimension, only one taskaspect: attend ignore)

    # mean decision normal vector with averaging over SA and EA:
    c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=2)  # use entire stimulus time          # this is a comparison group by neuron by   stats matrix
    norms_c_db = np.linalg.norm(c_db_means,axis=0)



    dbprojections = [0,1]         # only use the first two, which is the visual attend and ignore


    # ingredients:
    # activity projected onto the basis:
    # projections[cx,sx,gx]      # {visual,context,choice} x {pre,onstim} x {attend go, attend nogo, ignore go, ignore nogo}
    # decision normal vectors as basis:
    # c_db_means[dbprojections[0]][:,0]      # visual db normal

    
    # create orthonormal basis from the list of db normals above
    c_db_means_matrix = np.vstack( [c_db_means[dbprojections[i]][:,0] for i in range(len(dbprojections))] ).T
    #  use only the mean-> [:,0] from the distribution

    Q = c_db_means_matrix



    # transform the activity vectors to the new orthonormal basis
    projections = []  #  {attend,ignore} x {pre,onstim} x {go,nogo}
    for cx,comparisonid in enumerate(dbprojections):        # go over the desired projection coordinates: visual/audio, choice, context
        projected_prestim = []
        projected_onstim = []
        for clx in [0,1]:   # go over the classes (45 and 135 degrees in this instance)
            aux = []
            aux.append(  np.dot(   np.array(onstimactivity[0][cx][clx]),Q[:,cx] ) )
            projected_onstim.append(aux)
        projections.append( [projected_prestim,projected_onstim] )
#    projections = np.array(projections)
    

    # save projections
    if 0:
        pickle.dump(projections,open(cacheprefix+'continuous/withincontext,projectionstoignore_%s-%s.pck'%(dn,continuous_method),'wb'))



    return








def decoder_acc_relevantirrelevant(dn,block):
    # decode the relevant and the irrelevant stimuli, and congruent and conflicting stimuli trials separately, and for the two contexts


    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot

    
    n_neurons = block.segments[0].analogsignals[0].shape[1]

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
    taskaspects = ['visual,av','audio,av','visual,aa','audio,aa',\
                   'gonogo-congruent,av','gonogo-conflict,av','gonogo-congruent,aa','gonogo-conflict,aa']
    
    width = 50*pq.ms



    # assess performance

    perftask = preprocess.get_conditioned_behaviour(dn,comparisongroups,taskaspects)
    if 0:
        print(perftask.iloc[4:,:])
        return




    acrossdecoders = []
    for cx,comparison in enumerate(taskaspects):
        # if cx<4: continue
        print(comparison)
        if recalculate:
            acrossresponses = preprocess.collect_stimulusspecificresponses(block, comparisongroups[cx])
        
            acrossdecoder = nedi.get_responsedecoding(acrossresponses, width=width)

            pickle.dump(acrossdecoder,open(cacheprefix+'acc/stimulus,relevant-irrelevant_%s_%s-%s.pck'%(comparison,dn,continuous_method),'wb'))
        else:
            acrossdecoder = pickle.load(open(cacheprefix+'acc/stimulus,relevant-irrelevant_%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))

        acrossdecoders.append(acrossdecoder)
                    



    if doplot:

        fig,ax = plt.subplots(2,2,figsize=(2*12,2*8))

        colors = [  [['navy','lightgreen'],['dodgerblue','darkgreen']],\
                    [['rebeccapurple','deeppink'],['seagreen','coral']]   ]

        times = acrossdecoders[0][1].times


        for qx in [0,1]:    # questions: relevant-irrelevant  and  congruent-conflict
            for cx in [0,1]:    # attend visual, attend audio

                axs = ax[qx,cx]

                for sx in [0,1]:      # visual, audio

                    ad = acrossdecoders[qx*4+cx*2+sx][1]  # test timecourse

                    axs.plot(times, ad[:,0], color=colors[qx][cx][sx], lw=2, label=taskaspects[qx*4+cx*2+sx])
                    axs.fill_between(times, (ad[:,0]-ad[:,2]).squeeze(), (ad[:,0]+ad[:,2]).squeeze(), color=colors[qx][cx][sx], alpha=0.4)
                
                axs.legend(frameon=False)
                
                figs.setxt(axs)
                axs.set_ylim(0.45,1.05)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.5)

                if qx==0: axs.set_title('attend %s'%['visual','audio'][cx],fontsize=20)
                if cx==0: axs.set_ylabel('%s\ncv accuracy'%['relevant-irrelevant','congruent-conflict'][qx])
                


        fig.suptitle(dn+', %d neurons, relevant/irrelevant and congruent/conflicting trials'%n_neurons)



        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'relevantirrelevant+congruentconflict_%s_%s-%dms'%(dn,continuous_method,T['dt'])+ext)



    return







def decoder_acc_singletrialneuralbehaviour(dn,block):
     # decode the relevant and the irrelevant stimuli, and congruent and conflicting stimuli trials separately, and for the two contexts

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot

    
    n_neurons = block.segments[0].analogsignals[0].shape[1]

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
    taskaspects = ['visual,av','audio,av','visual,aa','audio,aa',\
                   'gonogo-congruent,av','gonogo-conflict,av','gonogo-congruent,aa','gonogo-conflict,aa']
    
    width = 50*pq.ms



    # assess performance

    perftask = preprocess.get_conditioned_behaviour_singletrial(dn,comparisongroups,taskaspects)

    # perftask_2c = pd.concat(perftask,axis=1)
    # perftask_unite = [ pd.concat(perftask[cx],axis=0) for cx in range(len(perftask)) ]

    # perftask_unite = [ pd.concat(  perftask[cx], axis=0, join='outer'  ) for cx in range(len(taskaspects)) ]


    # block boundaries:
    b = preprocess.getblocklengths(dn).cumsum()
    print(b)
    if blv[1]>bla[1]:
        b = b[ np.array([2,3,0,1],dtype=np.int16)  ]
    



    # for i in range(len(taskaspects)):
    #     if i>=1: continue
    #     print(perftask[i])
    #     print(pd.notna(perftask[i]['nogo']))
    #     print(pd.notna(perftask[i].iloc[:,1]))
    #     print(perftask[i].loc[pd.notna(perftask[i]['nogo']),'nogo'])
    #     print(pd.notna(perftask[i]['nogo']) & perftask[i]['nogo']==1)
    #     print(pd.notna(perftask[i]['nogo']) & perftask[i]['nogo']==0)

    #     perftask[i]['nogocorrect'] = perftask[i]['nogo']==1
    #     perftask[i]['nogoerror'] = perftask[i]['nogo']==0
    #     print(perftask[i])
    # 
    # return






    probasignals_tasks = []
    for cx,comparison in enumerate(taskaspects):
        # if cx<4: continue
        print(comparison)
        if recalculate:
            acrossresponses = preprocess.collect_stimulusspecificresponses(block, comparisongroups[cx])
        
            # acrossdecoder = nedi.get_responsedecoding(acrossresponses, width=width)
            probasignals, acc = nedi.get_pointwisedecoderprobability(acrossresponses, width=width)

            pickle.dump((probasignals, acc),open(cacheprefix+'acc/relevantirrelevant+congruentconflict,singletrialneuralbehaviour_%s_%s-%s.pck'%(comparison,dn,continuous_method),'wb'))
        else:
            probasignals, acc = pickle.load(open(cacheprefix+'acc/relevantirrelevant+congruentconflict,singletrialneuralbehaviour_%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))

        probasignals_tasks.append(probasignals)
                    



    # for i in range(len(taskaspects)):
    #     print(perftask[i].shape,probasignals_tasks[i].shape)




    if doplot:

       
        # fig,ax = plt.subplots(4,4,figsize=(4*8,4*8))

        # for tx in [0,1]:                 # tasks
        #     for cx in [0,1]:             # context
        #         for mx in [0,1]:         # modality/congconfl
        #             for clx in [0,1]:    # classes
        #                 axs = ax[tx*2+mx,cx*2+clx]
        #                 # x = perftask[tx*4+cx*2+mx,clx,b[cx*2]:b[cx*2+1]]
        #                 x = perftask[tx*4+cx*2+mx][clx]
        #                 axs.bar(x=range(x.shape[0]),height=x)
        #                 axs.set_title(taskaspects[tx*4+cx*2+mx]+' '+['go','nogo'][clx],fontsize=20)
        



        # save = 0 or globalsave
        # if save:
        #     fig.savefig(resultpath+'relevantirrelevant+congruentconflict,singletrialneuralbehaviour,nonly_%s_%s-%dms'%(dn,continuous_method,T['dt'])+ext)

    
    



        fig,ax = plt.subplots(4,2,figsize=(3*12,4*12))
        # we want to compare correct against error trials' class probability predictions
        # we have to do it for go and nogo trials separately (as they are the classes)

        for tx in [0,1]:                 # tasks
            for cx in [0,1]:             # context
                for mx in [0,1]:         # modality/congconfl
                    ix = tx*4+cx*2+mx
                    axs = ax[tx*2+mx,cx]
                    # x = perftask[ix]
                    y_early = probasignals_tasks[ix].magnitude[:,150:175].mean(axis=1)
                    y_late = probasignals_tasks[ix].magnitude[:,300:450].mean(axis=1)
                    
                    for clx,signal in enumerate(['go','nogo']):    # choose go and nogo signals respectively
                        # x = pd.notna(perftask[ix].iloc[:,clx])
                        # s = ['correct','error'][clx]
                        s = ['go','nogo'][clx]

                        mask_correct = perftask[ix][signal]==1     # correct trials
                        mask_error = perftask[ix][signal]==0       # error trials

                        # axs.boxplot(x=[y_early[mask_correct],y_early[mask_error]],positions=[0+clx*2,6+clx*2],labels=['%s\nearly\ncorrect'%s,'%s\nearly\nerror'%s])
                        # axs.boxplot(x=[y_late[mask_correct],y_late[mask_error]], positions=[1+clx*2,7+clx*2],labels=['%s\nlate\ncorrect'%s,'%s\nlate\nerror'%s])
                        axs.boxplot(x=[y_early[mask_correct],y_early[mask_error]],positions=[0+clx*2.5,6+clx*2.5],labels=['%s\nearly\ncorrect'%s,'%s\nearly\nerror'%s], notch=True)
                        axs.boxplot(x=[y_late[mask_correct],y_late[mask_error]], positions=[1+clx*2.5,7+clx*2.5],labels=['%s\nlate\ncorrect'%s,'%s\nlate\nerror'%s], notch=True)

                        # axs.boxplot(x=[y_early[x==0],y_early[x==1]],positions=[0+clx*2,6+clx*2],labels=['go\nearly\n%s'%s,'nogo\nearly\n%s'%s])
                        # axs.boxplot(x=[y_late[x==0],y_late[x==1]], positions=[1+clx*2,7+clx*2],labels=['go\nlate\n%s'%s,'nogo\nlate\n%s'%s])

                    axs.set_yticks([0,1])
                    axs.set_yticklabels(['P(go)=1','P(nogo)=1'])
                    figs.plottoaxis_chancelevel(axs,0.5)

                    if cx==0: axs.set_ylabel('class prob.')


                    axs.set_title(taskaspects[tx*4+cx*2+mx],fontsize=20)

        fig.suptitle(dn+', probability of predicted stimuli representations in two contexts in correct and error trials')
        
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'relevantirrelevant+congruentconflict,singletrialneuralbehaviour_%s_%s-%dms'%(dn,continuous_method,T['dt'])+ext)


















def decoder_behaviour(dn,block):
    # task performance (correct vs. erroneous choices)     in single trial neural representation (decoder predictive)
    width = 150*pq.ms
    blv,bla = preprocess.getorderattended(dn)
    bl_ls = np.array(preprocess.getblocklengths(dn))
    bl_ls = bl_ls[  np.hstack((blv,bla)).astype('int16')-1  ]
    cc = bl_ls[0]+bl_ls[1] # context change index
    
    
#    comparisongroups  =   [ [   [ [blv[1]],[],    [] ],   [ [bla[1]],[],     [] ]   ] ]
    comparisongroups  =   [ [   [ blv,[],    [] ],   [ bla,[],     [] ]   ] ]
    # first class is the visual, second is the audio context; only multimodal blocks? Not necessarily
    
    taskaspects = ['context']

    hit,miss,corrrej,fal = preprocess.assessbehaviouralperformance(dn,multimodalonly=False)
    behperflabels = ['hit','miss','corrrej','fal']
    
    for cx,comparison in enumerate(taskaspects):
        responses = preprocess.collect_stimulusspecificresponses(block, comparisongroups[cx])
        
            
        if 0:
            decoderprobapre = nedi.get_decoderprobability(responses,T,width=width,preon='pre')
            decoderprobaon = nedi.get_decoderprobability(responses,T,width=width,preon='on')
            pickle.dump((decoderprobapre,decoderprobaon),open(cacheprefix+'continuous/responsedecoderproba,%s_%s-%s.pck'%(comparison,dn,continuous_method),'wb'))
        else:
            decoderprobapre,decoderprobaon = pickle.load(open(cacheprefix+'continuous/responsedecoderproba,%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))
                    
    
    decoderprobapre[cc:,0] = 1-decoderprobapre[cc:,0]
    decoderprobaon[cc:,0] = 1-decoderprobaon[cc:,0]
    
    
    
    fig, ax = plt.subplots(2,2,figsize=(32,14))
    
    for sx in range(2):
        if sx==0:
            D = decoderprobapre
        else:
            D = decoderprobaon
#        A = -np.log(D)
    
        axs = ax[0,sx]
        axs.plot(hit,D[hit,0],        'o',color='green')
        axs.plot(miss,D[miss,0],      'o',color='red')
        axs.plot(corrrej,D[corrrej,0],'o',color='darkorange')
        axs.plot(fal,D[fal,0],        'o',color='purple')
        for ls in np.cumsum(bl_ls)[:-1]:
            axs.plot([ls,ls],[-0.2,1.2],color='darkgrey',alpha=0.5)
        axs.set_ylim([-0.05,1.05])
        axs.legend(behperflabels)
        axs.set_xlim(-1,np.sum(bl_ls)+1)
        axs.set_yticks([0,0.5,1])
        figs.plottoaxis_chancelevel(axs,0.5)
        axs.set_ylabel('probability of correct context')
        axs.set_xlabel('trial #')
        axs.set_title('single trial probabilities of correct context from %s '%['sa','ea'][sx])
        
        axs = ax[1,sx]
        axs.boxplot(    [D[hit[hit<cc],0],D[miss[miss<cc],0],\
                         D[corrrej[corrrej<cc],0],D[fal[fal<cc],0]], notch=True,positions=np.arange(4),labels=behperflabels       )
        axs.boxplot(    [D[hit[hit>cc],1],D[miss[miss>cc],1],\
                         D[corrrej[corrrej>cc],1],D[fal[fal>cc],1]], notch=True,positions=np.arange(5,9),labels=behperflabels       )
        axs.set_xlim(-1,9)
        axs.set_ylim([-0.05,1.05])
        axs.set_yticks([0,0.5,1])
        figs.plottoaxis_chancelevel(axs,0.5)
        axs.set_xlabel('visual context                                                  audio context')
        axs.set_ylabel('probability of correct context')


#
#        axs = ax[2,sx]
#        axs.boxplot(    [A[hit[hit<cc],0],A[miss[miss<cc],0],\
#                         A[corrrej[corrrej<cc],0],A[fal[fal<cc],0]], notch=True,positions=np.arange(4),labels=behperflabels       )
#        axs.boxplot(    [A[hit[hit>cc],1],A[miss[miss>cc],1],\
#                         A[corrrej[corrrej>cc],1],A[fal[fal>cc],1]], notch=True,positions=np.arange(5,9),labels=behperflabels       )
#        axs.set_xlim(-1,9)
#        axs.set_yscale('log')
##        axs.set_ylim(0.001,1-0.001)
##        figs.plottoaxis_chancelevel(axs,0.5)
#        axs.set_xlabel('visual context                                                  audio context')
#        axs.set_ylabel('odds of correct context')




    
    

    fig.suptitle(dn+' behavioral perforance vs. single trial probabilities of correct context\nblocks are visual,audio ordered')

    save = 0
    if save:
        fig.savefig(resultpath+'behavior,context,singletrialprobability_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
        











def exportdecoderaccuracies(dn):
    # a helper routine to export decoder accuracies for multiple task variables
    print('exporting',dn)
    taskaspects = ['visual','audio','context','choice']
    
    n_samples = 10
    
    sampling_rate = (T['dt']*n_samples).magnitude   # in ms
    
    # data is (taskaspects,{train,test},timepoints,{mean,2sem}) numpy array, with sampling rate
    data = []
    for cx,comparison in enumerate(taskaspects):
        acrossdecoder = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
    
        acrossdecoder = acrossdecoder[:2]
        acrossdecoder[0] = neph.downsample(acrossdecoder[0][:590,:],n_samples)
        acrossdecoder[1] = neph.downsample(acrossdecoder[1][:590,:],n_samples)
        data.append(acrossdecoder)
        

    # and the reward
    if dn in ['DT019','DT020','DT021','DT022']:
        data.append(np.zeros( ( len(acrossdecoder), acrossdecoder[0].shape[0], acrossdecoder[1].shape[1]) )    )
    else:
        acrossdecoder = pickle.load(open(cacheprefix+'continuous/rewarddecodes,syncstim,angles-%s_%s.pck'%(dn,continuous_method),'rb'))
        acrossdecoder = acrossdecoder[:2]
        acrossdecoder[0] = neph.downsample(acrossdecoder[0][:590,:],n_samples)
        acrossdecoder[1] = neph.downsample(acrossdecoder[1][:590,:],n_samples)
        data.append(acrossdecoder)
    
    
        
    data = np.array(data)[:,:,:,:]
    
    
    
    # data is ready to export
    print(data.shape, 'sampling rate %s ms'%sampling_rate)
    
    np.save('../notes/data/decoderaccuracies-samplingrate,%dms-%s.npy'%(sampling_rate,dn),data)
    
    
    # plt.plot(data[:,1,:,0].T)
    
    # data = np.load('../notes/data/decoderaccuracies-samplingrate,%dms-%s.npy'%(sampling_rate,dn))


    
    return









# **************************************************************************************
#                                SUBSPACES









def subspaces_examples(dn,block):
    # testing visualization of subspace projections
    
    from mpl_toolkits.mplot3d import Axes3D

    
    n_features = 2
    T = 2000
    T2 = 1000
    s = 0.05
    c = s*0.25*0
    mu1=np.array([0.5,0,0])
    mu2=np.array([-0.5,0,0])
    
    
    Xa = ( np.random.multivariate_normal( mean=mu1,cov=np.eye(3)*(s-c)+np.ones((3,3))*c, size=T2 ) )
    Xb = ( np.random.multivariate_normal( mean=mu2,cov=np.eye(3)*(s-c)+np.ones((3,3))*c, size=T2 ) )
    
    X = np.vstack(  [Xa,Xb]  )
            
    classes = np.hstack(   [ np.zeros((T2,),dtype='int16'),np.ones((T2,),dtype='int16')          ]       )
    
    
    X = sp.stats.zscore(X)

    
    _,dec = nedi.pointwiseclassifier(X,classes,lastclassone=T2,rr=1)

    decisionnormal = (dec.coef_/np.linalg.norm(dec.coef_))[0]





    # for first figure

    
    v0 = (mu1+mu2)/2. ;# v0 = np.zeros(2)
    d_mu = mu2 - mu1
    


    
    
    U,S,V = np.linalg.svd(X)
    
    eta = np.vstack((v0,v0+d_mu))
    n = np.vstack((v0,v0+decisionnormal))
    pc = [  np.stack((v0,v0+V[0])),  np.stack((v0,v0+V[1]))]
    
    
    # angles
    a_orthosvd = nedi.innerproductangle( V[0,:], V[1,:])
    a_normalsvd = nedi.innerproductangle( V[0,:], decisionnormal)








    # projections for second figure

    tpc_coords = [0,1]
    cpc_coords = [0,1]


    Ua,Sa,Va_h = np.linalg.svd(Xa)

    Ub,Sb,Vb_h = np.linalg.svd(Xb)


    
    U,S,V_h = np.linalg.svd( X )
    V = V_h.T
    
    
    B = V[:,tpc_coords].dot(V[:,tpc_coords].T)

    PXa = np.dot(  Xa.dot(Va_h), B)
    PXb = np.dot(  Xb.dot(Vb_h), B)
    PsXa = np.dot(Xa,V)
    PsXb = np.dot(Xb,V)

    PVa = np.dot(B,Va_h.T)
    PVb = np.dot(B,Vb_h.T)
    Pcoeffs = np.dot(B,decisionnormal)












    
    
#    print(a_normalsvd)
#    
#    
#
#    print(dec.intercept_)
#    print(dec.coef_,decisionnormal)
#    print(S,V)
#    
#    print('angles')
#    print(a_orthosvd,a_normalsvd)

    
    
    
    
    fig = plt.figure(num=1,figsize=(28,14))
    ax = fig.add_subplot(1,2,1,projection='3d')
    ax.plot( X[:T2,0], X[:T2,1], X[:T2,2],  'o', color='mediumvioletred', alpha=0.2)
    ax.plot( X[T2:,0], X[T2:,1], X[T2:,2],  'o', color='darkturquoise', alpha=0.2)


    ax.plot(eta[:,0],eta[:,1],eta[:,2],color='fuchsia',linewidth=3,alpha=0.9)
    ax.plot(n[:,0],n[:,1],n[:,2],color='green',linewidth=3,alpha=0.9)
    ax.plot(pc[0][:,0],pc[0][:,1],pc[0][:,2],color='darkgoldenrod',linewidth=3,alpha=0.9)
    ax.plot(pc[1][:,0],pc[1][:,1],pc[1][:,2],color='gold',linewidth=3,alpha=0.9)
    
#    ax.set_xlim([-2,2]);ax.set_ylim([-2,2])
#    ax.set_xlim([0,2]);ax.set_ylim([0,2]);ax.set_zlim([0,2])





    ax = fig.add_subplot(1,2,2,projection='3d')

    ax.plot(  PXa[:,0], PXa[:,1], PXa[:,2], 'o', linewidth=2,alpha=0.2,color='mediumvioletred', label='class 1' )
    ax.plot(  PXb[:,0], PXb[:,1], PXa[:,2], 'o', linewidth=2,alpha=0.2,color='darkturquoise', label='class 2' )

#    ax.plot(  PsXa[:,0], PsXa[:,1], PXa[:,2], 'o', linewidth=2,alpha=0.2,color='darkgoldenrod', label='coords class 1' )
#    ax.plot(  PsXb[:,0], PsXb[:,1], PXa[:,2], 'o', linewidth=2,alpha=0.2,color='gold', label='coords class 2' )

    
    ax.plot(  [0, PVa[0,cpc_coords[0]]],  [0, PVa[1,cpc_coords[0]]], [0, PVa[2,cpc_coords[0]]], linewidth=3,alpha=0.6,color='darkmagenta',label='class 1 pc 1' )
    ax.plot(  [0, PVa[0,cpc_coords[1]]],  [0, PVa[1,cpc_coords[1]]], [0, PVa[2,cpc_coords[1]]], linewidth=3,alpha=0.6,color='orchid',label='class 1 pc 2' )

    ax.plot(  [0, PVb[0,cpc_coords[0]]],  [0, PVb[1,cpc_coords[0]]], [0, PVa[2,cpc_coords[0]]], linewidth=3,alpha=0.6,color='darkblue',label='class 2 pc 1' )
    ax.plot(  [0, PVb[0,cpc_coords[1]]],  [0, PVb[1,cpc_coords[1]]], [0, PVa[2,cpc_coords[1]]], linewidth=3,alpha=0.6,color='royalblue',label='class 2 pc 2' )



    ax.plot(  [0, Pcoeffs[0]],  [0, Pcoeffs[1]], [0, Pcoeffs[2]], linewidth=3,alpha=1.0,color='green',label='decision' )

    ax.set_xlabel('PC 1'); ax.set_ylabel('PC 2'); ax.set_zlabel('PC 3'); 


    
    return










def subspaces(dn,block,examine='allexpcond'):
    # creates all subspace projections
    # saves decision vectors and projections for other routines
    # any other subspace routine needs to be preceded by this routine with recalculate = True
    
    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot
    
    blv,bla = preprocess.getorderattended(dn)
    if examine=='allexpcond':
        comparisongroups  = [ \
                              [  [[ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                                 [[ [blv[1]], [45],    [5000] ], [ [blv[1]],[135],     [5000] ]  ] ],\
                              [  [[ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                                 [[ [bla[1]],  [45],[5000] ], [ [bla[1]],   [45],[10000] ]    ] ],\
                              [  [[  blv,  [],[] ],   [ bla,   [], [] ] ],\
                                 [[ [blv[1]],  [45],[5000] ], [ [bla[1]],   [45],[5000] ]   ] ], \
                              [  [[  [2,4],[], []  ],[[2,4],[], []  ] ], \
                                 [[  [2,4],[], []  ],[[2,4],[], []  ] ]   ],
                              [  [[ [blv[1]],[45],    [] ], [ [blv[1]],[135],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [bla[1]],[45],    [] ], [ [bla[1]],[135],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [bla[1]],[],    [5000] ], [ [bla[1]],[],     [10000] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [blv[1]],[],    [5000] ], [ [blv[1]],[],     [10000] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[  [blv[1]],[], []  ],[[blv[1]],[], []  ] ], \
                                 [[  [blv[1]],[], []  ],[[blv[1]],[], []  ] ]   ],
                              [  [[  [bla[1]],[], []  ],[[bla[1]],[], []  ] ], \
                                 [[  [bla[1]],[], []  ],[[bla[1]],[], []  ] ]   ],
                              [  [[  [5],  [45],[] ],   [ [5],   [135], []   ]],\
                                            []], \
                              [  [[  [5],  [90],[] ],   [ [5],   [180], []   ]],\
                                            []], \
                              [  [[  [5],  [270],[] ],   [ [5],   [180], [] ]],\
                                            []], \
                              [  [[  [5],  [60],[] ],   [ [5],   [150], []   ]],\
                                            []], \
                              [  [[  [5],  [120],[] ],   [ [5],   [210], []   ]],\
                                            []], \
                              [  [[ [blv[1]],[],    [] ], [ [blv[1]],[],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [bla[1]],[],    [] ], [ [bla[1]],[],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                            ]
        taskaspects = ['visual','audio','context','choice','visualattend','visualignore',\
                       'audioattend','audioignore','choiceattendvisual','choiceattendaudio',\
                       'char,45-135','char,90-180','char,270-180','char,60-150','char,120-210',\
                       'attend visual','ignore visual']
        logical = 'and'
        trialconditions = ['all']#,'identical,go']
    elif examine=='attendignore':
        comparisongroups  = [ \
                              [  [ [ [blv[1]], [45],        [] ],  [ [blv[1]],[135],        [] ] ],\
                                 [ [ [blv[1]], [45],    [5000] ],  [ [blv[1]],[135],    [5000] ]  ] ],\
                              [  [ [ [bla[1]], [45],        [] ],  [ [bla[1]],[135],        [] ] ],\
                                 [ [ [bla[1]], [45],    [5000] ],  [ [bla[1]],[135],    [5000] ] ]  ],\
                              [  [ [ [bla[1]],   [],    [5000] ],  [ [bla[1]],[],      [10000] ] ],\
                                 [ [ [bla[1]], [45],    [5000] ],  [ [bla[1]],[45],    [10000] ] ]  ],\
                              [  [ [ [blv[1]],   [],    [5000] ],  [ [blv[1]],[],      [10000] ] ],\
                                 [ [ [blv[1]], [45],    [5000] ],  [ [blv[1]],[45],    [10000] ] ]  ]\
                            ]
        taskaspects = ['attend visual','ignore visual','attend audio','ignore audio']
        logical = 'and'
        trialconditions = ['all']#,'identical,go']
        

    lmask = np.arange(block.segments[0].analogsignals[0].shape[1])
    layeredfilename = 'all'
    layeredtitle = 'all'
    
    n_neuron = block.segments[0].analogsignals[0].shape[1]

    width = 50*pq.ms
    if continuous_method=='count': width = 1000*pq.ms


    if not recalculate:    #   this section draws the angles between aspects in "taskaspects", recalculate first in the following section
        c_db = []        # coefficients of decision boundaries
        
        stimulusspecificity = 'all'
        for cx,comparison in enumerate(taskaspects):
            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            wx = int((len(acrossdecoder)-7)/n_neuron)
            c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
        
        c_db = np.array(c_db)

        A = np.empty(  (len(taskaspects),len(taskaspects)  ), dtype=neo.AnalogSignal ) #, c_db[0].shape[1])  ) 

        if doplot: fig,ax = plt.subplots(len(comparisongroups),len(comparisongroups),figsize=(7*len(comparisongroups),7*len(comparisongroups)))

        for cx1,comparison1 in enumerate(taskaspects):
            for cx2,comparison2 in enumerate(taskaspects):
                if doplot: axs = ax[cx1,cx2]
                if cx2<=cx1:
                    if doplot: axs.set_visible(False)
                    continue
                else:
                    aux = np.zeros(c_db.shape[2])
                    for t in range(c_db[cx1].shape[1]):
                        aux[t] = nedi.innerproductangle(  c_db[cx1][:,t,0],  c_db[cx2][:,t,0]  )

                    A[cx1,cx2] = neo.AnalogSignal(aux, name='%s,%s'%(comparison1,comparison2),\
                             t_start=T['starttime'], sampling_period=T['dt'], units=pq.dimensionless)
                    A[cx2,cx1] = A[cx1,cx2]      # make it a symmetric matrix


                    if doplot:
                        axs.plot(A[cx1,cx2].times, A[cx1,cx2],color='purple',linewidth=2)
                        axs.set_title(comparison2)
                        axs.set_ylabel(comparison1)
                    
                        figs.plottoaxis_chancelevel(axs,np.pi/2,'orthogonal')
                        figs.plottoaxis_chancelevel(axs,0.0,'parallel')
                        axs.set_ylim([0-0.1,np.pi/2+0.1])
                        figs.plottoaxis_stimulusoverlay(axs,T)
                        axs.set_yticks([0,np.pi/2]); axs.set_yticklabels(['$0$','$\pi/2$'])
                        axs.set_xlim(-1500,4500)
                    

        if doplot:
            fig.suptitle(dn+', decision normal angles in neural vectorspace, between aspects')
            save=0
            if save:
                fig.savefig(resultpath+'aspects,decoders,subsp,character-%s_%s_%s-%dms%dms_%s'%(examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)

        if 1: # whether to save A into a vector
            pickle.dump(A,open(cacheprefix+'subspaces/angles,aspects,subspaces,character-%s_%s-%s-%dms%dms_%s.pck'%(examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn),'wb'))

        return      # no need to recaluculate or redraw, so return







    if recalculate:   # this section calculates (retrains) decoders and optionally displays coefficients in time
        if doplot:
            fig,ax = plt.subplots(9,2*len(taskaspects),  figsize=(11*2*len(taskaspects),9*4*2))
#            fig = plt.figure( figsize=(10,10)); ax = fig.gca()
    
    for cx,comparison in enumerate(taskaspects):
        # DT008 does not have a 45-135 in the 5th block, so treat as empty list element
        # if dn=='DT008' and cx==10: continue        # only needed for 5th block stimuli
        for sx,stimulusspecificity in enumerate(trialconditions):
            
            stimulusIDgroups = comparisongroups[cx][sx]

            print(cx,sx,comparison,stimulusspecificity,comparisongroups[cx][sx],stimulusIDgroups[0][0])


    
            # collect neural responses
            if not comparison[:6]=='choice': # visual, audio, context:
                acrossresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups,logical)
            else:  # choice:
                acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn,onlyinblock=stimulusIDgroups[0][0])

    
    
            if recalculate:  # whether to retrain the decoders
                # if dn=='DT008' and cx in [9]: print('YES')

                acrossdecoder = nedi.get_responsedecoding(acrossresponses, width=width)
                print('saving file','subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity))
                pickle.dump(acrossdecoder,open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'wb'))
            else:
                print('loading projection subspaces from saved file','subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity))
                acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            
            

            if doplot:   # display optionally the coefficients in time
                axs = ax[0,cx*2+sx]
                figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs)       # plot the performance
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.5,'chance level')
                if cx==0 and sx==0: axs.legend(); axs.set_ylabel('roc auc at each timepoint')
                axs.set_title(comparison+', '+stimulusspecificity)
                
                
                # decision coefficients
                
                wx = int((len(acrossdecoder)-7)/n_neuron)
                
                axs = ax[1,cx*2+sx]
                figs.plottoaxis_decodercoeffs(acrossdecoder[7::wx],axs)       # plot the subspace directions
                figs.plottoaxis_chancelevel(axs,0.0,'invariant')
    #            axs.set_ylim([0-0.1,np.pi/2+0.1])
                figs.plottoaxis_stimulusoverlay(axs,T)
    #            axs.set_yticks([0,np.pi/2]); axs.set_yticklabels(['$0$','$\pi/2$'])
                if cx==0 and sx==0: axs.legend(); axs.set_ylabel('neurons, first timepoint\ndecoder decision boundary\nnormal coefficients')
                axs.set_title(comparison+', '+stimulusspecificity)

            
                axs = ax[2,cx*2+sx]
                figs.plottoaxis_decodercoeffs(acrossdecoder[7:7+wx],axs)       # plot the subspace directions
                figs.plottoaxis_chancelevel(axs,0.0,'invariant')
                figs.plottoaxis_stimulusoverlay(axs,T)
                if cx==0 and sx==0: axs.legend(); axs.set_ylabel('timepoints, neuron 1\ndecoder decision boundary\nnormal coefficients')
                axs.set_xlabel('time from trial onset [ms]')
            
            
                axs = ax[3,cx*2+sx]
                figs.plottoaxis_decodercoeffs(acrossdecoder[7+wx:7+2*wx],axs)       # plot the subspace directions
                figs.plottoaxis_chancelevel(axs,0.0,'invariant')
                figs.plottoaxis_stimulusoverlay(axs,T)
                if cx==0 and sx==0: axs.legend(); axs.set_ylabel('timepoints, neuron 2\ndecoder decision boundary\nnormal coefficients')
                axs.set_xlabel('time from trial onset [ms]')

            
                axs = ax[4,cx*2+sx]
                figs.plottoaxis_decodercoeffs(acrossdecoder[7+2*wx:7+3*wx],axs)       # plot the subspace directions
                figs.plottoaxis_chancelevel(axs,0.0,'invariant')
                figs.plottoaxis_stimulusoverlay(axs,T)
                if cx==0 and sx==0: axs.legend(); axs.set_ylabel('timepoints, neuron 3\ndecoder decision boundary\nnormal coefficients')
                axs.set_xlabel('time from trial onset [ms]')




                # subspaces
                
                # trajectories over the timecourse
                axs = ax[5,cx*2+sx]
                figs.plottoaxis_subspacestrajectory(acrossresponses,acrossdecoder[7::wx],[0,1],axs)
                if cx==0 and sx==0: axs.legend();
                axs.set_xlabel('PC 1'); axs.set_ylabel('%sPC 2'%(['','trajectories\n'][cx+sx==0]))
                axs.set_title(comparison+', '+stimulusspecificity)



                # subspace projections at checkpoint
                axs = ax[6,cx*2+sx]
                tpc_coords = [0,1]
                figs.plottoaxis_neuralandboundary(acrossresponses,acrossdecoder[7::wx],tpc_coords=tpc_coords,ax=axs)
                figs.plottoaxis_crosshair(axs)
                if cx==0 and sx==0: axs.legend();
                axs.set_xlabel('PC %d'%(tpc_coords[0]+1)); axs.set_ylabel('%sPC %d'%(['','projections\n'][cx+sx==0],tpc_coords[1]+1))
                axs.set_title(comparison+', '+stimulusspecificity)
                
                axs = ax[7,cx*2+sx]
                tpc_coords = [0,2]
                figs.plottoaxis_neuralandboundary(acrossresponses,acrossdecoder[7::wx],tpc_coords=tpc_coords,ax=axs)
                figs.plottoaxis_crosshair(axs)
                axs.set_xlabel('PC %d'%(tpc_coords[0]+1)); axs.set_ylabel('PC %d'%(tpc_coords[1]+1))
                
                axs = ax[8,cx*2+sx]
                tpc_coords = [1,2]
                figs.plottoaxis_neuralandboundary(acrossresponses,acrossdecoder[7::wx],tpc_coords=tpc_coords,ax=axs)
                figs.plottoaxis_crosshair(axs)
                axs.set_xlabel('PC %d'%(tpc_coords[0]+1)); axs.set_ylabel('PC %d'%(tpc_coords[1]+1))
                
                
                
    if doplot:
        save = 1
        if save:
            fig.savefig(resultpath+'decoder+svdpc,traj,subsp-%s_%s_%s-%dms%dms_%s'%(examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)





            
    return







    
    
    
    
    
    
    
    


def subspaces_projections(dn,block,examine='allexpcond',play='visual',preon=1):
    # project activities to subspaces defined by span of a combination of task variable decoder decision boundary normal vectors

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot
    
    layeredfilename = 'all'
    
    n_neuron = block.segments[0].analogsignals[0].shape[1]
    
    stimulusspecificity = 'all'
    blv,bla = preprocess.getorderattended(dn)
    
    # comparion groups and taskaspects below hold the main trial groupings for activities and decoders
    # later variables will index into these
    comparisongroups  = [ \
                              [  [[ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                                 [[ [blv[1]], [45],    [5000] ], [ [blv[1]],[135],     [5000] ]  ] ],\
                              [  [[ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                                 [[ [bla[1]],  [45],[5000] ], [ [bla[1]],   [45],[10000] ]    ] ],\
                              [  [[ [blv[1]],  [45],[] ],   [ [bla[1]],   [45], [] ] ],\
                                 [[  blv,  [],[] ],     [ bla,   [],[] ]   ] ], \
                              [  [[ [blv[1]],  [135],[] ],   [ [bla[1]],   [135], [] ] ],\
                                 [[  blv,  [],[] ],     [ bla,   [],[] ]   ] ], \
                              [[2,4],[2,4]], \
                              [  [[ [blv[1]],[45],    [] ], [ [blv[1]],[135],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [bla[1]],[45],    [] ], [ [bla[1]],[135],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [bla[1]],[],    [5000] ], [ [bla[1]],[],     [10000] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [blv[1]],[],    [5000] ], [ [blv[1]],[],     [10000] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [[blv[1]],[]], \
                              [[bla[1]],[]], \
                              [  [[  [5],  [45],[] ],   [ [5],   [135], []   ]],\
                                            []], \
                              [  [[  [5],  [90],[] ],   [ [5],   [180], []   ]],\
                                            []], \
                              [  [[  [5],  [270],[] ],   [ [5],   [180], [] ]],\
                                            []], \
                              [  [[  [5],  [60],[] ],   [ [5],   [150], []   ]],\
                                            []], \
                              [  [[  [5],  [120],[] ],   [ [5],   [210], []   ]],\
                                            []], \
                              [  [[ [blv[1]],[],    [] ], [ [blv[1]],[],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [bla[1]],[],    [] ], [ [bla[1]],[],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ blv,[],    [] ], [   bla,[],     [] ] ],\
                                 [[  [],  [],[] ],     [ [],   [],[] ]   ] ], \
                              [  [[ [bla[1]],  [],[5000] ],   [ [blv[1]],   [], [5000] ] ],\
                                 [[  blv,  [],[] ],     [ bla,   [],[] ]   ] ], \
                              [  [[ [bla[1]],  [],[10000] ],   [ [blv[1]],   [], [10000] ] ],\
                                 [[  blv,  [],[] ],     [ bla,   [],[] ]   ] ], \
                              [  [[ [blv[1]],  [45],[] ],   [ [blv[1]],   [135], [] ] ],\
                                 [[  blv,  [],[] ],     [ bla,   [],[] ]   ] ], \
                              [  [[ [bla[1]],  [],[5000] ],   [ [blv[1]],   [], [10000] ] ],\
                                 [[  blv,  [],[] ],     [ bla,   [],[] ]   ] ], \
                               
                         ]
    taskaspects = ['visual','audio','context','context','choice','visualattend','visualignore',\
                       'audioattend','audioignore','choiceattendvisual','choiceattendaudio',\
                       'char,45-135','char,90-180','char,270-180','char,60-150','char,120-210',\
                       'attend visual','ignore visual','context','context','context','context','context']
    #                   0       1         2     3          4           5             6
    #                   7                 8                 9                 10
    #                   11                 12           13               14          15
    #                   16                   17             18      19        20        21        22
    
    # basisdict = dict(zip(taskaspects,np.arange(len(taskaspects),dtype=np.int16)))
        
    logical = 'and'
    trialconditions = ['all']   #,'identical,go']
    # if dn=='DT008': taskaspects[11] = []        # to avoid char45-135 for DT008 those trials do not exist
    # play = 'visual'
    # play = 'audio'
    # preon = 0       # 0: pre, 1: onstim
    
    preonlabels = ['pre','on']
    # play='audio'
    
    # triallistindex is a booster to quickly select trials for specific questions (called plays here) from the function argument
    if play=='visual': triallistindex = [2,3]           # use this for visual: which task aspect to use the classes for as data points (2,3 = context visual go and nogo)
    elif play=='audio': triallistindex = [19,20]         # use this for audio
    elif play=='choice': triallistindex = [2,3,19,20]        # and this for choice
    

    # triallistindex1D = [12,13]
    # triallistindex1D = [12,13]


    # dbprojections = [0,1]        # audio-visual
    # dbprojections = [0,1,2,4]    # audio-visual-context-choice
    # dbprojections = [0,4]        # visual-choice


# 3D:
    
    # dbprojections is for the coordinates to project the activities onto
    # in 3D plot later there will be a basisvectors for the dbnv directions, that are indexing! dbprojections
    if play=='visual': dbprojections = [0,2,4,5,6,9,10]            # visual-context-choice for 3D, +visual attend, visual ignore
    elif play=='audio': dbprojections = [1,2,4,7,8,9,10]           # audio-context-choice for 3D, +audio attend, audio ignore +choicevisual,choiceaudio
    elif play=='choice': dbprojections = [0,1,4,3,5,6,7,8,9,10,2]  #[0,1,4,5,6,7,8,9,10,2]        # visual-audio-choice for 3D, +visual attend, audio attend +choicevisual,choiceaudio
    elif play=='context': dbprojections = [2,0,4,5,6,9,10]            # context-visual-choice for 3D, +visual attend, visual ignore

    # dbprojections = [0,7,8,9,10,11]      # visual and visual in rf characterization block

    depth = 1000                  # 1000, or 3000
    depth_idx = int((depth*pq.ms/T['dt']).magnitude)





    # collect responses
    prestimactivity = []       # this will have for each comparisongroup a [cx][classes][trials,neurons] list of matrices
    onstimactivity = []
    for cx,comparison in enumerate(taskaspects):
        for sx,stimulusspecificity in enumerate(trialconditions):
            stimulusIDgroups = comparisongroups[cx][sx]
            # collect neural responses
            if not comparison[:6]=='choice': # visual, audio, context:
                acrossresponses = preprocess.collect_stimulusspecificresponses(block,stimulusIDgroups,logical)
            else:  # choice:
                acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn,onlyinblock=stimulusIDgroups)
            # acrossresponses is [class][trials][trajectory,neurons]

            aux1 = []
            aux2 = []
            for clx in range(2):
                if dn=='DT008' and cx==11: # doesn't have characterization
                    aux1.append(np.zeros((1,1)))
                    aux2.append(np.zeros((1,1)))
                else:
                    aux1.append(np.array( acrossresponses[clx] )[:,:T['stimstart_idx'],:].mean(axis=1).squeeze())
                    # aux2.append(np.array( acrossresponses[clx] )[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1).squeeze())
                    aux2.append(np.array( acrossresponses[clx] )[:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,:].mean(axis=1).squeeze())
            prestimactivity.append(aux1)
            onstimactivity.append(aux2)
            print(cx,comparison,onstimactivity[cx][0].shape,onstimactivity[cx][1].shape)

    # collect decision boundary normal vectors

    c_db = []          # vector components (coefficients) of the decision normal in the activity space
    
    for cx,comparison in enumerate(taskaspects):
        if dn=='DT008' and cx==11:      # doesn't have characterization
            c_db.append(np.zeros((10,596,3)))
        else:
            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
            wx = int((len(acrossdecoder)-7)/n_neuron)
            c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
        
    c_db = np.array(c_db) # [comparisongroup,neurons,trajectory,stats]

    A = np.empty(  (len(taskaspects),len(taskaspects)  ), dtype=neo.AnalogSignal ) #, c_db[0].shape[1])  ) 

    # mean context decision normal vector with averaging over SA:
#    c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=2)            # this is a comparison group by neuron by   stats matrix
    c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,:].mean(axis=2)  # only first 1 sec          # this is a comparison group by neuron by   stats matrix
    norms_c_db = np.linalg.norm(c_db_means,axis=0)
    
    
#    pca_features = []; expl_var = []; pca_components = []
#    for cx,comparison in enumerate(taskaspects):
#        for sx,stimulusspecificity in enumerate(trialconditions):
#            pca_features,expl_var,pca_components = nedi.pca(,n_components=3)
    
    



    #   3 D   with scalar products
    if 1: # projection to stimulus, choice and context decision normal bases
        ixv45,ixv135,ixa45,ixa135,ixv5000,ixv10000,ixa5000,ixa10000 = preprocess.getstimulusidents(dn,block,multimodalonly=True)
        ixhit,ixmiss,ixcorrrej,ixfal = preprocess.assessbehaviouralperformance(dn)
        
        
#        colors = ['navy','darkturquoise','darkred','fuchsia']
        if play=='visual':
            colors = ['b','c','m','r']
            basiscolors = ['navy','blue','dodgerblue','fuchsia','sienna','gold','darkorange']
            basisvectors = [0,3,4,1,2,5,6]      # indexes of dbprojections, dbprojections was: [0,2,4,5,6,9,10]
            labels = ['visual decision normal','attended visual dec. normal','ignored visual dec. normal',\
                      'context decision normal','choice decision normal','attended visual choice','ignored visual choice',\
                      'ignore visual,  45','ignore visual, 135','attend visual,  45, hit','attend visual,  45, miss','attend visual, 135, corr.rej.','attend visual, 135, false alarm']
            _,ixgohit,_ = np.intersect1d(ixv45,ixhit,return_indices=True)
            _,ixnogocorrrej,_ = np.intersect1d(ixv135,ixcorrrej,return_indices=True)
            _,ixgomiss,_ = np.intersect1d(ixv45,ixmiss,return_indices=True)
            _,ixnogofal,_ = np.intersect1d(ixv135,ixfal,return_indices=True)
        elif play=='audio':
            colors = ['darkgreen','lime','sienna','darkorange']
            basiscolors = ['darkgreen','seagreen','limegreen','fuchsia','sienna','gold','darkorange']
            basisvectors = [0,3,4,1,2,5,6]
            labels = ['audio decision normal','attended audio dec. normal','ignored audio dec. normal',\
                      'context decision normal','choice decision normal','attended audio choice','ignored audio choice',\
                      'ignore audio,  5k','ignore audio, 10k','attend audio,  5k, hit','attend audio,  5k, miss','attend audio, 10k, corr.rej.','attend audio, 10k, false alarm']
            _,ixgohit,_ = np.intersect1d(ixa5000,ixhit,return_indices=True)
            _,ixnogocorrrej,_ = np.intersect1d(ixa10000,ixcorrrej,return_indices=True)
            _,ixgomiss,_ = np.intersect1d(ixa5000,ixmiss,return_indices=True)
            _,ixnogofal,_ = np.intersect1d(ixa10000,ixfal,return_indices=True)
        elif play=='choice':
            colors = ['b','c','darkgreen','lime']
            basiscolors = ['navy','darkgreen','sienna','blue','dodgerblue','seagreen','limegreen','gold','orange','fuchsia']
            basisvectors = [0,1,4,5,6,7,8,9,10,2]        # 0,1,3,4,5,6,7,8,9,10,2
                                                         # 0,4,5,1,6,7,10,3,8,9
            labels = ['visual decision normal','audio decision normal','choice decision normal','attended visual dec. normal','ignored visual dec. normal',\
                      'attended audio dec. normal','ignored audio dec. normal',\
                      'attended visual choice','attended audio choice','context decision normal',\
                      'attend visual,  45, hit','attend visual,  45, miss','attend visual, 135, corr.rej.','attend visual, 135, false alarm',\
                      'attend audio,  5k, hit','attend audio,  5k, miss','attend audio, 10k, corr.rej.','attend audio, 10k, false alarm']
            _,ixgohit,_ = np.intersect1d(ixv45,ixhit,return_indices=True)
            _,ixnogocorrrej,_ = np.intersect1d(ixv135,ixcorrrej,return_indices=True)
            _,ixgomiss,_ = np.intersect1d(ixv45,ixmiss,return_indices=True)
            _,ixnogofal,_ = np.intersect1d(ixv135,ixfal,return_indices=True)
            _,ixgohit2,_ = np.intersect1d(ixa5000,ixhit,return_indices=True)
            _,ixnogocorrrej2,_ = np.intersect1d(ixa10000,ixcorrrej,return_indices=True)
            _,ixgomiss2,_ = np.intersect1d(ixa5000,ixmiss,return_indices=True)
            _,ixnogofal2,_ = np.intersect1d(ixa10000,ixfal,return_indices=True)
        if play=='context':
            colors = ['b','c','m','r']
            basiscolors = ['fuchsia','navy','darkorange']
            basisvectors = [0,1,2,]      # indexes of dbprojections, dbprojections was: [0,2,4,5,6,9,10]
            labels = ['context decision normal','visual decision normal','choice decision normal',\
                      'attended visual dec. normal','ignored visual dec. normal',\
                      'attended visual choice','ignored visual choice',\
                      'ignore visual,  45','ignore visual, 135','attend visual,  45, hit','attend visual,  45, miss','attend visual, 135, corr.rej.','attend visual, 135, false alarm']
            _,ixgohit,_ = np.intersect1d(ixv45,ixhit,return_indices=True)
            _,ixnogocorrrej,_ = np.intersect1d(ixv135,ixcorrrej,return_indices=True)
            _,ixgomiss,_ = np.intersect1d(ixv45,ixmiss,return_indices=True)
            _,ixnogofal,_ = np.intersect1d(ixv135,ixfal,return_indices=True)
            

        # ingredients:
        # activity projected onto the basis:
        # projections[cx,sx,gx]      # {visual,context,choice} x {pre,onstim} x {attend go, attend nogo, ignore go, ignore nogo}
        # decision normal vectors as basis:
        # c_db_means[dbprojections[0]][:,0]      # visual db normal
        # c_db_means[dbprojections[1]][:,0]      # context db normal
        # c_db_means[dbprojections[2]][:,0]      # choice db normal       etc. (depending on what is the index within dbprojections at positions 0,1,2, here example dbprojections = [0,2,3,...])

        
        # create a concatenated list of bases (together forming a basis) from the list of db normals above
        c_db_means_matrix = np.vstack( [c_db_means[dbprojections[i]][:,0] for i in range(len(dbprojections))] ).T
        # columns are choice, context, visual dbnvs; use only the mean-> [:,0] from the distribution



        # find an orthonormal basis representation via the Gram-Schmidt process
        # use the orthonormalorder list below to choose the first three coordinates that will be drawn to as close as possible:
        # choose the first in orthonormalorder as the first basis, then find the perpendicular to it in the direction
        # of the second basis, then a new one that is orthogonal to both 1st and 2nd, to and so on

        # the followings should all be done with play='visual' with dbprojections=[0,2,4,5,6,9,10] 
        # orthonormalorder = [0,1,2]; displayorder = [0,1,2]; orthoorderlabel = 'default'
        # orthonormalorder = [2,1,0]; displayorder = [2,1,0]; orthoorderlabel = 'chcxvi'    # this is when choice is kept as coord. z, context as next best orthogonal, and visual/audio as third best
        # orthonormalorder = [1,0,2]; displayorder = [1,0,2]; orthoorderlabel = 'cxvich'    # this is for representing context and visual for Fig3
        # orthonormalorder = [2,0,1]; displayorder = [1,2,0]; orthoorderlabel = 'chvicx'    # this is for representing choice and visual for Fig5
        # orthonormalorder = [2,3,1]; displayorder = [0,1,2]; orthoorderlabel = 'chvia'    # this is for representing attend visual and choice for Fig6 subspaces
        # orthonormalorder = [2,4,1]; displayorder = [0,1,2]; orthoorderlabel = 'chvii'    # this is for representing ignore visual and choice for Fig6 subspaces

        # this is only for play='choice'
        orthonormalorder = [0,1,2]; displayorder = [0,1,2]; orthoorderlabel = 'chviau'    # this is for representing choice around both stimuli modality
        # orthonormalorder = [2,1,2]; displayorder = [0,1,2]; orthoorderlabel = 'chviau'    # this is for representing choice around both stimuli modality


        # do the Gram-Schmidt process to find the QR decomposition for the 3 first most important selected quasi-directions
        Q,R = np.linalg.qr(c_db_means_matrix[:,orthonormalorder])     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix



        # now we need to rearrange the coordinates to a desired order for the first three coordinates, for display;
        # this keeps the orthogonality found in QR decomposition, it just changes the order for the horizontal, deep and vertical axes of the coordinate system
        # the order we want here is not necessarily reversing that in orthonormalorder, as dbprojections does not necessarily
        # contain the optimal 3D display order, but rather what makes sense as grouping them for display order in legend
        
        # the following section has been revised and replaced by displayorder variable at orthonormalorder lines and after these lines
        ## Q = Q[:,::-1]               # change back to visual/audio,context,choice;  don't use for fig3,5 and 6, only for original "kept as coord. z line
        ## Q[:,np.array([0,1,2],dtype=np.int)] = Q[:,np.array([1,0,2],dtype=np.int)]   # this is to order back for Fig3
        ## Q[:,np.array([0,1,2],dtype=np.int)] = Q[:,np.array([1,2,0],dtype=np.int)]   # this is to order back for Fig6
        ## Q[:,np.array([0,1,2],dtype=np.int)] = Q[:,np.array([0,1,2],dtype=np.int)]   # this is to order back for Fig6
        ## Q[:,np.array([0,1,2],dtype=np.int)] = Q[:,np.array([0,1,2],dtype=np.int)]   # this is to order for for choice visual audio, with visual and audio as the first 2

        # change back the new coordinates of the original basis vectors, with the proper new order for display
        Q = Q[:,np.array([0,1,2],dtype=np.int)] = Q[:,np.array(displayorder,dtype=np.int)]
        R = np.dot(Q.T,c_db_means_matrix)        # this is just for the record, R won't be used, as Q will be used directly on everything that needs transformation
        


        # project decision boundary normal basis vectors in their skewed coordinate system to the new orthonormal coordinate system
        dbnbasiscomps = np.dot(Q.T,c_db_means_matrix) / np.linalg.norm(c_db_means_matrix,axis=0,keepdims=True)         # orthonormbasis x dbbasistypes
        
        print('stats',c_db_means[0].shape)
        print('skewed basis',c_db_means_matrix.shape)
        print('orthonormal basis',Q.shape)
        print('transformed',R.shape)
        print((100*R).astype('int16'))
        print('db components in onb:\n',dbnbasiscomps,'\n',dbnbasiscomps.shape,'\nbasis vectors',[taskaspects[b] for b in basisvectors],'\n\n')


        # however the coordinates are to be flipped to the ++ quadrant, with the np.sign operator
        # if below is set to true, it will save the dbn bases and activity projections with the flipped display (subsequent QR-found orthogonal) coordinates
        flip = 0

        # transform the activity vectors to the new orthonormal basis
        projections = []  #{visual,context,choice} x {pre,onstim} x {attend go, attend nogo, ignore go, ignore nogo}
        for coordx,cx in enumerate(dbprojections[:3]):        # go over the desired projection coordinates: visual/audio, choice, context
            print('projection index, dbnv number',coordx,cx,'from dbprojections',dbprojections[:3])
            projected_prestim = []
            projected_onstim = []
            for clx in [0,1]:   # go over the classes (attend and ignore in e.g. in the context instance)
                for trlx in triallistindex: # go over the activity groups: here context attend, and context attend other
                    # project pre- and onstimulus activities separately
                    # print(clx,trlx)
                    # print('pre',len(prestimactivity), 'on',len(onstimactivity),)
                    # print('pre',len(prestimactivity[trlx]),len(prestimactivity[trlx][clx]),'on',len(onstimactivity[trlx]),len(onstimactivity[trlx][clx]))
                    projected_prestim.append(  np.dot(   prestimactivity[trlx][clx],Q[:,coordx]  ) )
                    projected_onstim.append(  np.dot(   onstimactivity[trlx][clx],Q[:,coordx] ) )
                    if flip:
                        projected_prestim[-1] *= np.sign(dbnbasiscomps[coordx,cx])
                        projected_onstim[-1] *= np.sign(dbnbasiscomps[coordx,cx])
            projections.append( [projected_prestim,projected_onstim] )
        projections = np.array(projections)

        # we also need to flip the dbn vectors similarly
        if flip:
            dbnbasiscomps *= np.sign(dbnbasiscomps)




        # calculate real orthogonal 1D projections for the actual normals for visual and context (choice will be in the original third coordinate (swapped from QR first))
#        Q_vi,R = np.linalg.qr(c_db_means_matrix[:,[0,1,2]])     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
#        Q_cx,R = np.linalg.qr(c_db_means_matrix[:,[1,0,2]])     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
        
#        projections_vi = []
#        projected_prestim = []
#        projected_onstim = []
#        for clx in [0,1]:   # go over the classes (attend and ignore in this instance)
#            for trlx in triallistindex: # go over the activity groups: here context attend, and context attend other
#                # project pre- and onstimulus activities separately
##                    print(clx,trlx,len(onstimactivity),len(onstimactivity[trlx]),len(onstimactivity[trlx][clx]))
#                projected_prestim.append(  np.dot(   prestimactivity[trlx][clx],Q_vi[:,0]  ) )
#                projected_onstim.append(  np.dot(   onstimactivity[trlx][clx],Q_vi[:,0] ) )
#        projections_vi.append( [projected_prestim,projected_onstim] )
#        projections_vi = np.array(projections_vi)
#
#        projections_cx = []
#        projected_prestim = []
#        projected_onstim = []
#        for clx in [0,1]:   # go over the classes (attend and ignore in this instance)
#            for trlx in triallistindex: # go over the activity groups: here context attend, and context attend other
#                # project pre- and onstimulus activities separately
##                    print(clx,trlx,len(onstimactivity),len(onstimactivity[trlx]),len(onstimactivity[trlx][clx]))
#                projected_prestim.append(  np.dot(   prestimactivity[trlx][clx],Q_cx[:,0]  ) )
#                projected_onstim.append(  np.dot(   onstimactivity[trlx][clx],Q_cx[:,0] ) )
#        projections_cx.append( [projected_prestim,projected_onstim] )
#        projections_cx = np.array(projections_cx)
#



        # save projections, and decide if figure should be shown
        if recalculate:
            pickle.dump((dbnbasiscomps,projections,colors,basiscolors,basisvectors,(ixgohit,ixnogocorrrej,ixgomiss,ixnogofal)),\
                        open(cacheprefix+'subspaces/dbprojections+projections,%s-%s_%s_%s-%s-%dms%dms_%s.pck'%(play,orthoorderlabel,examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn),'wb'))
        if not doplot: return



        # PLOT

        
        # figure scale:
        L = np.max(  np.hstack(  [ np.abs(projections[j,preon,k]) for j in range(3) for k in range(4)  ]   ) ) * 1.
#        L=8

        
        
        # plot the projections in the orthonormal coordinates system
        fig = plt.figure(figsize=(32,32))        #num=1,


        ax = fig.add_subplot(2,2,1,projection='3d')

        axs = ax
        # the axes will be: x visual, y context, z choice
        
        # plot the basis:
        x0 = 0; y0 = 0; z0 = 0    #-L*0.9*c1*s2
        # plot visual, visualattend, visualignore, context anc choice; only the || of the direction is interesting
#        for r,color in zip([0,3,4,1,2],['navy','dodgerblue','purple','fuchsia','gold']):
        for r,color in zip(basisvectors,basiscolors):
            axs.plot( [x0,x0+L/3*(dbnbasiscomps[0,r])],\
                      [y0,y0+L/3*(dbnbasiscomps[1,r])],\
                      [z0,z0+L/3*(dbnbasiscomps[2,r])],\
                  color=color,linewidth=4,alpha=0.9 )


        # plot the activities
            
        # plot the ignore
        if not play=='choice':
            for k in [2,3]:       # the last two projections: ignore conditions
#                print(projections[0,preon,k],projections[1,preon,k],projections[2,preon,k])
                axs.plot(projections[0,preon,k], projections[1,preon,k], projections[2,preon,k], 'o',color=colors[k],alpha=0.8)
        
        # plot the attended with performance
        axs.plot(projections[0,preon,0][ixgohit],      projections[1,preon,0][ixgohit],      projections[2,preon,0][ixgohit],      'o',color=colors[0],alpha=0.8)
        axs.plot(projections[0,preon,0][ixgomiss],     projections[1,preon,0][ixgomiss],     projections[2,preon,0][ixgomiss],     'x',markersize=10,color=colors[0],alpha=0.9)
        axs.plot(projections[0,preon,1][ixnogocorrrej], projections[1,preon,1][ixnogocorrrej], projections[2,preon,1][ixnogocorrrej], 'o',color=colors[1],alpha=0.8)
        axs.plot(projections[0,preon,1][ixnogofal],     projections[1,preon,1][ixnogofal],     projections[2,preon,1][ixnogofal],     'x',markersize=10,color=colors[1],alpha=0.9)
        if play=='choice':
            axs.plot(projections[0,preon,2][ixgohit2],      projections[1,preon,2][ixgohit2],      projections[2,preon,2][ixgohit2],      'o',color=colors[2],alpha=0.8)
            axs.plot(projections[0,preon,2][ixgomiss2],     projections[1,preon,2][ixgomiss2],     projections[2,preon,2][ixgomiss2],     'x',markersize=10,color=colors[2],alpha=0.9)
            axs.plot(projections[0,preon,3][ixnogocorrrej2], projections[1,preon,3][ixnogocorrrej2], projections[2,preon,3][ixnogocorrrej2], 'o',color=colors[3],alpha=0.8)
            axs.plot(projections[0,preon,3][ixnogofal2],     projections[1,preon,3][ixnogofal2],     projections[2,preon,3][ixnogofal2],     'x',markersize=10,color=colors[3],alpha=0.9)





        # we need 3D ellipsoids here
#        for k in range(4):
#            figs.confidence_ellipse(projections[0,1,k]*np.cos(awayangle), projections[1,1,k], ax, n_std=2.0, facecolor=colors[k],alpha=0.15)


        axs.set_xlim(-L,L)
        axs.set_ylim(-L,L)
        axs.set_zlim(-L,L)
#        axs.axis('off')
#        axs.set_title(dn+' evoked')

        axs.legend(labels,frameon=False)#,ncol=2)#loc='upper left')

        labels = [taskaspects[dbprojections[0]],taskaspects[dbprojections[1]]+' orth. to vi comp.',taskaspects[dbprojections[2]]+' orth. to vi&cx comp.']
        labels = [taskaspects[dbprojections[0]]+' orth. to ch&cx comp.',taskaspects[dbprojections[1]]+' orth. to ch comp.',taskaspects[dbprojections[2]]]
        axs.set_xlabel(labels[0])
        axs.set_ylabel(labels[1])
        axs.set_zlabel(labels[2])
        


        bp = np.array([[0,1],[0,2],[1,2]]) # basis pairs for the 2D unfolding
        for b in [0,1,2]:      # coordinate axes
            ax = fig.add_subplot(2,2,b+2)
            axs = ax
            
            # plot skewed basis:
            for r,color in zip(basisvectors,basiscolors):
                axs.plot( [x0,x0+L/2*(dbnbasiscomps[bp[b,0],r])],\
                          [y0,y0+L/2*(dbnbasiscomps[bp[b,1],r])],\
                          color=color,linewidth=4,alpha=0.9 )


            # plot the ignore
            if not play=='choice':
                for k in [2,3]:
#                    print(projections[0,preon,k],projections[1,preon,k],projections[2,preon,k])
                    axs.plot(projections[bp[b,0],preon,k], projections[bp[b,1],preon,k], 'o',color=colors[k],alpha=0.8)
            
            # plot the attended with performance
            axs.plot(projections[bp[b,0],preon,0][ixgohit],      projections[bp[b,1],preon,0][ixgohit],      'o',color=colors[0],alpha=0.8)
            axs.plot(projections[bp[b,0],preon,0][ixgomiss],     projections[bp[b,1],preon,0][ixgomiss],     'x',markersize=10,color=colors[0],alpha=0.9)
            axs.plot(projections[bp[b,0],preon,1][ixnogocorrrej], projections[bp[b,1],preon,1][ixnogocorrrej], 'o',color=colors[1],alpha=0.8)
            axs.plot(projections[bp[b,0],preon,1][ixnogofal],     projections[bp[b,1],preon,1][ixnogofal],     'x',markersize=10,color=colors[1],alpha=0.9)
            if play=='choice':
                axs.plot(projections[bp[b,0],preon,2][ixgohit2],      projections[bp[b,1],preon,2][ixgohit2],      'o',color=colors[2],alpha=0.8)
                axs.plot(projections[bp[b,0],preon,2][ixgomiss2],     projections[bp[b,1],preon,2][ixgomiss2],     'x',markersize=10,color=colors[2],alpha=0.9)
                axs.plot(projections[bp[b,0],preon,3][ixnogocorrrej2], projections[bp[b,1],preon,3][ixnogocorrrej2], 'o',color=colors[3],alpha=0.8)
                axs.plot(projections[bp[b,0],preon,3][ixnogofal2],     projections[bp[b,1],preon,3][ixnogofal2],     'x',markersize=10,color=colors[3],alpha=0.9)
            
            # plot the cov ellipsoids
            for k in range(4):
                figs.confidence_ellipse(projections[bp[b,0],preon,k], projections[bp[b,1],preon,k], ax, n_std=2.0, facecolor=colors[k],alpha=0.15)

            axs.set_xlabel(labels[bp[b,0]])
            axs.set_ylabel(labels[bp[b,1]])






        
        fig.suptitle(dn +' (%s) %s, context and choice bases, %s, %d ms, %s - V:%d, A:%d'%(orthoorderlabel, play, preonlabels[preon], depth, continuous_method,blv[1],bla[1]))
        # fig.suptitle(dn +' visual, audio and %s bases, %s, %d ms, %s - V:%d, A:%d'%(play, preonlabels[preon], depth, continuous_method,blv[1],bla[1]))
        
        
        
        save = 0
        if save:
            fig.savefig(resultpath+'contextclasses,projectonto,%sviau,%s,d%dms,3D+3x2D-%s_%s_%s-%dms%dms_%s'%(play[:2],preonlabels[preon],depth,examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)













def subspaces_decode_nullspace1visual1(dn,block):
    # take the DV, and find its nullspace, then decode from its nullspace
    # take the projection onto visual DV and decode from there


    recalculate = 0 or globalrecalculate
    doplot = 1


    blv,bla = preprocess.getorderattended(dn)
    comparisongroups = [ \
                          [        [ [2,4],  [45],[] ],       [ [2,4],   [135], [] ]         ],\
                          [        [ [2,4],  [],[5000] ],       [ [2,4],   [], [10000] ]         ],\
                          [        [ blv,  [],[] ],       [ bla,   [], [] ]         ],\
                          [ [], [] ],
                       ]
    taskaspects = ['visual','audio','context','choice']            # this is only done for the context for now

    width = 50*pq.ms
    
    
    
    
    print(dn)
    
    if recalculate:

        for cx,comparison in enumerate(taskaspects):
            print(comparison)
    
    
    
    
            # get the DV
    
            n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
            celltypes = block.annotations['celltypes']
            
            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
    
            c_db_matrix = np.zeros((n_neurons,n_trajectory))       # prepare the dbnv vector holder
    
            wx = int((len(acrossdecoder)-7)/n_neurons)            # we are only concerned with the coefficients coming after the train test (2) and PCs (5).
            width = wx * T['dt']    # we will retroactively do a similar feature with decoder
            
            coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
    
            # we do not take the mean of the stimulus over timecourse, because we want to remove the actual nullspace at each timepoint
            c_db_matrix = coeff[:,:,0]      # neurons by trajectory
            c_db_matrix /= np.linalg.norm(c_db_matrix,axis=0)
            


            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'visual','all'),'rb'))
    
            c_db_matrix_visual = np.zeros((n_neurons,n_trajectory))       # prepare the dbnv vector holder
    
            wx = int((len(acrossdecoder)-7)/n_neurons)            # we are only concerned with the coefficients coming after the train test (2) and PCs (5).
            width = wx * T['dt']    # we will retroactively do a similar feature with decoder
            
            coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
    
            # we do not take the mean of the stimulus over timecourse, because we want to remove the actual nullspace at each timepoint
            c_db_matrix_visual = coeff[:,:,0]      # neurons by trajectory
            c_db_matrix_visual /= np.linalg.norm(c_db_matrix_visual,axis=0)






    
    
            print('projections at all timepoints', c_db_matrix.shape)
    
    
    
            # remove any activity that is projected on to the DV, so that
            # any activity that remains is in the nullspace of DV
    
            if not comparison=='choice': # visual, audio, context:
                acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx])
            else:  # choice:
                acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
    
            
            acrossresponses_nullspace = [[],[]]
            acrossresponses_visualspace = [[],[]]
            for aix,responses in enumerate(acrossresponses):    # get either class of response collect
                n_trials_toshow = 8
                trials_toshow = np.random.permutation(len(responses))[:n_trials_toshow]
                trx = 0 # counts trials displayed
                
                if doplot: fig,ax = plt.subplots(n_trials_toshow,n_neurons,figsize=(n_neurons*12,n_trials_toshow*8))
                
                
                for tr,response_trial in enumerate(responses):                # get the analogsignal that is at this trial, (neurons,trajectory)
                    # print(type(response_trial),response_trial.shape)
    
                    # find the DV projection with the original neural activity space coordinate system
                    
                    nullspace_trial = response_trial.copy()
                    projection_trial = response_trial.copy()
                    visualspace_trial = response_trial.copy()
                    DV_projections = np.zeros((n_neurons,n_neurons,n_trajectory-wx))
                    DV_visual_projections = np.zeros((n_neurons,n_neurons,n_trajectory-wx))
                    
                    for t in range(n_trajectory-wx):
                        # shorthand vertical vector for the DV
                        v = c_db_matrix[:,t][:,np.newaxis]
                        DV_projections[:,:,t] = v @ np.linalg.inv(v.T @ v) @ v.T
                        projection_trial[t,:] = DV_projections[:,:,t] @ response_trial[t,:].magnitude
                        nullspace_trial[t,:] = response_trial[t,:] - DV_projections[:,:,t] @ response_trial[t,:].magnitude        # will result in (neurons,trajectory)


                        v = c_db_matrix_visual[:,t][:,np.newaxis]
                        DV_visual_projections[:,:,t] = v @ np.linalg.inv(v.T @ v) @ v.T
                        visualspace_trial[t,:] = DV_visual_projections[:,:,t] @ response_trial[t,:].magnitude        # will result in (neurons,trajectory)
                        
                        # test if programming is ok (linear algebra rank)
                        # print('rank context projection', np.linalg.matrix_rank(DV_projections[:,:,t]),\
                        #       'rank nullspace', np.linalg.matrix_rank(np.eye(DV_projections.shape[0])-DV_projections[:,:,t]),\
                        #       'rank visual projection', np.linalg.matrix_rank(DV_visual_projections[:,:,t]))
                        
                    
                    acrossresponses_nullspace[aix].append(nullspace_trial)
                    acrossresponses_visualspace[aix].append(visualspace_trial)


                    
                    if doplot:
                        if tr in trials_toshow:
                            for n in range(n_neurons):
                                axs = ax[trx,n]
                                axs.plot(nullspace_trial.times,response_trial[:,n],color='navy',lw=3,alpha=0.5,label='original activity')
                                axs.plot(nullspace_trial.times,projection_trial[:,n],color='fuchsia',lw=3,alpha=0.5,label='decision vector projection')
                                axs.plot(nullspace_trial.times,nullspace_trial[:,n],color='darkred',lw=3,alpha=0.5,label='nullspace projection')
                                axs.plot(nullspace_trial.times,visualspace_trial[:,n],color='dodgerblue',lw=3,alpha=0.5,label='visual space projection')
                                figs.plottoaxis_chancelevel(axs)
                                figs.setxt(axs)
                                figs.plottoaxis_stimulusoverlay(axs,T)
                                if trx==0: axs.set_title('neuron %d'%(n+1))
                                if n==0: axs.set_ylabel('trial %d'%(tr+1))
                                axs.legend(frameon=False)
                            trx += 1
    
                if doplot: fig.suptitle('%s, %s'%(dn,comparison))
    

            return

        # perform the decoding in the nullspace of the DV
        # and over the visual DV subspace

            acrossdecoder_nullspace = nedi.get_responsedecoding(acrossresponses_nullspace,width=width)
            pickle.dump(acrossdecoder_nullspace,open(cacheprefix+'subspaces/nullspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'wb'))
            acrossdecoder_visualspace = nedi.get_responsedecoding(acrossresponses_visualspace,width=width)
            pickle.dump(acrossdecoder_visualspace,open(cacheprefix+'subspaces/visualspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'wb'))
        else:
            acrossdecoder_nullspace = pickle.load(open(cacheprefix+'subspaces/nullspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
            acrossdecoder_visualspace = pickle.load(open(cacheprefix+'subspaces/visualspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))



    if doplot:
            
        fig,ax = plt.subplots(1,len(taskaspects),figsize=(len(taskaspects)*12,1*8))
        
        for cx,comparison in enumerate(taskaspects):
            if len(taskaspects)==1: axs = ax
            else: axs = ax[cx]
            # times = acrossdecoder_nullspace[1].analogsignals[0].times
            
            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
            acrossdecoder_nullspace = pickle.load(open(cacheprefix+'subspaces/nullspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
            acrossdecoder_visualspace = pickle.load(open(cacheprefix+'subspaces/visualspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))

            figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs)       # plot the performance
            figs.plottoaxis_decoderrocauc(acrossdecoder_nullspace[:2],axs,colorlist=['navy','darkred'])       # plot the performance
            figs.plottoaxis_decoderrocauc(acrossdecoder_visualspace[:2],axs,colorlist=['gold','rebeccapurple'])       # plot the performance
            figs.plottoaxis_stimulusoverlay(axs,T)
            figs.plottoaxis_chancelevel(axs,0.5)
            
            axs.set_title(comparison)


        fig.suptitle('%s - original (orange), sameDVnullspace (red), visualDVspace (purple)'%(dn))



        save = 0
        if save:
            fig.savefig(resultpath+'decoders,over,dv-nullspace,dv-visual,projectedactivities-%s-%s'%(dn,continuous_method)+ext)










def subspaces_decode_nullspacerecurrent(dn,block):
    # take the DV, and find its nullspace, then decode from its nullspace
    # repeat until the 0 subspace is all that left
    print('calculating recurrent projection nullspaces...')

    recalculate = 0 or globalrecalculate
    doplottrials = 0
    doplot = 1


    blv,bla = preprocess.getorderattended(dn)
    comparisongroups = [ \
                          [        [ [2,4],  [45],[] ],       [ [2,4],   [135], [] ]         ],\
                          [        [ [2,4],  [],[5000] ],       [ [2,4],   [], [10000] ]         ],\
                          [        [ blv,  [],[] ],       [ bla,   [], [] ]         ],\
                          [ [], [] ],
                       ]
    taskaspects = ['visual','audio','context','choice']            # this is only done for the context for now

    width = 50*pq.ms
    
    
    n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
    

    print(dn)
    
    if recalculate:
        for cx,comparison in enumerate(taskaspects):
            print(comparison)
    
            # iteratively project onto smaller and smaller subspaces:
            ranks = np.zeros((n_neurons-1,n_trajectory-5))
            previous_rank = n_neurons

            acrossdecoder_nullspaces = []
            for px in range(n_neurons-1):
                print('going for %d/%d orthogonal complement subspace'%(px+1,n_neurons))
    
        
        
                # get the previous DV. At first its the original DV:
                if px==0:
                    acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
                else:  # later it will be the DV over the previously calculated orthogonal complement subspace
                    acrossdecoder = acrossdecoder_nullspaces[-1]

                


                c_db_matrix = np.zeros((previous_rank,n_trajectory))       # prepare the dbnv vector holder
        
                wx = int((len(acrossdecoder)-7)/previous_rank)            # we are only concerned with the coefficients coming after the train test (2) and PCs (5).
                width = wx * T['dt']    # we will retroactively do a similar feature with decoder
                
                coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,previous_rank,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    # mean over the feature width of 50 ms forward of the decoder so only num neurons are feautres remaining
        
                # we do not take the mean of the stimulus over timecourse, because we want to remove the actual nullspace at each timepoint
                c_db_matrix = coeff[:,:,0]      # subspace neural basis by trajectory
                c_db_matrix /= np.linalg.norm(c_db_matrix,axis=0)
                



                print('projections at all timepoints', c_db_matrix.shape)
        
                # calculate projection matrixes for all time t
                DV_projections = []# np.zeros((n_neurons,n_neurons,n_trajectory-wx))
                for t in range(n_trajectory-wx):
                    # shorthand vertical vector for the DV
                    v = c_db_matrix[:,t][:,np.newaxis]
                    # project with 1-DV onto the original neural space coordinate system DV-orthogonal complement
                    A = np.eye(previous_rank) - v @ np.linalg.inv(v.T @ v) @ v.T
                    Q,R = np.linalg.qr(A)
                    r = np.linalg.matrix_rank(R)
                    DV_projections.append( R[:r,:] )        # reduce the projection to include only subspaces with orthogonal complement
                    # ranks[px,:,t] = np.linalg.matrix_rank(DV_projections[:,:,t]), np.linalg.matrix_rank(np.eye(DV_projections.shape[0])-DV_projections[:,:,t])
                    ranks[px,t] = r
                print('    rank time average: ', ranks[px,:].mean())
                previous_rank = np.round(ranks[px,:].mean()).astype(np.int16)




                # remove any activity that is projected on to the DV, so that
                # any activity that remains is in the nullspace of the DV projector /in the orthogonal complement of the DV
        
                if px==0: # if first iteration, do the original activity
                    if not comparison=='choice': # visual, audio, context:
                        acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx])
                    else:  # choice:
                        acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
                else: # after that go for each subsequent orthogonal complement starting from the previous
                    acrossresponses = acrossresponses_nullspace
                
                # go over both trial groups (classes) of activity
                acrossresponses_nullspace = [[],[]]
                for aix,responses in enumerate(acrossresponses):    # get either class of response collect
                    n_trials_toshow = 8
                    trials_toshow = np.random.permutation(len(responses))[:n_trials_toshow]
                    trx = 0 # counts trials displayed
                    
                    if doplottrials: fig,ax = plt.subplots(n_trials_toshow,n_neurons,figsize=(n_neurons*12,n_trials_toshow*8))
                    



                    for tr,response_trial in enumerate(responses):                # get the analogsignal that is at this trial, (neurons,trajectory)
                        # print(type(response_trial),response_trial.shape)
        
                        # find the DV projection with the original neural activity space coordinate system
                        
                        nullspace_trial = response_trial[:,:previous_rank].copy()
                        
                        for t in range(n_trajectory-wx):
                            # remove the DV space, i.e. project onto the orthogonal complement
                            nullspace_trial[t,:] = DV_projections[t] @ response_trial[t,:].magnitude

                        acrossresponses_nullspace[aix].append(nullspace_trial)


                        
                        if doplottrials:
                            if tr in trials_toshow:
                                for n in range(previous_rank):
                                    axs = ax[trx,n]
                                    axs.plot(nullspace_trial.times,response_trial[:,n],color='navy',lw=3,alpha=0.5,label='original activity')
                                    axs.plot(nullspace_trial.times,nullspace_trial[:,n],color='darkred',lw=3,alpha=0.5,label='nullspace projection')
                                    figs.plottoaxis_chancelevel(axs)
                                    figs.setxt(axs)
                                    figs.plottoaxis_stimulusoverlay(axs,T)
                                    if trx==0: axs.set_title('neuron %d'%(n+1))
                                    if n==0: axs.set_ylabel('trial %d'%(tr+1))
                                    axs.legend(frameon=False)
                                trx += 1
        
                    if doplottrials: fig.suptitle('%s, %s'%(dn,comparison))


        # perform the decoding in the nullspace of the DV

                acrossdecoder_nullspace = nedi.get_responsedecoding(acrossresponses_nullspace,width=width)
                acrossdecoder_nullspaces.append(acrossdecoder_nullspace)
            
            # display ranks for crosscheck
            print(ranks.mean(axis=1).ravel())

            # save the collection of recurrent decoders of all four variable over the current orthogonal complement subspace of the variable
            pickle.dump((acrossdecoder_nullspaces,ranks),open(cacheprefix+'subspaces/nullspace,recurrent,decodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'wb'))
        else:
            acrossdecoder_nullspaces,ranks = pickle.load(open(cacheprefix+'subspaces/nullspace,recurrent,decodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))



    if doplot:
        colorlists = [ ['dodgerblue','darkorange'], ['navy','darkred'] ]
        # ranks = np.zeros((n_neurons+1,2,596))
        fig,ax = plt.subplots(1,len(taskaspects),figsize=(len(taskaspects)*12,1*8))
        
        for cx,comparison in enumerate(taskaspects):
            # if cx>1: continue
            if len(taskaspects)==1: axs = ax
            else: axs = ax[cx]
            # times = acrossdecoder_nullspace[1].analogsignals[0].times

            acrossdecoder_nullspaces,ranks = pickle.load(open(cacheprefix+'subspaces/nullspace,recurrent,decodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
            print('ranks:',ranks.mean(axis=1).ravel())
            print(comparison,len(acrossdecoder_nullspaces),len(acrossdecoder_nullspaces[0]),len(acrossdecoder_nullspaces[0][0]))
            for px in range(n_neurons):
                if px==0:
                    acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
                else:
                    acrossdecoder = acrossdecoder_nullspaces[px-1]
                    colorlists[1][1] = np.array([1,0,0.5])*(px/n_neurons/2)+np.array([0.25,0,0.12])
                figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs,colorlist=colorlists[px>0],plottrain=False,onlysem=True,\
                    smooth=[1],label='r %d'%(n_neurons-px),lw=[4,1][px>0])       # plot the performance
                # axs.plot(acrossdecoder[1][:,0],lw=0.5,color=colorlists[px>0][1],alpha=0.5)       # plot the performance
                # axs.legend(frameon=False)
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.5)
            
            axs.set_title(comparison)


        fig.suptitle('%s - original (orange), recurrent nullspaces (red)'%(dn))



        save = 1
        if save:
            fig.savefig(resultpath+'decoders,over,dv-nullspace,recurrent,projectedactivities-%s-%s'%(dn,continuous_method)+ext)












def subspaces_neuralparticipation(dn,block):
    # explores activity subspaces labelling each dimension with cell types and layer depths
    examine = 'allexpcond'
    layeredfilename = 'all'
    stimulusspecificity = 'all'
    
    
    blv,bla = preprocess.getorderattended(dn)
    
    comparisongroups  = [ \
                              [  [[ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                                 [[ [blv[1]], [45],    [5000] ], [ [blv[1]],[135],     [5000] ]  ] ],\
                              [  [[ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                                 [[ [bla[1]],  [45],[5000] ], [ [bla[1]],   [45],[10000] ]    ] ],\
                              [  [[ [blv[1]],  [],[] ],   [ [bla[1]],   [], [] ] ],\
                                 [[  blv,  [],[] ],     [ bla,   [],[] ]   ] ], \
                              [[],[]] \
                         ]
    taskaspects = ['visual','audio','context','choice']
    
    
    

    ixct = preprocess.getcelltypes(block)
    ixnarrow,ixbroad,ixsuperf,ixinput,ixdeep5,ixdeep6 = ixct
    n_neuron = np.count_nonzero(np.array(ixnarrow)+1)+np.count_nonzero(np.array(ixbroad)+1)
    ctcolors = ['deepskyblue','lightcoral','darkviolet','deeppink','blue','red','navy','darkred']
    ctlabels = ['superf\nbroad','superf\nnarrow','input\nbroad','input\nnarrow','deep5\nbroad','deep5\nnarrow','deep6\nbroad','deep6\nnarrow']
    



    # collect decision boundary normal vectors

    c_db = []          # vector components (coefficients) of the decision normal in the activity space
    
    for cx,comparison in enumerate(taskaspects):
        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%(examine,layeredfilename,dn,continuous_method,comparison,stimulusspecificity),'rb'))
        wx = int((len(acrossdecoder)-7)/n_neuron)
        c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
        
    c_db = np.array(c_db) # [comparisongroup,neurons,trajectory,stats]

    print(dn,c_db.shape)

    forwardtime = 1000*pq.ms          # only first 1 sec

    # this is a comparison group by neuron by   stats matrix, use only mean from stats
    c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimstart_idx']+int((forwardtime/T['dt']).magnitude),0].mean(axis=2)
    
    print('c_db_means',c_db_means.shape)
    
    if 1:
        pickle.dump(c_db_means,open(cacheprefix+'subspaces/subspacecontrib,coeffs-%s-%s'%(continuous_method,dn),'wb'))           # save the coefficients
        return

#    norms_c_db = np.linalg.norm(c_db_means,axis=0)

    limit = np.max(np.abs(c_db_means))
    orderall = np.argsort(np.abs(c_db_means[2,:]))[::-1]
    
    ctgroups = []
    for i in range(4):
        ctgroups.append(np.intersect1d(ixct[1],ixct[i+2]).astype('int16'))  # add broad and then narrow for each group
        ctgroups.append(np.intersect1d(ixct[0],ixct[i+2]).astype('int16'))
    
    print(ctgroups)
    
    ctmeans = np.zeros((len(taskaspects),8))
    orderct = np.zeros((len(taskaspects),8),dtype='int16')
    for cx in range(len(taskaspects)):
        for i in range(8):
            if ctgroups[i].size>0:
                ctmeans[cx,i] = np.mean(np.abs(c_db_means[cx,ctgroups[i]]))

        orderct[cx,:] = np.argsort(ctmeans[cx,:])[::-1]
    
#    orderct
    print(orderct)
    
    colors=['navy','darkgreen','mediumvioletred','darkorange']

    # plot all cells
    fig,ax = plt.subplots(len(taskaspects),2,figsize=(32,32))
    for cx,comparison in enumerate(taskaspects):
        axs = ax[cx,0]
        axs.bar(np.arange(n_neuron),c_db_means[cx,orderall],color=colors[cx])
        axs.set_ylabel(taskaspects[cx])
        axs.set_ylim(-limit*1.1,limit*1.1)
        if cx==0: axs.set_title('neurons ordered by %s'%taskaspects[2])

        axs = ax[cx,1]
        for io in range(8):
            i = orderct[cx,io]
            if ctgroups[i].size>0:
                axs.bar(io*20+np.arange(len(ctgroups[i]))/2, c_db_means[cx,ctgroups[i]], width=0.3, color=ctcolors[i])
                axs.text(io*20,limit*0.8,ctlabels[i], color=ctcolors[i],fontsize=14)
        if cx==0: axs.set_title('cell types, ordered by group participation')
        axs.set_ylim(-limit*1.1,limit*1.1)
        axs.set_xticks([])


    fig.suptitle(dn+' neural participation in decoders, %d ms'%(forwardtime))

    save = 0
    if save:
        fig.savefig(resultpath+'neuralparticipation-%s_%s_%s-%dms%dms_%s'%(examine,layeredfilename,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
    
    
    
    
    
    
    
    
    
    
    
    
    
def subspaces_orthogonalcontrol_shuffle(dn,block):
    
    # do a gradual decrease of trials corresponding to the task variable, and 
    # do several label shuffling and a random decoder
    # get mean and std vector of noise distributions
    # compare to mean std of task variable
    
    recalculate = 0
    doplot = 1
    
    # comparisongroups = [ [ [2,4], [], [] ]              ]#      fully random labels
    # allresponses = preprocess.collect_stimulusspecificresponses(block, comparisongroups)

    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [   [ [ [2,4], [45],    [] ], [ [2,4],[135],     [] ]  ], \
                            [ [ [2,4],  [],[5000] ],  [ [2,4],   [],[10000] ]    ],    \
                            [ [ [blv[1]],  [],[] ],      [ [bla[1]],   [],[] ]    ],     \
                            [[],[]] \
                        ]

    taskaspects = ['visual','audio','context','choice']

    width = 50*pq.ms

    n_neurons = block.segments[0].analogsignals[0].shape[1]

    n_resample = 10
    grads = [0,1,2,3,4,5]
    n_grad = len(grads)-1
    
    if recalculate:
        for cx,comparison in enumerate(taskaspects):
    
            # collect neural responses
            if not comparison=='choice': # visual, audio, context:
                acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx])
            else:  # choice:
                acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)

            # n_grad = len(acrossresponses[0])+len(acrossresponses[1]) # this will allow for the first grads trials one by one, not proportions
            acrossdecoders = []
            for gx,g in enumerate(grads):   # we either do 1 resample run with gradually increasing number of trials to shuffle, or just all shuffled but more runs for averaging
                if n_resample<=1:
                    print(comparison,'%d/%d'%(g,n_grad))
                    acrossdecoders.append( nedi.get_responsedecoding(acrossresponses,width=width,cv=True,shuffle=g/n_grad) )
                elif g==grads[-1]:
                    for rx in range(n_resample):
                        print('%s, resample %d/%d'%(comparison,rx+1,n_resample))
                        acrossdecoders.append( nedi.get_responsedecoding(acrossresponses,width=width,cv=True,shuffle=g/n_grad) )
                
            if n_resample<=1:
                pickle.dump(acrossdecoders,open(cacheprefix+'subspaces/shuffledecoders,gradual-%s-%d-%s_%s.pck'%(comparison,n_grad,continuous_method,dn),'wb'))
            else:
                pickle.dump(acrossdecoders,open(cacheprefix+'subspaces/shuffledecoders,resampled-%s-full,r%d-%s_%s.pck'%(comparison,n_resample,continuous_method,dn),'wb'))

    # else:
    #     acrossdecoders = pickle.load(open(cacheprefix+'subspaces/shuffledecoders,gradual-%s-%d-%s_%s.pck'%(comparison,n_grad,continuous_method,dn),'rb'))
            





    

    if doplot:
        if 1:       # fully randomized n_resample rolls
            taskcolors = ['navy','darkgreen','mediumvioletred','darkorange']

            acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[2])    # get all multimodal trials involved
            n_trials = len(acrossresponses[0]) + len(acrossresponses[1])
            mul = np.sqrt(10/n_trials)   # CV to leave one out s.e.m. of runs

            acrossdecoders = [] # get together all random runs
            for cx,comparison in enumerate(taskaspects):
                acrossdecoders_taskaspect = pickle.load(open(cacheprefix+'subspaces/shuffledecoders,resampled-%s-full,r%d-%s_%s.pck'%(comparison,n_resample,continuous_method,dn),'rb'))
                acrossdecoders.extend(acrossdecoders_taskaspect)



            fig,ax = plt.subplots(2,len(taskaspects), figsize=(len(taskaspects)*9,2*9))


            colors=['dodgerblue','dodgerblue','red','gold']
            
            
            D = np.array([acrossdecoders[rx][1] for rx in range(n_resample*len(taskaspects))])
            m = D.mean(axis=0)                        #  mean across resample randomized trial labels
            e = 2*D.std(axis=0)/np.sqrt(n_resample*len(taskaspects))   # 2 sem across resample randomized trial labels

            axs = ax[0,0]
            # mean of cv-means over resamples
            axs.plot(acrossdecoders[0][1].times,m[:,0], lw=3, color=colors[0],\
                label='randomization mean and 2 sem over\nresamples of cv-means')
            # with 2 sem around mean of cv-means over resamples
            axs.fill_between(acrossdecoders[0][1].times, m[:,0]-e[:,0], m[:,0]+e[:,0], color=colors[1], alpha=0.7)
            axs.legend(frameon=False)
            axs.set_ylabel('randomization')

            axs = ax[0,1]
            # mean of cv
            axs.plot(acrossdecoders[0][1].times,m[:,0]+m[:,2], lw=2, color=colors[2],\
                label='randomization mean and 2 sem over\nresamples of 10 cv-means+cv-2sems')
            # 2sem cv
            axs.fill_between(acrossdecoders[0][1].times, m[:,0]+m[:,2]-e[:,0]+e[:,2], m[:,0]+m[:,2]+e[:,0]+e[:,2], color=colors[2], alpha=0.5)
            axs.legend(frameon=False)


            axs = ax[0,2]
            # mean of cv
            axs.plot(acrossdecoders[0][1].times,m[:,0]+m[:,2]*mul, lw=2, color=colors[3],\
                label='randomization mean and 2 sem over\nresamples of leave one out cv-means+cv-2sems')
            # 2sem cv
            axs.fill_between(acrossdecoders[0][1].times, m[:,0]+m[:,2]*mul-e[:,0]+e[:,2]*mul, m[:,0]+m[:,2]*mul+e[:,0]+e[:,2]*mul,\
                color=colors[3], alpha=0.5)
            axs.legend(frameon=False)




            for cx,comparison in enumerate(taskaspects):
                # cv shade:
                # axs.fill_between(acrossdecoders[0][1].times, m[:,0]-m[:,2], m[:,0]+m[:,2], color=taskcolors[cx], alpha=0.2)


                axs = ax[1,cx]
                acrossdecoders_full = pickle.load(open(cacheprefix+'subspaces/shuffledecoders,gradual-%s-%d-%s_%s.pck'%(comparison,n_grad,continuous_method,dn),'rb'))

                axs.plot(acrossdecoders_full[0][1].times,acrossdecoders_full[0][1][:,0], lw=2, color=taskcolors[cx])
                # 2std cv
                axs.fill_between(acrossdecoders_full[0][1].times, (acrossdecoders_full[0][1][:,0]-acrossdecoders_full[0][1][:,2]).squeeze(),\
                                                                  (acrossdecoders_full[0][1][:,0]+acrossdecoders_full[0][1][:,2]).squeeze(),\
                                    color=taskcolors[cx], alpha=0.5)
                
                
                
                # the last item uses leave one out multipler for theoretical best s.e.m. estimation
                ys = np.array([m[:,0].mean(axis=0), (m[:,0]+e[:,0]).mean(axis=0), (m[:,0]+m[:,2]+e[:,0]+e[:,2]).mean(axis=0),\
                     (m[:,0]+m[:,2]*mul+e[:,0]+e[:,2]).mean(axis=0)])
                # print(ys)
                # print(     np.concatenate([m[:,0],m[:,2],e[:,0],e[:,2]]).mean(axis=0)       )
                lbs = ['random mean','random 2 s.e.m.','random mean+2sem of 10 CVs mean+sem','random mean+2sem of leave one out CVs mean+2sem']
                colors=['dodgerblue','dodgerblue','red','gold']
                lws = [2,1,2,4]
                for k in range(len(ys)):
                    axs.plot([T['starttime'],T['endtime']],[ys[k],ys[k]],'--',color=colors[k],lw=lws[k],label=lbs[k])
                axs.legend(frameon=False)

                if cx==0: axs.set_ylabel('orig. decoders +\ntime-averages of shuffle')
                axs.set_title(comparison)



                for k in [0,1]:
                    axs = ax[k,cx]
                    axs.set_ylim([0.45,1.05])
                    figs.setxt(axs)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    figs.plottoaxis_chancelevel(axs,0.5)


            ax[0,3].remove()


            fig.suptitle('%s, %d neurons baseline decoder accuracies:\nrandomized labels, mean of %d rolls'%(dn,n_neurons, n_resample*len(taskaspects)))




            save = 0
            if save:
                fig.savefig(resultpath+'shuffle,x%d,noise-%s-%dms,%dms_%s'%(n_resample,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)





        if 0:     # gradual plot


            fig,ax = plt.subplots(2,len(taskaspects), figsize=(len(taskaspects)*9,2*9))

            for cx,comparison in enumerate(taskaspects):

                acrossdecoders = pickle.load(open(cacheprefix+'subspaces/shuffledecoders,gradual-%s-%d-%s_%s.pck'%(comparison,n_grad,continuous_method,dn),'rb'))
                print(len(acrossdecoders),len(acrossdecoders[0]),len(acrossdecoders[0][0]))
                for gx,g in enumerate(grads):


                    # s = neph.smooth(acrossdecoders[gx][1][:,0],mode='same')
                    # axs.plot(acrossdecoders[gx][1].times,s,linewidth=3,label='%d/%d'%(g,n_grad))
                    axs = ax[0,cx]
                    axs.plot(acrossdecoders[gx][1].times,acrossdecoders[gx][1][:,0],linewidth=0.5,label='%d/%d'%(g,n_grad))
                    
                    # axs.fill_between(acrossdecoders[gx][1].times, acrossdecoders[gx][1][:,0].squeeze()-acrossdecoders[gx][1][:,2].squeeze(),\
                    #                                               acrossdecoders[gx][1][:,0].squeeze()+acrossdecoders[gx][1][:,2].squeeze(), alpha=0.3)

                    axs = ax[1,cx]
                    axs.plot(acrossdecoders[gx][1].times,neph.smooth(acrossdecoders[gx][1][:,0],20,'same'),linewidth=0.5,label='%d/%d'%(g,n_grad))

                for k in [0,1]:
                    axs = ax[k,cx]
                    axs.legend(frameon=False)
                    axs.set_ylim([0.45,1.05])
                    # axs.set_yticks([0.5,1.0])
                    figs.setxt(axs)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    figs.plottoaxis_chancelevel(axs,0.5)

            fig.suptitle(dn)

            save = 0
            if save:
                fig.savefig(resultpath+'randomlabels,noiseceiling-%s-%dms,%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)



        if 0:
            n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
            n_trajectory -= (width/T['dt']).magnitude.astype(np.int16)
            
            c_db_matrix_shuffle = np.zeros((n_neurons,n_trajectory,n_shuffle))       # prepare the dbnv vector holder
            for i in range(n_shuffle):
                shuffledecoder = shuffledecoders[i]
            
            
                wx = int((len(shuffledecoder)-7)/n_neurons)            # we are only concerned with the coefficients coming after the train test (2) and PCs (5).
                width = wx * T['dt']    # we will retroactively do a similar feature with decoder
                
                coeff = np.reshape(np.array(shuffledecoder[7:]), (wx,n_neurons,shuffledecoder[7].shape[0],shuffledecoder[7].shape[1]) ).mean(axis=0)
            
                c_db_matrix_shuffle[:,:,i] = coeff[:,:,0]      # neurons by trajectory
            ln = np.linalg.norm(c_db_matrix_shuffle,axis=0)
            shuffle_ms = np.array([ ln.mean(axis=1), ln.std(axis=1), ln.std(axis=1)/np.sqrt(n_shuffle) ])


            times = shuffledecoders[0][1].times



            # now load multiple real decoders and compare DV vector mean and std

            taskaspects = ['visual','audio','context','choice']
            taskcolors = ['navy','darkgreen','mediumvioletred','darkorange']

            fig,ax = plt.subplots(2,len(taskaspects),figsize=(12*len(taskaspects),8*2))

            for cx,comparison in enumerate(taskaspects):
                acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))

                c_db_matrix = np.zeros((n_neurons,n_trajectory))       # prepare the dbnv vector holder

                wx = int((len(acrossdecoder)-7)/n_neurons)            # we are only concerned with the coefficients coming after the train test (2) and PCs (5).
                width = wx * T['dt']    # we will retroactively do a similar feature with decoder
                
                coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)

                c_db_matrix = coeff[:,:,0]      # neurons by trajectory
                ln = np.linalg.norm(c_db_matrix,axis=0)
                decoder_ms = ln

            
                
                axs = ax[0,cx]
                
                axs.plot(times, shuffle_ms[0],color='red')
                axs.fill_between(times, shuffle_ms[0]-shuffle_ms[1], shuffle_ms[0]+shuffle_ms[1], color='red',alpha=0.3,  label='shuffled')
            
                axs.plot(times, decoder_ms, color=taskcolors[cx])
                
                if cx==0: axs.legend(frameon=False)
            
                if cx==0: axs.set_ylabel('norm of decision vector (decoder coeffs.)')
                axs.set_title(comparison)
                
                axs.set_ylim([0,1.5])
                figs.setxt(axs)
                figs.plottoaxis_stimulusoverlay(axs,T)


                axs = ax[1,cx]

                figs.plottoaxis_decoderrocauc(shuffledecoders[0],axs,colorlist=['purple','red'],plottrain=True)
                figs.setxt(axs)
                figs.plottoaxis_chancelevel(axs)
                axs.set_ylim([0.35,1.05])
                axs.set_yticks([0.5,0.75,1.])
                figs.plottoaxis_stimulusoverlay(axs,T)



            
            fig.suptitle(dn)

            save = 0
            if save:
                fig.savefig(resultpath+'shuffle,noiseceiling-%s-%dms,%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                
        


    print('shuffle done')



    return













def subspaces_leaveoneout(dn,block):
    # do a predicted class probability of leave one neuron out 

    recalculate_task = 0
    recalculate_shuffle = 0

    doplot = 1


    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [   [ [ [2,4], [45],    [] ], [ [2,4],[135],     [] ]  ], \
                            [ [ [2,4],  [],[5000] ],  [ [2,4],   [],[10000] ]    ],    \
                            [ [ [blv[1]],  [],[] ],      [ [bla[1]],   [],[] ]    ],     \
                            [[],[]] \
                        ]

    taskaspects = ['visual','audio','context','choice']

    width = 50*pq.ms
    wx = 5

    n_neurons = block.segments[0].analogsignals[0].shape[1]

    n_resample = 40
    

    if recalculate_task:
        acrossdecoders = []
        for cx,comparison in enumerate(taskaspects):
    
            # collect neural responses
            if not comparison=='choice': # visual, audio, context:
                acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx])
            else:  # choice:
                acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
            n_trials = len(acrossresponses[0]) + len(acrossresponses[1])

            # rr = n_trials: leave one out CV! (takes a very long time, about 8-12x more than 10 fold CV!!!)
            # acrossdecoder = nedi.get_responsedecoding(acrossresponses,width=width,cv=True,rr=n_trials)
            # pickle.dump(acrossdecoder,open(cacheprefix+'subspaces/leaveoneout_%s-%s_%s.pck'%(comparison,continuous_method,dn),'wb'))


            class_probabilities,accuracies = nedi.get_pointwisedecoderprobability(acrossresponses,width=width,rr=n_trials)
            acrossdecoder_proba = class_probabilities,accuracies
            pickle.dump(acrossdecoder_proba,open(cacheprefix+'subspaces/leaveoneout,proba_%s-%s_%s.pck'%(comparison,continuous_method,dn),'wb'))



    if recalculate_shuffle:
        # random shuffle control:
        acrossdecoders = []
        for rx in range(n_resample):
            print('%s, resample %d/%d'%(comparison,rx+1,n_resample))
            acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[2])   # choose all in blocks 2 4 multimodal
            acrossdecoders.append( nedi.get_responsedecoding(acrossresponses,width=width,rr=cv,shuffle=1.) )
        pickle.dump(acrossdecoders,open(cacheprefix+'subspaces/shuffledecoders,resampled-full,r%d-%s_%s.pck'%(n_resample,continuous_method,dn),'wb'))





    if doplot:
        
        if 1:   # with probabilities
            taskcolors = ['navy','darkgreen','mediumvioletred','darkorange']


            fig,ax = plt.subplots(2,len(taskaspects),figsize=(len(taskaspects)*9,2*9))

            for cvtx in [0,1]:
                for cx,comparison in enumerate(taskaspects):
                    print(dn,comparison)

                    if cvtx==0:
                        class_probabilities,accuracies = pickle.load(open(cacheprefix+'subspaces/leaveoneout,proba_%s-%s_%s.pck'%(comparison,continuous_method,dn),'rb'))
                        times = block.segments[0].analogsignals[0].times[:-wx]
                        accuracies[:,0] = 1-accuracies[:,0]; print('adjusting for flipped class probabilities...')
                    elif cvtx==1:
                        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
                        # print(len(acrossdecoder),len(acrossdecoder[0]),len(acrossdecoder[0][0]))

                    axs = ax[cvtx,cx]

                    if cvtx==0:
                        axs.plot(times,accuracies[:,0],color=taskcolors[cx])
                        axs.fill_between(times,accuracies[:,0]-accuracies[:,2],accuracies[:,0]+accuracies[:,2],\
                                color=taskcolors[cx],alpha=0.5)
                        axs.set_xlim([times[0],times[-1]])
                        axs.set_ylim([0.45,1.01])
                    if cvtx==1:
                        figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs,colorlist=['black',taskcolors[cx]], onlysem=True,lw=1)
                    
                    figs.setxt(axs)
                    figs.plottoaxis_chancelevel(axs,0.5)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    
                    if cvtx==0: axs.set_title(comparison)
                    if cvtx==0 and cx==0: axs.set_ylabel('leave one out')
                    elif cx==0: axs.set_ylabel('10 fold CV')

            fig.suptitle(dn+' leave one out and 10 fold CV comparison')

            save = 1
            if save:
                fig.savefig(resultpath+'leaveoneout,bandtest-%s-%dms,%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                




        if 0:   # with rounded predictions - gives huge error bands
            taskcolors = ['navy','darkgreen','mediumvioletred','darkorange']


            fig,ax = plt.subplots(2,len(taskaspects),figsize=(len(taskaspects)*9,2*9))

            for cvtx in [0,1]:
                for cx,comparison in enumerate(taskaspects):

                    if cvtx==0:
                        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/leaveoneout_%s-%s_%s.pck'%(comparison,continuous_method,dn),'rb'))
                    elif cvtx==1:
                        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))

                    print(dn,comparison,len(acrossdecoder),len(acrossdecoder[0]),len(acrossdecoder[0][0]))

                    axs = ax[cvtx,cx]
                    # axs.plot(acrossdecoder[1].times,acrossdecoder[1][:,0],color=taskcolors[cx])
                    # axs.fill_between(acrossdecoder[1].times,(acrossdecoder[1][:,0]-acrossdecoder[1][:,2]).squeeze(),(acrossdecoder[1][:,0]+acrossdecoder[1][:,2]).squeeze(),\
                    #            color=taskcolors[cx],alpha=0.5)
                    figs.plottoaxis_decoderrocauc(acrossdecoder[:2],axs,colorlist=['black',taskcolors[cx]])
                    
                    figs.setxt(axs)
                    figs.plottoaxis_chancelevel(axs,0.5)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    
                    if cvtx==0: axs.set_title(comparison)
                    if cvtx==0 and cx==0: axs.set_ylabel('leave one out')
                    elif cx==0: axs.set_ylabel('10 folc CV')

            fig.suptitle(dn+' leave one out and 10 fold CV comparison')



        















def subspaces_behaviour_ACC(dn,block):
    # learn stimulus DV
    # project activity to DV
    # show correct and error trials; are error trials closer to discrimination threshold than correct?

    doplot = 0 or globaldoplot

    n_neuron = block.segments[0].analogsignals[0].shape[1]

    blv,bla = preprocess.getorderattended(dn)
    # comparisongroups  = [\
    #                         [ [ [blv[1]], [45],      [] ], [ [blv[1]],   [135],      [] ]  ], \
    #                         [ [ [blv[1]],   [],  [5000] ], [ [blv[1]],      [], [10000] ]  ], \
    #                         [ [ [bla[1]], [45],      [] ], [ [bla[1]],   [135],      [] ]  ], \
    #                         [ [ [bla[1]],   [],  [5000] ], [ [bla[1]],      [], [10000] ]  ], \

    #                         [ [ [blv[1]], [45],  [5000] ], [ [blv[1]],   [135], [10000] ]  ], \
    #                         [ [ [blv[1]], [45], [10000] ], [ [blv[1]],   [135],  [5000] ]  ], \
    #                         [ [ [bla[1]], [45],  [5000] ], [ [bla[1]],   [135], [10000] ]  ], \
    #                         [ [ [bla[1]], [135], [5000] ], [ [bla[1]],    [45], [10000] ]  ], \
    #                     ]

    # taskaspects = ['visual,av','audio,av','visual,aa','audio,aa',\
    #                'gonogo-congruent,av','gonogo-conflict,av','gonogo-congruent,aa','gonogo-conflict,aa']



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


    width = 50*pq.ms
    wx = int(width/T['dt'])




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


    cmap = figs.getcorrcolormap('correlation')




    if doplot:

        if 0:

            fig,ax = plt.subplots(2,len(taskaspects)//2, figsize=(len(taskaspects)//2*9,2*8))
            
            colors = [ ['darkgreen','red'],   ['dodgerblue','orange'] ]
            label_behav = ['correct','error']
            label_action = ['go','nogo']

            for cx,comparison in enumerate(taskaspects):
                
                
                axs = ax[cx//4,cx%4]
            
                # project activity onto the DV

                
                for clx,signal in enumerate(label_action):    # choose go and nogo signals respectively

                    X = np.array(activity[cx][clx])                        # [trials,trajectory,neurons]
                    Y = np.array( [  X[:,t,:] @ W[cx,:,t]    for t in range(W.shape[2])     ] ).T         # [trials,trajectory]
                    Y[:,:T['stimstart_idx']] = np.nan        # no valid DV projection possible before stimulus

                    P = perftask[cx][signal]
                    P = P[P.notna()]
                    mask_correct = P==1     # correct trials
                    mask_error = P==0       # error trials

                    # print('X',X.shape, 'Y',Y.shape, sum(mask_correct), sum(mask_error))

                    for ex,mask in enumerate([mask_correct,mask_error]):

                        m = Y[mask,:].mean(axis=0)
                        s = Y[mask,:].std(axis=0)
                        e = s/np.sqrt(Y[mask,:].shape[0])

                        axs.plot(times,m,color=colors[clx][ex], label='%s %s %d'%(label_action[clx],label_behav[ex],sum(mask)))
                        axs.fill_between(times,m-e,m+e,color=colors[clx][ex], alpha=0.3)
                        axs.fill_between(times,m-s,m+s,color=colors[clx][ex], alpha=0.1)

                axs.set_ylim(-1.24,1.24)
                axs.set_xlim(times[0],times[-1])
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0)

                axs.legend(frameon=False)
                if cx==0 or cx==4:
                    axs.set_ylabel(['relevant, irrelevant\nDV projection','congruent, conflicting\nDV projection'][cx//4])
                axs.set_title(comparison,fontsize=20)
                if cx>3: axs.set_xlabel('time from stimulus onset [ms]')
                    

            fig.suptitle(dn+' %d neurons, activity projected onto DV at each timepoint, rel-irrel,cong-confl*contexts, shades: sem, std'%n_neuron)


            save = 0 or globalsave
            if save:
                fig.savefig(resultpath+'ricc-dvprojection,timeresolved_%s_%s-%dms'%(dn,continuous_method,T['dt'].magnitude)+ext)









        if 1:

            fig,ax = plt.subplots(1,2, figsize=(2*9,1*8))
            
            colors = [ 'seagreen','fuchsia' ]
            colors = [ ['navy','red'],   ['darkgreen','red'] ]
            # colors = [ ['darkgreen','red'],   ['dodgerblue','orange'] ]

            label_behav = ['correct','error']
            label_action = ['go','nogo']



            for cix,cx in enumerate([0,3]):            # take only the relevants

                Y = np.zeros((0,len(times)))      # grow this array
                mask_correct = np.zeros((0),dtype=np.int16)
                mask_error = np.zeros((0),dtype=np.int16)
                
            
                # project activity onto the DV

                
                for clx,signal in enumerate(label_action):    # choose go and nogo signals respectively

                    X = np.array(activity[cx][clx])                        # [trials,trajectory,neurons]
                    Y_aux = np.array( [  X[:,t,:] @ W[cx,:,t]    for t in range(W.shape[2])     ] ).T         # [trials,trajectory]
                    Y_aux *= 2*(clx-0.5)        # symmetrise around the decision boundary
                    Y = np.vstack((Y,Y_aux))
                    
                    P = perftask[cx][signal]
                    P = P[P.notna()]
                    mask_correct_aux = P==1     # correct trials
                    mask_error_aux = P==0       # error trials

                    mask_correct = np.hstack((mask_correct,mask_correct_aux))
                    mask_error = np.hstack((mask_error,mask_error_aux))

                    # print(cx,clx,'    ',Y_aux.shape,Y.shape,mask_correct.shape,mask_error.shape)
            
                axs = ax[cix]
                Y[:,:T['stimstart_idx']] = np.nan        # no valid DV projection possible before stimulus


                # print('X',X.shape, 'Y',Y.shape, sum(mask_correct), sum(mask_error))

                for ex,mask in enumerate([mask_correct,mask_error]):

                    m = Y[mask,:].mean(axis=0)
                    s = Y[mask,:].std(axis=0)
                    e = s/np.sqrt(Y[mask,:].shape[0])

                    axs.plot(times,m,color=colors[cix][ex], label='%s %d'%(label_behav[ex],sum(mask)))
                    axs.fill_between(times,m-e,m+e,color=colors[cix][ex], alpha=0.3)
                    axs.fill_between(times,m-s,m+s,color=colors[cix][ex], alpha=0.1)

                axs.set_ylim(-1.5,1.5)
                axs.set_xlim(times[0],times[-1])
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0)

                axs.legend(frameon=False)
                if cx==0: axs.set_ylabel('go-nogo-symmetrised DV projection')
                axs.set_title('relevant, %s'%taskaspects[cx],fontsize=20)
                axs.set_xlabel('time from stimulus onset [ms]')
                    

            fig.suptitle(dn+' %d neurons, activity projected onto DV at each timepoint, symmetric for go-nogo, shades: sem, std'%n_neuron)



            save = 0 or globalsave
            if save:
                fig.savefig(resultpath+'ricc-dvprojection,timeresolved,symm+avaa_%s_%s-%dms'%(dn,continuous_method,T['dt'].magnitude)+ext)





    return













def subspaces_timeresolved(dn,block):
    # show DV geometric configuration changes in time
    # absolute direction is PCA 1-3
    # DV direction is decision boundary normal
    # magnitude is unit vector, times separation strength: d'

    recalculate = 0 or globalrecalculate
    doplot = 1 or globaldoplot


    n_neuron = block.segments[0].analogsignals[0].shape[1]


    blv,bla = preprocess.getorderattended(dn)

    comparisongroups  = [\
                            [ [ [blv[1]], [45],      [] ], [ [blv[1]],   [135],      [] ]  ], \
                            [ [ [blv[1]],   [],  [5000] ], [ [blv[1]],      [], [10000] ]  ], \
                            [ [blv[1]], [] ],\
                            [ [ [bla[1]], [45],      [] ], [ [bla[1]],   [135],      [] ]  ], \
                            [ [ [bla[1]],   [],  [5000] ], [ [bla[1]],      [], [10000] ]  ], \
                            [ [bla[1]], [] ] \
                        ]

    taskaspects = ['visual,av','audio,av','choice,av','visual,aa','audio,aa','choice,aa']

    width = 50*pq.ms
    wx = int((width/T['dt']).magnitude)

    acrossdecoders = []
    W = []
    for cx,comparison in enumerate(taskaspects):
        # collect activity
        if not comparison[:6]=='choice': # visual, audio, context:
            acrossresponses = preprocess.collect_stimulusspecificresponses(block, comparisongroups[cx])
        else:  # choice:
            acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block, dn, onlyinblock=comparisongroups[cx][0])

        print(cx,comparison)
        # activity.append(   )
        
        # collect decoders, already calculated, in this function we just display a different crossection:
        # acrossdecoder = pickle.load(open(cacheprefix+'acc/stimulus,relevant-irrelevant_%s_%s-%s.pck'%(comparison,dn,continuous_method),'rb'))
        
        if recalculate:
            acrossdecoder = nedi.get_responsedecoding(acrossresponses)
            pickle.dump(acrossdecoder,open(cacheprefix+'subspaces/subspaces,timeresolved-%s-%s_%s.pck'%(comparison,dn,continuous_method),'wb'))
        else:
            acrossdecoder = pickle.load(open(cacheprefix+'subspaces/subspaces,timeresolved-%s-%s_%s.pck'%(comparison,dn,continuous_method),'rb'))
        acrossdecoders.append(acrossdecoder)
        w,w_mean = nedi.get_decoder_dv(acrossdecoder, n_neuron, n_preweightelements=7, timeaverage_start_dx=150, timeaverage_end_dx=150+150)
        W.append(w)



    W = np.array(W) # [taskaspects,neurons,trajectory]

    # print(len(acrossdecoder),len(acrossdecoder[0]))

    times = block.segments[0].analogsignals[0].times[:-wx]

    # find the orthogonal axes
    Q,R = np.linalg.qr(W[:,:,T['stimstart_idx']:T['stimend_idx']].mean(axis=2))     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix

    # Z = np.einsum('ij,jkl',Q.T,W)
    # Z = np.dot(Q.T,W) 
    # print(Z.shape)






    if doplot:


        colors_dv = ['navy','darkgreen','darkorange']
        labels = ['visual context','audio context']
        # angle = np.convolve(60 + np.cumsum(np.cumsum(0.2*np.random.randn(len(times)))),np.ones(15)/15,'same')
        angle = 60 + np.arange(len(times))/10






        fig = plt.figure(figsize=(2*8,2*8))

        axs1 = fig.add_subplot(2,2,1)
        axs2 = fig.add_subplot(2,2,2)

        axs3 = fig.add_subplot(2,2,3,projection='3d')
        axs4 = fig.add_subplot(2,2,4,projection='3d')



        def animate(t):
            vline = []
            for aix,axs in enumerate([axs1,axs2]):
                axs.clear()

                # if t==0:
                for zx in [0,1,2]:
                    axs.plot(times,acrossdecoders[aix*3+zx][1],color=colors_dv[zx])
            
                figs.setxt(axs)
                axs.set_ylim(0.45,1.05)
                axs.set_xlim(times[0],times[-1])
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.plottoaxis_chancelevel(axs,0.5)


                axs.set_title(labels[aix])
                axs.set_xlabel('time from stimulus onset [ms]')
                if aix==0: axs.set_ylabel('cv accuracy')

                # if t>0: vline.pop(0).remove()
                # print(t,times[t],vline)
                vline = axs.axvline(times[t],0,1.5,lw=2,alpha=0.7,color='black')


            for aix,axs in enumerate([axs3,axs4]):
                Qt,_ = np.linalg.qr(W[aix*3:(aix+1)*3,:,t])     # Q holds the orthonormal basis vectors as columns
                Qth = []
                for thx in np.arange(t,max(0,t-10),-1):
                    Qtha,_ = np.linalg.qr(W[aix*3:(aix+1)*3,:,thx])
                    Qth.append(Qtha)
                print(t,times[t])
                
                # chiral symmetry:
                # Qt = np.abs(Qt)
                # Qt = Qt * np.sign(Qt[0,0])

                # DVs of visual, audio, choice:
                axs.clear()
                for zx,color in enumerate(colors_dv):
                    if t==0: label=['visual DV','audio DV','choice DV'][zx]
                    else: label=None
                    a = axs.plot([0,Qt[0,zx]],[0,Qt[1,zx]],[0,Qt[2,zx]], lw=2, color=color, alpha=0.7,label=label)

                    # history
                    for ithx in range(len(Qth)):
                        af = axs.plot([0,Qth[ithx][0,zx]],[0,Qth[ithx][1,zx]],[0,Qth[ithx][2,zx]], lw=2, color=color, alpha=0.3)
                axs.set_xlim(-1,1)
                axs.set_ylim(-1,1)
                axs.set_zlim(-1,1)
                axs.view_init(elev = None, azim=angle[t])
                axs.legend(frameon=False)

        plt.ioff()

        # frames = np.arange(200,220)
        frames = np.arange(0,len(times),10)

        anime = FuncAnimation(fig,animate,frames=frames,interval=100)

        # f = resultpath+'animation.gif'
        f = resultpath+'animation-%s.gif'%dn
        writergif = animation.PillowWriter(fps=30) 
        anime.save(f, writer=writergif)




    return























# perform PCA


def subspaces_PCA(dn,block,returnall=False,preon=0):
    # iterate over variables
    # get activity grouped for classification
    # iterate over PCA components
    # create X from concatenating trials
    # do PCA on r components, return accuracy
    # plot number of principal components vs. accuracy

    recalculate = 0 or globalrecalculate
    doplot = 0 or globaldoplot
    
    n_neurons = block.segments[0].analogsignals[0].shape[1]


    # pcadimlist_all = [60]
    pcadimlist_all = [1,2,3,4,5,7,10,15,30,60]
    # pcadimlist_all = [15,50]
    # pcadimlist_all = [1,4]

    # avoid doing more pcs than neurons
    arg_max_r = np.where(np.array(pcadimlist_all)<n_neurons)[0][-1]
    pcadimlist = pcadimlist_all[:arg_max_r+1]
    pcadimlist.append(n_neurons)


    taskaspects = ['visual','audio','context','choice',\
                   'task,45-135','char,45-135','char,90-180','char,270-180']
    
    taskaspects = taskaspects[0:4]  # 0:4 = variables,
    # taskaspects = taskaspects[4:8]  # 4:8 = orientation characterization
    # taskaspects = ['visual','context']
    # taskaspects = ['task,45-135','char,45-135']
    # taskaspects = ['char,45-135','char,90-180','char,270-180']
    
    # where to run the PCA on:
    if preon==0:
        pcamean = 'sa'; pcameanlabel = 'spontaneous'
    else:
        pcamean = 'ea'; pcameanlabel = 'evoked'
    # other possibilities:
    # pcamean = 'allt'
    # pcamean = 'sa-perclass'

    # number of randomized orthogonal projection runs to average accuracies over
    n_ortho = 10




    
    blv,bla = preprocess.getorderattended(dn)
    pcagroups = [      [[2,4],[],[]]   ]
    comparisongroups  = [ \
                          [ [ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                          [ [ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ] ],\
                          [ [  blv,  [],[] ],   [ bla,   [], [] ] ],\
                          [[],[]], \
                          [ [ [blv[1]], [45],      [] ],  [   [blv[1]],[135],    [] ] ] ,\
                          [ [ [5],    [45],        [] ],  [        [5],[135],    [] ] ] ,\
                          [ [ [5],    [90],        [] ],  [        [5],[180],    [] ] ] ,\
                          [ [ [5],   [270],        [] ],  [        [5],[180],    [] ] ] \
                        ]


    if doplot:
        fig,ax = plt.subplots(len(taskaspects),len(pcadimlist)+2,figsize=(6*(len(pcadimlist)+2),6*len(taskaspects)))

    if returnall:
        data = []
    for cx,comparison in enumerate(taskaspects):

        # if cx in [4,5,6]: continue
        
        print(comparison)
        
        # choose sliding crossvalidation for the few trial sets
        if comparison in ['task,45-135','char,45-135','char,90-180','char,270-180']: cv = True
        else: cv = False
        
        # collect neural responses
        if not comparison=='choice': # visual, audio, context:
            acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],'and')
        else:  # choice:
            acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
#        print(np.array(acrossresponses[0]).shape)


        acrossdecoders = []
        acrossrandomdecoders = []
        if recalculate:

            for r in pcadimlist[:-1]: # run over pca dims, except the last (ne neurons)


                # projection matrices, V and V_r, for given r dimension
                if pcamean=='sa': # class invariant PCA
                    responses = preprocess.collect_stimulusspecificresponses(block,pcagroups,'and')
                    _,E,V = nedi.pca(np.array(responses[0])[:,T['start_idx']:T['stimstart_idx'],:].mean(axis=1),n_components=r)
                    print(dn,r,E)
                if pcamean=='ea': # class invariant PCA
                    responses = preprocess.collect_stimulusspecificresponses(block,pcagroups,'and')
                    _,E,V = nedi.pca(np.array(responses[0])[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1),n_components=r)
                    print(dn,r,E)


                # create a random projection as well, for baseline; for all variations of PCA, use the V_r matrix set up
                aux_random_acc = []
                for ox in range(n_ortho):
                    V_r = sp.stats.ortho_group.rvs(dim=responses[0][0].shape[1])[:r,:]
                    acrossrandompcaresponses = []
                    for aix in range(len(acrossresponses)):
                        sh = np.array(acrossresponses[aix]).shape
                        X_r = np.zeros( (sh[0],sh[1],r) )           # this is for a random projection
                        aux_X = np.array(acrossresponses[aix])
                        for t in range(len(acrossresponses[aix][0])):
                            X_r[:,t,:] = np.dot(aux_X[:,t,:],V_r.T)
                        randompcaresponse = []
                        for trx in range(sh[0]):  # get through the trials
                            randompcaresponse.append( neo.AnalogSignal( X_r[trx,:,:],\
                                                                 name=acrossresponses[0][0].name+'randomorthonormal%d'%r,\
                                                                 t_start=acrossresponses[0][0].t_start,\
                                                                 sampling_period=acrossresponses[0][0].sampling_period, units=acrossresponses[0][0].units   ) )
                        acrossrandompcaresponses.append(randompcaresponse)
                        
                    aux_random_acc.append(nedi.get_responsedecoding(acrossrandompcaresponses,cv=cv)) # perform classification on the orthonormal rotation
                
                acrossrandomdecoder = aux_random_acc[0]
                for ox in range(n_ortho): # make an average over random orthogonal rotations
                    for k in range(len(aux_random_acc[ox])):
                        if ox==0: acrossrandomdecoder[k] = aux_random_acc[ox][k]/n_ortho
                        else:     acrossrandomdecoder[k] += aux_random_acc[ox][k]/n_ortho
                
                acrossrandomdecoders.append(acrossrandomdecoder)    # collect for plot




                # and here do the PCA transforms
                acrosspcaresponses = []
                for aix in range(len(acrossresponses)):
                    sh = np.array(acrossresponses[aix]).shape
                    X = np.zeros( (sh[0],sh[1],r) )
                    X_a = np.zeros( (sh[0],sh[1],r) )
                    X_r = np.zeros( (sh[0],sh[1],r) )           # this is for a random projection
                    aux_X = np.array(acrossresponses[aix])
                    


                    # classinvariant pca:
                    if pcamean=='sa' or pcamean=='ea':      # Here V comes from class invariant pca from above
                        for t in range(len(acrossresponses[aix][0])): 
                            X[:,t,:] = np.dot(aux_X[:,t,:],V.T)


#                    # aggregated by class:
                    if pcamean=='sa-perclass':
                        X_tr,E,V = nedi.pca(aux_X[:,T['start_idx']:T['stimstart_idx'],:].mean(axis=1),n_components=r)
    #                    X_tr,E,V = nedi.pca(aux_X[:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=1),n_components=r)
                        for t in range(len(acrossresponses[aix][0])):
                            X[:,t,:] = np.dot(aux_X[:,t,:],V.T)
                    
                    # per each timepoint:
                    elif pcamean=='allt':
                        for t in range(sh[1]):
                            X_a[:,t,:],E,V = nedi.pca(aux_X[:,t,:],n_components=r)
                            X[:,t,:] = np.dot(aux_X[:,t,:],V.T)
                        



                    
#                    fig,ax = plt.subplots(1,3,figsize=(18,6))
#                    axs = ax[0]
#                    axs.plot(acrossresponses[0][0].times,np.array(acrossresponses[0]).mean(axis=(0,2)),color='navy')
#                    figs.plottoaxis_chancelevel(axs)
#                    figs.plottoaxis_stimulusoverlay(axs,T)
#
#                    axs = ax[1]
#                    for k in range(r):
#                        axs.plot(acrossresponses[0][0].times,X[:,:,k].mean(axis=(0)))
#                    figs.plottoaxis_chancelevel(axs)
#                    figs.plottoaxis_stimulusoverlay(axs,T)
#
#                    axs = ax[2]
#                    for k in range(r):                   # uncomment X_a above!!!!!
#                        axs.plot(acrossresponses[0][0].times,X_a[:,:,k].mean(axis=(0)))
#                    figs.plottoaxis_chancelevel(axs)
#                    figs.plottoaxis_stimulusoverlay(axs,T)
#
#                    return

                    
                    pcaresponse = []
                    for trx in range(sh[0]):  # get through the trials
                        pcaresponse.append( neo.AnalogSignal( X[trx,:,:],\
                                                             name=acrossresponses[0][0].name+'pca%d'%r,\
                                                             t_start=acrossresponses[0][0].t_start,\
                                                             sampling_period=acrossresponses[0][0].sampling_period, units=acrossresponses[0][0].units   ) )

                    acrosspcaresponses.append(pcaresponse)
                    
                acrossdecoder = nedi.get_responsedecoding(acrosspcaresponses,cv=cv) # perform classification on the selected principal components
                acrossdecoders.append(acrossdecoder)    # collect for plot



                # save for each pca dim and task variable
                pickle.dump((acrossdecoder,acrossrandomdecoder),open(cacheprefix+'pca/responsedecodes+O,pcadim%d,%s-%s-%s-%s.pck'%(r,pcamean,dn,continuous_method,comparison),'wb'))
            
            
            
            # now do this for the real number of neurons, without PCA, but also use random orthogonal matrix rotation of it for the random-baseline
            acrossdecoder = nedi.get_responsedecoding(acrossresponses,cv=cv)
            acrossdecoders.append(acrossdecoder)    # collect for plot
            
            acrossrandomresponses = []
            print(dn,n_neurons,'all neurons')
            V_r = sp.stats.ortho_group.rvs(dim=responses[0][0].shape[1])
            for aix in range(len(acrossresponses)):  # both classes
                sh = np.array(acrossresponses[aix]).shape
                X_r = np.zeros( (sh[0],sh[1],n_neurons) )           # this is for a random projection
                aux_X = np.array(acrossresponses[aix])
                for t in range(len(acrossresponses[aix][0])): 
                    X_r[:,t,:] = np.dot(aux_X[:,t,:],V_r.T)
                randomresponse = []
                for trx in range(sh[0]):  # get through the trials
                    randomresponse.append( neo.AnalogSignal( X_r[trx,:,:],\
                                                         name=acrossresponses[0][0].name+'randomorthonormal%d'%r,\
                                                         t_start=acrossresponses[0][0].t_start,\
                                                         sampling_period=acrossresponses[0][0].sampling_period, units=acrossresponses[0][0].units   ) )

                acrossrandomresponses.append(randomresponse)
            acrossrandomdecoder = nedi.get_responsedecoding(acrossrandomresponses,cv=cv) 
            acrossrandomdecoders.append(acrossrandomdecoder)    # collect for plot
            pickle.dump((acrossdecoder,acrossrandomdecoder),open(cacheprefix+'pca/responsedecodes+O,pcadim%d,%s-%s-%s-%s.pck'%(n_neurons,pcamean,dn,continuous_method,comparison),'wb'))
            
        else:
            for r in pcadimlist: # go for all pca dims and the last one as the actual number of neurons
                acrossdecoder,acrossrandomdecoder = pickle.load(open(cacheprefix+'pca/responsedecodes+O,pcadim%d,%s-%s-%s-%s.pck'%(r,pcamean,dn,continuous_method,comparison),'rb'))
                acrossdecoders.append(acrossdecoder)       # collect for plot
                acrossrandomdecoders.append(acrossrandomdecoder)    # collect for plot
                
    
    
    
            
        p_sa = []
        p_ea = []
        p_sa_r = []
        p_ea_r = []
        for rx in range(  len(pcadimlist)   ):
            p_sa.append(acrossdecoders[rx][1][T['start_idx']:T['stimstart_idx'],:].mean(axis=0))
            p_ea.append(acrossdecoders[rx][1][T['stimstart_idx']:T['stimend_idx'],:].mean(axis=0))
            p_sa_r.append(acrossrandomdecoders[rx][1][T['start_idx']:T['stimstart_idx'],:].mean(axis=0))
            p_ea_r.append(acrossrandomdecoders[rx][1][T['stimstart_idx']:T['stimend_idx'],:].mean(axis=0))



        if returnall:
            data.append( [ pcadimlist,p_sa,p_sa_r,p_ea,p_ea_r ] )# save for each variable






        if doplot:
            for rx in range(  len(pcadimlist)   ):
                
                axs = ax[cx,rx+2]
                colors_random_projection=['slategray','firebrick','']

                figs.plottoaxis_decoderrocauc(acrossrandomdecoders[rx],axs,colorlist=colors_random_projection)
                figs.plottoaxis_decoderrocauc(acrossdecoders[rx],axs)
                figs.plottoaxis_chancelevel(axs,0.5)
                figs.plottoaxis_stimulusoverlay(axs,T)
                if cx==0 and rx<len(pcadimlist)-1: axs.set_title('# pca dims: %d'%pcadimlist[rx])
                elif cx==0 and rx==len(pcadimlist)-1: axs.set_title('using all %d neurons'%pcadimlist[rx])
                if cx==len(taskaspects)-1: axs.set_xlabel('[ms]')


            axs = ax[cx,1] # use only the test trajectories:

            axs.plot( pcadimlist, p_sa, 'o-', color='seagreen', label='<SA>'  )
            axs.plot( pcadimlist, p_sa_r, 'o-', color='seagreen', alpha=0.3  )
            axs.plot( pcadimlist, p_ea, 'o-', color='purple', label='<EA>' )
            axs.plot( pcadimlist, p_ea_r, 'o-', color='purple',alpha=0.3 )
            axs.set_ylabel(taskaspects[cx])
            axs.set_xlim(0,10*int(np.ceil(n_neurons/10))+1)
            axs.set_ylim(0.45,1.01)
            figs.plottoaxis_chancelevel(axs,0.5)
            axs.legend()
            if cx==0: axs.set_title('pre+on mean acc.')
            if cx==len(taskaspects)-1: axs.set_xlabel('# pca dims')




            color = ['navy','darkgreen','mediumvioletred','orange','navy','darkblue','dodgerblue','skyblue'][cx]
            axs = ax[0,0] # SA
            axs.plot( pcadimlist, p_sa, 'o-', color=color, label=taskaspects[cx])
            axs.plot( pcadimlist, p_sa_r, 'o-', color=color, alpha=0.3)
            axs.set_ylabel('acc.')
            axs.set_xlim(0,10*int(np.ceil(n_neurons/10))+1)
            axs.set_ylim(0.45,1.01)
            figs.plottoaxis_chancelevel(axs,0.5)
            axs.legend();
            axs.set_title('pre mean acc. (SA)')

            axs = ax[1,0] # EA
            axs.plot( pcadimlist, p_ea, 'o-', color=color, label=taskaspects[cx])
            axs.plot( pcadimlist, p_ea_r, 'o-', color=color, alpha=0.3)
            axs.set_ylabel('acc.')
            axs.set_xlim(0,10*int(np.ceil(n_neurons/10))+1)
            axs.set_ylim(0.45,1.01)
            figs.plottoaxis_chancelevel(axs,0.5)
            axs.set_title('on mean acc. (EA)')
            axs.set_xlabel('# pca dims')

            









    
    if doplot:
        fig.suptitle(dn+' accuracies of decoders using only the best principal components\npca performed on the mean FR in %s activity; +random orthonormal baseline'%pcameanlabel)
        save = 1
        if save:
            fig.savefig(resultpath+'pcareduceddecoders+O,%s-%s-%dms%dms_%s'%(pcamean,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)


    if returnall:
        return data
















def subspacedynamics_PCA(dn,block):
    # perform PCA on neural activities in the multimodal blocks concatenated as a whole
    # draw neural activity projections on the PCA axes throughout the trials


    recalculate = 0 or globalrecalculate
    dump = 1
    doplot = 0 or globaldoplot

    taskaspects = ['visual','audio','context','choice']
    taskcolors = [ ['navy','darkgreen','purple','saddlebrown'],\
                   ['steelblue','c','mediumvioletred','orange'] ]

    blv,bla = preprocess.getorderattended(dn)
    pcagroups = [      [[2,4],[],[]]   ]         # use the multimodal blocks as activities
    comparisongroups  = [ \
                          [ [ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                          [ [ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ] ],\
                          [ [  blv,  [],[] ],   [ bla,   [], [] ] ],\
                          [[],[]]   \
                        ]

    classnames = [['45°','135°'],['5kHz','10kHz'],['attend visual','attend audio'],['lick','withhold lick']]
    
    
    n_pc = 4
    
    print('PCA to total concatenated activity ',dn)



    # projection matrices, V and V_r; use the entire trial from pre -1500 activity to post 3000+1500 ms
    responses = preprocess.collect_stimulusspecificresponses(block,pcagroups,'and')
    n_trials = len(responses[0])
    n_trajectory,n_neurons = responses[0][0].shape
    X = np.array(responses[0])[:,T['start_idx']:T['end_idx'],:].reshape(-1,n_neurons)
    _,E,V = nedi.pca(X,n_components=n_neurons)  # do all the possible PCA-s, but later use only the first n_pc ones
    print('PCA eigenvalues',E)
    print('neurons %d, PCs %d,'%(n_neurons,2),'    PCA transform shape',V.shape)
    print('explained variance % and cum %:')
    totalvariance = np.trace(np.cov(X.T))
    print(dn,  E/totalvariance)
    print('     ',np.cumsum(E/totalvariance))



    projected_dynamics = []    # this will be (taskaspects,classes,[trials,2])
    for cx,comparison in enumerate(taskaspects):
        # now project each activity onto the common PCs vectorspace

        # if cx in [4,5,6]: continue
        print(comparison)
        
        # collect neural responses
        if not comparison=='choice': # visual, audio, context:
            acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],'and')
        else:  # choice:
            acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
#        print(np.array(acrossresponses[0]).shape)



        projected_dynamic = []
        for aix,classresponse in enumerate(acrossresponses):
            # project the entire dynamics trajectory onto the pcas axes for each trial
            X = np.dot(np.array(classresponse),V.T)  #    (trials,trajectory,pc) = (trials,trajectory,neurons) * (,,pc,neurons).T
            projected_dynamic.append(X)
        
        projected_dynamics.append(projected_dynamic)



    # calculate distance metric in the PCs of the rate subpaces
        
    distances = []     # (numpc),(task),(mean,se),(trajectory)
    for px in range(n_pc):
        # distance in multiple dimensions with errors:
        distance = []
        for cx,comparison in enumerate(taskaspects):
                
                
                # get the trial averages and their 2 standard errors
                # these each hold (trajectory,pcs) analogsignal like activities
                trial = [ np.array(projected_dynamics[cx][aix]).mean(axis=0) for aix in [0,1] ]
                trial_e = [ np.array(projected_dynamics[cx][aix]).std(axis=0)*2/np.sqrt(len(projected_dynamics[cx][aix])) for aix in [0,1] ]

                # calculate vector distance for each point, in the pc subspace
                # take the difference
                stderror_components = trial_e[0]+trial_e[1]
                distance_trajectory_components = np.abs(trial[0]-trial[1])#-stderror_components
                # apply distance postivity condition, due to standard error
                distance_trajectory_components = np.max( np.stack( [distance_trajectory_components , np.zeros(trial[0].shape)],axis=2), axis=2  )

                d_m = np.sqrt(np.sum(distance_trajectory_components[:,:px+1]**2,axis=1))
                d_e = np.sqrt(np.sum(stderror_components[:,:px+1]**2,axis=1))
                
                distance.append([d_m,d_e])
        
        distances.append(distance)


    if dump:    # save for display in publications
        pickle.dump((projected_dynamics,distances,responses[0][0].times),open(cacheprefix+'subspaces/subspacedynamics,3D-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'wb'))
        
        
        
        
        



    if doplot:

        pcmaxs = 0*np.ones(n_pc)
        pcmins = 0*np.ones(n_pc)

        alpha = 0.8
        # time parametrization
        skip_idx = 20
        t_all = responses[0][0].times[skip_idx:-skip_idx]
        t_0 = responses[0][0].times[skip_idx:T['stimstart_idx']+1]
        t_1 = responses[0][0].times[T['stimstart_idx']:T['stimend_idx']+1]
        t_oo = responses[0][0].times[T['stimend_idx']:-skip_idx]

        
        fig,ax = plt.subplots(n_pc+3,len(taskaspects),figsize=(len(taskaspects)*8,(n_pc+3)*8))
        for cx,comparison in enumerate(taskaspects):
            
            
            for aix,classresponse in enumerate(projected_dynamics[cx]): # gp through classes
                
                # single trials:
                # classresponse_selected = classresponse[trialset]
                # for trialx,trial in enumerate(classresponse_selected):    # go through trials
                # axs = ax[trialx,cx]
                
                # trial average:
                
                trial = np.array(classresponse).mean(axis=0)
                trial_e = np.array(classresponse).std(axis=0)*2/np.sqrt(len(classresponse))
                
                trial_pre = trial[skip_idx:T['stimstart_idx']+1,:]
                trial_on = trial[T['stimstart_idx']:T['stimend_idx']+1,:]
                trial_post = trial[T['stimend_idx']:-skip_idx,:]
                
                
                
                
                # 3d 3 PC
                if aix==0:
                    ax[1,cx].remove()
                    ax[1,cx] = fig.add_subplot(6,len(taskaspects),len(taskaspects)+cx+1,projection='3d')
                axs = ax[1,cx]
                
                axs.plot( trial_pre[:,0], trial_pre[:,1], trial_pre[:,2], lw=1,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( trial_on[:,0], trial_on[:,1], trial_on[:,2], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( trial_post[:,0], trial_post[:,1], trial_post[:,2], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( [trial_on[0,0]], [trial_on[0,1]], [trial_on[0,2]], 'o',color=taskcolors[aix][cx],alpha=alpha )
                
                axs.set_xlabel('PC 1')
                axs.set_ylabel('PC 2')
                axs.set_zlabel('PC 3')
                
                
                
                # 2D 2 PC
                axs = ax[2,cx]
                
                axs.plot( trial_pre[:,0], trial_pre[:,1], lw=1,color=taskcolors[aix][cx],alpha=alpha)
                axs.plot( trial_on[:,0], trial_on[:,1], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( trial_post[:,0], trial_post[:,1], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( [trial_on[0,0]], [trial_on[0,1]], 'o',color=taskcolors[aix][cx],alpha=alpha )
                
                
                axs.set_xlabel('PC 1')
                if cx==0: axs.set_ylabel('PC 2')
                
                
                # 1D 1 PC
                for px in range(n_pc):

                    # gather for identical axes limits
                    pcmaxs[px] = np.max(( pcmaxs[px], np.max(trial[skip_idx:-skip_idx,px])   ))
                    pcmins[px] = np.min(( pcmins[px], np.min(trial[skip_idx:-skip_idx,px])   ))


                    axs = ax[px+3,cx]
                    
                    axs.plot( t_0, trial_pre[:,px], lw=1,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( t_1, trial_on[:,px], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( t_oo, trial_post[:,px], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( [t_1[0]], [trial_on[0,px]], 'o',color=taskcolors[aix][cx],alpha=alpha,label=classnames[cx][aix] )
                    
                    axs.fill_between(t_all, trial[skip_idx:-skip_idx,px]-trial_e[skip_idx:-skip_idx,px], \
                                            trial[skip_idx:-skip_idx,px]+trial_e[skip_idx:-skip_idx,px], \
                                     color=taskcolors[aix][cx],alpha=alpha/4.5)#,label=label)

                    figs.setxt(axs)
                    axs.set_xlim(skip_idx*T['dt']+T['starttime'], -skip_idx*T['dt']+T['endtime'])
                    

                    if aix==1:
                        if px==0: axs.legend(frameon=False)
                    if cx==0: axs.set_ylabel('PC %d'%(px+1))
                    
                    
                    
                    
                    # now plot distance metrics:
                    if aix==0:
                        axs = ax[0,cx]

                        axs.plot(t_all,distances[px][cx][0][skip_idx:-skip_idx],\
                                 color=taskcolors[0][cx],lw=2,alpha=0.8-px*1/n_pc/2,label='#PCs$\leq$%d'%(px+1))
                        
                        s_e_m = np.max( np.stack( (distances[px][cx][0][skip_idx:-skip_idx]-distances[px][cx][1][skip_idx:-skip_idx],\
                                                   np.zeros(len(t_all))),axis=1     ), axis=1 )
                        axs.fill_between(t_all,s_e_m,\
                                               distances[px][cx][0][skip_idx:-skip_idx],\
                                         color=taskcolors[0][cx],lw=2,alpha=(0.8-px*1/n_pc/2)/8)
                        
                        
                        if px==n_pc-1:
                            axs.legend(frameon=False)
                            figs.setxt(axs)
                            figs.plottoaxis_stimulusoverlay(axs,T)
                            figs.plottoaxis_chancelevel(axs)
                            axs.set_title(comparison)
                            if cx==0: axs.set_ylabel('distance in PC space')

                    

                axs.set_xlabel('trial time from stimulus onset [ms]')


        pcmins *=1.1
        pcmaxs *=1.1
        for cx in range(len(taskaspects)):
            ax[1,cx].set_xlim(pcmins[0],pcmaxs[0])
            ax[1,cx].set_ylim(pcmins[1],pcmaxs[1])
            ax[1,cx].set_zlim(pcmins[2],pcmaxs[2])

            ax[2,cx].set_xlim(pcmins[0],pcmaxs[0])
            ax[2,cx].set_ylim(pcmins[1],pcmaxs[1])

            for px in range(n_pc):
                axs = ax[3+px,cx]
                axs.set_ylim(pcmins[px],pcmaxs[px])
                figs.plottoaxis_chancelevel(axs)
                figs.plottoaxis_stimulusoverlay(axs,T)
                
                
        fig.suptitle(dn+' trial averaged projected neuron dynamics onto %d PCs from trial+trajectory-concatenated activity'%n_pc)
        save = 0
        if save:
            fig.savefig(resultpath+'trials,over,totalconcatpca,avg+dist-%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)


    
    return








def subspacedynamics_PCA_context(dn,block):
    # perform PCA on neural activities in the multimodal blocks concatenated as a whole
    # draw neural activity projections on the PCA axes throughout the trials


    recalculate = 1 or globalrecalculate
    dump = 0
    doplot = 1 or globaldoplot

    pcawhere = ['full','ea']
    taskaspects = ['attend visual','ignore visual','45°','135°']
    taskcolors = [ ['darkgreen','green','darkgreen','darkred'],\
                   ['darkred','orangered','green','orangered'] ]

    blv,bla = preprocess.getorderattended(dn)
    pcagroups = [      [[2,4],[],[]]   ]         # use the multimodal blocks as activities
    
    comparisongroups  = [ \
                          [ [  [blv[1]],  [45],[] ],   [  [blv[1]],  [135],[] ] ],\
                          [ [  [bla[1]],  [45],[] ],   [  [bla[1]],  [135],[] ] ], \
                          [ [  [blv[1]],  [45],[] ],   [  [bla[1]],  [45],[] ] ], \
                          [ [  [blv[1]],  [135],[] ],  [  [bla[1]],  [135],[] ] ], \
                        ]

    classnames = [['45°','135°'],['45°','135°'],['attend','ignore'],['attend','ignore']]
    
    
    n_pc = 4
    
    print('PCA to total concatenated activity ',dn)



    # projection matrices, V and V_r; use the entire trial from pre -1500 activity to post 3000+1500 ms
    responses = preprocess.collect_stimulusspecificresponses(block,pcagroups,'and')
    n_trials = len(responses[0])
    n_trajectory,n_neurons = responses[0][0].shape

    for pcabx,pcabase in enumerate(pcawhere):    # go over full and only EA PCA basis

        if pcabase=='full':
            concatenatedactivity = np.array(responses[0])[:,T['start_idx']:T['end_idx'],:]
        elif pcabase=='ea':
            concatenatedactivity = np.array(responses[0])[:,T['stimstart_idx']:T['stimend_idx'],:]
        _,E,V = nedi.pca(concatenatedactivity.reshape(-1,n_neurons),n_components=n_pc)
        print('PCA eigenvalues',E)
        print('neurons %d, PCs %d,'%(n_neurons,2),'    PCA transform shape',V.shape)



        
        projected_dynamics = []    # this will be (tasks,classes,[trials,2])
        for cx,comparison in enumerate(taskaspects):
            # now project each activity onto the common PCs vectorspace
    
            # if cx in [4,5,6]: continue
            print(comparison)
            
            # collect neural responses
            if not comparison=='choice': # visual, audio, context:
                acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],'and')
            else:  # choice:
                acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
    #        print(np.array(acrossresponses[0]).shape)
    
    
    
            projected_dynamic = []
            for aix,classresponse in enumerate(acrossresponses):
                # project the entire dynamics trajectory onto the PCs axes for each trial
                X = np.dot(np.array(classresponse),V.T)  #    (trials,trajectory,pc) = (trials,trajectory,neurons) * (,,pc,neurons).T
                projected_dynamic.append(X)
            
            projected_dynamics.append(projected_dynamic)
    
    
    
        # calculate distance metric in the PCs of the rate subpaces
            
        distances = []     # (numpc),(task),(mean,se),(trajectory)
        for px in range(n_pc):
            # distance in multiple dimensions with errors:
            distance = []
            for cx,comparison in enumerate(taskaspects):
                    
                    
                    # get the trial averages and their 2 standard errors
                    # these each hold (trajectory,pcs) analogsignal like activities
                    trial = [ np.array(projected_dynamics[cx][aix]).mean(axis=0) for aix in [0,1] ]
                    trial_e = [ np.array(projected_dynamics[cx][aix]).std(axis=0)*2/np.sqrt(len(projected_dynamics[cx][aix])) for aix in [0,1] ]
    
                    # calculate vector distance for each point, in the pc subspace
                    # take the difference
                    stderror_components = trial_e[0]+trial_e[1]
                    distance_trajectory_components = np.abs(trial[0]-trial[1])#-stderror_components
                    # apply distance postivity condition, due to standard error
                    distance_trajectory_components = np.max( np.stack( [distance_trajectory_components , np.zeros(trial[0].shape)],axis=2), axis=2  )
    
                    d_m = np.sqrt(np.sum(distance_trajectory_components[:,:px+1]**2,axis=1))
                    d_e = np.sqrt(np.sum(stderror_components[:,:px+1]**2,axis=1))
                    
                    distance.append([d_m,d_e])
            
            distances.append(distance)
    
    
        if dump:    # save for display in publications
            pickle.dump((projected_dynamics,distances,responses[0][0].times),open(cacheprefix+'subspaces/subspacedynamics,3D,attendignore-%s_%s-%dms_%s.pck'%(dn,continuous_method,T['dt'].magnitude,pcabase),'wb'))
            
            
            
            
            
    
    
    
        if doplot:
    
            pcmaxs = 0*np.ones(n_pc)
            pcmins = 0*np.ones(n_pc)
    
            alpha = 0.8
            # time parametrization
            skip_idx = 20
            t_all = responses[0][0].times[skip_idx:-skip_idx]
            t_0 = responses[0][0].times[skip_idx:T['stimstart_idx']+1]
            t_1 = responses[0][0].times[T['stimstart_idx']:T['stimend_idx']+1]
            t_oo = responses[0][0].times[T['stimend_idx']:-skip_idx]
    
            if pcabx==0:
                fig,ax = plt.subplots(n_pc+3,len(taskaspects)*len(pcawhere)+1,\
                                      figsize=((len(taskaspects)*len(pcawhere)+1)*8,(n_pc+3)*8))
            for cx,comparison in enumerate(taskaspects):
                
                fx = cx+(pcabx*5) # the columns for taskaspects and pairing full and ea only pca
                
                for aix,classresponse in enumerate(projected_dynamics[cx]): # gp through classes
                    
                    # single trials:
                    # classresponse_selected = classresponse[trialset]
                    # for trialx,trial in enumerate(classresponse_selected):    # go through trials
                    # axs = ax[trialx,cx]
                    
                    # trial average:
                    
                    trial = np.array(classresponse).mean(axis=0)
                    trial_e = np.array(classresponse).std(axis=0)*2/np.sqrt(len(classresponse))
                    
                    trial_pre = trial[skip_idx:T['stimstart_idx']+1,:]
                    trial_on = trial[T['stimstart_idx']:T['stimend_idx']+1,:]
                    trial_post = trial[T['stimend_idx']:-skip_idx,:]
                    
                    
                    
                    
                    # 3d 3 PC
                    if aix==0:
                        ax[1,fx].remove()
                        ax[1,fx] = fig.add_subplot(n_pc+3,len(taskaspects)*len(pcawhere)+1,len(taskaspects)*len(pcawhere)+1+fx+1,projection='3d')
                    axs = ax[1,fx]
                    
                    axs.plot( trial_pre[:,0], trial_pre[:,1], trial_pre[:,2], lw=1,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( trial_on[:,0], trial_on[:,1], trial_on[:,2], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( trial_post[:,0], trial_post[:,1], trial_post[:,2], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( [trial_on[0,0]], [trial_on[0,1]], [trial_on[0,2]], 'o',color=taskcolors[aix][cx],alpha=alpha )
                    
                    axs.set_xlabel('PC 1')
                    axs.set_ylabel('PC 2')
                    axs.set_zlabel('PC 3')
                    
                    
                    
                    # 2D 2 PC
                    axs = ax[2,fx]
                    
                    axs.plot( trial_pre[:,0], trial_pre[:,1], lw=1,color=taskcolors[aix][cx],alpha=alpha)
                    axs.plot( trial_on[:,0], trial_on[:,1], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( trial_post[:,0], trial_post[:,1], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( [trial_on[0,0]], [trial_on[0,1]], 'o',color=taskcolors[aix][cx],alpha=alpha )
                    
                    
                    axs.set_xlabel('PC 1')
                    if cx==0: axs.set_ylabel('PC 2')
                    
                    
                    # 1D 1 PC
                    for px in range(n_pc):
    
                        # gather for identical axes limits
                        pcmaxs[px] = np.max(( pcmaxs[px], np.max(trial[skip_idx:-skip_idx,px])   ))
                        pcmins[px] = np.min(( pcmins[px], np.min(trial[skip_idx:-skip_idx,px])   ))
    
    
                        axs = ax[px+3,fx]
                        
                        axs.plot( t_0, trial_pre[:,px], lw=1,color=taskcolors[aix][cx],alpha=alpha )
                        axs.plot( t_1, trial_on[:,px], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                        axs.plot( t_oo, trial_post[:,px], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                        axs.plot( [t_1[0]], [trial_on[0,px]], 'o',color=taskcolors[aix][cx],alpha=alpha,label=classnames[cx][aix] )
                        
                        axs.fill_between(t_all, trial[skip_idx:-skip_idx,px]-trial_e[skip_idx:-skip_idx,px], \
                                                trial[skip_idx:-skip_idx,px]+trial_e[skip_idx:-skip_idx,px], \
                                         color=taskcolors[aix][cx],alpha=alpha/4.5)#,label=label)
    
                        figs.setxt(axs)
                        axs.set_xlim(skip_idx*T['dt']+T['starttime'], -skip_idx*T['dt']+T['endtime'])
                        
    
                        if aix==1:
                            if px==0: axs.legend(frameon=False)
                        if cx==0: axs.set_ylabel('PC %d'%(px+1))
                        
                        
                        
                        
                        # now plot distance metrics:
                        if aix==0:
                            axs = ax[0,fx]
    
                            axs.plot(t_all,distances[px][cx][0][skip_idx:-skip_idx],\
                                     color=taskcolors[0][cx],lw=2,alpha=0.8-px*1/n_pc/2,label='#PCs$\leq$%d'%(px+1))
                            
                            s_e_m = np.max( np.stack( (distances[px][cx][0][skip_idx:-skip_idx]-distances[px][cx][1][skip_idx:-skip_idx],\
                                                       np.zeros(len(t_all))),axis=1     ), axis=1 )
                            axs.fill_between(t_all,s_e_m,\
                                                   distances[px][cx][0][skip_idx:-skip_idx],\
                                             color=taskcolors[0][cx],lw=2,alpha=(0.8-px*1/n_pc/2)/8)
                            
                            
                            if px==n_pc-1:
                                axs.legend(frameon=False)
                                figs.setxt(axs)
                                figs.plottoaxis_stimulusoverlay(axs,T)
                                figs.plottoaxis_chancelevel(axs)
                                axs.set_title('PCA found on %s\n'%pcabase+comparison)
                                if cx==0: axs.set_ylabel('distance in PC space')
    
    
                    axs.set_xlabel('trial time from stimulus onset [ms]')
    
    
            pcmins *=1.1
            pcmaxs *=1.1
            for cx in range(len(taskaspects)):
                fx = cx+(pcabx*5)
                ax[1,fx].set_xlim(pcmins[0],pcmaxs[0])
                ax[1,fx].set_ylim(pcmins[1],pcmaxs[1])
                ax[1,fx].set_zlim(pcmins[2],pcmaxs[2])
    
                ax[2,fx].set_xlim(pcmins[0],pcmaxs[0])
                ax[2,fx].set_ylim(pcmins[1],pcmaxs[1])
    
                for px in range(n_pc):
                    axs = ax[3+px,fx]
                    axs.set_ylim(pcmins[px],pcmaxs[px])
                    figs.plottoaxis_chancelevel(axs)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    
        for k in range(n_pc+3):
            ax[k,4].axis('off')
                
        fig.suptitle(dn+' trial averaged projected neuron dynamics onto %d PCs from trial+trajectory-concatenated activity\nPCA found on SA+EA (left) and EA only (right)'%n_pc)
        save = 1
        if save:
            fig.savefig(resultpath+'trials,over,total+eaonlyconcatpca,avg+dist-%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)


    
    return







def subspacedynamics_PCA_signalnoisecorrelation(dn,block):
    # total variance = signal variance + noise variance
    # Cov yi,yj = Ex[ Cov yi yj | x ] + Covx( E[yi|x]E[yj|x] )
    # that is
    # 1) we take a task variable
    # 2) condition on one value of it
    # 3) calculate the Covyiyj and EyiEyj
    # 4) we average and calc. covariance over the values of x:  Ex Covyiyj  and Covx EyiEyj


    doplot = 1 or globaldoplot
    dump = 1

    pcagroups = [      [[2,4],[],[]]   ]         # use the multimodal blocks as activities
    
    comparisongroups = [  ]
    
    taskaspects = ['visual','audio','context','choice']
    taskcolors = ['navy','darkgreen','mediumvioletred','darkorange']

    # calculate the PCA projector matrix, V, which will project to the first k principal component: X @ V[:,:k]
    

    n_pc_display = 4
    # responses will be n_trials, n_trajectory, n_neurons
    responses_pca = np.array(preprocess.collect_stimulusspecificresponses(block,pcagroups,'and')[0])
    n_trials_full,n_trajectory,n_neurons = responses_pca.shape
    X = responses_pca[:,T['start_idx']:T['end_idx'],:].reshape(-1,n_neurons)
    Y,E,V = nedi.pca(X,n_components=n_neurons)
    
    # load projected activities
    # projected_dynamics will be (taskaspects)(taskclass)(trials,trajectory,pc) = (trials,trajectory,neurons) * (,,pc,neurons).T
    projected_dynamics,_,times = pickle.load(open(cacheprefix+'subspaces/subspacedynamics,3D-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'rb'))
    n_tasks = len(projected_dynamics)       # 4, should equal len(taskaspects)
    n_classes = len(projected_dynamics[0])     #  2
    n_trajectory = projected_dynamics[0][0].shape[1]
    
    # projections_alltrials = (responses_pca @ V.T)[:,:,:n_pc]     # only the first few pc
    projections_alltrials = (responses_pca @ V.T)[:,:,:]           # project activities to all available pcs (# pc = # neurons)

    total_variance_concatpca = np.repeat(E[np.newaxis,:],n_trajectory,axis=0)

    total_covariance_pca = np.array([ np.cov(projections_alltrials[:,t,:].T) for t in range(n_trajectory) ])
    total_variance_pca = total_covariance_pca.diagonal(axis1=1,axis2=2)     # just count the variance, no need for covariance in the PCA space
    print(total_variance_concatpca.shape)
    

    # results collector:
    correlations = np.zeros(( n_tasks,n_trajectory,60,5 ))  # avg. noise variance. and signal-average variance, and totals

    for cx in range(n_tasks):
        print(taskaspects[cx])
        
        noise_covariance = []
        signal_average = []
        for aix in range(n_classes):
            # print(cx+1,aix)

            n_trials,n_trajectory,n_pc = projected_dynamics[cx][aix].shape
            X = projected_dynamics[cx][aix]
            # print('    original data shape:',X.shape)
            # choose either:
            if 0:
                X = X.reshape(-1,n_neurons)       # reshape to reflect total trial and time variance
                # X is    observations  by  principal component projections
                n_observations = n_trials*n_trajectory
            else:
                X = np.moveaxis(X,1,0)  # move the trials to average over, and keep the trajectory in the first untouched dimension
                n_observations = n_trials
            
            noise_covariance.append(  np.array([ np.cov(X[t,:,:].T) for t in range(n_trajectory) ])   )

            signal_average.append(   X.mean(axis=-2)  )   # assume the observation dimension and average over it (e.g. trials)


        noise_covariance = np.array(noise_covariance)  #   this holds  n_class,n_trajectory,n_trials,n_pcs
        signal_average = np.array(signal_average) #   this holds  n_class,n_trajectory,n_pcs
        print('    per class shape',noise_covariance.shape, signal_average.shape)
        
        average_noise_covariance = np.mean(noise_covariance,axis=0)       # expected value of the noise_cov over condition values
        
        
        # signal_average_variance = signal_average.transpose(1,2,0).var(axis=2)
        # M = signal_average.transpose(2,1,0)
        signal_average_covariance = np.einsum('...i,...k->...ik',signal_average,signal_average)      # create the outer product for the combinations of average activities
        # signal_average_variance = np.trace(signal_covariance,axis1=2,axis2=3) / n_classes 
        signal_average_variance = signal_average.var(axis=0)    # variance over the values of the condition
        

        print('    inside  shape',noise_covariance.shape, signal_average.shape)
        print('    outside shape',average_noise_covariance.shape, signal_average_variance.shape)

        # print('timecourse average noise:')
        # print(average_noise_covariance.mean(axis=0))
        # print('timecourse average signal:')
        # print(signal_average_variance.mean(axis=0))

        # total_covariance = average_noise_covariance + average_signal_covariance

        # Now we would cumulate the first couple of PCs, and register their cumulative variance, i.e. we are
        # only concerned with the expected values and the diagonals of the covariance matrix.
        # Since the contribution of later PCs are small, we stick to not cumulative, so the plots are not identical
        
        correlations[cx,:,:n_pc,0] = average_noise_covariance.diagonal(axis1=1,axis2=2)
        correlations[cx,:,:n_pc,1] = signal_average_variance
        # correlations[cx,:,:,0] = average_noise_covariance.diagonal(axis1=1,axis2=2).cumsum(axis=1)
        # correlations[cx,:,:,1] = signal_average_variance.diagonal(axis1=1,axis2=2).cumsum(axis=1)


        correlations[cx,:,:n_pc,2] = np.repeat(total_variance_pca[:,:n_pc].sum(axis=1)[:,np.newaxis],n_pc,axis=1)
        correlations[cx,:,:n_pc,3] = np.repeat(total_variance_pca[:,:].sum(axis=1)[:,np.newaxis],n_pc,axis=1)
        correlations[cx,:,:n_pc,4] = np.repeat(total_variance_concatpca[:,:].sum(axis=1)[:,np.newaxis],n_pc,axis=1)
        

        # correlations[cx,:,:,3] = total_variance_neuron

    if dump: # save the trajectory=time averaged values
        average_covariance = correlations.mean(axis=1)
        pickle.dump((average_covariance,n_neurons),open(cacheprefix+'subspaces/noise+signal-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'wb'))


    if doplot:
        
        fig,ax = plt.subplots(n_pc_display,n_tasks,figsize=(n_tasks*7,n_pc_display*7))

        for cx,comparison in enumerate(taskaspects):
            for j in range(n_pc_display):
                axs = ax[j,cx]
                
                axs.plot(times,correlations[cx,:,j,2],'grey',lw=2,label='total variance 4 pc')
                axs.plot(times,correlations[cx,:,j,3],'black',lw=2,label='total variance all %d pc'%n_neurons)
                axs.plot(times,correlations[cx,:,j,4],'--',color='grey',lw=1,label='total concat. variance all %d pc'%n_neurons)
                
                axs.plot(times,correlations[cx,:,j,0],'--',color=taskcolors[cx],lw=2,alpha=0.5,label='average noise variance')
                axs.plot(times,correlations[cx,:,j,1],color=taskcolors[cx],lw=3,alpha=0.8,label='signal average variance')
                
                # axs.plot(times,correlations[j,cx,:,1],lw=2,color=taskcolors[cx],alpha=0.8,\
                #          label='average signal var.')
                
                axs.set_xlim(-1000,4000)
                
                figs.plottoaxis_stimulusoverlay(axs,T)
                figs.setxt(axs)
                figs.plottoaxis_chancelevel(axs)

                if cx==0 and j==0: axs.legend(frameon=False,ncol=1)

                if cx==0: axs.set_ylabel('PC %d\nvariance [SU]'%(j+1))
                if j==0: axs.set_title(comparison)
                    
                    
                    
            # fig.suptitle(dn+' '+comparison+'\nanc = average noise covariance (over pc),\nsc kl = average signal variance (over classes)')
            fig.suptitle(dn+'        Total Variance (black) =\nAverage Noise Variance (color thin dash) + Variance of Signal Averages (color thick solid')

            save = 0
            if save:
                # fig.savefig(resultpath+'covariance-noise,signal-%s_%s-%dms%dms_%s'%(comparison,continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)
                fig.savefig(resultpath+'covariance-noise,signal-perpcnc_%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)




    return










def subspacedynamics_decisionvector_basis(dn,block):
    # perform DBNV on neural activities in the multimodal blocks concatenated as a whole
    # draw neural activity projections on the DBNV axes throughout the trials


    recalculate = 0 or globalrecalculate
    dump = 1
    doplot = 0 or globaldoplot

    taskaspects = ['visual','audio','context','choice']
    taskcolors = [ ['navy','darkgreen','purple','saddlebrown'],\
                   ['steelblue','c','mediumvioletred','orange'] ]




    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [ \
                          [ [ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                          [ [ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ] ],\
                          [ [  blv,  [],[] ],   [ bla,   [], [] ] ],\
                          [[],[]],   \
                          [[],[]],   \
                          [[],[]]   \
                        ]

    classnames = [['45°','135°'],['5kHz','10kHz'],['attend visual','attend audio'],['lick','withhold lick']]
    
    
    print(dn, 'Activity trajectory projections onto DBNVs')






    # projection basis vectors:
        
    # choose the order of the orthogonalization:
    # basisaspects = ['visual','context','choice']; basisaspectcolors = ['dodgerblue','fuchsia','gold']             # for context to be related to visual
    # basisaspects = ['choice','context','visual']; basisaspectcolors = ['gold','fuchsia','dodgerblue']             # for context to be related to choice
    # run this and the next separately for choice and context!
    # basisaspects = ['context','choiceattendvisual','visual']; basisaspectcolors = ['fuchsia','gold','dodgerblue']             # for context to be related to choice
    basisaspects = ['context','choiceattendaudio','visual']; basisaspectcolors = ['fuchsia','gold','dodgerblue']             # for context to be related to choice


    # do the entire projection to the above basis:
    filenamesuffix = ''.join([ s[:2] for s in basisaspects ]).replace('o','x')
    depth_idx = int((1500*pq.ms / T['dt']).magnitude)
    n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
    c_db_matrix = np.zeros((n_neurons,3))       # prepare the dbnv vector holder
    n_ba = len(basisaspects)
    for bx,aspect in enumerate(basisaspects):
        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,aspect,'all'),'rb'))
        wx = int((len(acrossdecoder)-7)/n_neurons)
        coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
        # [neurons,trajectory,stats]
        # mean of on stimulus
        c_db_matrix[:,bx] = coeff[:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,0].mean(axis=1)

    # c_db_matrix projects to the 3 dbnvs each row:   c_db_matrix @ neural_input
    Q,R = np.linalg.qr(c_db_matrix)     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
    V = Q                          # activity transformer matrix
    # print('c_db',c_db_matrix.shape,'Q',Q.shape,'R',R.shape)
    B = Q.T @ c_db_matrix / np.linalg.norm(c_db_matrix,axis=0,keepdims=True)        # transformed dbnvs (in columns) in orthonormal coordinates (in rows)
    # print( 'V',V.shape, 'B', B.shape )

    # to test, choose no orthogonalizatoin, just project onto the skewed basis:
    # V = c_db_matrix
    # B = np.eye(3)


    projected_dynamics = []    # this will be (taskaspects,classes,[trials,dbnv])
    for cx,comparison in enumerate(taskaspects):
        # now project each activity onto the common DBNVs vectorspace

        # if cx in [4,5,6]: continue
        print(comparison)
        
        # collect neural responses
        if not comparison[:6]=='choice': # visual, audio, context:
            acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],'and')
        elif comparison=='choice':  # choice different versions
            acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
        elif comparison=='choiceattendvisual':
            acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn,onlyinblock=[blv[1]])
        elif comparison=='choiceattendaudio':
            acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn,onlyinblock=[bla[1]])
                
#        print(np.array(acrossresponses[0]).shape)



        projected_dynamic = []
        for aix,classresponse in enumerate(acrossresponses):
            # project the entire dynamics trajectory onto the dbnvs axes for each trial
            X = np.dot(np.array(classresponse),V)  #    (trials,trajectory,dbnv) = (trials,trajectory,neurons) * (,,dbnv,neurons).T
            projected_dynamic.append(X)
        
        projected_dynamics.append(projected_dynamic)






    # now project all multimodal trials with grouping, so that they can be selected by index
    if filenamesuffix=='vicxch':   # 2 times double = 4 activitylist for the contexts.
        projectionlist = [[[blv[1]],[45],[]],[[blv[1]],[135],[]],[[bla[1]],[45],[]],[[bla[1]],[135],[]]]
        responses = preprocess.collect_stimulusspecificresponses(block,projectionlist)
    elif filenamesuffix=='cxchvi':  # choose the list so that it will be av hmcf or aa hmcf    8 length list of response lists; this will be run two times for projecting to each context's choice dbnv
        # projectionlist = [[[blv[1]],[45],[]],[[blv[1]],[135],[]],[[bla[1]],[45],[]],[[bla[1]],[135],[]]]
        responses = preprocess.collect_stimulusspecificresponses_choice_detailed(block,dn,onlyinblock=[blv[1]])
        responses.extend(preprocess.collect_stimulusspecificresponses_choice_detailed(block,dn,onlyinblock=[bla[1]]))

    # projections_av = [ np.dot(np.array(trial.analogsignals[0],V) for trial in block.segments if trial.annotations['block'] in [blv[1]]  ]
    #    (trials,trajectory,dbnv) = (trials,trajectory,neurons) * (,,dbnv,neurons).T
    projections_alltrials = [ [np.dot(np.array(responses_trial),V) for responses_trial in response ] for response in responses ]
    # print(len(projections_alltrials),len(projections_alltrials[0]),projections_alltrials[0][0].shape)






    # calculate distance metric in the space from bases of the task variables
        
    distances = []     # (numpc),(task),(mean,se),(trajectory)
    for bx in range(n_ba):
        # distance in multiple dimensions with errors:
        distance = []
        for cx,comparison in enumerate(taskaspects):
                
                
                # get the trial averages and their 2 standard errors
                # these each hold (trajectory,pcs) analogsignal like activities
                trial = [ np.array(projected_dynamics[cx][aix]).mean(axis=0) for aix in [0,1] ]
                trial_e = [ np.array(projected_dynamics[cx][aix]).std(axis=0)*2/np.sqrt(len(projected_dynamics[cx][aix])) for aix in [0,1] ]

                # calculate vector distance for each point, in the pc subspace
                # take the difference
                stderror_components = trial_e[0]+trial_e[1]
                distance_trajectory_components = np.abs(trial[0]-trial[1])#-stderror_components
                # apply distance postivity condition, due to standard error
                distance_trajectory_components = np.max( np.stack( [distance_trajectory_components , np.zeros(trial[0].shape)],axis=2), axis=2  )

                d_m = np.sqrt(np.sum(distance_trajectory_components[:,:bx+1]**2,axis=1))
                d_e = np.sqrt(np.sum(stderror_components[:,:bx+1]**2,axis=1))
                
                distance.append([d_m,d_e])
        
        distances.append(distance)


    if dump:    # save for display in publications
        # pickle.dump((projected_dynamics,distances,responses[0][0].times),open(cacheprefix+'subspaces/subspacedynamics,3D-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'wb'))
        choiceattendtype = basisaspects[1][12:13]
        pickle.dump((projections_alltrials,B),open(cacheprefix+'subspaces/subspacedynamics,projected+dbnv,%s,%s-%s_%s-%dms.pck'%(filenamesuffix,'a'+choiceattendtype,dn,continuous_method,T['dt'].magnitude),'wb'))
        
        
        
        
        



    if doplot:

        maxs = 0*np.ones(n_ba)
        mins = 0*np.ones(n_ba)

        alpha = 0.8
        # time parametrization
        skip_idx = 20
        t_all = block.segments[0].analogsignals[0].times[skip_idx:-skip_idx]
        t_0 = block.segments[0].analogsignals[0].times[skip_idx:T['stimstart_idx']+1]
        t_1 = block.segments[0].analogsignals[0].times[T['stimstart_idx']:T['stimend_idx']+1]
        t_oo = block.segments[0].analogsignals[0].times[T['stimend_idx']:-skip_idx]

        
        fig,ax = plt.subplots(n_ba+3,len(taskaspects),figsize=(len(taskaspects)*8,(n_ba+3)*8))
        for cx,comparison in enumerate(taskaspects[:4]):
            
            
            for aix,classresponse in enumerate(projected_dynamics[cx]): # gp through classes
                
                # single trials:
                # classresponse_selected = classresponse[trialset]
                # for trialx,trial in enumerate(classresponse_selected):    # go through trials
                # axs = ax[trialx,cx]
                
                # trial average:
                
                trial = np.array(classresponse).mean(axis=0)
                trial_e = np.array(classresponse).std(axis=0)*2/np.sqrt(len(classresponse))
                
                trial_pre = trial[skip_idx:T['stimstart_idx']+1,:]
                trial_on = trial[T['stimstart_idx']:T['stimend_idx']+1,:]
                trial_post = trial[T['stimend_idx']:-skip_idx,:]
                
                

                
                # 3d 3 dbnv
                if aix==0:
                    ax[1,cx].remove()
                    ax[1,cx] = fig.add_subplot(6,len(taskaspects),len(taskaspects)+cx+1,projection='3d')
                axs = ax[1,cx]
                
                # activities
                axs.plot( trial_pre[:,0], trial_pre[:,1], trial_pre[:,2], lw=1,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( trial_on[:,0], trial_on[:,1], trial_on[:,2], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( trial_post[:,0], trial_post[:,1], trial_post[:,2], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( [trial_on[0,0]], [trial_on[0,1]], [trial_on[0,2]], 'o',color=taskcolors[aix][cx],alpha=alpha )
                
                # basis vectors
                if aix==0:
                    for bx,basisaspect in enumerate(basisaspects):
                        # axs.plot([0,B[bx,0]],[0,B[bx,1]],[0,B[bx,2]],lw=3,color=basisaspectcolors[bx])
                        axs.plot([0,B[0,bx]],[0,B[1,bx]],[0,B[2,bx]],lw=3,color=basisaspectcolors[bx])
                
                axs.set_xlabel('visual dbnv')
                axs.set_ylabel('context dbnv orth.n.')
                axs.set_zlabel('choice dbnv orth.n.')
                
                
                
                # 2D 2 dbnv
                axs = ax[2,cx]
                
                # activities
                axs.plot( trial_pre[:,0], trial_pre[:,1], lw=1,color=taskcolors[aix][cx],alpha=alpha)
                axs.plot( trial_on[:,0], trial_on[:,1], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( trial_post[:,0], trial_post[:,1], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                axs.plot( [trial_on[0,0]], [trial_on[0,1]], 'o',color=taskcolors[aix][cx],alpha=alpha )
                
                # basis vectors
                if aix==0:
                    for bx,basisaspect in enumerate(basisaspects):
                        # axs.plot([0,B[bx,0]],[0,B[bx,1]],lw=3,color=basisaspectcolors[bx])
                        axs.plot([0,B[0,bx]],[0,B[1,bx]],lw=3,color=basisaspectcolors[bx])

                
                axs.set_xlabel('visual dbnv')
                if cx==0: axs.set_ylabel('context dbnv orthonormalized')
                
                
                # 1D 1 DBNV
                for bx in range(n_ba):

                    # gather for identical axes limits
                    maxs[bx] = np.max(( maxs[bx], np.max(trial[skip_idx:-skip_idx,bx])   ))
                    mins[bx] = np.min(( mins[bx], np.min(trial[skip_idx:-skip_idx,bx])   ))


                    axs = ax[bx+3,cx]
                    
                    axs.plot( t_0, trial_pre[:,bx], lw=1,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( t_1, trial_on[:,bx], lw=3,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( t_oo, trial_post[:,bx], '--',lw=1,color=taskcolors[aix][cx],alpha=alpha )
                    axs.plot( [t_1[0]], [trial_on[0,bx]], 'o',color=taskcolors[aix][cx],alpha=alpha,label=classnames[cx][aix] )
                    
                    axs.fill_between(t_all, trial[skip_idx:-skip_idx,bx]-trial_e[skip_idx:-skip_idx,bx], \
                                            trial[skip_idx:-skip_idx,bx]+trial_e[skip_idx:-skip_idx,bx], \
                                     color=taskcolors[aix][cx],alpha=alpha/4.5)#,label=label)

                    figs.setxt(axs)
                    axs.set_xlim(skip_idx*T['dt']+T['starttime'], -skip_idx*T['dt']+T['endtime'])
                    

                    if aix==1:
                        if bx==0: axs.legend(frameon=False)
                    if cx==0: axs.set_ylabel(basisaspects[bx]+' dbnv'+[' orthonormalized',''][bx==0])
                    
                    
                    
                    
                    # now plot distance metrics:
                    if aix==0 and bx==2:
                        axs = ax[0,cx]

                        axs.plot(t_all,distances[bx][cx][0][skip_idx:-skip_idx],\
                                 color=taskcolors[0][cx],lw=2,alpha=0.8-bx*1/n_ba/2)
                        
                        s_e_m = np.max( np.stack( (distances[bx][cx][0][skip_idx:-skip_idx]-distances[bx][cx][1][skip_idx:-skip_idx],\
                                                   np.zeros(len(t_all))),axis=1     ), axis=1 )
                        axs.fill_between(t_all,s_e_m,\
                                               distances[bx][cx][0][skip_idx:-skip_idx],\
                                         color=taskcolors[0][cx],lw=2,alpha=(0.8-bx*1/n_ba/2)/8)
                        
                        
                        if bx==n_ba-1:
                            # axs.legend(frameon=False)
                            figs.setxt(axs)
                            figs.plottoaxis_stimulusoverlay(axs,T)
                            figs.plottoaxis_chancelevel(axs)
                            axs.set_title(comparison+' trial grouping')
                            if cx==0: axs.set_ylabel('distance in task variable dbnv\northonormalized basis')

                    

                axs.set_xlabel('trial time from stimulus onset [ms]')


        mins *=1.1
        maxs *=1.1
        for cx in range(len(taskaspects)):
            ax[1,cx].set_xlim(mins[0],maxs[0])
            ax[1,cx].set_ylim(mins[1],maxs[1])
            ax[1,cx].set_zlim(mins[2],maxs[2])

            ax[2,cx].set_xlim(mins[0],maxs[0])
            ax[2,cx].set_ylim(mins[1],maxs[1])

            for bx in range(n_ba):
                axs = ax[3+bx,cx]
                axs.set_ylim(mins[bx],maxs[bx])
                figs.plottoaxis_chancelevel(axs)
                figs.plottoaxis_stimulusoverlay(axs,T)
                
                
        fig.suptitle(dn+' trial averaged projected neuron dynamics onto dbnvs')
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'trials,over,dbnv,avg-%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)







    return
















def subspacedynamics_decisionvector_basis_initialconditions(dn,block):

    # show starting positions of activities for trials in different contexts
    # project everything onto   (vis,cx-o.n.) basis
    # show with different colors: GO hit,  and         GO miss
    #                             NOGO corr rej. and   NOGO false alarm.

    recalculate = 0 or globalrecalculate
    dump = 0
    doplot = 1 or globaldoplot

    taskaspects = ['attend visual','ignore visual']
    taskcolors = [ ['green','orange','lightseagreen','red'],\
                   ['darkgreen','darkorange','darkcyan','darkred'] ]
    
    travcolor = 'black'
        
    blv,bla = preprocess.getorderattended(dn)
    
    comparisongroups  = [ \
                          [ [  [blv[1]],  [45],[] ],   [  [blv[1]],  [135],[] ] ],\
                          [ [  [bla[1]],  [45],[] ],   [  [bla[1]],  [135],[] ] ], \
                        ]

    classnames = [['45°','135°'],['45°','135°'],['attend','ignore'],['attend','ignore']]
    
    perfnames = [ ['45° hit','45° miss','135° correct rejection','135° false alarm'], \
                  ['45° & audio hit','45° & audio miss','135° & audio correct rejection','135° & audio false alarm'] ]


    # projection basis vectors:
    depth_idx = int((1500*pq.ms / T['dt']).magnitude)
    n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
    basisaspects = ['visual','context']
    n_ba = len(basisaspects)
    basisaspectcolors = ['dodgerblue','fuchsia']   # np.array(taskcolors[0])[np.array([0,2,3],dtype=np.int16)]
    c_db_matrix = np.zeros((n_neurons,n_ba))       # prepare the dbnv vector holder
    for bx,aspect in enumerate(basisaspects):
        acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,aspect,'all'),'rb'))
        wx = int((len(acrossdecoder)-7)/n_neurons)
        coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
        # mean of on stimulus
        c_db_matrix[:,bx] = coeff[:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,0].mean(axis=1)

    # c_db_matrix projects to the 3 dbnvs each row:   c_db_matrix @ neural_input
    Q,R = np.linalg.qr(c_db_matrix)     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
    V = Q                          # activity transformer matrix
    # print('c_db',c_db_matrix.shape,'Q',Q.shape,'R',R.shape)
    B = Q.T @ c_db_matrix / np.linalg.norm(c_db_matrix,axis=0,keepdims=True)        # transformed dbnvs (in columns) in orthonormal coordinates (in rows)
    # print( 'V',V.shape, 'B', B.shape )




        
    projected_dynamics = []    # this will be (tasks,classes,[trials,2])
    for cx,comparison in enumerate(taskaspects):
        # now project each activity onto the dbnv basis vectorspace

        # collect neural responses
        acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],correctonly=1)
        acrossresponses.extend( preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],erroronly=1) )
        acrossresponses = [acrossresponses[i] for i in [0,2,1,3]]  # reorder for h,m,cr,fa

        projected_dynamic = []
        for aix,classresponse in enumerate(acrossresponses):
            # project the entire dynamics trajectory onto the PCs axes for each trial
            if len(classresponse)>0:
                X = np.dot(np.array(classresponse),V)  #    (trials,trajectory,dbnvs) = (trials,trajectory,neurons) * (,,neurons,dbnv)
            else:
                X = np.empty((0,n_trajectory,n_ba))
            projected_dynamic.append(X)
        
        projected_dynamics.append(projected_dynamic)



    if dump:
        pickle.dump( (projected_dynamics,V,B), open(cacheprefix+'subspaces/projections,dbnv-visual,context_%s.pck'%(dn),'wb'))


    if doplot:

        maxs = 0*np.ones(n_ba)
        mins = 0*np.ones(n_ba)

        alpha = 0.8
        alfac = 2
        # time parametrization
        skip_idx = 20
        t_all = block.segments[0].analogsignals[0].times[skip_idx:-skip_idx]
        t_pre = block.segments[0].analogsignals[0].times[skip_idx:T['stimstart_idx']+1]
        t_on = block.segments[0].analogsignals[0].times[T['stimstart_idx']:T['stimend_idx']+1]
        t_post = block.segments[0].analogsignals[0].times[T['stimend_idx']:-skip_idx]

        
        fig,ax = plt.subplots(2*len(taskaspects),8,figsize=((8)*8, (2*len(taskaspects))*8) )
        for cx,comparison in enumerate(taskaspects):
            
            
            for aix,classresponse in enumerate(projected_dynamics[cx]): # gp through classes
                
                if len(classresponse)==0: continue
            
                # trial average:
                
                trials = np.array(classresponse)            # all single trials
                trial_av = np.array(classresponse).mean(axis=0)        # trial averages
                trial_e = np.array(classresponse).std(axis=0)*2/np.sqrt(len(classresponse))
                
                trials_pre = trials[:,skip_idx:T['stimstart_idx']+1,:]
                trials_on = trials[:,T['stimstart_idx']:T['stimend_idx']+1,:]
                trials_post = trials[:,T['stimend_idx']:-skip_idx,:]

                trial_av_pre = trial_av[skip_idx:T['stimstart_idx']+1,:]
                trial_av_on = trial_av[T['stimstart_idx']:T['stimend_idx']+1,:]
                trial_av_post = trial_av[T['stimend_idx']:-skip_idx,:]


                for preon in [0,1]:
    
                    # 2D 2 dbnv
                    axs = ax[2*cx+preon,aix]
                    
                    axs.set_title(perfnames[cx][aix])
                    
                    # single trials:
                    # classresponse_selected = classresponse[trialset]
                    # for trialx,trial in enumerate(classresponse_selected):    # go through trials
                    # axs = ax[trialx,cx]
                    
                    
                    # print('{',cx,aix,comparison,'}',trials.shape,trials_pre.shape)
                    
                    
                    # activities (each single trial, criss crossing all over the place)
                    if preon==0:
                        axs.plot( trials_pre[:,:,0].T, trials_pre[:,:,1].T, lw=0.5,color=taskcolors[cx][aix],alpha=alpha/alfac)
                        axs.plot( [trials_pre[:,0,0]], [trials_pre[:,0,1]], 'd',color=taskcolors[cx][aix],alpha=alpha/alfac )
                    else:
                        axs.plot( trials_on[:,:,0].T, trials_on[:,:,1].T, lw=0.5,color=taskcolors[cx][aix],alpha=alpha/alfac )
                        # axs.plot( trials_post[:,:,0].T, trials_post[:,:,1].T, '--',lw=0.25,color=taskcolors[cx][aix],alpha=alpha/alfac )
                        axs.plot( [trials_on[:,0,0]], [trials_on[:,0,1]], 'o',color=taskcolors[cx][aix],alpha=alpha/alfac )
    
                    # trial averaged activities
                    axs.plot( [trial_av_pre[0,0]], [trial_av_pre[0,1]], 'd',markersize=15,color=travcolor,alpha=alpha,label='-1500 ms' )
                    axs.plot( trial_av_pre[:,0], trial_av_pre[:,1], lw=2,color=travcolor,alpha=alpha,label='trial av. pre' )
                    axs.plot( [trial_av_on[0,0]], [trial_av_on[0,1]], 'o',markersize=15,color=travcolor,alpha=alpha, label='0 ms')
                    axs.plot( trial_av_on[:,0], trial_av_on[:,1], lw=5,color=travcolor,alpha=alpha,label='trial av. on' )
                    # axs.plot( trial_av_post[:,0], trial_av_post[:,1], '--',lw=2,color=taskcolors[cx][aix],alpha=alpha )
                    
                    
                    
                    # basis vectors
                    if aix==0 and cx==0 and preon==0:
                        axs.legend(frameon=False)
                        for bx,basisaspect in enumerate(basisaspects):
                            xoff = -np.sign(B[0,bx]) * 1.9
                            yoff = -np.sign(B[1,bx]) * 1.9
                            axs.plot([xoff,B[0,bx]+xoff],[yoff,B[1,bx]+yoff],lw=3,color=basisaspectcolors[bx])
    
                    
                    if cx==1 and preon==1: axs.set_xlabel('visual dbnv')
                    if aix==0: axs.set_ylabel('%s\nsingle trials: %s\ncontext dbnv orthonormalized'%(comparison,['pre','on'][preon]))
                    



                for bx,basisaspect in enumerate(basisaspects):
                
                # for aix,classresponse in enumerate(projected_dynamics[cx]): # gp through classes
    
                    # plot 1D projections only
                    axs = ax[cx*2+bx,4+aix]

                    # single trial activities
                    axs.plot( t_pre, trials_pre[:,:,bx].T, lw=0.5,color=taskcolors[cx][aix],alpha=alpha/alfac)
                    for l in range(len(trials_pre)):
                        axs.plot( t_pre[0], trials_pre[l,0,bx], 'd',color=taskcolors[cx][aix],alpha=alpha/alfac )
                        axs.plot( t_on[0],trials_on[l,0,bx], 'o',color=taskcolors[cx][aix],alpha=alpha/alfac )
                    axs.plot( t_on, trials_on[:,:,bx].T, lw=0.5,color=taskcolors[cx][aix],alpha=alpha/alfac )
                    axs.plot( t_post, trials_post[:,:,bx].T, '--',lw=0.5,color=taskcolors[cx][aix],alpha=alpha/alfac )
    
                    # trial averaged activities
                    axs.plot( t_pre[0], trial_av_pre[0,bx],  'd',markersize=15,color=travcolor,alpha=alpha)
                    axs.plot( t_pre, trial_av_pre[:,bx], lw=2,color=travcolor,alpha=alpha)
                    axs.plot( [t_on[0]], [trial_av_on[0,bx]],  'o',markersize=15,color=travcolor,alpha=alpha)
                    axs.plot( t_on, trial_av_on[:,bx], lw=5,color=travcolor,alpha=alpha)
                    axs.plot( t_post, trial_av_post[:,bx],'--',lw=2, color=travcolor,alpha=alpha)
    
                    
    
                    axs.set_ylim(-2,2)
                    figs.setxt(axs)
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    figs.plottoaxis_chancelevel(axs)
    
    
                    if aix==0: axs.set_ylabel(taskaspects[cx]+'\n'+basisaspect+' dbnv')
                    if bx==1: axs.set_xlabel('trial time from stimulus onset [ms]')

                    axs.set_title(perfnames[cx][aix])


        # mins *=1.1
        # maxs *=1.1
        for aix in [0,1,2,3]:
            for cx in [0,1,2,3]:
                # ax[cx,aix].set_xlim(mins[0],maxs[0])
                # ax[cx,aix].set_ylim(mins[1],maxs[1])
                ax[cx,aix].set_xlim(-2,2)
                ax[cx,aix].set_ylim(-2,2)

                
        fig.suptitle(dn+' n=%d, context and choice dependent single trial neural dynamics projected onto dbnvs'%n_neurons)
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'trials,over,dbnv,ics-%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)







    return





def subspacedynamics_decisionvector_basis_initialconditions_highres(dn,block):

    # show starting positions of activities for trials in different contexts
    # project everything onto   (vis,cx-o.n.) basis      (cx-o.n. means closest vector to context, but orthogonal to visual)
    # show with different colors: GO hit,  and         GO miss
    #                             NOGO corr rej. and   NOGO false alarm.

    recalculate = 1 or globalrecalculate
    dump = 0
    doplot = 1 or globaldoplot

    taskaspects = ['attend visual','ignore visual']
    taskcolors = [ ['green','orange','lightseagreen','red'],\
                   ['darkgreen','darkorange','darkcyan','darkred'] ]
    
    travcolor = 'black'
        
    blv,bla = preprocess.getorderattended(dn)
    
    comparisongroups  = [ \
                          [ [  [blv[1]],  [45],[] ],   [  [blv[1]],  [135],[] ] ],\
                          [ [  [bla[1]],  [45],[] ],   [  [bla[1]],  [135],[] ] ], \
                        ]

    classnames = [['45°','135°'],['45°','135°'],['attend','ignore'],['attend','ignore']]
    
    perfnames = [ ['total','45° hit','45° miss','135° correct rejection','135° false alarm'], \
                  ['total','45° & audio hit','45° & audio miss','135° & audio correct rejection','135° & audio false alarm'] ]


    # projection basis vectors:
    depth_idx = int((1500*pq.ms / T['dt']).magnitude)
    n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
    basisaspects = ['context']
    n_ba = len(basisaspects)
    basisaspectcolors = ['fuchsia']   # np.array(taskcolors[0])[np.array([0,2,3],dtype=np.int16)]
    c_db_matrix = np.zeros((n_neurons,n_trajectory))       # prepare the dbnv vector holder

    basisaspect = basisaspects[0]
    acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,basisaspect,'all'),'rb'))
    wx = int((len(acrossdecoder)-7)/n_neurons)
    coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
    # mean of on stimulus
    c_db_matrix = coeff[:,:,0].T      # neurons times trajectory
    

    # c_db_matrix projects to the 3 dbnvs each row:   c_db_matrix @ neural_input
    Q,R = np.linalg.qr(c_db_matrix)     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
    print('c_db_matrix',c_db_matrix.shape,'Q',Q.shape,'R',R.shape)
    V = Q                          # activity transformer matrix
    B = (Q.T @ c_db_matrix) / np.linalg.norm(c_db_matrix,axis=0,keepdims=True)        # transformed dbnvs (in columns) in orthonormal coordinates (in rows)
    print( 'V',V.shape, 'B', B.shape )




        
    projected_dynamics = []    # this will be (tasks,classes,[trials,2])
    for cx,comparison in enumerate(taskaspects):
        # now project each activity onto the dbnv basis vectorspace

        # collect neural responses
        acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],correctonly=1)
        acrossresponses.extend( preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx],erroronly=1) )
        
        # create
        acrossresponses_toconcat = [acrossresponse for acrossresponse in acrossresponses if len(np.array(acrossresponse).shape)==3]
        acrossresponses.append( np.concatenate(  acrossresponses_toconcat  ) )
        # reorder for all, and usual behaviour order
        acrossresponses = [acrossresponses[i] for i in [4,0,2,1,3]]  # reorder for h,m,cr,fa

        projected_dynamic = []
        for aix,classresponse in enumerate(acrossresponses):
            # project the entire dynamics trajectory onto the PCs axes for each trial
            if len(classresponse)>0:
                A = np.array(classresponse)#.swapaxes(1,2)
                X = np.dot(A,V.T)      #    (trials,trajectory,dbnv_projection_times) = (trials,trajectory,neurons) * (,,neurons,dbnv_projection_times)
                X = X / np.linalg.norm(V.T,axis=0,keepdims=True)       # norm down the projection operator
            else:
                X = np.empty((0,n_trajectory,n_trajectory-wx))
            projected_dynamic.append(X)

        
        projected_dynamics.append(projected_dynamic)



    if doplot:

        maxs = 0*np.ones(n_ba)
        mins = 0*np.ones(n_ba)

        alpha = 0.8
        alfac = 2
        # time parametrization
        skip_idx = 20
        t_all = block.segments[0].analogsignals[0].times[skip_idx:-skip_idx]
        t_pre = block.segments[0].analogsignals[0].times[skip_idx:T['stimstart_idx']+1]
        t_on = block.segments[0].analogsignals[0].times[T['stimstart_idx']:T['stimend_idx']+1]
        t_post = block.segments[0].analogsignals[0].times[T['stimend_idx']:-skip_idx]

        
        fig,ax = plt.subplots(len(taskaspects),5,figsize=(5*8, len(taskaspects)*8) )
        
        for cx,comparison in enumerate(taskaspects):
            
            for aix,classresponse in enumerate(projected_dynamics[cx]): # gp through classes
                
                axs = ax[cx,aix]
                axs.set_title(perfnames[cx][aix])
                if aix==0: axs.set_ylabel(comparison + '\nprojection to dbnv at')
                if cx==1: axs.set_xlabel('activity projection from')


                if len(classresponse)==0: continue
            
                # trial average:
                
                trials = np.array(classresponse)            # all single trials
                trial_av = np.array(classresponse).mean(axis=0)        # trial averages;  [projection from activity here, projection to dbnv here]
                trial_e = np.array(classresponse).std(axis=0)*2/np.sqrt(len(classresponse))
                
                
                cmap = figs.getcorrcolormap('correlation')
                cf = axs.pcolormesh(trial_av[:,:].T, vmin=-2,vmax=2, cmap=cmap)        
                
                axs.set_aspect('equal')
                
                ticks=[150,450]
                ticklabels=['0 ms','3000 ms']
                axs.set_xticks(ticks)
                axs.set_xticklabels(ticklabels)                
                axs.set_yticks(ticks)
                axs.set_yticklabels(ticklabels,rotation=90)                
                

                # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                fig.colorbar(cf,ax=axs)#,ticks=np.arange(-0.3,0.31,0.1))



                
        fig.suptitle(dn+' n=%d, context and choice dependent neural dynamics along timecourse (horizontal)\nprojected onto context dbnvs (color) found along timecourse (vertical)'%n_neurons)
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'trials,over,dbnv-highres,ics-%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)







    return




















def subspacedynamics_DVPCA(dn,block):
    # show dynamics of the context DV projected onto a couple of principal component axis.

    dump = 0
    doplot = 1 or globaldoplot

    

    n_components = 3
    n_clusters = 8

    # projection basis vectors:

    n_trajectory,n_neurons = block.segments[0].analogsignals[0].shape
    basisaspects = ['context']
    n_ba = len(basisaspects)
    basisaspectcolors = ['mediumvioletred']   # np.array(taskcolors[0])[np.array([0,2,3],dtype=np.int16)]
    # c_db_matrix = np.zeros((n_neurons,n_trajectory))       # prepare the dbnv vector holder

    basisaspect = basisaspects[0]
    acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,basisaspect,'all'),'rb'))
    wx = int((len(acrossdecoder)-7)/n_neurons)
    coeff = np.reshape(np.array(acrossdecoder[7:]), (wx,n_neurons,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)
    # mean of on stimulus
    DV = coeff[:,:,0].T      # (observations in time, neurons)
    print(DV.shape)
    

    # # c_db_matrix projects to the 3 dbnvs each row:   c_db_matrix @ neural_input
    # Q,R = np.linalg.qr(c_db_matrix)     # Q holds the orthonormal basis vectors as columns, R is the transformed c_db_means_matrix
    # print('c_db_matrix',c_db_matrix.shape,'Q',Q.shape,'R',R.shape)
    # V = Q                          # activity transformer matrix
    # B = (Q.T @ c_db_matrix) / np.linalg.norm(c_db_matrix,axis=0,keepdims=True)        # transformed dbnvs (in columns) in orthonormal coordinates (in rows)
    # print( 'V',V.shape, 'B', B.shape )



    times = block.segments[0].analogsignals[0].times[:-wx]
    X = np.zeros((DV.shape[0],n_components))
    
    S,E,V = nedi.pca(DV,n_components=n_neurons,ret='singularvalue')
    print( DV[:,0].shape, DV[0,:].shape, V.shape, X.shape )

    # print(V.T)
    S = np.cumsum(S)/np.sum(S)*100
    for t in range(DV.shape[0]):
        X[t,:] = DV[t,:] @ V[:n_components,:].T
    print('eigenpectra cumulative % ', S)
    print('explained variance cumulative % ', np.cumsum(E)*100)


    # print(DV[:15,:])
    
    
    
    
    # cluster the findings
    from sklearn.cluster import AgglomerativeClustering
    
    
    connectivity = sp.sparse.diags( [ np.arange(X.shape[0]-1), np.arange(X.shape[0]-1) ], [-1,1]  )
    # connectivity = sp.sparse.diags( [ np.arange(X.shape[0]-1)], [1]  )
    clustering = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward',connectivity=connectivity)
    C = clustering.fit_predict(X)
    
    
    
    tl = [[T['start_idx'],T['stimstart_idx']],      [T['stimstart_idx'],T['stimend_idx']],\
                      [T['stimend_idx'],T['end_idx']]]
    # colors = ['navy','darkgreen','mediumvioletred','darkorange']
    alphas = [0.5,1,0.75]
    alphas = [0.8,0.8,0.8]
    lws = [1,1,1]
    colors = ['black','mediumvioletred','grey']
    labels = ['pre','on','post']
    
    
    
    
    
    if 1:
        # fig,ax = plt.subplots(1,1,figsize=(1*8,1*8))
        # fig = plt.figure(figsize=(2*12,2*12))
        fig, ax = plt.subplots(3,6,figsize=(6*8,3*8))

        pairs = [[0,1],[0,2],[1,2]]

        for px in [0,1,2]:
            axs = ax[0,px]
                
            for tx in [0,1,2]:
                axs.plot(X[ tl[tx][0]:tl[tx][1], pairs[px][0]], X[ tl[tx][0]:tl[tx][1], pairs[px][1]],\
                         'o-',lw=lws[tx], alpha=alphas[tx], color=colors[tx],label=labels[tx],\
                             markersize=3)
            axs.set_xlabel('PC %d'%(pairs[px][0]+1),labelpad=-20)
            axs.set_ylabel('PC %d'%(pairs[px][1]+1),labelpad=-20)
            xl = axs.get_xlim(); yl = axs.get_ylim()
            axs.set_xticks([-0.5,-0.1,0,0.1,0.5]); axs.set_yticks([-0.5,-0.1,0,0.1,0.5])
            axs.set_xlim(xl); axs.set_ylim(yl); 
            

        ax[0,3].remove()
        ax[0,3] = fig.add_subplot(3,6,4,projection='3d')
        axs = ax[0,3]
        for tx in [0,1,2]:
            axs.plot(X[ tl[tx][0]:tl[tx][1], 0], X[ tl[tx][0]:tl[tx][1], 1],  X [tl[tx][0]:tl[tx][1], 2],\
                     'o-',lw=lws[tx], alpha=alphas[tx], color=colors[tx],label=labels[tx],\
                         markersize=3)
        axs.set_xlabel('PC 1'); axs.set_ylabel('PC 2'); axs.set_zlabel('PC 3')
        xl = axs.get_xlim(); yl = axs.get_ylim(); zl = axs.get_zlim(); 
        axs.set_xticks([-0.5,-0.1,0,0.1,0.5]); axs.set_yticks([-0.5,-0.1,0,0.1,0.5]); axs.set_zticks([-0.5,-0.1,0,0.1,0.5]); 
        axs.set_xlim(xl); axs.set_ylim(yl); axs.set_zlim(zl);
        axs.legend(frameon=False)
        # fig.tight_layout()
        fig.suptitle(dn+' N=%d\ncontext DV projections onto PCs over trial time-course variability'%n_neurons+\
                     '\n cumulative eigen %%: %4.1f %4.1f %4.1f'%(S[0],S[1],S[2]))




        axs = ax[0,4]

        for tx in [0,1,2]:
            axs.plot(times[tl[tx][0]:tl[tx][1]], np.array(acrossdecoder[1])[tl[tx][0]:tl[tx][1],0],\
                     color=colors[tx])
        axs.set_ylim(0.45,1.05)
        axs.set_yticks([0.5,1.0])
        figs.plottoaxis_chancelevel(axs,0.5)
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.set_ylabel('context accuracy',labelpad=-30)
        axs.set_title('pre-on-post')

        figs.invisibleaxes(axs)





        for px in [0,1,2]:
            axs = ax[1,px]
                
            for clidx in range(n_clusters):
                axs.plot( X[ C==clidx, pairs[px][0]], X[ C==clidx, pairs[px][1] ],\
                         'o-',markersize=3)
            axs.set_xlabel('PC %d'%(pairs[px][0]+1),labelpad=-20)
            axs.set_ylabel('PC %d'%(pairs[px][1]+1),labelpad=-20)
            xl = axs.get_xlim(); yl = axs.get_ylim()
            axs.set_xticks([-0.5,-0.1,0,0.1,0.5]); axs.set_yticks([-0.5,-0.1,0,0.1,0.5])
            axs.set_xlim(xl); axs.set_ylim(yl); 
            

        ax[1,3].remove()
        ax[1,3] = fig.add_subplot(3,6,10,projection='3d')
        axs = ax[1,3]
        for clidx in range(n_clusters):
            axs.plot( X[ C==clidx, pairs[px][0]], X[ C==clidx, pairs[px][1] ], X[ C==clidx, pairs[px][1] ],\
                         'o-',markersize=3 )
        axs.set_xlabel('PC 1'); axs.set_ylabel('PC 2'); axs.set_zlabel('PC 3')
        xl = axs.get_xlim(); yl = axs.get_ylim(); zl = axs.get_zlim(); 
        axs.set_xticks([-0.5,-0.1,0,0.1,0.5]); axs.set_yticks([-0.5,-0.1,0,0.1,0.5]); axs.set_zticks([-0.5,-0.1,0,0.1,0.5]); 
        axs.set_xlim(xl); axs.set_ylim(yl); axs.set_zlim(zl);
        # axs.legend(frameon=False)
        





        axs = ax[1,4]
        for clidx in range(n_clusters):
            mask = C==clidx
            
            axs.plot(times[mask], np.array(acrossdecoder[1])[mask,0])
        axs.set_ylim(0.45,1.05)
        axs.set_yticks([0.5,1.0])
        figs.plottoaxis_chancelevel(axs,0.5)
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.set_ylabel('context accuracy',labelpad=-30)
        axs.set_title('unsup. clusters')

        figs.invisibleaxes(axs)







        # distance based clusters


        
        axs = ax[2,5]
        kds = [5,5,5] # euclidean distance
        labels = ['eu.dist.','valleys','peaks']
        D = np.zeros(X.shape[0])
        
        for kx,kd in enumerate(kds):
            for t in range(len(D)):
                D[t] = np.sqrt(  np.mean(  X[ max(0,t-kd)  :   min(len(D),t+kd), :   ] @ X[t,:] )   )
            
            # if kx==0:
            #     mnix, mxix = neph.localminimamaxima(D,sm=5)
            #     print(mnix,mxix)
            #     axs.plot(times[mnix],D[mnix],'ob',markersize=15)
            #     axs.plot(times[mxix],D[mxix],'or',markersize=15)

            # if kx==1: D = neph.smooth(D,3,'same'); print(D.shape)
            # if kx==2: D = neph.smooth(D,5,'same'); D = np.convolve(D,np.array([-0.5,5,-0.5])/6,'same'); print(D.shape)
            # if kx==3: D = np.r_[0,0,np.diff(neph.smooth(D,10,'same'),2)]*100; clustm = np.where(np.abs(D-np.mean(D))<0.04)[0]
            if kx==1:
                D = neph.smooth(D,3,'same') - neph.smooth(D,60,'same'); D = D - neph.smooth(D,60,'same')
                clustm = D<np.mean(D)
                D[D>np.mean(D)] = np.nan
            if kx==2:
                D = neph.smooth(D,3,'same') - neph.smooth(D,60,'same'); D = D - neph.smooth(D,60,'same')
                D[D<=np.mean(D)] = np.nan

            axs.plot(times,D,lw=2,label=labels[kx])
        
        axs.legend(frameon=False)
        axs.set_title('mov.av. neighbours\' eucl. dist.')
        figs.setxt(axs)
        axs.set_ylim(axs.get_ylim())
        figs.plottoaxis_stimulusoverlay(axs,T)

        figs.invisibleaxes(axs)





        C = np.zeros(clustm.shape)
        counter = 0
        current = clustm[0]
        for tx in range(len(clustm)):
            # print(tx,clustm[tx],current,counter)
            if clustm[tx] ^ current:
                current = clustm[tx]
                counter += 1
            C[tx] = counter
        # print(C)
        
        n_clusters = counter



        for px in [0,1,2]:
            axs = ax[2,px]
                
            for clidx in range(n_clusters):
                axs.plot( X[ C==clidx, pairs[px][0]], X[ C==clidx, pairs[px][1] ],\
                         'o-',markersize=3)
            axs.set_xlabel('PC %d'%(pairs[px][0]+1),labelpad=-20)
            axs.set_ylabel('PC %d'%(pairs[px][1]+1),labelpad=-20)
            xl = axs.get_xlim(); yl = axs.get_ylim()
            axs.set_xticks([-0.5,-0.1,0,0.1,0.5]); axs.set_yticks([-0.5,-0.1,0,0.1,0.5])
            axs.set_xlim(xl); axs.set_ylim(yl); 


        ax[2,3].remove()
        ax[2,3] = fig.add_subplot(3,6,16,projection='3d')
        axs = ax[2,3]
        for clidx in range(n_clusters):
            axs.plot( X[ C==clidx, pairs[px][0]], X[ C==clidx, pairs[px][1] ], X[ C==clidx, pairs[px][1] ],\
                         'o-',markersize=3 )
        axs.set_xlabel('PC 1'); axs.set_ylabel('PC 2'); axs.set_zlabel('PC 3')
        xl = axs.get_xlim(); yl = axs.get_ylim(); zl = axs.get_zlim(); 
        axs.set_xticks([-0.5,-0.1,0,0.1,0.5]); axs.set_yticks([-0.5,-0.1,0,0.1,0.5]); axs.set_zticks([-0.5,-0.1,0,0.1,0.5]); 
        axs.set_xlim(xl); axs.set_ylim(yl); axs.set_zlim(zl);
        # axs.legend(frameon=False)




        axs = ax[2,4]
        for clidx in range(n_clusters):
            mask = C==clidx
            
            axs.plot(times[mask], np.array(acrossdecoder[1])[mask,0])
        axs.set_ylim(0.45,1.05)
        axs.set_yticks([0.5,1.0])
        figs.plottoaxis_chancelevel(axs,0.5)
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.set_ylabel('context accuracy',labelpad=-30)
        axs.set_title('distance based clusters')

        figs.invisibleaxes(axs)




        ax[0,4].axis('off')
        ax[0,5].axis('off')
        ax[1,5].axis('off')






        save = 1 or globalsave
        if save:
            fig.savefig(resultpath+'context,DV,timecoursePCA-%s-%dms%dms_%s'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn)+ext)















