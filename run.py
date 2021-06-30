# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:17:21 2019

@author: mahajnal
"""


from config import *


# import preprocess

from physiology import *
from decodersubspace import *
from unsupervised import *

from glmpredictneural import *

from memorymodels import *
from spikemodels import *

from aggregateallmice import *

    






def tests(dn,block):
    
    fig,ax = plt.subplots(1,3,figsize=(32,8))
    
    for k in [0,1,2]:
        axs = ax[k]
        axs.plot(block.segments[k+30].analogsignals[0][1500:1600,:])
        print(block.segments[k+30].analogsignals[0][1500:1510,:])
    
    return

    trialnum = 30+5
    neuron = 8
    bins = 100
    spiketrain = np.array(block.segments[trialnum].spiketrains[neuron])
    
    spikecounts = np.array([ np.sum((spiketrain>t0) & (spiketrain<t0+bins)) for t0 in np.arange(0,3000,bins) ])
    
    print(dn,'neuron %d, trial: %d'%(neuron,trialnum))
    print(spiketrain)
    print(spikecounts)


# ***************************************************************************







def select_data():
    
    
    # datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']

    
    
    
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT022','DT030','DT031','DT032']
    # DT022 'list index out of range' is a very few trial mouse!!!!!! So the followings do not contain them:         # it's not the case anymore since ks2 double check!
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    
    # ME108 has technically zero spikes, the whole recordings is noise continuously since ks2 sorting!!!
    # ME112 has zero good nondrifting single activity units since ks2 sorting!!!
    
    # DT008 is not possible with all subspaces of receptive field characterization, as it does not contain 45 and 135 degrees!
    # datanames = ['ME108','ME110','ME112','ME113','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']

    
    # furthur problem is audio noise contamination into LFP in DT030,31,32; below don't contain them:
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021']
    # not that audio noise is eliminated with the ks2 sorting now

    # datanames = ['DT008','DT019']
    
    # runspeed:
    # datanames = ['DT009','ME108','ME110','ME112','ME113','DT030','DT031','DT032']
    # datanames = ['DT009','DT014','DT017','DT018','DT019','DT021']

    # nucleus thalami lateralis posterioris
    # datanames = ['LP001','LP003','LP004','LP005']
    # datanames = ['LP001']
    # datanames = ['LP003']
    # datanames = ['LP004']
    # datanames = ['LP005']


    # cortex cingulatis anterioris
    # datanames = ['AC001','AC003','AC004','AC006','AC007','AC008','AC009']
    # datanames = ['AC001','AC003','AC006','AC007']
    # datanames = ['AC001','AC003']
    # datanames = ['AC006','AC007']
    # datanames = ['AC001']
    # datanames = ['AC003']
    # datanames = ['AC006']
    # datanames = ['AC007']



    # individuals
    # datanames = ['ME103']
    # datanames = ['ME110']
    # datanames = ['ME113']
    # datanames = ['DT008']
    # datanames = ['DT009']
    # datanames = ['DT014']
    # datanames = ['DT017']
    # datanames = ['DT018']
    # datanames = ['DT019']
    # datanames = ['DT020']
    # datanames = ['DT021']
    # datanames = ['DT022']
    # datanames = ['DT030']
    # datanames = ['DT031']
    # datanames = ['DT032']
    
    # datanames = ['AC001']
    # datanames = ['AC003']
    # datanames = ['AC004']
    # datanames = ['AC006']
    datanames = ['AC007']
    
    # # datanames = ['AC008']
    # # datanames = ['AC009']
    # datanames = ['AC006','AC007']
    # datanames = ['AC001','AC003','AC004','AC006','AC007']

    # datanames = ['DT008','DT009','DT014']
    # datanames = ['DT009','DT017','DT018','DT022']
    # datanames = ['DT009','DT014','DT017','DT019']
    # datanames = ['DT017','DT019']
    
    
    # all
    # datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']
    
    # correct runspeed:
    # datanames = ['ME110','ME113','DT009','DT014','DT017','DT030','DT031','DT032']
    
    # split into 2
    # datanames = ['DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022']
    # datanames = ['ME103','ME110','ME113','DT030','DT031','DT032']
    
    # split into 4
    # datanames = []
    # datanames.append(['ME103','ME110','ME113','DT008'])
    # datanames.append(['DT009','DT014','DT017','DT018'])
    # # datanames.append(['DT014','DT017','DT018'])
    # datanames.append(['DT019','DT020','DT021','DT022'])
    # datanames.append(['DT030','DT031','DT032'])

    # split into 5
    # datanames = []
    # datanames.append(['ME103','ME110','ME113'])
    # datanames.append(['DT008','DT009','DT014'])
    # datanames.append(['DT017','DT018','DT019'])
    # datanames.append(['DT020','DT021','DT022'])
    # datanames.append(['DT030','DT031','DT032'])
    


    # concatenate them again for row
    # datanames = np.concatenate(datanames)


    # longer timewindow
    # datanames = ['ME110','DT008','DT030','DT032']
    # datanames = ['DT030','DT032']

    
    # spike sorting
    # datanames = ['DT008','DT009']
    # datanames = ['ME108','ME110']
    # datanames = ['ME113']
    # datanames = ['ME108','ME112','DT017','DT018']
    # datanames = ['DT014','DT032']
    # datanames = ['DT030','DT031','DT032']
    # datanames = ['DT008','DT009','DT014','DT017']

    # datanames = ['ME103','ME113']#,'DT022']
    
    # spike sort grouping
    # datanames = ['ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    # datanames = ['DT009','DT014','DT017','DT030','DT031','DT032']        # for subspace
    # datanames = ['DT021','DT022']
    # datanames = ['DT020','DT021','DT022']
    # datanames = ['DT018','DT019','DT021']
    
    # datanames = ['ME110','ME113','DT008']
    # datanames = ['DT009','DT014','DT017']
    # datanames = ['DT018','DT019','DT021']
    # datanames = ['DT030','DT031','DT032']



    # datanames = ['GT102']         # this is the double block experiment: v,av,a,aa,v,av,a,aa,rf
    
    
    
    
    


    # smart group:
#    datanames = ['DT008','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    
    # silly group:
#    datanames = ['ME108','ME110','ME112','ME113']


    # post July 5th group
#    datanames = ['ME103','DT013','DT015','DT020']
#    datanames = ['DT013']


    # predictors:
#    datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    

    
    
    return datanames
    







def analysis(datanames, multiple=False):


    # choose analyis type
    # multiple = True
    
    # aggregating over multiple mice routines
    if multiple:
        
        # full routines:
            


        # aggregatemice_relevantirrelevant_behaviour(datanames)
        # aggregatemice_relevantirrelevant_singletrialneuralbehaviour(datanames)
        aggregatemice_subspaces_timeresolved_behaviour_ACC(datanames)



        # cluster_celltypes(datanames)

        # aggregatemice_decoderaccuracyaverages(datanames)
        
        
        # decoderstobehaviouralperformance(datanames);
        # decodersbehaviourwithtraining(datanames);
        # preprocess.loadtrainingbehaviouraldata(dn)
        
        # aggregatemice_representationmeasures_featuresvector(datanames)
        
        # aggregatemice_subspaceangles(datanames)
        # aggregatemice_spontaneousdecoder(datanames,'allexpcond')
        # aggregatemice_spontaneousdecoder(datanames,'attendignore')
        # aggregatemice_crossgradientdecoder(datanames)
        # aggregatemice_cueperiod_neuraltobehaviour(datanames)
    
        # aggregatemice_layercelltypecontributions(datanames)
        # aggregatemice_layertypecontributions_firingratestatistics(datanames)
        # aggregatemice_layertypecontributions_firingratestatistics(datanames,'allexpcond')
        # aggregatemice_layertypecontributions_firingratestatistics(datanames,'context')
        # aggregatemice_layertypecontributions_predictiveglm(datanames)


        # aggregatemice_behaviourvscontextprobability(datanames)        # change to d'


        # aggregate_spontaneousprojections(datanames)
        # aggregatemice_attendignore(datanames)
        
        # aggregatemice_subspaces_PCA(datanames)
        # aggregatemice_subspaces_PCA_saveonly(datanames)
        # aggregatemice_subspaces_PCA_explainedvariance(datanames)
        # aggregatemice_subspaces_PCA_signalvariance(datanames)

        # aggregatemice_subspaces_behaviour(datanames)



        # aggregatemice_subspace_projectiontodbnvs_behaviour(datanames)

        # aggregatemice_subspace_anglesbetweendbnvs(datanames)
        # aggregatemice_subspace_contextrotationdecay(datanames)
        # aggregatemice_calculatedecayconstant(datanames)
        # aggregatemice_contextpreonsubpop(datanames)


        # aggregatemice_chancelevel(datanames)




        
        return               # this is to avoid the following preloaded type aggregate routines


        # preload blocks routines:
        blocks = []
        for dn in datanames:
            block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=False)
            blocks.append(  block  )
    
        # nebay.aggregate_analogsignals(datanames, blocks, T)
        aggregatemice_crossorentation(datanames, blocks, T)







    # single mouse routines
    else:

        # for dn in datanames:
            dn = datanames
            # preprocess.convert_trialstocsv(dn); return            # use this once to get the trials_good_start.csv for the mouse

            # use this to export spike times in a simple format
            # preprocess.loaddatamouse(dn,T,continuous_method,recalculate=True,exportspiketrains=True)


            # loads or creates spike counts, needed for all processes below:
            block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=True,recalculate=False)
            
            # export trial cut spike trains from the converted neo format:
            # preprocess.exportspiketrainsfromneo(dn,block)
            
            


            
            # tests(dn,block)

            # decoder_driftsearch(dn,block)


            # print(dn,'neurons:',block.segments[0].analogsignals[0].shape[1])
            # preprocess.exporttrialtimes(dn)
            # preprocess.exportdecoderaccuracies(dn)
            # preprocess.exporthdf5(block,dn)
            

            # orientationselectivity(dn,block)




            # compareacrossconditions_maxnormdiffdecoder(dn,block,examine='allexpcond')
            # compareacrossconditions_maxnormdiffdecoder(dn,block,examine='gonogo')
            # compareacrossconditions_maxnormdiffdecoder(dn,block,examine='attendignore')
            # compareacrossconditions_maxnormdiffdecoder(dn,block,examine='character')




            # decoder_celltypes(dn,block)    # test context decoder on only putative excitatory and inhibitory cells, but only with high number of units




            # comparisonacrossconditions_choice(dn,block)
            # angles_between_decoders(dn,block)

            # decoder_crosstest(dn,block)
            # decoder_crosstest_highres(dn,block)

            
            # decoder_crosstestgradient(dn,block)
            # decoder_crosscontexttest(dn,block)
            # decoder_withincontext_saveprojections(dn,block)

            # for trainidx in [0,1]:
            #     for testidx in [0,1]:
            #         decoder_crosscontexttest(dn,block,trainidx=trainidx,testidx=testidx)
            # decoder_crosscontexttest(dn,block,trainidx=0,testidx=1)
            # decoder_crosscontexttest(dn,block,trainidx=1,testidx=0)
            # decoder_crosscontexttest(dn,block)
            
            


            # decoder_acc_relevantirrelevant(dn,block)
            # decoder_acc_singletrialneuralbehaviour(dn,block)


            
            
            # decoder_behaviour(dn,block)
            # latent_behaviour(dn)
            
            
            
    
            # joint_differentiabilityrunspeed(dn,block)
            
            # runspeed_regression(dn,block)
            # runspeed_taskconditions(dn,block)

            # corr_runspeed_dbnv(dn,block)
            # decoder_equalizedrunspeeddistribution(dn,block)


            # subspaces(dn,block,examine='attendignore')
            # subspaces(dn,block,examine='allexpcond')
            
            # subspaces_projections(dn,block,examine='allexpcond',play='visual',preon=0)
            # subspaces_projections(dn,block,examine='allexpcond',play='visual',preon=1)
            # subspaces_projections(dn,block,examine='allexpcond',play='audio',preon=0)
            # subspaces_projections(dn,block,examine='allexpcond',play='audio',preon=1)
            # subspaces_projections(dn,block,examine='allexpcond',play='choice',preon=1)



            # subspaces_behaviour_ACC(dn,block)
            subspaces_timeresolved(dn,block)



            # subspaces_decode_nullspace1visual1(dn,block)
            # subspaces_decode_nullspacerecurrent(dn,block)

            # subspaces_orthogonalcontrol_shuffle(dn,block)
            # subspaces_leaveoneout(dn,block)

            # subspaces_neuralparticipation(dn,block)
            # subspaces_residualvariance(dn,block,examine='allexpcond')
            # subspaces_examples(None,None)
            
            # subspaces_PCA(dn,block,preon=0)
            # subspaces_PCA(dn,block,preon=1)
            
            # subspacedynamics_PCA(dn,block)
            # subspacedynamics_PCA_context(dn,block)
            # subspacedynamics_PCA_signalnoisecorrelation(dn,block)
            
            # subspacedynamics_DVPCA(dn,block)
            


            # subspacedynamics_decisionvector_basis(dn,block)
            # subspacedynamics_decisionvector_basis_initialconditions(dn,block)
            # subspacedynamics_decisionvector_basis_initialconditions_highres(dn,block)


            # predictive_glm(dn,block)
            
            

            # reward(dn,block)




            # tcafactors(dn,block)


            # variationallatentgaussianprocess(dn,block)





def main():
    usedask = 0
    multiple = 0

    datanames = select_data()
    if usedask:
        daskbag = db.from_sequence(datanames,npartitions=2)
        daskbag.map(analysis).compute()
    else:
        # parallel groups
        if multiple:
            analysis(datanames, multiple=multiple)
        else:
            if type(datanames[0])==list:
                select = 3
                for dn in datanames[select]:
                    analysis(dn)
            # single group
            else:
                for dn in datanames:
                    # if dn=='DT008': continue
                    analysis(dn)



if __name__ == '__main__':
    main()
    plt.show()