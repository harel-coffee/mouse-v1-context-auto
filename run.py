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

from aggregateallmice import *

    






def select_data():
    # all behaviourally symmetric V1
    datanames = ['ME110','ME113','DT009','DT014','DT017','DT021','DT022','MT020_2']
    # datanames = [['ME110','ME113','DT009'],['DT014','DT017','DT021'],['DT021','DT022','MT020_2']]




    
    return datanames
    







def analysis(datanames, multiple=False):


    # choose analyis type
    # multiple = True
    
    # aggregating over multiple mice routines
    if multiple:
        
        # full routines:
            
        # aggregatemice_numberofhighlowperformancetrials(datanames)
        # aggregatemice_modelLLs(datanames)

        # aggregatemice_relevantirrelevant_behaviour(datanames)
        # aggregatemice_singletrialneuralbehaviour_relevantirrelevantcongruentconflicting(datanames)
        # aggregatemice_singletrialneuralbehaviour_context(datanames)
        # aggregatemice_subspaces_timeresolved_behaviour_ACC(datanames)




        # cluster_celltypes(datanames)

        # aggregatemice_decoderaccuracyaverages(datanames)
        
        
        # decoderstobehaviouralperformance(datanames);
        # decodersbehaviourwithtraining(datanames);
        
        # aggregatemice_representationmeasures_featuresvector(datanames)
        
        # aggregatemice_subspaceangles(datanames)
        # aggregatemice_spontaneousdecoder(datanames,'allexpcond')
        # aggregatemice_spontaneousdecoder(datanames,'attendignore')
        # aggregatemice_acrosscontextcomparison(datanames)
        # aggregatemice_crosscontextdecoder(datanames)
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
        # aggregatemice_lowdimensionnullspacecontext(datanames)
        # aggregatemice_contextpreonsubpop(datanames)

        # aggregatemice_behavioursymmetrycontext(datanames)


        # aggregatemice_chancelevel(datanames)
        # aggregatemice_chancelevel(datanames, reducedtrials=True)




        
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
            # # preprocess.convert_trialstocsv_stimmatrix_early(dn); return # use this once to get the trials_good_start.csv for the mouse (stimmatrix multi-shift version)
            # preprocess.convert_trialstocsv_stimmatrix(dn); return # use this once to get the trials_good_start.csv for the mouse (stimmatrix multi-shift version)

            # EXPORT spike times in a simple format without trial segmentation
            # preprocess.loaddatamouse(dn,T,continuous_method,recalculate=True,exportspiketrains='save')
            # block = preprocess.loaddatamouse(dn,T,continuous_method,recalculate=True,exportspiketrains='load')


            # LOADS or CREATES spike counts, needed for all processes below
            # normalize True for decoding, False for firing rate and raw PCA
            # recalculate True to create neo format cache data files with smoothed inst.fr, False to reload
            block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=True,recalculate=False)
            
            # EXPORT trial cut spike trains from the converted neo format:
            # preprocess.exportspiketrainsfromneo(dn,block)
            # EXPORT firing rates from neo
            # preprocess.exportfiringratesfromneo(dn,block); return



            

            # preprocess.loadtrainingbehaviouraldata(dn, recalculate=True, plot=True)
            # behaviour_generate_movingaveragestatistics()
            # behaviour_likelihood_simplemodels(dn)
            # behaviour_likelihood_idealobserver(dn)
            # behaviour_likelihood_idealobserver(dn, onlybehaviourma=True)             # just the behaviour plot, no likelihood needed
            # behaviour_likelihood_sigmoidlinearmodels(dn)

            # behaviour_symmetry_highperformance(dn)
            # behaviour_symmetry(dn,block,displaystatsonly=True)
            # behaviour_symmetry(dn,block)
            # behaviour_symmetry_context(dn,block)
            # behaviour_symmetry_context(dn,block,equalize=True,n_bootstrap=20)

            
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
            # angles_between_decoders_motion(dn,block)

            # decoder_crosstest(dn,block)
            # decoder_crosstest_highres(dn,block,method='singlecontext')
            # decoder_crosstest_highres(dn,block,method='allcontexts')
            # decoder_crosstest_highres(dn,block,method='correctnogocontext')
            # decoder_crosstest_highres(dn,block,method='conditionedcontext')
            # decoder_crosstest_highres(dn,block,method='congruency')
            # decoder_crosstest_highres(dn,block,method='sparsecontext')
            # if not dn=='DT014': decoder_crosstest_highres(dn,block,method='sparsecontext-small')
            # decoder_crosstest_highres(dn,block,method='visualnullspacecontext')
            # decoder_crosstest_highres(dn,block,method='recurrentnullspacecontext')
            

            
            # decoder_crosstestgradient(dn,block)
            # decoder_acrosscontextcomparison(dn,block)
            # decoder_crosscontexttest(dn,block)
            # decoder_withincontext_saveprojections(dn,block)

            # for trainidx in [0,1]:
            #     for testidx in [0,1]:
            #         decoder_crosscontexttest(dn,block,trainidx=trainidx,testidx=testidx)
            # decoder_crosscontexttest(dn,block,trainidx=0,testidx=1)
            # decoder_crosscontexttest(dn,block,trainidx=1,testidx=0)
            # decoder_crosscontexttest(dn,block)
            
            # decoder_context_clevercluelessperiods(dn,block)
            


            # decoder_acc_relevantirrelevant(dn,block)
            # decoder_acc_relevantirrelevant(dn,block,area='acc/comparisonV1/')       # compare with V1
            # decoder_acc_singletrialneuralbehaviour_relevantirrelevantcongruentconflicting(dn,block)
            # decoder_acc_singletrialneuralbehaviour_context(dn,block)           # for now this requires ks2ifr!


            
            
            # decoder_behaviour(dn,block)
            # latent_behaviour(dn)
            
            
            
    
            # joint_differentiabilityrunspeed(dn,block)
            
            # runspeed_regression(dn,block)
            # runspeed_taskconditions(dn,block)

            # corr_runspeed_dbnv(dn,block)
            # decoder_equalizedrunspeeddistribution(dn,block)
            
            # display_movementpca(dn)
            # decode_movementpca_singleconcat(dn,block)
            # decode_movementpca(dn,block,calculuslevel='posture')
            # decode_movementpca(dn,block,calculuslevel='motion')
            # decode_movementbodyparts(dn,block)
            # decode_movementpca_tasks(dn)
            # subspace_decode_motiontoneuron_reducedrank(dn, block, calculuslevel='posture')
            # subspace_decode_motiontoneuron_reducedrank(dn, block, calculuslevel='motion')
            # comparestationarycontext(dn, block)
            # comparemovementdistributioncontext(dn, block,roi='total')
            # comparemovementdistributioncontext(dn, block,roi='bodyparts')
            # comparecorrectincongruentgonogocontext(dn,block)


            # subspaces(dn,block,examine='allexpcond')
            # subspaces(dn,block,examine='attendignore')
            
            # subspaces_projections(dn,block,examine='allexpcond',play='visual',preon=0)
            # subspaces_projections(dn,block,examine='allexpcond',play='visual',preon=1)
            # subspaces_projections(dn,block,examine='allexpcond',play='audio',preon=0)
            # subspaces_projections(dn,block,examine='allexpcond',play='audio',preon=1)
            # subspaces_projections(dn,block,examine='allexpcond',play='choice',preon=1)



            # subspaces_behaviour_ACC(dn,block)        # ACC
            # subspaces_behaviour_ACC(dn,block,area='acc/comparisonV1/')       # compare with V1
            # subspaces_timeresolved(dn,block)



            # subspaces_decode_nullspacespanmotion(dn,block)
            # subspaces_decode_nullspace1visual1(dn,block)
            # subspaces_decode_extendedspacevisualcontext(dn,block)
            # subspaces_decode_extendedspacevisualcontext_cvsplit(dn,block)
            # subspaces_decode_contextnullspacevisual(dn,block)
            # subspaces_decode_nullspacerecurrent(dn,block)
            # subspaces_decode_nullspacerecurrent(dn,block,onlyfirst=True)

            # subspaces_orthogonalcontrol_shuffle(dn,block)
            # subspaces_orthogonalcontrol_shuffle(dn,block,reducedtrials=True)
            # subspaces_leaveoneout(dn,block)

            # subspaces_neuralparticipation(dn,block)
            # subspaces_residualvariance(dn,block,examine='allexpcond')
            # subspaces_examples(None,None)
            
            # subspaces_PCA(dn,block,preon=0)
            # subspaces_PCA(dn,block,preon=1)
            
            # subspacedynamics_PCA(dn,block,normalize=True)
            # subspacedynamics_PCA_context(dn,block)
            # subspacedynamics_PCA_signalnoisecorrelation(dn,block)
            
            # subspacedynamics_DVPCA(dn,block)
            


            # subspacedynamics_decisionvector_basis(dn,block)
            # subspacedynamics_decisionvector_basis_initialconditions(dn,block)
            # subspacedynamics_decisionvector_basis_initialconditions_highres(dn,block)




            # variancemethods_dPCA(dn,block)



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
                select = 0
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