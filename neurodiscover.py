# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:00:23 2019

@author: mahajnal
"""

import multiprocessing
from os import cpu_count
n_cpu = cpu_count()

import numpy as np
import scipy as sp
from scipy.stats import ortho_group
import matplotlib.pyplot as plt

import quantities as pq
import neo

import sklearn.dummy as skdu
import sklearn.linear_model as sklm
import sklearn.model_selection as skms
import sklearn.metrics as skme
import sklearn.decomposition as skdc
import sklearn.feature_selection as skfs
import sklearn.neural_network as sknn
import sklearn.gaussian_process as gp
import sklearn.mixture as gm
import sklearn

import tensorly as tl
from tensorly.decomposition import tucker, parafac, non_negative_tucker


import neurophysiology as neph

import neurobayesian as neba

import pickle

import warnings

# calculators


def get_maxnormdifference(rd):
    # input is response variable: time   x  channels    x   stats
    d0 = np.abs(rd[0])
    d2 = np.abs(rd[2])
    dx = np.argmax(d0,axis=1)
    dm = np.array([ d0[k,dx[k]] for k in range(len(dx)) ])
    daux = np.array([ d2[k,dx[k]] for k in range(len(dx)) ])
#        de = np.max((np.zeros(dm.shape),dm-2.*de),axis=0)       # clip below zero
    de = dm-2.*daux                         # use the whole band
    
    dm = neo.AnalogSignal(dm, name='mean difference '+rd[0].name,t_start=rd[0].t_start,sampling_period=rd[0].sampling_period,units=rd[0].units)
    de = neo.AnalogSignal(de, name='s.e.m. of diff. '+rd[0].name,t_start=rd[0].t_start,sampling_period=rd[0].sampling_period,units=rd[0].units)
    
    return dm,de




def residuals(x,y,s_,i_):
    return np.mean((s_*x+i_ - y)**2)




def ll_bernoulli(k, p):
    # k is a vector of 1s and 0s (data)
    # p is a vector of probabilities (estimates from model)
    # return bernoulli log likelihood for each data point
    ll = (     k * np.log(  p  )   +   (1-k) * np.log( (1-p) )     )
    return ll







def partialcorrelation(C):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in C, controlling 
    for the remaining variables in C.
    Parameters
    ----------
    C : array-like, shape (n, p)
        Array with the different variables. Each column of C is taken as a variable
    Returns
    -------
    P : array-like, shape (p, p)
        P[i, j] contains the partial correlation of C[:, i] and C[:, j] controlling
        for the remaining variables in C.

    Partial Correlation in Python (clone of Matlab's partialcorr)
    This uses the linear regression approach to compute the partial 
    correlation (might be slow for a huge number of variables). The 
    algorithm is detailed here:
        http://en.wikipedia.org/wiki/Partial_correlation#Using_linear_regression
    Taking X and Y two variables of interest and Z the matrix with all the variable minus {X, Y},
    the algorithm can be summarized as
        1) perform a normal linear least-squares regression with X as the target and Z as the predictor
        2) calculate the residuals in Step #1
        3) perform a normal linear least-squares regression with Y as the target and Z as the predictor
        4) calculate the residuals in Step #3
        5) calculate the correlation coefficient between the residuals from Steps #2 and #4; 
        The result is the partial correlation between X and Y while controlling for the effect of Z
    Date: Nov 2014
    Author: Fabian Pedregosa-Izquierdo, f@bianp.net
    Testing: Valentina Borghesani, valentinaborghesani@gmail.com
    """
    
    C = np.asarray(C)
    p = C.shape[1]
    n = C.shape[0]
    P_corr = np.zeros((p, p,2), dtype=np.float)
    R_lsq = np.zeros((p,p,2,n), dtype=np.float)
    for i in range(p):
        P_corr[i, i] = 1
        for j in range(i+1, p):
            idx = np.ones(p, dtype=np.bool)
            idx[i] = False
            idx[j] = False
            beta_i = sp.linalg.lstsq(C[:, idx], C[:, j])[0]
            beta_j = sp.linalg.lstsq(C[:, idx], C[:, i])[0]

            res_j = C[:, j] - C[:, idx].dot( beta_i)
            res_i = C[:, i] - C[:, idx].dot(beta_j)
            R_lsq[i,j,0,:] = res_i
            R_lsq[i,j,1,:] = res_j
            
            # corr = sp.stats.pearsonr(res_i, res_j)     # rho and p
            corr = sp.stats.linregress(res_i, res_j)[2:4]     # rho and p
            P_corr[i, j,:] = corr
            P_corr[j, i,:] = corr
        
    return P_corr,R_lsq









    

def innerproductangle(a,b,degree=False):
    # a,b are vectors
    
    gamma =  np.arccos(   np.dot(a,b)  /   (  np.linalg.norm(a)  *  np.linalg.norm(b)  )    )

    # since these are not vectors but generally infinite lines to both directions, 
    # take the continuation and reduce to the angle of these lines (always < 90 degrees)
    if gamma>np.pi/2:
        gamma = np.abs(np.pi-gamma)
    
    if degree:
        gamma *= 180 / np.pi

    return gamma


def rotationmatrix3d(ax=0,ay=0,az=0,angle=True):
    
    if angle:
        ax = ax * np.pi/180
        ay = ay * np.pi/180
        az = az * np.pi/180
    
    Mx = np.array([[1,0,0],[0,np.cos(ax),-np.sin(ax)],[0,np.sin(ax),np.cos(ax)]])
    My = np.array([[np.cos(ay),0,np.sin(ay)],[0,1,0],[-np.sin(ay),0,np.cos(ay)]])
    Mz = np.array([[np.cos(az),-np.sin(az),0],[np.sin(az),np.cos(az),0],[0,0,1]])
    
    M = Mz @ My @ Mx
    
    return M



#def tensorreconstructnorm(U):
#    N = np.zeros(4)
#    X_ = np.zeros([U[lx] for lx,l in range(len(U))])
#    for r in range(U[0].shape[1]):
#        N






# information criteria
    
def datazeroermatrix(nc=4):
    halfnc = 1
    mc = []   # model coefficients are in the columns, each row is a model combination
    # first each variable alone:
    for cx in range(nc):
        mc.append(np.zeros(nc+halfnc,dtype='int16'))
        mc[cx][cx]=1
    # then all four variables, except the one in question (subtracted model)
    for cxi,cx in enumerate(np.arange(nc,nc*2)):
        mc.append(np.ones(nc+halfnc,dtype='int16'))
        mc[cx][cxi]=0
    mc = np.array(mc)
    mc = np.r_[mc,mc]        # make a double list, first without, the second part will be with runspeed
    mc[:2*nc,nc+halfnc-1] = 0
    # make the next: runspeed only
    mc[2*nc:,nc+halfnc-1] = 1
    # now make the next two the full model, with all task variables, without and with runspeed
    mc = np.r_[mc,np.zeros((1,nc+halfnc),dtype='int16'),np.ones((2,nc+halfnc),dtype='int16')]
    mc[-3,nc+halfnc-1] = 1
    mc[-2,nc+halfnc-1] = 0
    # insert a last with only a single runspeed:
    mc = np.r_[mc, np.zeros((1,nc+halfnc),dtype='int16')]
    mc[-1,-1] = 1
    return mc
    

def aic(y,y_,k):
    n = len(y)         # sample size
    rss = np.sum( (y-y_)**2 )          # residual sum of squares (total error; will be averaged with 1/n in the next line)
    aic = 2*k + n * ( np.log(  2*np.pi * rss / n )  + 1)
    # Then the quantity exp((AIC_min − AIC_i)/2) can be interpreted as being proportional to the probability that the ith model minimizes the (estimated) information loss
    return aic



def dprime_normal(x1,x2):
    m = np.zeros(2)
    v = np.zeros(2)
    for k in [0,1]:
        m[k]=np.mean([x1,x2][k])
        v[k]=np.var([x1,x2][k])
    dprime = (m[0]-m[1]) / np.sqrt((v.sum())/2)
    return dprime





# data structure discovery
    

def pca(X,n_components=2,ret='transformed'):
    # X must be n_samples x n_features

    PCAinstance = skdc.PCA(n_components=n_components)
    X_transformed = PCAinstance.fit_transform(X)
    
    S = PCAinstance.singular_values_
    E = PCAinstance.explained_variance_
    V = PCAinstance.components_

    if ret=='singularvalue':
        return S,E,V
    else:
        return X_transformed,E,V






def fa(X,n_components=2):
    # X must be n_samples x n_features

    FactorAnalysisInstance = skdc.FactorAnalysis(n_components=n_components)
    X_transformed = FactorAnalysisInstance.fit_transform(X)
    
    ll = FactorAnalysisInstance.loglike_
    m = FactorAnalysisInstance.mean_
    e = FactorAnalysisInstance.noise_variance_
    transformer = FactorAnalysisInstance.components_

    return X_transformed,transformer,ll,m,e,









def bootstrapequalizesmallest(labels,labellist,indices_f,N_cat):
    # equalize data so that the fewest members is not penalized
    y = labels[indices_f]
    h, be = np.histogram(y, bins=N_cat)
    bootstrap_n_sample = np.min(h)
    indices = []
    for l in labellist:       # l runs through all labels
        matching_indices = np.where(y==l)[0]
        indices.extend( indices_f[matching_indices[:bootstrap_n_sample]] )
    return indices






# random orthomatrices
def get_randomorthocumulativeprojections(d,n):
    # return a list of n orthogonal matrices with dimension d
    Ms = []
    for i in range(n):
        M = ortho_group.rvs(dim=d)
        Ms.append(M)
    return Ms








def rgb2gray(I):
    c = np.array([0.2989, 0.5870, 0.1140])     # convert to grayscale weights of RGB
    G = np.dot(I[:,:,:3],c)
    if I.shape[2]>3:    # use transparency mask the 4th dimension
        G = G*I[:,:,3]
    return G








def pointwiseclassifier(X,y, lastclassone=1, cv=False, rr=10, metric='roc', maxpc=5, solver='liblinear',shuffle=None):
    
    roc_tr = []
    roc_te = []
    a_ns = []       # angle between decision boundary normal and principal components
    c_db = []       # decision boundary vector (decoder coefficients)

    if cv==False:

        for r in range(rr):
            X1 = X[:lastclassone,:]
            X2 = X[lastclassone:,:]
            y1 = y[:lastclassone]
            y2 = y[lastclassone:]
            X1_tr, X1_te, y1_tr, y1_te = skms.train_test_split(X1, y1, test_size=0.25)
            X2_tr, X2_te, y2_tr, y2_te = skms.train_test_split(X2, y2, test_size=0.25)
            X_tr = np.vstack([X1_tr, X2_tr])
            X_te = np.vstack([X1_te, X2_te])
            y_tr = np.hstack([y1_tr, y2_tr])
            y_te = np.hstack([y1_te, y2_te])
            
            dec = sklm.LogisticRegression(solver=solver,penalty='l2')
            
            if not shuffle==None:     # partially randomize the indices (partial shuffle), given by the percentage/proportion shuffle.
                sh_N = len(y_tr)
                print('shuffle train %4.2f, %d/%d'%(shuffle,shuffle*sh_N,sh_N))
                indices = np.random.permutation(sh_N)[ :int(shuffle*sh_N) ]
                y_tr[indices] = round(np.random.rand())      # 1-y_tr[indices], y_tr[indices] = y_tr[indices]
                sh_N = len(y_te)
                print('shuffle  test %4.2f, %d/%d'%(shuffle,shuffle*sh_N,sh_N))
                indices = np.random.permutation(sh_N)[ :int(shuffle*sh_N) ]
                y_te[indices] = round(np.random.rand())        # 1-y_te[indices]
                
            
            dec.fit(X_tr, y_tr)

            y_tr_ = dec.predict(X_tr)
            y_te_ = dec.predict(X_te)
            y1_te_ = dec.predict(X1_te)
            y2_te_ = dec.predict(X2_te)


            if metric=='roc':
                roc_tr.append(skme.roc_auc_score(y_tr,y_tr_))
                roc_te.append(skme.roc_auc_score(y_te,y_te_))
            elif metric=='classacc':  # this case no training, first class and second class
                roc_tr.append(skme.accuracy_score(y1_te,y1_te_))
                roc_te.append(skme.accuracy_score(y2_te,y2_te_))
                
                
            # get decision boundary vector coordinates in vectorspace basis of SUA activities
            c_db.append(dec.coef_)
    
            # get pca components to decision normal angles        
            decisionnormal = (dec.coef_/np.linalg.norm(dec.coef_))[0]
            U,S,V = np.linalg.svd(X)
            angle_normaltosvd = []
    
            for comp in range(min(maxpc,V.shape[0])):
                angle_normaltosvd.append(innerproductangle( V[comp,:], decisionnormal))
            a_ns.append(angle_normaltosvd)
            
        n_folds = rr
        
    elif cv==True:
        if type(rr)==int:
            rr = np.min(  ((y==0).sum(), (y==1).sum(), rr) )
        
#        print((y==0).sum(),(y==1).sum(),rr)



    
        cv_results = skms.cross_validate(sklm.LogisticRegression(solver=solver,penalty='l2'), \
                                    X, y, cv=rr,return_train_score=True, return_estimator=True)
        roc_tr = cv_results['train_score']
        roc_te = cv_results['test_score']
        
        n_folds = len(cv_results['estimator'])
        for r in range(n_folds):
            c_db.append(cv_results['estimator'][r].coef_)

            decisionnormal = (cv_results['estimator'][r].coef_/  \
                               np.linalg.norm(cv_results['estimator'][r].coef_))[0]
            U,S,V = np.linalg.svd(X)
            angle_normaltosvd = []
            for comp in range(min(maxpc,V.shape[0])):
                angle_normaltosvd.append(innerproductangle( V[comp,:], decisionnormal))
            a_ns.append(angle_normaltosvd)
    


    roc_dist = np.array([roc_tr, roc_te])        # te,tr    x     1,2,3,4... reruns
    roc = np.stack(   (np.mean(roc_dist,axis=1),2*np.std(roc_dist,axis=1),2*np.std(roc_dist,axis=1)/np.sqrt(n_folds)), axis=1)

    if n_folds==1: return roc,dec                    # return for subspaces_examples with n_folds==1


    a_ns = np.array(a_ns).T                 # pca 1, pca2...     x  1,2,3,4.... reruns
    a = np.stack(   (np.mean(a_ns,axis=1),2*np.std(a_ns,axis=1),2*np.std(a_ns,axis=1)/np.sqrt(n_folds)),    axis=1)
    
    c_db = np.array(c_db).squeeze().swapaxes(0,1)      #   1     x     1,2,3,4... reruns      x      dimensions
    coeffs = np.stack(   [np.mean(c_db,axis=1),2*np.std(c_db,axis=1),2*np.std(c_db,axis=1)/np.sqrt(n_folds)],    axis=1)
    
    # now combine train, test and pc angles
    pack = np.vstack((roc,a,coeffs))
    
    
    return pack






def pointwiseclassifierreturnestimator(X,y,n_cvfolds=10,cvtype=skms.StratifiedKFold,solver='liblinear'):
    cv_results = skms.cross_validate(sklm.LogisticRegression(solver=solver,penalty='l2'), \
                                X, y, cv=n_cvfolds,return_train_score=True, return_estimator=True)
    roc_tr = [ cv_results['train_score'].mean(), cv_results['train_score'].std()/np.sqrt(n_cvfolds) ]
    roc_te = [ cv_results['test_score'].mean(), cv_results['test_score'].std()/np.sqrt(n_cvfolds) ]

    # collect single trial test probabilities
    proba_te = np.zeros((len(y),2))
    S = skms.StratifiedKFold(n_splits=n_cvfolds)
    for ix, (i_tr, i_te) in enumerate(S.split(X,y)):
        proba_te[i_te,:] = cv_results['estimator'][ix].predict_proba(X[i_te,:])
        # print(X[i_te,:].shape,p.shape, p[0], p[0]>0.5, y[i_te[0]])



    return roc_tr, roc_te, proba_te, cv_results['estimator']





def pointwiseclassifiercrosstest(X,y,Xcross,ycross,n_cvfolds=5,cvtype=skms.StratifiedKFold,solver='liblinear'):
    cv_results = skms.cross_validate(sklm.LogisticRegression(solver=solver,penalty='l2'), \
                                X, y, cv=n_cvfolds,return_train_score=True, return_estimator=True)
    roc_tr = [ cv_results['train_score'].mean(), cv_results['train_score'].std()/np.sqrt(n_cvfolds) ]
    roc_te = [ cv_results['test_score'].mean(), cv_results['test_score'].std()/np.sqrt(n_cvfolds) ]

    # collect cross test predictions
    roc_te_cross = np.zeros(len(ycross))
    S = skms.StratifiedKFold(n_splits=n_cvfolds)
    for ix, (i_tr, i_te) in enumerate(S.split(Xcross,ycross)):
        roc_te_cross[i_te,:] = cv_results['estimator'][ix].predict(Xcross[i_te,:])

    return roc_tr, roc_te, roc_te_cross, cv_results['estimator']








def pointwiseclassifierproba(X,y,lastclassone=1,rr=5,cv=skms.StratifiedKFold):
    proba = skms.cross_val_predict(sklm.LogisticRegression(solver='liblinear'),X,y,cv=cv(n_splits=rr),method='predict_proba')

    return proba
    







def crossblock_pointwiseclassifier(X_t,X_p,y_t,y_p,cv=None,rr=10,maxpc=5):
    
    roc_t_tr = []
    roc_t_te = []
    roc_p_te = []
    a_ns = []       # angle between decision boundary normal and principal components
    c_db = []       # decision boundary vector (decoder coefficients)

    for r in range(rr):
        X_t_tr, X_t_te, y_t_tr, y_t_te = skms.train_test_split(X_t, y_t, stratify=y_t, test_size=0.25)
        X_p_tr, X_p_te, y_p_tr, y_p_te = skms.train_test_split(X_p, y_p, stratify=y_p, test_size=0.25)
        
        if cv==None:
            dec = sklm.LogisticRegression(solver='liblinear',penalty='l2')
        else:
            dec = sklm.LogisticRegressionCV(cv=cv,solver='liblinear',penalty='l2')
        
        dec.fit(X_t_tr, y_t_tr)
    
        y_t_tr_ = dec.predict(X_t_tr)
        y_t_te_ = dec.predict(X_t_te)
        y_p_te_ = dec.predict(X_p_te)
    
        roc_t_tr.append(skme.roc_auc_score(y_t_tr,y_t_tr_))
        roc_t_te.append(skme.roc_auc_score(y_t_te,y_t_te_))
        roc_p_te.append(skme.roc_auc_score(y_p_te,y_p_te_))
            
            
        # get decision boundary vector coordinates in vectorspace basis of SUA activities
        c_db.append(dec.coef_)

        # get pca components to decision normal angles        
        decisionnormal = (dec.coef_/np.linalg.norm(dec.coef_))[0]
        U,S,V = np.linalg.svd(X_t)
        angle_normaltosvd = []

        for comp in range(min(maxpc,V.shape[0])):
            angle_normaltosvd.append(innerproductangle( V[comp,:], decisionnormal))
        a_ns.append(angle_normaltosvd)
        

    roc_dist = np.array([roc_t_tr, roc_t_te, roc_p_te])        # te,tr    x     1,2,3,4... reruns
    roc = np.stack(   (np.mean(roc_dist,axis=1),2*np.std(roc_dist,axis=1),2*np.std(roc_dist,axis=1)/np.sqrt(rr)), axis=1)

    # beware, that for crossblock, pca and coefficients start from index #3 = the 4th place
    a_ns = np.array(a_ns).T                 # pca 1, pca2...     x  1,2,3,4.... reruns
    a = np.stack(   (np.mean(a_ns,axis=1),2*np.std(a_ns,axis=1),2*np.std(a_ns,axis=1)/np.sqrt(rr)),    axis=1)
    
    c_db = np.array(c_db).squeeze().swapaxes(0,1)      #   1     x     1,2,3,4... reruns      x      dimensions
    coeffs = np.stack(   [np.mean(c_db,axis=1),2*np.std(c_db,axis=1),2*np.std(c_db,axis=1)/np.sqrt(rr)],    axis=1)
    
    # now combine train, test and pc angles
    pack = np.vstack((roc,a,coeffs))
    
    
    return pack






    






def offclassifier(X_train,X_test,y_all, lastclassone=1, cv=None, rr=2):
    # the purpose of this method is to separate the timecourse into 3 sections:
    # 1st stimulus condition: train and validation
    # 2nd stimulus condition: test
    # so here test and train is not just split between trials, but over time bunches
    roc_tr = []
    roc_te = []
    
    for r in range(rr):
#        try:
            if cv==None:
                dec = sklm.LogisticRegression(solver='liblinear',penalty='l2')
            else:
                dec = sklm.LogisticRegressionCV(cv=cv,solver='liblinear',penalty='l2')
            

            X1 = X_train[:lastclassone,:]
            X2 = X_train[lastclassone:,:]
            X1_test = X_test[:lastclassone,:]
            X2_test = X_test[lastclassone:,:]
            y1 = y_all[:lastclassone]
            y2 = y_all[lastclassone:]
            X1_tr, _, y1_tr, _ = skms.train_test_split(X1, y1, test_size=0.25)
            X2_tr, _, y2_tr, _ = skms.train_test_split(X2, y2, test_size=0.25)
            X1_te, _, y1_te, _ = skms.train_test_split(X1_test, y1, test_size=0.25)
            X2_te, _, y2_te, _ = skms.train_test_split(X2_test, y2, test_size=0.25)
            X_tr = np.vstack([X1_tr, X2_tr])
            X_te = np.vstack([X1_te, X2_te])
            y_tr = np.hstack([y1_tr, y2_tr])
            y_te = np.hstack([y1_te, y2_te])



            dec.fit(X_tr, y_tr)
        
            y_tr_ = dec.predict(X_tr)
            y_te_ = dec.predict(X_te)
        
            roc_tr.append(skme.roc_auc_score(y_tr,y_tr_))
            roc_te.append(skme.roc_auc_score(y_te,y_te_))
#        except ValueError:
#            print('resample')
            
        
    roc_dist = np.array([roc_tr, roc_te])        # te,tr    x     1,2,3,4... reruns
    roc = np.stack(   (np.mean(roc_dist,axis=1),2*np.std(roc_dist,axis=1),2*np.std(roc_dist,axis=1)/np.sqrt(rr)), axis=1)

    return roc











def offclassifier_multitest(X_train,X_test_multi,y_all, lastclassone=1, cv=None, rr=5):
    # the purpose of this method is to separate the timecourse into 3 sections:
    # 1st stimulus condition: train and validation
    # 2nd stimulus condition: test
    # so here test and train is not just split between trials, but over time bunches
    # this version allows for a single trained decoder to test several Xs
    roc_tr = []
    roc_te_all = []


    
    for r in range(rr):
#        try:
            if cv==None:
                dec = sklm.LogisticRegression(solver='liblinear',penalty='l2')
            else:
                dec = sklm.LogisticRegressionCV(cv=cv,solver='liblinear',penalty='l2')
            

            X1 = X_train[:lastclassone,:]
            X2 = X_train[lastclassone:,:]
            y1 = y_all[:lastclassone]
            y2 = y_all[lastclassone:]
            X1_tr, _, y1_tr, _ = skms.train_test_split(X1, y1, test_size=0.25)
            X2_tr, _, y2_tr, _ = skms.train_test_split(X2, y2, test_size=0.25)
            X_tr = np.vstack([X1_tr, X2_tr])
            y_tr = np.hstack([y1_tr, y2_tr])



            dec.fit(X_tr, y_tr)
        
            y_tr_ = dec.predict(X_tr)

            roc_tr.append(skme.roc_auc_score(y_tr,y_tr_))

            
            roc_te = []     # hold test accuracies for all test timepoints in a list
            for tx,X_test in enumerate(X_test_multi):

                X1_test = X_test[:lastclassone,:]
                X2_test = X_test[lastclassone:,:]
                X1_te, _, y1_te, _ = skms.train_test_split(X1_test, y1, test_size=0.25)
                X2_te, _, y2_te, _ = skms.train_test_split(X2_test, y2, test_size=0.25)
                X_te = np.vstack([X1_te, X2_te])
                y_te = np.hstack([y1_te, y2_te])


                y_te_ = dec.predict(X_te)
                roc_te.append( skme.roc_auc_score(y_te,y_te_) )


            roc_te_all.append( roc_te )      # this will be   (cv)(testtimepoints)



#        except ValueError:
#            print('resample')
            

    # the followings are not needed to separate train and test, as we do not have train for different timepoints        
    # roc_dist = np.array([roc_tr, roc_te])        # te,tr    x     1,2,3,4... reruns
    # roc = np.stack(   (np.mean(roc_dist,axis=1),2*np.std(roc_dist,axis=1),2*np.std(roc_dist,axis=1)/np.sqrt(rr)), axis=1)


    roc_dist = np.array(roc_te_all)
    roc = np.stack(   (np.mean(roc_dist,axis=0),2*np.std(roc_dist,axis=0),2*np.std(roc_dist,axis=0)/np.sqrt(rr)), axis=1)

    # roc is  (timepoints)(mean,std,2sem)

    return roc











def offgradientclassifier(X_tr,X_te,y_tr,y_te,cv=None,rr=2):
    # takes the training set from the dual modality
    # tests alongside all single modalities (or whatever, given in te_lists)
    # returns training_roc, test_roc_list

    pry_tr = []
    pry_te = []

    for r in range(rr):
        if cv==None:
            dec = sklm.LogisticRegression(solver='liblinear',penalty='l2')
        else:
            dec = sklm.LogisticRegressionCV(cv=cv,solver='liblinear',penalty='l2')

#        if cv==None:
#            dec = sklm.LogisticRegression(solver='liblinear',penalty='l2')
#        else:
#            dec = sklm.LogisticRegressionCV(cv=cv,solver='liblinear',penalty='l2')

        dec.fit(X_tr, y_tr)
    
        proba_y_tr_ = dec.predict_proba(X_tr)
        pry_tr.append(proba_y_tr_[:,y_tr].mean())
        
        proba_y_te_ = dec.predict_proba(X_te)
        pry_te.append(proba_y_te_[:,int(y_te)].mean())
    
    pry_dist = np.array([pry_tr, pry_te])        # te,tr    x     1,2,3,4... reruns
    pry = np.stack(   (np.mean(pry_dist,axis=1),2*np.std(pry_dist,axis=1),2*np.std(pry_dist,axis=1)/np.sqrt(rr)), axis=1)

    return pry






def multiclassifier(responses,labels,T,rr=5,dummy=False,uniform=False,bootstrap=True):
    
    N = len(responses)            # number of categories
    N_cat = len(np.unique(labels))
    # construct train and test set:
    
    solver = 'newton-cg'      # this is needed to ensure convergence in multinomial classifier fit
    
#    for n in range(N):
    f1 = np.zeros((rr,N_cat))
    for r in range(rr):
        X_tr, X_te, y_tr, y_te = skms.train_test_split(responses, labels, stratify=labels, test_size=0.25)

        # reduce the samples into equal sizes independent on class type
        if bootstrap:    # use only if shuffled
            ix = 0
            for y,X in zip([y_tr,y_te],[X_tr,X_te]):
                # find the maximum number of samples to bootstrap for
                labellist = np.unique(y)       # this should be sorted, check if issues
                h, be = np.histogram(y, bins=len(labellist))
                bootstrap_n_sample = np.min(h)
                indices = []
                for l in labellist:       # l runs through all labels
                    matching_indices = np.where(y==l)[0]
                    indices.extend( matching_indices[:bootstrap_n_sample] )
                if ix==0:
                    X_tr = X[indices,:]
                    y_tr = y[indices]
                else:
                    X_te = X[indices,:]
                    y_te = y[indices]
                ix += 1
        
    
    
        if uniform:
            if dummy:
                dec = skdu.DummyClassifier(strategy='uniform')     #stratified
            else: # class_weight='balanced'
                dec = sklm.LogisticRegression(multi_class='multinomial',class_weight=None,solver=solver,penalty='l2')
        else:
            if dummy:
                dec = skdu.DummyClassifier(strategy='stratified')     #
            else: # class_weight='balanced'
                dec = sklm.LogisticRegression(multi_class='multinomial',class_weight='balanced',solver=solver,penalty='l2')
        
        with warnings.catch_warnings():
            try:
                dec.fit(X_tr, y_tr)
            except Warning:
                r = r-1
                continue
        
        y_tr_ = dec.predict(X_tr)
        y_te_ = dec.predict(X_te)
    
#        print(skme.classification_report(y_te, y_te_))
    
#        cm_tr = skme.multilabel_confusion_matrix(y_tr,y_tr_)
#        cm_te = skme.multilabel_confusion_matrix(y_te,y_te_)
    
        _,_,f1_tr,_ = skme.precision_recall_fscore_support(y_tr,y_tr_)
        _,_,f1_te,_ = skme.precision_recall_fscore_support(y_te,y_te_)
        
        f1[r,:] = f1_te
    
    f1_mean = np.nanmean(f1,axis=0)
    f1_sem = np.nanstd(f1,axis=0)/np.sqrt(rr)


    return f1_mean,f1_sem









def multiclassifierproba(responses,labels,T,rr=5,dummy=False,uniform=False):
    N = len(responses)            
    labellist = np.unique(labels)
    N_cat = len(labellist)       # number of categories
    # construct train and test set:
    
    solver = 'newton-cg'      # this is needed to ensure convergence in multinomial classifier fit
    
    
    
    
    cv = skms.StratifiedShuffleSplit(n_splits=3)       # =rr
    dec = sklm.LogisticRegression(multi_class='multinomial',class_weight=None,solver=solver,penalty='l2')
    
    
    for sx, (i_tr_f, i_te_f) in enumerate(cv.split(responses, labels)):
       
        i_tr = bootstrapequalizesmallest(labels,labellist,i_tr_f,N_cat)
        i_te = bootstrapequalizesmallest(labels,labellist,i_te_f,N_cat)


    print(labels[i_tr])
    print(labels[i_te])
    dec.fit(responses[i_tr], labels[i_tr])
    y_tr_ = dec.predict_proba(responses[i_tr])
    y_te_ = dec.predict_proba(responses[i_te])
    




    cross_val_proba


    return

#    for n in range(N):
    f1 = np.zeros((rr,N_cat))

    for r in range(rr):
        X_tr, X_te, y_tr, y_te = skms.train_test_split(responses, labels, stratify=labels, test_size=0.25)

        # reduce the samples into equal sizes independent on class type
        ix = 0
        for y,X in zip([y_tr,y_te],[X_tr,X_te]):
            # find the maximum number of samples to bootstrap for
            labellist = np.unique(y)       # this should be sorted, check if issues
            h, be = np.histogram(y, bins=len(labellist))
            bootstrap_n_sample = np.min(h)
            indices = []
            for l in labellist:       # l runs through all labels
                matching_indices = np.where(y==l)[0]
                indices.extend( matching_indices[:bootstrap_n_sample] )
            if ix==0:
                X_tr = X[indices,:]
                y_tr = y[indices]
            else:
                X_te = X[indices,:]
                y_te = y[indices]
            ix += 1
        
    
    
        if uniform:
            if dummy:
                dec = skdu.DummyClassifier(strategy='uniform')     #stratified
            else: # class_weight='balanced'
                dec = sklm.LogisticRegression(multi_class='multinomial',class_weight=None,solver=solver,penalty='l2')
        else:
            if dummy:
                dec = skdu.DummyClassifier(strategy='stratified')     #
            else: # class_weight='balanced'
                dec = sklm.LogisticRegression(multi_class='multinomial',class_weight='balanced',solver=solver,penalty='l2')
        
        with warnings.catch_warnings():
            try:
                dec.fit(X_tr, y_tr)
            except Warning:
                r = r-1
                continue
        
        y_tr_ = dec.predict(X_tr)
        y_te_ = dec.predict(X_te)
    
        _,_,f1_tr,_ = skme.precision_recall_fscore_support(y_tr,y_tr_)
        _,_,f1_te,_ = skme.precision_recall_fscore_support(y_te,y_te_)
        
        f1[r,:] = f1_te
    
    f1_mean = np.nanmean(f1,axis=0)
    f1_sem = np.nanstd(f1,axis=0)/np.sqrt(rr)


    return f1_mean,f1_sem











    









def pointwiselinearregression(X,y,wx=1,rr=5,maxpc=5):
#    print('pointwise linear regression',X.shape,y.shape)
    
    n_trials = X.shape[0]
    
    cv_results = skms.cross_validate(sklm.Ridge(), \
                                X, y, cv=rr,return_train_score=True, return_estimator=True)
    sc_tr = cv_results['train_score']
    sc_te = cv_results['test_score']

    E = cv_results['estimator']
#    y_tr_ = [ E[i].get_params()  ]
#    print(y_tr_)
#    .reshape(-1,wx)

#    y_tr_ = skms.cross_val_predict(sklm.Ridge(), \
#                                X, y, cv=rr)

    splits = []            # this is needed for prediction output
    for tr_ix,te_ix in skms.KFold(rr).split(X,y):
        splits.append([tr_ix,te_ix])
    

    a_ns = []
    coeffs = []
    intercepts = []
    
#    R = np.zeros((len(y),1,2,2))         # trials, trajectories placeholder, train-test, mean+std
    
    y_tr_ = np.zeros((n_trials,rr))
    y_te_ = np.zeros((n_trials,rr))
    
    for r in range(rr):
        coeffs.append(cv_results['estimator'][r].coef_)
        
        U,S,V = np.linalg.svd(X)

        angle_weightstosvd = []
        for comp in range(min(maxpc,V.shape[0])):
            angle_weightstosvd.append(innerproductangle( V[comp,:], coeffs[r].squeeze()))
        a_ns.append(angle_weightstosvd)

        # predict y:
        y_tr_[splits[r][0],r] += cv_results['estimator'][r].predict(X[splits[r][0],:]) / (rr-1)
        y_te_[splits[r][1],r]                         = cv_results['estimator'][r].predict(X[splits[r][1],:])
        

    coeffs = np.array(coeffs).swapaxes(0,1)      #   1     x     1,2,3,4... reruns      x      dimensions
    
    c = np.stack(   [np.mean(coeffs,axis=1),2*np.std(coeffs,axis=1),2*np.std(coeffs,axis=1)/np.sqrt(rr)],    axis=1)
    
    a_ns = np.array(a_ns)#.squeeze()[:,np.newaxis]                 # pca 1, pca2...     x  1,2,3,4.... reruns
    a = np.stack(   (np.mean(a_ns,axis=0),2*np.std(a_ns,axis=0),2*np.std(a_ns,axis=0)/np.sqrt(rr)),    axis=1)

    sc_dist = np.array([sc_tr, sc_te])        # te,tr    x     1,2,3,4... reruns
    sc = np.stack(   (np.mean(sc_dist,axis=1),2*np.std(sc_dist,axis=1),2*np.std(sc_dist,axis=1)/np.sqrt(rr)), axis=1)

#    print(a_ns.shape,'\n',sc.shape,a.shape,c.shape)

    # now combine train, test and pc angles
    pack = np.vstack((sc,a,c))
#    print(X.shape,y.shape)
#    print('pack',pack.shape,sc.shape,a.shape,c.shape)

    
    
    
    #run is: (trials, trajectory, traintest, stats)
    run = np.array([   [y_tr_.mean(axis=1),y_tr_.std(axis=1)],   [y_te_.mean(axis=1),y_te_.std(axis=1)]   ])
    
#    print('runshape',run.shape)

    return pack, run
    










# imagerecognition models
    

def imagetracedecoder(X,y):
    
    X_tr, X_te, y_tr, y_te = skms.train_test_split(X, y, test_size=0.25)
    
    # print('imagetracedecoder')
    # print(X_tr.shape,X_te.shape,y_tr.shape,y_te.shape)
    
    # F = sknn.MLPRegressor()

    n_features = y_tr.shape[1]
    y_tr_ = np.zeros(y_tr.shape)
    y_te_ = np.zeros(y_te.shape)
    
    
    # train to predict each feature weights of transformed image from neural activity
    for n in range(n_features):
        F = sklm.Ridge()
        
        F = F.fit(X_tr, y_tr[:,n])
    
        y_tr_[:,n] = F.predict(X_tr)
        y_te_[:,n] = F.predict(X_te)
        # y_tr_ = F.predict(np.random.randn(1,X_tr.shape[1])).squeeze()
        # y_tr_ = F.predict(X_tr)[np.random.randint(X_tr.shape[0]),:]
        # y_te_ = F.predict(X_te)[np.random.randint(X_te.shape[0]),:]
        
        # print(y_tr_.shape, y_te_.shape)
        
        # and the error is training error and testing error

    # calculate scores and total errors
    s_tr_ = np.array( np.mean( (y_tr-y_tr_)**2,axis=0 ) )
    s_te_ = np.array( np.mean( (y_tr-y_tr_)**2,axis=0 ) )
    e = np.array( [ np.mean( (y_tr-y_tr_)**2 ), np.mean( (y_te-y_te_)**2 ) ] )
    
#    return F.coefs_[0].mean(axis=1), y_te_, e
    return s_tr_, s_te_, e







def glm_spikesfromimagefeatures(X,Y):
    # X is image features, Y is neurons
    # number of trials is shape[0]
    
    n_neurons = Y.shape[1]
    n_features = X.shape[1]
    
    X_tr, X_te, Y_tr, Y_te = skms.train_test_split(X, Y, test_size=0.25)

    Y_tr_ = np.zeros(Y_tr.shape)
    Y_te_ = np.zeros(Y_te.shape)
    mll = np.zeros(n_neurons)

    for r in range(rr):
        glms = []
        for n in range(n_neurons):
            glms.append( sklm.BayesianRidge(compute_score=True) )        # return the MLL as well
            glms[n].fit(X_tr,Y_tr[:,n])
        
        
            Y_tr_[:,n] = ( glms[n].predict(X_tr) )
            Y_te_[:,n] = ( glms[n].predict(X_te) )
            # W[r,n,:] = glms[n].coef_
            mll[n] = glms[n].scores_[-1]          # get the marginal log likelihood of the data from the last iteration
        
        
    # fig,ax = plt.subplots(1,2,figsize=(24,16))
    # if 1:
    #     ax[0].plot(idxs_tr,Y_tr,'o',color='navy',alpha=0.2)
    #     ax[0].plot(idxs_tr,Y_tr_,'o',color='firebrick',alpha=0.2)
    #     ax[0].plot(idxs_te,Y_te,'o',color='dodgerblue',alpha=0.2)
    #     ax[0].plot(idxs_te,Y_te_,'o',color='fuchsia',alpha=0.2)
        
    #     ax[1].plot(idxs_tr,Y_tr_-Y_tr,'o',color='darkred',alpha=0.2)
    #     ax[1].plot(idxs_te,Y_te_-Y_te,'o',color='red',alpha=0.2)
            
            
    #        print(Y_te.shape,Y_te_.shape)        
        R2[rr] = skme.r2_score(Y_te.ravel(),Y_te_.ravel())
        
        A[rr] = aic(Y_te,Y_te_,n_features)
        
        M[rr] = mll

    
    if len(R2[R2>0])>0:
        R2_mean = np.nanmean(R2[R2>0],axis=0)
        R2_sem = np.nanstd(R2[R2>0],axis=0)/np.sqrt(np.sum(R2>0))
    else:
        R2_mean = 0
        R2_sem = np.inf
    

    W_mean = np.mean(W,axis=0)
    W_sem = np.std(W,axis=0)/np.sqrt(rr)
    
    AIC_mean = np.mean(A,axis=0)
    AIC_sem = np.std(A,axis=0)/np.sqrt(rr)
    
    L_mean = np.mean(L,axis=0)
    L_sem = np.std(L,axis=0)/np.sqrt(rr)
    
    # return M,R2,A,Y_tr,Y_tr_
    return W_mean, W_sem, R2_mean, R2_sem, AIC_mean, AIC_sem, L_mean, L_sem
    
    
    
    





# organizers routines



def get_mapping_imagefeaturesfromspikes(responses,image,width=50*pq.ms):
    # train a decoder: predictors are the neural data in responses. Dependent is the image feature map

    times = responses[0].times
    n_trials = len(responses)
    
    responses_concat = np.array(responses)
    print('neural response',responses_concat.shape)
#    responses_concat = np.random.randn(responses_concat.shape[0],responses_concat.shape[1],responses_concat.shape[2])
    
    # prepare the feature image to be present at each trial
    target = np.repeat(image[np.newaxis,:],n_trials,axis=0)
    print('target',target.shape)
    # try other inputs for test purposes
    # target = (0.2*np.random.randn(n_trials,len(image))[:,np.newaxis]+image[np.newaxis,:]).squeeze()
    # print('target',target.shape)

    wx = int(width/(times[1]-times[0]))
    
    Ir_tr_ = np.zeros( (len(times),len(image)) )     # training reconstruction
    Ir_te_ = np.zeros( (len(times),len(image)) )     # testing reconstruction
    E = np.zeros( (len(times),2) )                       # train and test error
    
    for tx,t in enumerate(times[:-wx]):
        if tx<40 or tx>111: continue              # skip non-stimulus time for now
        predictors = responses_concat[:,tx:tx+wx,:]
        predictors = predictors.reshape(predictors.shape[0],predictors.shape[1]*predictors.shape[2])
    
        Ir_tr_[tx,:], Ir_te_[tx,:], E[tx,:] = imagetracedecoder(predictors,target)
        
   
    return Ir_tr_, Ir_te_, E









def get_mapping_spikesfromimagefeatures(imagefeatures,neuralresponses, width=50*pq.ms):
    
    
    times = neuralresponses[0].times
    n_trials = len(neuralresponses)
    n_neurons = neuralresponses[0].shape[1]
    
    responses_concat = np.array(neuralresponses)
    print('neural response',responses_concat.shape)
#    responses_concat = np.random.randn(responses_concat.shape[0],responses_concat.shape[1],responses_concat.shape[2])
    
    # prepare the feature image to be present at each trial
    predictors = np.repeat(imagefeatures[np.newaxis,:],n_trials,axis=0)
    print('features as predictors',predictors.shape)
    

    wx = int(width/(times[1]-times[0]))
    
    # Y_tr = np.zeros( (len(times),n_neurons) )     # training data
    # Y_te = np.zeros( (len(times),n_neurons) )     # testing data
    # Y_tr_ = np.zeros( (len(times),n_neurons) )     # training reconstruction
    # Y_te_ = np.zeros( (len(times),n_neurons) )     # testing reconstruction
    R2 = np.zeros( len(times) )     # scores for overall all neurons
    A = np.zeros( len(times) )                       # model MLL of the data
    M = np.zeros( (len(times), n_neurons * wx) )                       # model AIC
    
    for tx,t in enumerate(times[:-wx]):
        # if tx<40 or tx>111: continue              # skip non-stimulus time for now
        targets = responses_concat[:,tx:tx+wx,:]
        targets = targets.reshape(targets.shape[0],targets.shape[1]*targets.shape[2])
    
        # Y_te[tx,:,:],Y_te_[tx,:,:],Y_tr[tx,:,:],Y_tr_[tx,:,:],R2[tx,:],A[tx,:] =\
        M[tx,:],R2[tx],A[tx],_,_ = glm_spikesfromimagefeatures(predictors,targets)
        
        
   
    return M,R2,A









def removechannelsfromresponses(acrossresponse,channels_to_remove):
    for lx,triallist in enumerate(acrossresponse):
        for tx,trial in enumerate(triallist):
            acrossresponse[lx][tx] = neph.remochannels(acrossresponse[lx][tx],channels_to_remove)





def get_responsedifference(across_responses):
    responsedifference = []
    responsedifference.append(neph.trialaveraging(across_responses[0])[0] - neph.trialaveraging(across_responses[1])[0])
    responsedifference.append(neph.trialaveraging(across_responses[0])[1] + neph.trialaveraging(across_responses[1])[1])
    responsedifference.append(neph.trialaveraging(across_responses[0])[2] + neph.trialaveraging(across_responses[1])[2])
#    responsedifference = np.stack(responsedifference,axis=2)
    return responsedifference




def get_responsedecoding(across_responses,width=50*pq.ms,cv=False,rr=10,metric='roc',solver='liblinear',shuffle=None):
    if rr=='loo': rr = skms.LeaveOneOut()
    responses = np.concatenate(across_responses,axis=0)
    times = across_responses[0][0].times
    print('time',times[1]-times[0])
    wx = int(width/(times[1]-times[0]))
    print(responses.shape,times[-1])
    print('decoder width',wx,width)
    print('class numbers:',len(across_responses[0]),len(across_responses[1]))
#    responses = responses[:,T['epochstart_idx']:T['epochend_idx']+1,:]
    targets = np.hstack((np.zeros(len(across_responses[0]),dtype='int16'),\
                         np.ones(len(across_responses[1]),dtype='int16')))

    if not shuffle==None:     # partially randomize the indices (partial shuffle), given by the percentage/proportion from argument shuffle.
        if shuffle>0:
            sh_N = len(targets)
            num_shuffle = int(shuffle*sh_N)
            print('shuffle %4.2f, %d/%d'%(shuffle,num_shuffle,sh_N))
            indices = np.random.permutation(sh_N)[ :num_shuffle ]
            # y[indices] = round(np.random.rand())        # this method will give uncontrolled number in each class
            # flipping method is better
            # ix_y_swap_0 = np.where([y[indices]==0])[0]
            # ix_y_swap_1 = np.where([y[indices]==1])[0]
            # y_save = y.copy()
            # y[indices][ix_y_swap_0] = np.round(np.random.rand(len(ix_y_swap_0)))
            # y[indices][ix_y_swap_1] = np.round(np.random.rand(len(ix_y_swap_1)))
            targets[indices] = np.round(np.random.rand(num_shuffle)).astype(np.int16)

            # print(indices)
            # print(ix_y_swap_0,ix_y_swap_1)
            # print(y_save)
            # print(y)
            # print(y_save-y)

            # print(len(indices), np.sum(np.abs(y_save-y)))


    roc=[]
    for tx,t in enumerate(times[:-wx]):
        predictors = responses[:,tx:tx+wx,:]
        predictors = predictors.reshape(predictors.shape[0],predictors.shape[1]*predictors.shape[2])
    
        if 1:   #use t==0 or an interval  to speed up computations during debug
            roc.append( pointwiseclassifier(predictors,targets,len(across_responses[0]),cv=cv,rr=rr,metric=metric,solver=solver,shuffle=shuffle) )
        else:
            roc.append(roc[tx-1])        # this is for speeding up computations for debug
            

    aux = np.array(np.swapaxes(roc,0,1))
    sig_t_start = [ trial.t_start for a_c in across_responses for trial in a_c ]
    print(len(sig_t_start))
    sig_t_start = np.mean(sig_t_start)*times.units
    print(sig_t_start)

    rocsignals = []
    for s in range(aux.shape[0]):
        if s<2:
            if metric=='roc': namestr = ['rocsignals,train ','rocsignals,test '][s]
            elif metric=='classacc': namestr = ['accsignals,cl1 ','accsignals,cl2 '][s]
        elif s<7:
            namestr = 'anglevsdtopc%d '%(s-2+1)     # only PC1-5 angles shown
        else:
            namestr = 'decisionboundarycoefficients%d '%(s-7+1)   # feature coefficients: width x neurons
        rocsignals.append( neo.AnalogSignal(  aux[s,:,:],\
           name=namestr+across_responses[0][0].name, t_start=sig_t_start,\
           sampling_period=across_responses[0][0].sampling_period, units=across_responses[0][0].units   ) )
    return rocsignals









def get_pointwisedecoderprobability(across_responses,width=50*pq.ms,rr=10,cv=skms.KFold):
    # returns timecourse of class CV-d prediction probabilities, and overall accuracy timecourse
    responses = np.concatenate(across_responses,axis=0)
    times = across_responses[0][0].times
    print('time',times[1]-times[0])
    wx = int(width/(times[1]-times[0]))
    print(responses.shape,times[-1])
    print('decoder width',wx,width)
    print('class numbers:',len(across_responses[0]),len(across_responses[1]))
#    responses = responses[:,T['epochstart_idx']:T['epochend_idx']+1,:]
    targets = np.hstack((np.zeros(len(across_responses[0]),dtype='int16'),\
                         np.ones(len(across_responses[1]),dtype='int16')))
    proba=[]
    for tx,t in enumerate(times[:-wx]):
        predictors = responses[:,tx:tx+wx,:]
        predictors = predictors.reshape(predictors.shape[0],predictors.shape[1]*predictors.shape[2])
    
        proba.append( pointwiseclassifierproba(predictors,targets,len(across_responses[0]),cv=cv,rr=rr) )
        # else:
        #     proba.append( np.zeros((len(targets),2)) )

    aux = np.array(np.swapaxes(proba,0,1))    # switch axes, so that the first dim becomes class, second time; third stays as proba class 1 and 2
    # probasignals is a (trials,times,classes); here we create only the probability of the first class, which has index 1, strange, eh?
    probasignals = neo.AnalogSignal(  aux[:,:,1],\
           name='probability'+across_responses[0][0].name, t_start=across_responses[0][0].t_start,\
           sampling_period=across_responses[0][0].sampling_period, units=across_responses[0][0].units   )

    # print(probasignals.shape,type(probasignals))
    predictions = (np.array(aux[:,:,0])<=0.5).astype(np.int16)
    # print(targets)
    # print(predictions[:,0])

    acc = np.zeros( (len(times[:-wx]), 3) )
    for tx,t in enumerate(times[:-wx]):
        acc[tx,0] = skme.accuracy_score(targets,predictions[:,tx])       # mean
        acc[tx,1] = np.sqrt( acc[tx,0]*(1-acc[tx,0]) )                                 # std
        acc[tx,2] = 2*acc[tx,1] / np.sqrt( len(targets) )                             # 2 sem

    # print(probasignals.shape,acc.shape)

    return probasignals, acc













    
    
def get_decoderprobability(across_responses,T=None,width=50*pq.ms,preon='on'):
    
    responses = np.concatenate(across_responses,axis=0)
    times = across_responses[0][0].times
    print('time',times[1]-times[0])
    wx = int(width/(times[1]-times[0]))
    print(responses.shape,times[-1])
    print('decoder width',wx,width)
    print('class numbers:',len(across_responses[0]),len(across_responses[1]))
    if preon=='pre':
        predictors = responses[:,T['stimstart_idx']-wx:T['stimstart_idx'],:].reshape((responses.shape[0],-1))
    elif preon=='on':
        predictors = responses[:,T['stimstart_idx']:T['stimstart_idx']+wx,:].reshape((responses.shape[0],-1))
    targets = np.hstack((np.zeros(len(across_responses[0]),dtype='int16'),\
                         np.ones(len(across_responses[1]),dtype='int16')))
    proba = pointwiseclassifierproba(predictors,targets,len(across_responses[0]),rr=10)
    
    return proba
    
    
    
    
    
    
    
def get_pointwisedecoder_tdistribution(across_responses,sample_analogsignal,metric='roc'):
    # here we only use responses from the masked array,
    # where invalid entries are ruled out due to runspeed equalization
    
    # across_responses have the usual (class)(trials,trajectory,neurons) shape, but masked, without Analogsignal structure
    # indexlist had (class,trajectory) shape, each array element is a list of integer trial index
    # we have an sample_analogsignal that stores times, starts and sampling rate information
    
    
    # times = across_responses[0][0].times
    times = sample_analogsignal['times']
    # print('time',times[1]-times[0])
    

    
    # wx = int(width/(times[1]-times[0]))
    # print('decoder width',wx,width)
    
    # print('full class numbers:',len(across_responses[0]),len(across_responses[1]))
    # print('mean equalized class numbers:',[len(indexlist[0,t]) for t in rangeindices.shape[1] ],len(across_responses[1,0]))
#    responses = responses[:,T['epochstart_idx']:T['epochend_idx']+1,:]
    roc=[]
    
    for tx,t in enumerate(times):
        mask = [ ~(across_responses[clx].mask[:,tx,0]) for clx in range(len(across_responses)) ]      # select valid entries (Falses in masked array)
        # print(tx,t,np.sum(~across_responses[0].mask[:,tx,0]),np.sum(~across_responses[1].mask[:,tx,0]))
        # print(tx,t,np.sum(mask[0]),np.sum(mask[1]))
        # targets = np.hstack( (np.zeros( np.sum(~across_responses[0].mask[:,tx,0]),dtype='int16'),\
        #                       np.ones( np.sum(~across_responses[1].mask[:,tx,0]),dtype='int16')) )
        targets = np.hstack( (np.zeros( np.sum(mask[0]),dtype='int16'),\
                              np.ones( np.sum(mask[1]),dtype='int16')) )

        predictors = [ across_responses[clx][mask[clx],tx,:] for clx in [0,1] ]
        predictors = np.concatenate( predictors, axis=0   )
        # print(tx,t,targets.shape,predictors.shape)
    
        if 1:# tx==0:   #use t==0  to speed up computations during debug
            roc.append( pointwiseclassifier(predictors,targets,len(across_responses[0]),cv=True,metric=metric) )
        else:
            roc.append(roc[tx-1])        # this is for speeding up computations for debug
    

    aux = np.array(np.swapaxes(roc,0,1))


    rocsignals = []
    for s in range(aux.shape[0]):
        if s<2:
            if metric=='roc': namestr = ['rocsignals,train ','rocsignals,test '][s]
            elif metric=='classacc': namestr = ['accsignals,cl1 ','accsignals,cl2 '][s]
        elif s<7:
            namestr = 'anglevsdtopc%d '%(s-2+1)     # only PC1-5 angles shown
        else:
            namestr = 'decisionboundarycoefficients%d '%(s-7+1)   # feature coefficients: width x neurons
        rocsignals.append( neo.AnalogSignal(  aux[s,:,:],\
           name=namestr, t_start=sample_analogsignal['t_start'],\
           sampling_period=sample_analogsignal['sampling_period'], units=sample_analogsignal['units']   ) )
    return rocsignals 
    
    
    






def get_offdecoding(off_responses,T,width=50*pq.ms,windowwidth=10*pq.ms,trainingpoint=-1500*pq.ms):
    responses = np.concatenate(off_responses,axis=0)
    starttimes = [trainingpoint,trainingpoint+width]
    traintx = int((trainingpoint/T['dt']).magnitude)-T['offset_idx']
    wx = int((width/T['dt']).magnitude)
    swx = int((windowwidth/T['dt']).magnitude)
    
    print(responses.shape,starttimes,wx)
    print('class numbers:',len(off_responses[0]),len(off_responses[1]))
#    responses = responses[:,T['epochstart_idx']:T['epochend_idx']+1,:]
    targets = np.hstack((np.zeros(len(off_responses[0]),dtype='int16'),\
                         np.ones(len(off_responses[1]),dtype='int16')))
    roc=[]
    trainpredictors = responses[:,traintx:traintx+wx,:]
    trainpredictors = trainpredictors.reshape(trainpredictors.shape[0],trainpredictors.shape[1]*trainpredictors.shape[2])
    
    
    # only start testing time points after training time point (lower triangle return matrix)
    # iterator = np.arange(traintx+wx,T['stimend_idx']-wx-T['offset_idx']+1,dtype='int16')
    
    # test the entire set (backward forward test, full return matrix)
    iterator = np.arange(T['start_idx'], T['end_idx']-wx-swx, swx, dtype='int16')


    # collect test points as list
    testpredictors = []
    for tx,t in enumerate(iterator):

        testpredictor = responses[:,t:t+wx,:]
        testpredictors.append( testpredictor.reshape(testpredictor.shape[0],testpredictor.shape[1]*testpredictor.shape[2]) )
        
    # print('iterating tx',tx,t,responses.shape,responses[:,t:t+wx,:].shape, trainpredictors.shape,testpredictors.shape)
    roc = offclassifier_multitest(trainpredictors,testpredictors,targets,len(off_responses[0]))
    # roc is (timepoints)(stats), we need analogsignals for first dimensions
    
    print(roc.shape)    

    # aux = np.array(np.swapaxes(roc,0,1))
    # print('shape',aux.shape)
    rocsignal = neo.AnalogSignal(  roc,\
       name='rocsignals'+off_responses[0][0].name, t_start=starttimes[0],\
       sampling_period=off_responses[0][0].sampling_period, units=off_responses[0][0].units   )
    
    # changed: t_start = starttimes[1] for test is changed to start for the whole set, [0] remained as train time
        
    return rocsignal










def get_offgradientdecoding(off_responses_train,off_responses_test,classes_test,T,width=100*pq.ms,trainingpoint=-1500*pq.ms,ma=10):
    responses_train = np.concatenate(off_responses_train,axis=0)
    starttimes = [trainingpoint,trainingpoint+width]

    traintx = int((trainingpoint/T['dt']).magnitude)-T['offset_idx']
    wx = int((width/T['dt']).magnitude)
    
    print(responses_train.shape,starttimes,wx)
    print('class numbers, train:',len(off_responses_train[0]),len(off_responses_train[1]))
#    responses_train = responses_train[:,T['epochstart_idx']:T['epochend_idx']+1,:]
    targets_train = np.hstack((np.zeros(len(off_responses_train[0]),dtype='int16'),\
                         np.ones(len(off_responses_train[1]),dtype='int16')))
    trainpredictors = responses_train[:,traintx:traintx+wx,:]
    trainpredictors = trainpredictors.reshape(trainpredictors.shape[0],trainpredictors.shape[1]*trainpredictors.shape[2])
    

    # collect test trials one by one

    responses_test = off_responses_test#np.concatenate(off_responses_test,axis=0)
    targets_test = classes_test#.astype('int16')

    
    pry_triallist = []

    for rx,r in enumerate(responses_test):
        # iterator = np.arange(traintx+wx,T['end_idx']-T['start_idx']-wx,dtype='int16')      # original, this gives negative index error
        # print('offsets',traintx,wx,T['offset_idx'],T['stimend_idx'],'results',traintx+wx-T['offset_idx'],T['stimend_idx']-wx-T['offset_idx']+1)
        iterator = np.arange(traintx+wx,T['stimend_idx']-wx-T['offset_idx']+1,dtype='int16')
        # print('start',traintx,traintx+wx,'end',T['stimend_idx']-wx-T['offset_idx'])
        # print('start',iterator[0],'end',iterator[-1])
        pry=[]
        for tx,t in enumerate(iterator):
            # if t in [iterator[0],iterator[-1]]: print(rx,r.shape,'traintx',traintx,'lim',iterator[0],iterator[-1],'ix',t,t+wx)
            testpredictors = r[t:t+wx,:]
            testpredictors = testpredictors.reshape(1,testpredictors.shape[0]*testpredictors.shape[1])
            pry.append(offgradientclassifier(trainpredictors, testpredictors, targets_train, np.array([targets_test[rx]])))
        pry_triallist.append(pry)     # triallist:      trials       x       timecourse       x   mean,std,s.e.m

    aux = np.array(pry_triallist)
    aux = np.swapaxes(aux,1,2)             # swap timecourses and statistics dimensions
    prysignals_triallist = []
    for rx,r in enumerate(responses_test):
        prysignals = []
        prysignals.append( neo.AnalogSignal(  aux[rx,0,0:wx,:],\
           name='prysignals'+off_responses_train[0][0].name, t_start=starttimes[0],\
           sampling_period=off_responses_train[0][0].sampling_period, units=off_responses_train[0][0].units   ) )
        prysignals.append( neo.AnalogSignal(  aux[rx,1,:,:],\
           name='prysignals'+off_responses_train[0][0].name, t_start=starttimes[1],\
           sampling_period=off_responses_train[0][0].sampling_period, units=off_responses_train[0][0].units   ) )
        prysignals_triallist.append(prysignals)
    
    return prysignals_triallist







def get_crosscontextdecoding(across_responses_trained, across_responses_predicted, width=50*pq.ms):
    # example when trained and predicted across_responses are different:
    # train responses:   multimodal Ignore   75%
    # test responses:    multimodal Attend   25%
    #                    multimodal Ignore   25%
    # responses will always contain: attended and ignore
    
    responses_trained = np.concatenate(across_responses_trained,axis=0)
    responses_predicted = np.concatenate(across_responses_predicted,axis=0)
    times = across_responses_trained[0][0].times
    print('time',times[1]-times[0])
    wx = int(width/(times[1]-times[0]))
    print(responses_trained.shape,responses_predicted.shape,times[-1])
    print('decoder width',wx,width)
    print('class numbers:',len(across_responses_trained[0]),len(across_responses_trained[1]),\
                           len(across_responses_predicted[0]),len(across_responses_predicted[1]))
#    responses = responses[:,T['epochstart_idx']:T['epochend_idx']+1,:]
    targets_trained = np.hstack((np.zeros(len(across_responses_trained[0]),dtype='int16'),\
                                np.ones(len(across_responses_trained[1]),dtype='int16')))
    targets_predicted = np.hstack((np.zeros(len(across_responses_predicted[0]),dtype='int16'),\
                                np.ones(len(across_responses_predicted[1]),dtype='int16')))
    
    roc=[]
    for tx,t in enumerate(times[:-wx]):
        predictors_trained = responses_trained[:,tx:tx+wx,:]
        predictors_trained = predictors_trained.reshape(predictors_trained.shape[0],predictors_trained.shape[1]*predictors_trained.shape[2])
        predictors_predicted = responses_predicted[:,tx:tx+wx,:]
        predictors_predicted = predictors_predicted.reshape(predictors_predicted.shape[0],predictors_predicted.shape[1]*predictors_predicted.shape[2])
    
        roc.append(   crossblock_pointwiseclassifier(predictors_trained,predictors_predicted,targets_trained,targets_predicted)   )

    aux = np.array(np.swapaxes(roc,0,1))
    rocsignals = []
    for s in range(aux.shape[0]):
        if s<3:
            namestr = ['rocsignals,train ','rocsignals,test ','rocsignals,crosstest'][s]
        elif s<8:
            namestr = 'anglevsdtopc%d '%(s-3+1)     # only PC1-5 angles shown
        else:
            namestr = 'decisionboundarycoefficients%d '%(s-7+1)   # feature coefficients: width x neurons
        rocsignals.append( neo.AnalogSignal(  aux[s,:,:],\
           name=namestr+across_responses_trained[0][0].name, t_start=across_responses_trained[0][0].t_start,\
           sampling_period=across_responses_trained[0][0].sampling_period, units=across_responses_trained[0][0].units   ) )
    return rocsignals






def get_crossgradient(ptrx,bl,dn,timeaverage_range,continuous_method='instfr'):
                
    offstimdecoderlist = []
    if ptrx==0:
        for trx in [0,1,2]:
            offstimdecoderlist.append(pickle.load(open(cacheprefix+'continuous/offstimgradientdecodes,allstim,b%d-%s-%s-t%d.pck'%(bl,dn,continuous_method,trx),'rb')))
    else:
        trx = 3
        offstimdecoderlist.append(pickle.load(open(cacheprefix+'continuous/offstimgradientdecodes,allstim,b%d-%s-%s-t%d.pck'%(bl,dn,continuous_method,trx),'rb')))
    n_cuetrials = len(offstimdecoderlist[0])
    
    d = []
    for n,offstimdecoders in enumerate(offstimdecoderlist):
        d_tp = []
        for j,dec in enumerate(offstimdecoders):
            d_tp.append(dec[1][:480,:]) # choose the test signals, and use the first 4800 ms (upgrade to T['dt'] normed index!), that covers all training point starts from -1500 ms to -100 ms
        d.append(np.array(d_tp))
    d = np.array(d)
    
#                dtm = d[:,timeaverage_range_off[0]:timeaverage_range_off[1],:].mean(axis=1).squeeze()
#                dtstd = d[:,timeaverage_range_off[0]:timeaverage_range_off[1],:].std(axis=1).squeeze()*2
#                dtsem = dtstd/np.sqrt(timeaverage_range_off[1]-timeaverage_range_off[0])


    # calculate statistics of the probability estimate: mean over several starting points and timeaverages
    # for s.e.m.  use the   square n  of n = # starting points  times   # timepoints
    dtm = d[:,:,timeaverage_range[0]:timeaverage_range[1],:].mean(axis=(0,2)).squeeze()
    dtstd = d[:,:,timeaverage_range[0]:timeaverage_range[1],:].std(axis=(0,2)).squeeze()*2
    dtsem = dtstd/np.sqrt( d.shape[0]  * (timeaverage_range[1]-timeaverage_range[0])  )
    
    return dtm,dtstd,dtsem,n_cuetrials





    
    




def get_linearregression(response_data,dependent_data,width=50*pq.ms,cv=5):
    responses = np.array(response_data)
    dependents = np.array(dependent_data)
    
    times = response_data[0].times
    print('time',times[1]-times[0])
    wx = int(width/(times[1]-times[0]))
    print(responses.shape,times[-1])
    print(dependents.shape,times[-1])
    print('regression width',wx,width)

    regr = []# = np.zeros((times[:-wx:wx].shape[0],3))
    runs = []
    for tx,t in enumerate(times[:-wx:wx]):
        predictors = responses[:,tx:tx+wx,:]
        
        # now these are to be trained as i.i.d. process: pack width trajectories along with trials
        predictors = predictors.reshape(predictors.shape[0]*predictors.shape[1],predictors.shape[2])
#        predictors = predictors.reshape(predictors.shape[0],predictors.shape[1]*predictors.shape[2])
        
#        target = dependents[:,tx:tx+wx].mean(axis=1)   # here we go for the mean of the dependent variable in the time-window
        target = dependents[:,tx:tx+wx].reshape(-1)   # here we go for repacking for each trajectory point and trial
    
        repack, R = pointwiselinearregression(predictors,target,wx=wx,rr=cv)
        
        # R trajectory, traintest, statistics
#        repack.append( aux.reshape(-1,wx) )
#        print('repack',repack.shape,repack[0,0])
#        regr.append(repack.reshape(15,wx,3))
        regr.append(repack)
        runs.append(R)

    print('bf sw',np.array(regr).shape)
    regr = np.swapaxes(np.array(regr),0,1)
    print(regr.shape)
    regrsignals = []
    for s in range(regr.shape[0]):
        if s<2:
            namestr = ['regrsignals,train ','regrsignals,test '][s]
        elif s<7:
            namestr = 'regrsignals,angletopc%d'%(s-2+1)
        else:
            namestr = 'regressioncoefficients%d '%(s-7+1)   # regression coefficients: window width  x  n_neurons
        regrsignals.append( neo.AnalogSignal(  regr[s,:,:],\
           name=namestr+response_data[0].name, t_start=response_data[0].t_start,\
           sampling_period=response_data[0].sampling_period*wx, units=(dependent_data[0].units)**2   ) )

#    print('regressionsignals',regrsignals[1].shape)
    
    runs = np.moveaxis(np.array(runs),3,0)
#    runs = np.array(runs)
#    print(runs.shape)

    
    runspeedestimates =     [ [ neo.AnalogSignal(  runs[trx,:,trtex,:],\
           name='runspeed '+response_data[0].name+' linreg', t_start=response_data[0].t_start,\
           sampling_period=response_data[0].sampling_period*wx, units=(dependent_data[0].units)   ) \
               for trtex in range(2)    ] \
               for trx in range(len(runs)) ]

    
    #return with score and runspeed estimates
    
    return [regrsignals,runspeedestimates]








# gaussian process kernel methods
def get_gaussianprocessregression(response_data,dependent_data,width=50*pq.ms,cv=3):

    n_neurons = response_data[0].shape[1]

    n_trials = len(response_data)
    n_timepoints = response_data[0].shape[0]-1

    responses = np.array(response_data)
    dependents = np.array(dependent_data)
    
    
    times = response_data[0].times
    print('time',times[1]-times[0])
    wx = int(width/(times[1]-times[0]))       # wx is the index step for downsampling, width is the timewindow for downsampling
    print('\nneural response',responses.shape,times[-1])
    print('running speed',dependents.shape,times[-1])
    print('gaussian process regression width',width,'index step',wx,'\n')
    
    # now resample for the widths window with mean value within window:
    responses_aux = responses[:,:-1:wx,:]    # auxiliary variable for value assign
    dependents = dependents[:,:-1].reshape(n_trials, wx, -1).mean(axis=1)
    for n in range(n_neurons): responses_aux[:,:,n] = responses[:,:-1,n].reshape(n_trials, wx,-1).mean(axis=1)
    responses = responses_aux; del responses_aux
    print('neural response shape',responses.shape)
    print('running speed shape',dependents.shape,'\n')
    
    # now create the trial by trial block matrix
#    X = responses.swapaxes(1,2).reshape(n_trials,-1)
    X = responses
    y = dependents
#    X = responses.reshape(n_trials*n_timepoints,n_neurons)
#    y = dependents.reshape(-1)
    print('concatenated blocks for input to GP',X.shape,y.shape)
    
    
    sc,R = neba.gpkron(X,y)
    
    print(sc.shape,R.shape)
    
#    sc,R = gaussianprocessregression(X,y,wx,rr=cv)
#    R = np.moveaxis(np.moveaxis(np.array(R),3,0), 3,0)
    
    # sc is the scores, R2: train-test, trajectories, mean-std-sem
    # R is the predictions:   trials, trajectories, train-test, mean-std
    

    signals = []
    
    signals.append(  [ neo.AnalogSignal(  sc[trtex],\
           name='R2 explained variance '+response_data[0].name+' gaussprocess', t_start=response_data[0].t_start,\
           sampling_period=response_data[0].sampling_period*wx, units=pq.dimensionless      ) \
                for trtex in range(2)    ]  )

    # signals runspeed will be:   [trials][train-test] AnalogSignal[trajectories, mean-std]
    signals.append([ [ neo.AnalogSignal(  R[trx,:,trtex,:],\
           name='runspeed '+response_data[0].name+' gaussprocess', t_start=response_data[0].t_start,\
           sampling_period=response_data[0].sampling_period*wx, units=(dependent_data[0].units)   ) \
               for trtex in range(2)    ] \
               for trx in range(len(R)) ] )
    
#    print('signals R2',type(signals[0]),np.array(signals[0]).shape,\
#          '\nsignals run',type(signals[1]),np.array(signals[1]).shape)
    
    return signals
    





def gaussianprocessregression(X,y,wx,rr=3):

    
    
    rr=2

    # dummy test data:
#    n_trials = 100
#    n_trajectory = 30
#
#    X = np.stack( [np.random.randn(n_trials,n_trajectory)-5,np.cumsum(np.random.randn(n_trials,n_trajectory)/2+3,axis=1)-5, np.random.randn(n_trials,n_trajectory)/4+4 ], axis=2)
#    y = np.sin( np.arange(n_trajectory) ) +np.random.randn(n_trials,n_trajectory)/20
#    
#    print('shapes',X.shape,y.shape)
#    X = X.mean(axis=1)
#    y = y.mean(axis=1)
#    
#    [ plt.plot(X[:,:,k].T) for k in range(3) ]
#    plt.plot(y.T)
#
#    X_tr,X_te,y_tr,y_te = skms.train_test_split(X, y, test_size=0.25)
#    
#    y_tr_,y_te_ = neba.gaussianprocessregression(X_tr,y_tr,X_te,y_te)
#    
#    
#    return



    n_trials, n_trajectory = y.shape[0], y.shape[1]



    
    # use RBF (squared exponential) kernel which is the default;
    # we want to use anisotropic kernel, one for each cell:
    # then we want an index kernel as a full cov matrix amongst cells; can be constructed by n x n constant kernels
    # and combine the two as product
    irregularities = gp.kernels.RationalQuadratic(length_scale=20/wx,alpha = 2)
    fastkernel = gp.kernels.RBF( length_scale=20/wx )
    slowkernel = gp.kernels.RBF( length_scale=100/wx )
    noise = gp.kernels.WhiteKernel(noise_level=1.)
    constantkernel = gp.kernels.ConstantKernel(constant_value=30)       # ,constant_value_bounds=(0,150)
    K = 1*constantkernel + 1*irregularities + 2*slowkernel + 5 *fastkernel + 0.5*noise
#    K = 1*constantkernel  + 2*slowkernel + 5 *fastkernel + 0.5*noise
#    K = 1*constantkernel + 0.1*noise
#    K = 1*constantkernel
#    K = 1*slowkernel
#    K = 1*fastkernel
#    K = 0.1*noise
    
    print('data shapes:',X.shape,y.shape)
#    print('kernel',K,K.diag(np.arange(30)))
    scorer = skme.make_scorer(skme.r2_score,multioutput='uniform_average')
    cv_results = skms.cross_validate(gp.GaussianProcessRegressor(kernel=K), \
                                X, y, cv=rr, scoring=scorer, \
                                return_train_score=True, return_estimator=True)
#    sc_tr = cv_results['train_score']
#    sc_te = cv_results['test_score']
    print(cv_results.keys())
    
    print(cv_results['estimator'][0].kernel_)
    print(cv_results['test_score'])
    
#    plt.figure(figsize=(24,12))
#    plt.plot(X,y,'.',alpha=0.2);
#    plt.figure(figsize=(24,12))
#    plt.plot(X, cv_results['estimator'][0].predict(X),'.',alpha=0.2)
#    return


    splits = []            # this is needed for prediction output
    for tr_ix,te_ix in skms.KFold(rr).split(X,y):
        splits.append([tr_ix,te_ix])
    
    y_tr_ = np.zeros((n_trials,n_trajectory,rr))
    y_te_ = np.zeros((n_trials,n_trajectory,rr))
    sc_tr = []
    sc_te = []
    
    for r in range(rr):
        # predict y:
#        y_tr_[splits[r][0],:,r] += cv_results['estimator'][r].predict(X[splits[r][0],:]) / (rr-1)
        y_tr_[splits[r][0],:,r] += cv_results['estimator'][r].y_train_  / (rr-1)
        y_te_[splits[r][1],:,r]  = cv_results['estimator'][r].predict(X[splits[r][1],:])

    for r in range(rr):
        # get error measures:
        sc_tr.append(skme.r2_score(y[splits[r][0],:], y_tr_[splits[r][0],:,r], multioutput='raw_values'))
        sc_te.append(skme.r2_score(y[splits[r][1],:], y_te_[splits[r][1],:,r], multioutput='raw_values'))
    
    
    run = np.array([   [y_tr_.mean(axis=2),y_tr_.std(axis=2), 2*y_tr_.std(axis=2)/np.sqrt(len(splits[0][0])*rr)],\
                        [y_te_.mean(axis=2),y_te_.std(axis=2), 2*y_te_.std(axis=2)/np.sqrt(len(splits[0][1])*rr)]   ])
    
    


    # single trial run speed estimates, this will be          trials   x  trajectory  x   mean, std, sem of CV groups
#    strse = []
#    splits = []
#    for tr_ix,te_ix in skms.KFold(rr).split(X,y):
#        splits.append([tr_ix,te_ix])
#
#    R = np.zeros((len(y),y.shape[1],2,2))         # trials, trajectories, train-test, mean-std
#    sc = np.zeros((2,y.shape[1],rr))            # train/test,trajectory,runs
#    
#    for r,ix in enumerate(splits):     # CV:    go over all splits
#        X_tr =  X[ix[0],:]
#        y_tr =  y[ix[0],:]
#        X_te =  X[ix[1],:]
#        y_te =  y[ix[1],:]
#        
#        
#        est = gp.GaussianProcessRegressor(kernel=K)
#        est.fit(X_tr, y_tr)
#        ym_tr_,ys_tr_  = est.predict(X_tr,return_std=True)
#        ym_te_,ys_te_ = est.predict(X_te,return_std=True)
#
#
#        R[ix[0],:,0,0] += ym_tr_ / (rr-1)          # means over the train split parts
#        R[ix[0],:,0,1] += ys_tr_[:,np.newaxis] / (rr-1)
#        
#        R[ix[1],:,1,0] = ym_te_      # just add up as tests only there for one split
#        R[ix[1],:,1,1] = ys_te_[:,np.newaxis]
#        
#        
#        sc[0,:,r] = skme.r2_score(y_tr,ym_tr_,multioutput='raw_values')
#        sc[1,:,r] = skme.r2_score(y_te,ym_te_,multioutput='raw_values')


#    for trx in range(len(X)):
#        R = np.zeros((len(y),y.shape[1],2,rr))   # trials, trajectory, train-test, mean-std, split ix
#        for r in range(rr): # go over stats
#            m,s = cv_results['estimator'][r].predict(X,return_std=True)
#            R[:,:,0,r] = m
#            R[:,:,1,r] = s[:,np.newaxis]
#    strse = np.stack([R.mean(axis=-1), 2*R.std(axis=-1), 2*R.std(axis=-1)/np.sqrt(rr)],axis=3)

    sc = np.array([sc_tr,sc_te])
    
#    print('multiple estimates R$^2$',sc[1,:,:].mean(axis=0))

    # scores will be of       te,tr         x          mean, 2*std,  2*s.e.m.
    sc_dist = sc        # te,tr    x     1,2,3,4... reruns
    sc = np.stack(   (np.mean(sc_dist,axis=1),2*np.std(sc_dist,axis=1),2*np.std(sc_dist,axis=1)/np.sqrt(rr)), axis=1)
    
    sc = np.moveaxis(sc,2,1)   # move trajectory into the middle before stats


    
    pack = (sc,run)
    

    return pack









def get_linearregression_singleprocess(X,Y,indices,method=sklm.Ridge(),cv_method=skms.StratifiedKFold(5)):
    # elementary per timepoint decoder
    n_targets = Y.shape[2]
    # a = np.frombuffer(sharedvars{'v'}).reshape(sharedvars{'sh_acc'})
    a = np.zeros((len(indices), n_targets, 2,3))
    c = np.zeros((len(indices), n_targets, X.shape[2], 3))
    for tg in range(n_targets):
        for kx,k in enumerate(indices):
            x = X[:,k,:]
            y = Y[:,k,tg]
            cv = cv_method.get_n_splits(x,y)
            cv_results = skms.cross_validate(method, \
                                    x, y, cv=cv,return_train_score=True, return_estimator=True)
            sc_tr = cv_results['train_score']
            sc_te = cv_results['test_score']
            a[kx,tg,:,:] = np.array( \
                [ [ sc_tr.mean(), sc_tr.std(), sc_tr.std()/np.sqrt(cv) ], \
                  [ sc_te.mean(), sc_te.std(), sc_te.std()/np.sqrt(cv) ] ] )
            cs = np.array([ cv_results['estimator'][r].coef_ for r in range(cv) ]).squeeze()
            c[kx,tg,:,:] = np.array([cs.mean(axis=0), cs.std(axis=0), cs.std(axis=0)/np.sqrt(cv)]).T
    return indices, a, c


def get_crosslinearregression_singleprocess(X,Y,Xcross,Ycross,indices,method=sklm.Ridge(),cv_method=skms.StratifiedKFold(5)):
    # elementary per timepoint decoder
    # compares two subsets of observations, training in one subset, and testing in another, both CV folds
    n_targets = Y.shape[2]
    # a = np.frombuffer(sharedvars{'v'}).reshape(sharedvars{'sh_acc'})
    a = np.zeros((len(indices), n_targets, 2+2,3))         # need +2 for cross tests
    c = np.zeros((len(indices), n_targets, X.shape[2], 3))
    for tg in range(n_targets):
        for kx,k in enumerate(indices):
            x = X[:,k,:]
            y = Y[:,k,tg]
            cv = cv_method.get_n_splits(x,y)#,groups=y) # this will be the same in both cross-groups as equalized
            cv_results = skms.cross_validate(method, \
                                    x, y, cv=cv_method,return_train_score=True, return_estimator=True)
            # print(cv_method, 'nsplits', cv, 'y',sum(y),sum(y==0), sum(Ycross[:,k,tg]), sum(Ycross[:,k,tg]==0), '#estim',len(cv_results['estimator']))
            sc_tr = cv_results['train_score']
            sc_te = cv_results['test_score']
            sc_tr_cross = np.zeros(Ycross.shape[0])
            sc_te_cross = np.zeros(Ycross.shape[0])
            S = cv_method
            for ix, (i_tr, i_te) in enumerate(S.split(Xcross[:,k,:],Ycross[:,k,tg])):      #,groups=Ycross[:,k,tg])):
                sc_tr_cross[i_tr] = cv_results['estimator'][ix].score(Xcross[:,k,:][i_tr,:], Ycross[:,k,tg][i_tr])
                sc_te_cross[i_te] = cv_results['estimator'][ix].score(Xcross[:,k,:][i_te,:], Ycross[:,k,tg][i_te])


            a[kx,tg,:,:] = np.array( \
                [ [ sc_tr.mean(), sc_tr.std(), sc_tr.std()/np.sqrt(cv) ], \
                  [ sc_te.mean(), sc_te.std(), sc_te.std()/np.sqrt(cv) ], \
                  [ sc_tr_cross.mean(), sc_tr_cross.std(), sc_tr_cross.std()/np.sqrt(cv) ], \
                  [ sc_te_cross.mean(), sc_te_cross.std(), sc_te_cross.std()/np.sqrt(cv) ] ] )
            
            cs = np.array([ cv_results['estimator'][r].coef_ for r in range(cv) ]).squeeze()
            c[kx,tg,:,:] = np.array([cs.mean(axis=0), cs.std(axis=0), cs.std(axis=0)/np.sqrt(cv)]).T
    return indices, a, c


def get_linearregressionmultivariate_singleprocess(X,Y,indices,method=sklm.Ridge(),cv_method=skms.StratifiedKFold(5)):
    # elementary per timepoint decoder
    n_targets = Y.shape[2]

    # a = np.frombuffer(sharedvars{'v'}).reshape(sharedvars{'sh_acc'})
    a = np.zeros((len(indices), 1, 2,3))
    c = np.zeros((len(indices), n_targets, X.shape[2], 3))
    for kx,k in enumerate(indices):
        x = X[:,k,:]
        y = Y[:,k,:]
        cv = cv_method.get_n_splits(x,y)
        cv_results = skms.cross_validate(method, \
                                x, y, cv=cv_method,return_train_score=True, return_estimator=True)
        sc_tr = cv_results['train_score']
        sc_te = cv_results['test_score']
        a[kx,0,:,:] = np.array( \
            [ [ sc_tr.mean(), sc_tr.std(), sc_tr.std()/np.sqrt(cv) ], \
                [ sc_te.mean(), sc_te.std(), sc_te.std()/np.sqrt(cv) ] ] )
        cs = np.array([ cv_results['estimator'][r].coef_ for r in range(cv) ])
        c[kx,:,:,:] = np.stack([cs.mean(axis=0), cs.std(axis=0), cs.std(axis=0)/np.sqrt(cv)],axis=2)
    return indices, a, c


def get_linearregressionmultivariate(X,Y,times,method='regression',reguralization='ridge',alpha=1.0, cv=5,
                                         rank=None,crossindicestrain=None, crossindicestest=None):
    # create multiple decoders for each output in Y for each timepoint
    # X and Y must be: (observations,times,features) dimensions
    # method is either 'regression', 'reducedrankregression', or 'classification', 'crossclassification'
    # reguralization is either 'ridge' (default) or 'lasso', where alpha is the regularization parameter (default 1.0)
    # ranks is for multivariate reduced regression
    # crossindices is train+cv test and crosstest observations indices available


    if method=='regression':
        if reguralization=='ridge':
            linmethod = sklm.Ridge(alpha=alpha)
        elif reguralization=='lasso':
            linmethod = sklm.Lasso(alpha=alpha)
    elif method=='reducedrankregression': linmethod = ReducedRankRidge(fit_intercept=False, rank=rank)    # ridge_solver='liblinear', 
    elif method=='classification' or \
         method=='crossclassification': linmethod = sklm.LogisticRegression(solver='liblinear', penalty='l2')
    else: print('not allowed method')

    if cv=='loo' or cv=='LOO': cv_method = skms.LeaveOneOut()
    else: cv_method = skms.StratifiedKFold(cv)


    n_observations, n_timestamps, n_features = X.shape
    n_targets = Y.shape[2]


    indexpartition = np.array_split(np.arange(n_timestamps), n_cpu)

    if method=='crossclassification':
        acc = np.zeros((n_timestamps,n_targets,2+2,3))   # last two dims is (train/test + train/test cross, stats)
    else:
        acc = np.zeros((n_timestamps,n_targets,2,3))   # last two dims is (train/test,stats)
    coef = np.zeros((n_timestamps,n_targets,n_features,3))   # last dim is (stats)
    print(method, "\n", linmethod)
    # with multiprocessing.Pool(processes=n_cpu, initializer=init_worker, initargs=(acc_p, sh_acc)) as pool:
    with multiprocessing.Pool(processes=n_cpu) as pool:
        for prc in range(n_cpu):
            if method=='reducedrankregression':
                r = pool.apply_async( get_linearregressionmultivariate_singleprocess, (X,Y,indexpartition[prc],linmethod,cv_method) )
            elif method=='crossclassification':
                r = pool.apply_async( get_crosslinearregression_singleprocess, (X[crossindicestrain,:,:],Y[crossindicestrain,:,:], \
                                                                                X[crossindicestest,:,:],Y[crossindicestest,:,:], \
                                                                                indexpartition[prc],linmethod,cv_method) )
            else:
                r = pool.apply_async( get_linearregression_singleprocess, (X,Y,indexpartition[prc],linmethod,cv_method) )
            ind, a, c = r.get()
            # print('process',prc,ind, a.shape, acc.shape)
            acc[ind,:,:,:] = a
            coef[ind,:,:,:] = c
    return acc, coef









def get_singleregression(x, y, method=sklm.Ridge(), cv_method=skms.KFold(10)):
    n_observations = x.shape[0]
    acc = np.zeros((2,3))   # last two dims is (train/test,stats)
    coef = np.zeros(3)   # last dim is (stats)
    
    cv = cv_method.get_n_splits(x,y)
    
    cv_results = skms.cross_validate(method, \
                            x, y, cv=cv_method, return_train_score=True, return_estimator=True)
    sc_tr = cv_results['train_score']
    sc_te = cv_results['test_score']

    
    # acc (2,3) (train/test,stats)
    acc = np.array( \
        [ [ sc_tr.mean(), sc_tr.std(), sc_tr.std()/np.sqrt(cv) ], \
            [ sc_te.mean(), sc_te.std(), sc_te.std()/np.sqrt(cv) ] ] )
    
    # coef (3) (stats)
    cs = np.array([ cv_results['estimator'][r].coef_ for r in range(cv) ])
    coef = np.array([cs.mean(axis=0), cs.std(axis=0), cs.std(axis=0)/np.sqrt(cv)]).squeeze()

    # cv predictions
    y_te = y.copy()
    for cvx,(is_tr, is_te) in enumerate(cv_method.split(x)):
        y_te[is_te] = cv_results['estimator'][cvx].predict(x[is_te])
    
    return acc,coef,y_te.squeeze()








# 1st and 2nd order statistics

def getorderstatistics(acrossresponses):
    
    n_bootstrap = 20
    
    n_trials = min(len(acrossresponses[0]),len(acrossresponses[1]))
    splitpoint = round(n_trials/2)      # split trials to compare fairly


    # meanstatistics = response1,2  x    (means, std, s.e.m.)   x   timecourse  x  channels
    # correlations = response1,2  x   timecourse  x  channel_pairs_extracted_offdiaguppertriangle
    
    meanstatistics = []
    correlations = []

    bootstrappedm = []
    bootstrappedc = []
    for boot in range(n_bootstrap):
        indexlist = np.random.permutation(n_trials)
        indexsubgroup1 = indexlist[:splitpoint]
        indexsubgroup2 = indexlist[splitpoint:]
        aux1 = (np.array(acrossresponses[0]))[indexsubgroup1,:,:]
        aux2 = (np.array(acrossresponses[1]))[indexsubgroup2,:,:]
#        print(aux1.shape)
        bootstrappedm.append(  np.array([[aux1.mean(axis=0),aux1.var(axis=0),aux1.std(axis=0)/np.sqrt(len(indexsubgroup1))],\
                                        [aux2.mean(axis=0),aux2.var(axis=0),aux2.std(axis=0)/np.sqrt(len(indexsubgroup2))]] )  )
#        print('corr:',np.corrcoef(aux1[:,15,:].squeeze(),rowvar=False)[np.triu_indices(aux1.shape[2],k=1)].shape)
        bootstrapping1 = []    
        bootstrapping2 = []    
        for t in range(aux1.shape[1]):
            bootstrapping1.append(np.corrcoef(aux1[:,t,:].squeeze(),rowvar=False)[np.triu_indices(n=aux1.shape[2],k=1)])
            bootstrapping2.append(np.corrcoef(aux2[:,t,:].squeeze(),rowvar=False)[np.triu_indices(n=aux2.shape[2],k=1)])

#        bootstrappedc.append(  np.array([[ np.corrcoef(aux1[:,t,:].squeeze(),rowvar=False)[np.triu_indices(n=aux1.shape[2],k=1)]\
#                                             for t in range(aux1.shape[1])],\
#                                         [ np.corrcoef(aux2[:,t,:].squeeze(),rowvar=False)[np.triu_indices(n=aux1.shape[2],k=1)]\
#                                             for t in range(aux2.shape[1])] ]  )   )
        bootstrappedc.append(  [bootstrapping1,bootstrapping2]  )
    #        print('end:',len(bootstrappedc[0][0]),bootstrappedc[0][0][0].shape)

    bootstrappedm_array = np.array(bootstrappedm)
    bootstrappedmeans = bootstrappedm_array.mean(axis=0)
    singlemeans = bootstrappedm_array[0,:,:,:,:].squeeze()    # this is a single instance of split trial groups to show an example
    bootstrappedmeandifferences = (bootstrappedm_array[:,1,:,:,:]-bootstrappedm_array[:,0,:,:,:]).squeeze().mean(axis=0)  # this will show the bootstrapped differences
    meanstatistics = [ [ neo.AnalogSignal(singlemeans[j,k,:,:],\
                    name='means, variances, half trials'+acrossresponses[0][0].name,\
                    t_start=acrossresponses[0][0].t_start,
                    sampling_period=acrossresponses[0][0].sampling_period,\
                    units=acrossresponses[0][0].units)\
                    for k in range(3) ] for j in range(2) ]
    meanstatistics_bsdiff = [ neo.AnalogSignal(bootstrappedmeandifferences[k,:,:],\
                    name='means, variances, bootstrapped differences '+acrossresponses[0][0].name,\
                    t_start=acrossresponses[0][0].t_start,
                    sampling_period=acrossresponses[0][0].sampling_period,\
                    units=acrossresponses[0][0].units)\
                    for k in range(3) ]

    


    bootstrappedc_array = np.array(bootstrappedc)
    bootstrappedcorrelations = bootstrappedc_array.mean(axis=0)
    singlecorrelations = bootstrappedc_array[0,:,:,:].squeeze()    # this is a single instance of split trial groups to show an example
    bootstrappedcorrelationdifferences = (bootstrappedc_array[:,1,:,:]-bootstrappedc_array[:,0,:,:]).squeeze().mean(axis=0)  # this will show the bootstrapped differences


#    print(bootstrappedcorrelations.shape)
    correlations =  [ neo.AnalogSignal(singlecorrelations[j,:,:],\
                    name='correlations, half trials '+acrossresponses[0][0].name,\
                    t_start=acrossresponses[0][0].t_start,
                    sampling_period=acrossresponses[0][0].sampling_period,\
                    units=acrossresponses[0][0].units)\
                    for j in range(2) ]

    correlations_bsdiff = neo.AnalogSignal(bootstrappedcorrelationdifferences,\
                    name='correlations, bootstrapped differences '+acrossresponses[0][0].name,\
                    t_start=acrossresponses[0][0].t_start,
                    sampling_period=acrossresponses[0][0].sampling_period,\
                    units=acrossresponses[0][0].units)
    
    
#    print(len(meanstatistics),len(meanstatistics[0]),meanstatistics[0][0].shape,type(meanstatistics[0][0]))
#    print(len(correlations),correlations[0].shape,type(correlations[0]))




    
    return meanstatistics,meanstatistics_bsdiff,bootstrappedmeans,    correlations,correlations_bsdiff,bootstrappedcorrelations













def get_orthogonalprojections_ontopcaxes(X,coeffs,tpc_coords,checkpoint = 50*pq.ms, vs = 5.0):
    # project all original neural activites (X) onto the subspace spanned by the first tpc_coords coordinates of PCA
    # choose a single "checkpoint" timepoint, and use that to display the subspaces, pointclouds and vectors
    
    # checkpoint variable is the checkpoint after stimulus presentation start
    t = checkpoint - X[0][0].times[0]           # get the index of the checkpoint, with pre-stimulus start times
    tx = int((t/(10*pq.ms)).magnitude)

    #vs is vectorscale from unit norm, so that the arrows ar emore visible



    # split the neural activity by class
    Xa = np.array(X[0])[:,tx,:].squeeze()          # class 1
    Xb = np.array(X[1])[:,tx,:].squeeze()          # class 2
    
    
    Xab = np.vstack((Xa,Xb))      # all classes together

    
    



    # first calculate per class principal components
    
    Ua,Sa,Va_h = np.linalg.svd(Xa)
#    PaXa = Ua[:,tpc_coords]*Sa[tpc_coords]

    Ub,Sb,Vb_h = np.linalg.svd(Xb)
#    PbXb = Ub[:,tpc_coords]*Sb[tpc_coords]


    # project both class onto the total principal component projection:
    
    U,S,V_h = np.linalg.svd( Xab )
    V = V_h.T          # get V from the horiz. row space form to the basis-in-columnvectors form
    
    
    B = V[:,tpc_coords].dot(V[:,tpc_coords].T)      # domain (row space) projection to base matrix; this is a symmetric matrix: B==B.T True

    
    # transform the per-class-pc coordinates of the datapoints onto the common PCA
    PXa = np.dot(  Xa.dot(Va_h), B)
    PXb = np.dot(  Xb.dot(Vb_h), B)

    # here we transform the original neural coordinates into the total pc axes to compare with
    # the coordinates projected from the per-class-pc
    PsXa = np.dot(Xa,V)[:,tpc_coords]
    PsXb = np.dot(Xb,V)[:,tpc_coords]


    
    # project per-class-pc directions onto the first 2 axes of the all-class-pc
    # the last slice shows:   [tpc_coords,...] first few pc-s in the total-pc basis (the plot axes)
    #                         [...,cpc_coords] and display only the first few pc-s from the class PCA
    PVa = np.dot(B,Va_h.T)[tpc_coords,:]*vs
    PVb = np.dot(B,Vb_h.T)[tpc_coords,:]*vs

#


    # and finally project the average decoder decision boundary normal onto the
    # principal components of the display axis i.e. total-pc; only a single vector to be transformed
    # We have to take the decision coeeficients also at the 'checkpoint' timepoint

    Pcoeffs = np.dot(B,np.array(coeffs)[:,tx,0])[tpc_coords]*vs


    return PsXa, PsXb, PXa, PXb, PVa, PVb, Pcoeffs
















def decompose_tca_tensortools(X,R=10,nn=False):
    # this only works with the old tensortools package
    if nn:
        H = tensortools.ncp_hals(X=X,rank=R)
    else:
        H = tensortools.cp_als(X=X,rank=R)
    return H



def decompose_tca(X,R=10,nn=False):
    X = tl.tensor(X)
    if nn:
        _, H = non_negative_tucker(tensor=X,rank=R)
    else:
        _, H = parafac(tensor=X, rank=5)
    return H








def classifierfortcalatents(X,Y):
    
    # X: features, neural data
    # Y: classes, stimulus or task values; this is to be learnt and predict

#    X_tr, X_te, Y_tr, Y_te = ml.train_test_split(X, Y, test_size=0.8)

    dec = sklm.LogisticRegressionCV(cv=3)
#    dec.fit(X_tr, Y_tr)
    dec.fit(X, Y)

#    Y_tr_ = dec.predict(X_tr)
#    Y_te_ = dec.predict(X_te)
    Y_ = dec.predict(X)

    roc = skme.roc_auc_score(Y,Y_)

    return roc











def multiglm(X,Y,rr=5,showfig=False):
    
    if showfig: fig,ax = plt.subplots(2,rr,figsize=(30,16))
    
    R2 = np.zeros(rr)
    W = np.zeros((rr,Y.shape[1],X.shape[1]))
    mll = np.zeros((rr,Y.shape[1]))
    A = np.zeros(rr)
    L = np.zeros(rr)
    
    n_components = X.shape[1]
    n_neurons = Y.shape[1]
    
    
    
    for r in range(rr):
        X_tr, X_te, Y_tr, Y_te = skms.train_test_split(X, Y, test_size=0.25)
    
        n_tts = len(X_tr)
        idxs_tr = np.arange(0,n_tts)
        idxs_te = np.arange(n_tts,len(X))

        glms = []
        for n in range(n_neurons):
            glms.append( sklm.BayesianRidge(compute_score=True) )        # return the MLL as well
            glms[n].fit(X_tr,Y_tr[:,n])
        
        Y_tr_ = np.zeros(Y_tr.shape)
        Y_te_ = np.zeros(Y_te.shape)
        for n in range(n_neurons):
            Y_tr_[:,n] = ( glms[n].predict(X_tr) )
            Y_te_[:,n] = ( glms[n].predict(X_te) )
            W[r,n,:] = glms[n].coef_
            mll[r,n] = glms[n].scores_[-1]          # get the marginal log likelihood of the data from the last iteration
        
        
        if showfig:
            ax[0,r].plot(idxs_tr,Y_tr,'o',color='navy',alpha=0.2)
            ax[0,r].plot(idxs_tr,Y_tr_,'o',color='firebrick',alpha=0.2)
            ax[0,r].plot(idxs_te,Y_te,'o',color='dodgerblue',alpha=0.2)
            ax[0,r].plot(idxs_te,Y_te_,'o',color='fuchsia',alpha=0.2)
            
            ax[1,r].plot(idxs_tr,Y_tr_-Y_tr,'o',color='darkred',alpha=0.2)
            ax[1,r].plot(idxs_te,Y_te_-Y_te,'o',color='red',alpha=0.2)
            
        
#        print(Y_te.shape,Y_te_.shape)        
        score = skme.r2_score(Y_te,Y_te_)
        R2[r] = score
        
        A[r] = aic(Y_te,Y_te_,n_components)
        
        L[r] = mll[r,:].mean()        # here   the /N in the mean will be the log 1/sqrt(pi*s^2)^N   which is proportional to 1/N
        
    
    if len(R2[R2>0])>0:
        R2_mean = np.nanmean(R2[R2>0],axis=0)
        R2_sem = np.nanstd(R2[R2>0],axis=0)/np.sqrt(np.sum(R2>0))
    else:
        R2_mean = 0
        R2_sem = np.inf
    

    W_mean = np.mean(W,axis=0)
    W_sem = np.std(W,axis=0)/np.sqrt(rr)
    
    AIC_mean = np.mean(A,axis=0)
    AIC_sem = np.std(A,axis=0)/np.sqrt(rr)
    
    L_mean = np.mean(L,axis=0)
    L_sem = np.std(L,axis=0)/np.sqrt(rr)
    
    return W_mean, W_sem, R2_mean, R2_sem, AIC_mean, AIC_sem, L_mean, L_sem





def subsetglm(X,Y):
    # this will go through a specific models-in-rows matrix, mc, and perform glm using components selected as 1 in the columns of each row
    # The selector only works for single variate response variable, so we have to perform for each neuron manually.
    # Original:
#    selector = skfs.RFECV(estimator=sklm.Ridge, step=1, min_features_to_select=4,\
#                     scoring=skme.make_scorer(aic))
#    selector.fit(X,Y[:,0])
#    return selector.grid_scores_, selector.support_, selector.ranking_


    n_predictors = X.shape[1]
    
    
    mc = datazeroermatrix(n_predictors-1)
    print(mc)

    n_models = len(mc)

    W = []        # neural coefficients
    R2 = []       # explained variance 
    A = []        # AIC
    D = []        # deviance as marginal log likelihood


    for mx,mselector in enumerate(mc):
        X_selected = X[:,np.where(mselector)[0]]
        wm, ws, mm,sm, aicm, aics, devm, devs = multiglm(X_selected,Y,rr=10,showfig=False)
        W.append(wm)
        R2.append(mm)
        A.append(aicm)
        D.append(devm)
    
    R2 = np.array(R2)
    A = np.array(A)
    D = np.array(D)
    
    print(R2[:16].reshape(-1,4))
    print(R2[16:])
    print(A[:16].reshape(-1,4))
    print(A[16:])
    print(D[:16].reshape(-1,4))
    print(D[16:])

    return mc, R2, A, D






def multinomialglm(X,y):
    
    X_tr, X_te, y_tr, y_te = skms.train_test_split(X, y, test_size=0.25)
    
#    n_tts = len(X_tr)
#    idxs_tr = np.arange(0,n_tts)
#    idxs_te = np.arange(n_tts,len(X))

    glm = sklm.LogisticRegression(multi_class='multinomial',solver='lbfgs')
    glm.fit(X_tr,y_tr)
    
    y_tr_ = glm.predict(X_tr)
    y_te_ = glm.predict(X_te)

    sc_tr = glm.score(X_tr,y_tr)
    sc_te = glm.score(X_te,y_te)
    
    r2_tr = skme.r2_score(y_tr,y_tr_)
    r2_te = skme.r2_score(y_te,y_te_)
    
    return r2_tr,r2_te,sc_tr,sc_te,y_tr_,y_te_






def fitgaussianmixture(data,n_components=2,means_init=None):
    # observations times features
    
    model = gm.GaussianMixture(n_components=n_components, covariance_type='full',means_init=means_init)
    labels = model.fit_predict(data)

    # Mean (n_components, n_features)
    # Cov (n_components, n_features, n_features)
    
    return labels,model.means_, model.covariances_,model












def get_decoder_dv(acrossdecoder, n_neuron, n_preweightelements=7, timeaverage_start_dx=150, timeaverage_end_dx=150+150):

    # collect DBNVs, to project the activities collected above onto
    wx = int((len(acrossdecoder)-n_preweightelements)/n_neuron)
    
    # take the mean across feature-expansion timepoints of each neuron: [features,trajectory,stats] -> [neurons,trajectory,stats]
    w = np.reshape(np.array(acrossdecoder[n_preweightelements:]),
                  (wx,n_neuron,acrossdecoder[n_preweightelements].shape[0],acrossdecoder[n_preweightelements].shape[1]) ).mean(axis=0)
        
    # take only the mean stat stored in the last dimension first element
    w = np.array(w[:,:,0]) # [neurons,trajectory]


    # mean decision vector with timecourse-averaging:        [neurons]
    w_mean = w[:,timeaverage_start_dx:timeaverage_end_dx].mean(axis=1)


    return w,w_mean











# reduced rank regression from: https://github.com/krey/rrpy/blob/main/tests/test_reduced_rank_ridge.py

def _fit_rrr_no_intercept_all_ranks(X: np.ndarray, Y: np.ndarray, alpha: float, solver: str):
    ridge = sklearn.linear_model.Ridge(alpha=alpha, fit_intercept=False, solver=solver)
    beta_ridge = ridge.fit(X, Y).coef_
    Lambda = np.eye(X.shape[1]) * np.sqrt(np.sqrt(alpha))
    X_star = np.concatenate((X, Lambda))
    Y_star = X_star @ beta_ridge.T
    _, _, Vt = np.linalg.svd(Y_star, full_matrices=False)
    return beta_ridge, Vt

def _fit_rrr_no_intercept(X: np.ndarray, Y: np.ndarray, alpha: float, rank: int, solver: str, memory=None):
    memory = sklearn.utils.validation.check_memory(memory)
    fit = memory.cache(_fit_rrr_no_intercept_all_ranks)
    beta_ridge, Vt = fit(X, Y, alpha, solver)
    return Vt[:rank, :].T @ (Vt[:rank, :] @ beta_ridge)

class ReducedRankRidge(sklearn.base.MultiOutputMixin, sklearn.base.RegressorMixin, sklearn.linear_model._base.LinearModel):
    def __init__(self, alpha=1.0, fit_intercept=True, rank=None, ridge_solver='auto', memory=None):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.rank = rank
        self.ridge_solver = ridge_solver
        self.memory = memory

    def fit(self, X, y):
        if self.fit_intercept:
            X_offset = np.average(X, axis=0)
            y_offset = np.average(y, axis=0)
            # doesn't modify inplace, unlike -=
            X = X - X_offset
            y = y - y_offset
        self.coef_ = _fit_rrr_no_intercept(X, y, self.alpha, self.rank, self.ridge_solver, self.memory)
        self.rank_ = np.linalg.matrix_rank(self.coef_)
        if self.fit_intercept:
            self.intercept_ = y_offset - X_offset @ self.coef_.T
        else:
            self.intercept_ = np.zeros(y.shape[1])
        return self


# def test_reduced_rank_regression_rank():
#     X, Y = sklearn.datasets.make_regression(n_samples=100, n_features=50, n_targets=50, random_state=0, n_informative=25)
#     estimator = ReducedRankRidge(fit_intercept=True, rank=20).fit(X, Y)
#     assert estimator.rank_ <= 20, f"{estimator.rank_} > 20"

















if __name__ == '__main__':
    print('analysis methods to discover structure in neurophysiology data')

