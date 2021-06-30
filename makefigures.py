#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 16:36:50 2019

@author: mahajnal
"""


from run import *


import numpy as np
import scipy as sp
import quantities as pq
import pandas as pd

import os
import pickle
import h5py

#import matplotlib.image as mimg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as clrs
import matplotlib.lines
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import mpl_toolkits.mplot3d.art3d as art3d

#import seaborn as sns

import neo

import neurophysiology as neph
import neurodiscover as nedi
import neurobayesian as neba

import preprocess
import figs


plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'legend.fontsize': 20})
plt.rcParams.update({'lines.linewidth': 1})

# info about figure absolute sizes and nice fonts: https://jwalton.info/Embed-Publication-Matplotlib-Latex/

globalsave = 0

# continuous_method = 'instfr'      # JRC spike sorting
continuous_method = 'ks2ifr'        # kilosort2 spike sorting, with no drift



# _____________
# publications:
conference = False
resultpath = '../../../publish/journals/journal2020spring/figures/'
ext = '.pdf'
ext = '.png'

# resultpathprefix = '../results_ks/'
# resultpathseries = 'tribe/'
# resultpath = resultpathprefix + resultpathseries




def drawschematics(ax,wix=None,what=''):
    whatlist = ['decoderscheme legend','decoderscheme intime','decoderscheme crosstime',\
                'decoderscheme crossblock','decoderscheme intime preonly',\
                'decoderscheme crosstime prepre','decoderscheme crosstime preon',\
                'decoderscheme crossblock aa+ai','decoderscheme withinblock aa+ii',\
                'decoderscheme crossblock ii+ia','decoderscheme PCA space']
    if wix!=None:
        what = whatlist[wix]

    # helper aliases:
    trp = [0,0,0,0]    # transparent color
    lw = 3             # linewidth
    


    
    #these lines help when designing the plots; comment them out for production
#    fig = plt.figure(figsize=(12,12))
#    ax = fig.gca()
#    what = 'decoderscheme crossblock'
#    ax.plot([0,0],[-0.5,0.5],'o')
#    ax.set_xlim(0,1)
#    ax.set_ylim(0,1)
#    ax.grid('on')

    
    
    # plot each schematics
    
    if what=='decoderscheme legend':         # 0
        ax.add_patch(plt.Rectangle((0.25,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.65,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.35,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.75,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.arrow(0.35,0.6,0,-0.18, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.75,0.6,0,-0.18, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.4,0.6,0.26,-0.22, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.7,0.6,-0.26,-0.22, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )


        ax.text(0.01,0.68, 'training\ndata', transform=ax.transAxes    )
        ax.text(0.01,0.28, 'testing\ndata', transform=ax.transAxes    )
        
        ax.text(0.35,0.85, 'condition 1', ha='center', transform=ax.transAxes    )
        ax.text(0.75,0.85, 'condition 2', ha='center', transform=ax.transAxes    )

#        ax.text(0.55,0.95, 'decoding paradigms' , ha='center', transform=ax.transAxes    )


    elif what=='decoderscheme intime':        # 1
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.7,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.8,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.arrow(0.2,0.6,0,-0.18, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.5,0.6,0,-0.18, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.8,0.6,0,-0.18, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        
        ax.text(0.2,0.85, 'PRE', ha='center', transform=ax.transAxes    )
        ax.text(0.5,0.85, 'ON', ha='center', transform=ax.transAxes    )
        ax.text(0.8,0.85, 'POST', ha='center', transform=ax.transAxes    )



    elif what=='decoderscheme intime preonly':    # 2 
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.arrow(0.2,0.6,0,-0.18, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        
        ax.text(0.2,0.85, 'PRE', ha='center', transform=ax.transAxes    )



    elif what=='decoderscheme crosstime':        # 3
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.7,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.8,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.arrow(0.2,0.6,0.215,-0.215, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.5,0.6,0.215,-0.215, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        
        ax.text(0.2,0.85, 'PRE', ha='center', transform=ax.transAxes    )
        ax.text(0.5,0.85, 'PRE', ha='center', transform=ax.transAxes    )
        ax.text(0.8,0.85, 'ON', ha='center', transform=ax.transAxes    )



    elif what=='decoderscheme crossblock':      # 4
        ax.add_patch(plt.Rectangle((0.2,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.6,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.3,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.7,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.arrow(0.65,0.6,-0.26,-0.22, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )


        
        ax.text(0.3,0.85, 'initital,\ntransition', ha='center', transform=ax.transAxes    )
        ax.text(0.7,0.85, 'multi-\nmodal', ha='center', transform=ax.transAxes    )


    elif what=='decoderscheme crosstime prepre':      #  5
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.arrow(0.2,0.6,0.215,-0.215, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        
        ax.text(0.2,0.85, 'PRE', ha='center', transform=ax.transAxes    )
        ax.text(0.5,0.85, 'PRE', ha='center', transform=ax.transAxes    )


    elif what=='decoderscheme crosstime preon':       # 6
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.arrow(0.2,0.6,0.215,-0.215, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        
        ax.text(0.2,0.85, 'PRE', ha='center', transform=ax.transAxes    )
        ax.text(0.5,0.85, 'ON', ha='center', transform=ax.transAxes    )


    elif what=='decoderscheme crossblock aa+ai': # 7
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        c1,c2 = ['dodgerblue','navy']
        ax.arrow(0.2,0.6,0,-0.18, head_width=0.02, ec=c1, fc=c1, lw=lw, length_includes_head=True, transform=ax.transAxes, zorder=0    )
        ax.arrow(0.2,0.6,0.215,-0.215, head_width=0.02, ec=c2, fc=c2, lw=lw, length_includes_head=True, transform=ax.transAxes, zorder=0    )
        
        ax.text(0.2,0.85, 'attend   ', ha='center', transform=ax.transAxes    )
        ax.text(0.5,0.85, '   ignore', ha='center', transform=ax.transAxes    )
        
        cc = 'slategrey'   # comparator color
        ax.plot([0.2,0.2],[0.19,0.15],color=cc,lw=lw,transform=ax.transAxes)
        ax.add_patch(patches.Arc( (0.3,0.15),0.2,0.2,0,180,270, color=cc, lw=lw, transform=ax.transAxes))
        ax.arrow(0.3,0.05,0.14,0, head_width=0.02, ec=cc, fc=cc, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.5,0.19,0,-0.08, head_width=0.02, ec=cc, fc=cc, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.text(0.5,0.05,'-',color=cc, fontsize='x-large', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.add_patch(plt.Circle((0.5,0.05),0.04, ec=cc, fc=trp, lw=int(lw/2), transform=ax.transAxes))
        

    elif what=='decoderscheme withinblock aa+ii': # 8
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        c1,c2 = ['dodgerblue','navy']
        ax.arrow(0.2,0.6,0,-0.18, head_width=0.02, ec=c1, fc=c1, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.5,0.6,0,-0.18, head_width=0.02, ec=c2, fc=c2, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        
        ax.text(0.2,0.85, 'attend   ', ha='center', transform=ax.transAxes    )
        ax.text(0.5,0.85, '   ignore', ha='center', transform=ax.transAxes    )

        cc = 'slategrey'   # comparator color
        ax.plot([0.2,0.2],[0.19,0.15],color=cc,lw=lw,transform=ax.transAxes)
        ax.add_patch(patches.Arc( (0.3,0.15),0.2,0.2,0,180,270, color=cc, lw=lw, transform=ax.transAxes))
        ax.arrow(0.3,0.05,0.14,0, head_width=0.02, ec=cc, fc=cc, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.5,0.19,0,-0.08, head_width=0.02, ec=cc, fc=cc, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.text(0.5,0.05,'-',color=cc, fontsize='x-large', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.add_patch(plt.Circle((0.5,0.05),0.04, ec=cc, fc=trp, lw=int(lw/2), transform=ax.transAxes))


    elif what=='decoderscheme crossblock ii+ia': # 9
        ax.add_patch(plt.Rectangle((0.1,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.2,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc=trp, lw=lw, transform=ax.transAxes))
        c1,c2 = ['dodgerblue','navy']
        ax.arrow(0.5,0.6,-0.215,-0.215, head_width=0.02, ec=c1, fc=c1, lw=lw, length_includes_head=True, transform=ax.transAxes, zorder=0    )
        ax.arrow(0.5,0.6,0,-0.18, head_width=0.02, ec=c2, fc=c2, lw=lw, length_includes_head=True, transform=ax.transAxes, zorder=0    )
        
        ax.text(0.2,0.85, 'attend   ', ha='center', transform=ax.transAxes    )
        ax.text(0.5,0.85, '   ignore', ha='center', transform=ax.transAxes    )
        
        cc = 'slategrey'   # comparator color
        ax.plot([0.2,0.2],[0.19,0.15],color=cc,lw=lw,transform=ax.transAxes)
        ax.add_patch(patches.Arc( (0.3,0.15),0.2,0.2,0,180,270, color=cc, lw=lw, transform=ax.transAxes))
        ax.arrow(0.3,0.05,0.14,0, head_width=0.02, ec=cc, fc=cc, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.arrow(0.5,0.19,0,-0.08, head_width=0.02, ec=cc, fc=cc, lw=lw, length_includes_head=True, transform=ax.transAxes    )
        ax.text(0.5,0.05,'-',color=cc, fontsize='x-large', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.add_patch(plt.Circle((0.5,0.05),0.04, ec=cc, fc=trp, lw=int(lw/2), transform=ax.transAxes))


    elif what=='decoderscheme PCA space':        # 10   !!! (changed after crosscontext attends put in)
        ax.add_patch(plt.Rectangle((0.4,0.6),0.2,0.2, ec='k', fc='lightgrey', lw=lw, transform=ax.transAxes))
        ax.add_patch(plt.Circle((0.5,0.3),0.1, ec='k', fc='lightgrey', lw=lw, transform=ax.transAxes))
        ax.arrow(0.5,0.6,0,-0.18, head_width=0.02, ec='k', fc='k', lw=lw, length_includes_head=True, transform=ax.transAxes    )
        
        ax.text(0.5,0.85, 'ON', ha='center', transform=ax.transAxes    )







def drawsubspaces(axs,wix=None,what=None):
#    fig,axs = plt.subplots(1,1,figsize=(12,12))

    
    whats = ['dbnv','2dsubspace','pcasubspace','pcafullcomparedbnv','2dsubspacechoice']
    
    if wix!=None:
        what = whats[wix]
    
    if what == 'dbnv':
        axs.view_init(elev=25,azim=45)
        
        cl = 2
        lm = 1
        
        axs.plot([0,cl],[0,0],[0,0],color='black',linewidth=3,alpha=0.9)
        axs.plot([0,0],[0,cl],[0,0],color='black',linewidth=3,alpha=0.9)
        axs.plot([0,0],[0,0],[0,cl],color='black',linewidth=3,alpha=0.9)
        
    
        # rotate onto position in 3D
        ax = -30   # -15
        az = +40   # +15
        M = nedi.rotationmatrix3d(ax,0,az)
    
        # initialize plane object in 2D    
        R = np.array([[-1,-1,1,1],[1,-1,-1,1],[0,0,0,0]])*1.3
        R = M @ R
        Rs = R
        # draw object
        r = [list(zip(R[0],R[1],R[2]))]
        p = art3d.Poly3DCollection(r,alpha=0.3,edgecolor='teal',facecolor='teal')
        axs.add_collection3d(p)
    
        
        # points
        N = 14
        pd = 1
        s = 0.1
        R1 = np.array([np.random.randn(N)*s, np.random.randn(N)*s, np.random.randn(N)*s-pd])
        R2 = np.array([np.random.randn(N)*s, np.random.randn(N)*s, np.random.randn(N)*s+pd])
        # colors = ['mediumturquoise','darkcyan']
        colors = ['dodgerblue','red']
        for rx,Raux in enumerate([R1, R2]):
            R = M @ Raux 
            axs.plot(R[0],R[1],R[2],'o',color=colors[rx],alpha=0.7)
            axs.text(R[0].mean()-0.15,R[1].mean()+0.2,R[2].mean()+0.3-(1-rx)*0.7,'class %d'%(1-rx+1),(-1,1,0),color=colors[rx])
        

        # arrow
        # R = np.array([[0],[0],[0.7071]])
        R = np.array([0,0,0.9])
        R = M @ R
        Ra = R
        axs.quiver(0,0,0,R[0],R[1],R[2], color='teal', linewidth=4)
        sp = 0.7
        axs.plot([0, -R[0]], [0, -R[1]], [0, -R[2]],'--k',alpha=0.4,lw=1)
        axs.plot([-R[0]*sp, -R[0]], [-R[1]*sp, -R[1]], [-R[2]*sp,-R[2]],'--',color='grey',lw=2)

        
        # annotations
        axs.text(2.2,0,-0.15,'neuron #1','x')
        axs.text(0.04,1,-0.45,'neuron #2','y')
        axs.text(0,0,1.2,'neuron #3','z')
        axs.text2D(0.25,1,'activity\nspace',transform=axs.transAxes)
        axs.text(Rs[0,2]+0.35,Rs[1,2]-0.3,Rs[2,2]+0.1,'decision\nboundary',Rs[:,1]-Rs[:,2],color='teal',verticalalignment='bottom')
        axs.text(-0.1,0.1,0.1,'DV',Ra,color='teal')

        
        
        axs.set_xlabel('x')
        axs.set_ylabel('y')
        axs.set_zlabel('z')
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_zticks([])
        axs.set_xlim(-lm,lm)    
        axs.set_ylim(-lm,lm)    
        axs.set_zlim(-lm,lm)    
        axs.axis('off')
    
    
    
    
    elif what=='2dsubspace' or what=='2dsubspacechoice':
        if what=='2dsubspace': cx_pool=[0,1]
        elif what=='2dsubspacechoice': cx_pool=[2,1]
        
        axs.view_init(elev=25,azim=45)
        
        cl = 2
        lm = 1
        
        axs.plot([0,cl],[0,0],[0,0],color='black',linewidth=3,alpha=0.9)
        axs.plot([0,0],[0,cl],[0,0],color='black',linewidth=3,alpha=0.9)
        axs.plot([0,0],[0,0],[0,cl],color='black',linewidth=3,alpha=0.9)
        
        
        colors = ['mediumvioletred','navy','darkorange']
        taskaspects = ['context','visual','choice']
        # rotate onto position in 3D
        ax = +45
        az = +105
        M = nedi.rotationmatrix3d(ax,0,az)
        
        for cx in [0,1]:
        
            # arrows
            ml = 1/0.7071
            R = np.array([[0.7071*(cx)],[0.7071*(1-cx)+0.7071/5*(1-cx)],[0]])
            L = np.array([[0.7071*(cx)],[0.7071*(1-cx)],[0]]) * ml
            R = M @ R
            L = M @ L
            # axs.plot([0,L[0]],[0,L[1]],[0,L[2]],'--',color='darkgrey')
            axs.quiver(0,0,0,R[0],R[1],R[2], color=colors[cx_pool[cx]], linewidth=4)
            axs.text(0,0.3+0.1*(1-cx),-0.2+0.7*(1-cx),'%s DV'%taskaspects[cx_pool[cx]],R[:,0],color=colors[cx_pool[cx]])
        
        # subspace plane
        R = np.array([[-1,-1,1,1],[1,-1,-1,1],[0,0,0,0]])
        R = M @ R
        Rs = R
        r = [list(zip(R[0],R[1],R[2]))]
        p = art3d.Poly3DCollection(r,alpha=0.4,edgecolor='grey',facecolor='grey')
        axs.add_collection3d(p)
        
        axs.text2D(0.25,1,'activity\nspace',transform=axs.transAxes)
        axs.text(Rs[0,1]+0.3,Rs[1,1],Rs[2,1]+0.1+0.2*(1-cx),'DVs\' subspace',Rs[:,0]-Rs[:,1],color='grey',verticalalignment='bottom')
        l = 0.3
        # axs.text((1-l)*Rs[0,2]+l*R[0,3],(1-l)*Rs[1,2]+l*R[1,3]+0.3,(1-l)*Rs[2,2]+l*R[2,3],\
        #          'orthogonal\nbasis',Rs[:,2]-Rs[:,1],color='grey',verticalalignment='bottom')
        
        
        axs.set_xlim(-lm,lm)
        axs.set_ylim(-lm,lm)
        axs.set_zlim(-lm,lm)
        axs.axis('off')
    
    
    
    
    elif what=='pcasubspaces':

        axs.set_xlim(-lm,lm)    
        axs.set_ylim(-lm,lm)    
        axs.set_zlim(-lm,lm)    
        axs.axis('off')




    elif what=='pcafullcomparedbnv':
        axs.view_init(elev=25,azim=45)
        
        cl = 2
        lm = 1
        
        axs.plot([0,cl],[0,0],[0,0],color='black',linewidth=3,alpha=0.9)
        axs.plot([0,0],[0,cl],[0,0],color='black',linewidth=3,alpha=0.9)
        axs.plot([0,0],[0,0],[0,cl],color='black',linewidth=3,alpha=0.9)
        
        
        colors = ['mediumvioletred','navy']
        color = colors[0]
        taskaspects = ['context','visual']
        
        
        # original decision boundary normal vector
        # rotate onto position in 3D
        ax = +45
        az = +105
        M = nedi.rotationmatrix3d(ax,0,az)
        
        ml = 1/0.7071
        R = np.array([[0],[0.7071],[0]])
        L = R * ml
        R = M @ R
        Ra = R
        L = M @ L
        axs.plot([0,L[0]],[0,L[1]],[0,L[2]],'--',color='darkgrey') # show the full subspace spanned
        axs.quiver(0,0,0,R[0],R[1],R[2], color=color, linewidth=4)
        
        
        
        # define pca subspace
        ax = +10
        ay = +0
        az = -20
        M_p = nedi.rotationmatrix3d(ax,ay,az)
        
        # get the plane
        R = np.array([[-1,-1,1,1],[1,-1,-1,1],[0,0,0,0]])
        R = M_p @ R
        Rs = R
        r = [list(zip(R[0],R[1],R[2]))]
        p = art3d.Poly3DCollection(r,alpha=0.4,edgecolor='grey',facecolor='grey')
        axs.add_collection3d(p)
        
        # get the ellipse:
        C = plt.Circle((0,0),0.7071).get_verts()
        C = np.c_[C, np.zeros(len(C))]
        C[:,1] = C[:,1]/2
        C = C.T
        R = M_p @ C
        r = [list(zip(R[0],R[1],R[2]))]
        p = art3d.Poly3DCollection(r,edgecolor='grey',facecolor=(0,0,0,0),linestyle='--')
        axs.add_collection3d(p)
        
        # get the PCA defined decision boundary normal vector
        Rx = np.array([[-0.7071],[-0.1],[0]])
        Rx = M_p @ Rx
        R = Rx
        axs.quiver(0,0,0,R[0],R[1],R[2], color=color, linewidth=1)

        # get the x and y coordinate of the pca subspace
        # we need to get the unit vectors in the columns of Rb
        Rb = np.array([[1,0],[0,1],[0,0]])
        Rb = M_p @ Rb

        # get the shadow of the original dv to the pca subspace
        # we need the 
#        Rp = np.dot(Rb.T,L)
#        axs.plot([0,Rp[0]],[0,Rp[1]],[0,Rp[2]],'--',color=color)



        axs.text2D(0.25,1,'activity\nspace',transform=axs.transAxes)

        axs.text2D(0,0.58,'PCA$_{ON}$ subspace',transform=axs.transAxes,color='grey')
        axs.text(Rs[0,3],Rs[1,3]+0.4,Rs[2,3],'pc #1',Rs[:,0]-Rs[:,3],color='grey',verticalalignment='bottom')
        axs.text(Rs[0,3]+0.66,Rs[1,3]-0.2,Rs[2,3],'pc #2',Rs[:,2]-Rs[:,3],color='grey',verticalalignment='bottom')

        dm = 1.25
        axs.text(Ra[0,0]*dm+0.1,Ra[1,0]*dm-0.4,Ra[2,0]*dm,'%s DV\nin activity space'%taskaspects[0],Ra[:,0],color=color)
        axs.text(Rx[0,0]*dm,Rx[1,0]*dm+0.1,Rx[2,0]*dm-0.1,'%s DV\nin PCA$_{ON}$ subspace'%taskaspects[0],Rx[:,0],color=color)


        
        axs.set_xlim(-lm,lm)    
        axs.set_ylim(-lm,lm)    
        axs.set_zlim(-lm,lm)    
        axs.axis('off')
    






    return










def statshelper():

    recalculate = 0        # calucate and save stats if 1, load if 0

    skip = 20

    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting 
    # datanames = ['DT008','DT009']
    n_mice = len(datanames)


    taskaspects = ['visual','audio','context','choice']

    times = np.arange(601)[skip:-skip]

    if recalculate:
        #  ( mice, taskaspects, trajectories, classes,{mean,s.e.m.} )
        trajectory_matrix = np.zeros( (n_mice,4,len(times),2,3) )
        #  ( taskaspects, all_neurons, classes,{mean,s.e.m.} )
        stats_matrix = np.zeros( (4,0,5) )
        stats_matrix_celltypes = [ np.zeros( (4,0,5) ), np.zeros( (4,0,5) ) ]
        
        # firing rate stats:
        for n,dn in enumerate(datanames):
            block = preprocess.loaddatamouse(dn,T,continuous_method=continuous_method,normalize=False)        # use raw firing rates: normalzie=False
            n_neuron = block.segments[0].analogsignals[0].shape[1]
            blv,bla = preprocess.getorderattended(dn)
            comparisongroups  = [ \
                                    [ [ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                                    [ [ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                                    [ [  blv,  [],[] ],   [ bla,   [], [] ] ],\
                                    [ [] ]         ]
            local_stats = np.zeros((4,n_neuron,5))  # mean, variance, mean and var of trial-to-trial variance amongst neurons during stimulus
            local_stats_celltypes = [ np.zeros((4,np.sum(0==block.annotations['celltypes']),5)),
                                      np.zeros((4,np.sum(1==block.annotations['celltypes']),5)) ]
            for cx,comparison in enumerate(taskaspects):
            # collect neural responses
                if not comparison=='choice': # visual, audio, context:
                    acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx])
                else:  # choice:
                    acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)
        
                # trajectory_matrix[n,cx,:,cidx,:,0] =  np.array(acrossresponses[cidx])[:,:,:].mean(axis=0)
                # trajectory_matrix[n,cx,:,cidx,:,1] =  2*np.array(acrossresponses[cidx])[:,:,:].std(axis=0)/\
                #                             np.sqrt(len(acrossresponses[cidx]))
                for cidx in range(2):       # two classes
                    trajectory_matrix[n,cx,:,cidx,0] =  np.array(acrossresponses[cidx])[:,skip:-skip,:].mean(axis=(0,2))
                    trajectory_matrix[n,cx,:,cidx,1] =  2*np.array(acrossresponses[cidx])[:,skip:-skip,:].var(axis=(0,2))
                    trajectory_matrix[n,cx,:,cidx,2] =  2*trajectory_matrix[n,cx,:,cidx,1]/\
                                            np.sqrt(len(acrossresponses[cidx])*n_neuron)
                    

                local_stats[cx,:,0] = np.concatenate((np.array(acrossresponses[0])[:,skip:-skip,:],\
                                                      np.array(acrossresponses[1])[:,skip:-skip,:]),axis=0).mean(axis=(0,1))
                local_stats[cx,:,1] = np.concatenate((np.array(acrossresponses[0])[:,skip:-skip,:],\
                                                      np.array(acrossresponses[1])[:,skip:-skip,:]),axis=0).var(axis=(0,1))
                # mean firing rate of neurons, neural variance; average over time and trials; pres stim and during stim
                local_stats[cx,:,2] = np.concatenate((np.array(acrossresponses[0])[:,:T['stimstart_idx'],:],\
                                                      np.array(acrossresponses[1])[:,:T['stimstart_idx'],:]),axis=(0)).mean(axis=(0,1))
                local_stats[cx,:,3] = np.concatenate((np.array(acrossresponses[0])[:,T['stimstart_idx']:T['stimend_idx'],:],\
                                                      np.array(acrossresponses[1])[:,T['stimstart_idx']:T['stimend_idx'],:]),axis=(0)).mean(axis=(0,1))
                

                    
                # trial to trial variance of all neurons individually during stimulus:
                local_stats[cx,:,4] = np.concatenate((np.array(acrossresponses[0])[:,T['stimstart_idx']:T['stimend_idx'],:],\
                                                      np.array(acrossresponses[1])[:,T['stimstart_idx']:T['stimend_idx'],:]),axis=0).mean(axis=1).var(axis=0)
        
        
        
                for ct in [0,1]:
                    mask = ct==block.annotations['celltypes']
                    local_stats_celltypes[ct][cx,:,0] = np.concatenate((np.array(acrossresponses[0])[:,skip:-skip,mask],\
                                                          np.array(acrossresponses[1])[:,skip:-skip,mask]),axis=0).mean(axis=(0,1))
                    local_stats_celltypes[ct][cx,:,1] = np.concatenate((np.array(acrossresponses[0])[:,skip:-skip,mask],\
                                                          np.array(acrossresponses[1])[:,skip:-skip,mask]),axis=0).var(axis=(0,1))
                    # mean firing rate of neurons, neural variance; average over time and trials; pres stim and during stim
                    local_stats_celltypes[ct][cx,:,2] = np.concatenate((np.array(acrossresponses[0])[:,:T['stimstart_idx'],mask],\
                                                          np.array(acrossresponses[1])[:,:T['stimstart_idx'],mask]),axis=(0)).mean(axis=(0,1))
                    local_stats_celltypes[ct][cx,:,3] = np.concatenate((np.array(acrossresponses[0])[:,T['stimstart_idx']:T['stimend_idx'],mask],\
                                                          np.array(acrossresponses[1])[:,T['stimstart_idx']:T['stimend_idx'],mask]),axis=(0)).mean(axis=(0,1))
                    
    
                        
                    # trial to trial variance of all neurons individually during stimulus:
                    local_stats_celltypes[ct][cx,:,4] = np.concatenate((np.array(acrossresponses[0])[:,T['stimstart_idx']:T['stimend_idx'],mask],\
                                                          np.array(acrossresponses[1])[:,T['stimstart_idx']:T['stimend_idx'],mask]),axis=0).mean(axis=1).var(axis=0)
        
        
        
                                            
            stats_matrix = np.concatenate( (stats_matrix, local_stats), axis=1)

            for ct in [0,1]:
                if local_stats_celltypes[ct].shape[1]>0:
                    stats_matrix_celltypes[ct] = np.concatenate( (stats_matrix_celltypes[ct], local_stats_celltypes[ct]), axis=1)
        
        
        pickle.dump((trajectory_matrix,stats_matrix,stats_matrix_celltypes),open('../cache/phys/stats,trajectory_matrix-%s.pck'%(continuous_method),'wb'))
    
    else:
        trajectory_matrix,stats_matrix,stats_matrix_celltypes = pickle.load(open('../cache/phys/stats,trajectory_matrix-%s.pck'%(continuous_method),'rb'))



    
    
    # STATS for publication:
    cx = 0      # visual
    n_all_neurons = stats_matrix.shape[1]
    print('visual stim. mean firing rate from %4.2f +/- %4.2f to %4.2f +/- %4.2f and trial to trial variance on stimulus: %4.2f +/- %4.2f,   n neurons %d'%(\
              stats_matrix[cx,:,2].mean(), stats_matrix[cx,:,2].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[cx,:,3].mean(), stats_matrix[cx,:,3].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[cx,:,4].mean(),stats_matrix[cx,:,4].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix.shape[1]) )

    for ct in [0,1]:
        n_all_neurons_type = stats_matrix_celltypes[ct].shape[1]
        print('visual stim. %s neurons, mean firing rate from %4.2f +/- %4.2f to %4.2f +/- %4.2f and trial to trial variance on stimulus: %4.2f +/- %4.2f,   n neurons %d'%(\
              ['inhibitory','excitatory'][ct],
              stats_matrix_celltypes[ct][cx,:,2].mean(), stats_matrix_celltypes[ct][cx,:,2].std()*1/np.sqrt(n_all_neurons_type),\
              stats_matrix_celltypes[ct][cx,:,3].mean(), stats_matrix_celltypes[ct][cx,:,3].std()*1/np.sqrt(n_all_neurons_type),\
              stats_matrix_celltypes[ct][cx,:,4].mean(),stats_matrix_celltypes[ct][cx,:,4].std()*1/np.sqrt(n_all_neurons_type),\
              stats_matrix_celltypes[ct].shape[1]) )


        
    mask = ((stats_matrix[cx,:,3]-stats_matrix[cx,:,2])>0) # choose cells that are with positive firing rate change
    print('+ cells FR: visual stim. mean firing rate from %4.2f +/- %4.2f to %4.2f +/- %4.2f and trial to trial variance on stimulus: %4.2f +/- %4.2f,   n neurons %d'%(\
              stats_matrix[cx,mask,2].mean(), stats_matrix[cx,mask,2].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[cx,mask,3].mean(), stats_matrix[cx,mask,3].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[cx,mask,4].mean(),stats_matrix[cx,mask,4].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[:,mask,:].shape[1]) )
    mask = np.logical_not(mask)
    print('- cells FR: visual stim. mean firing rate from %4.2f +/- %4.2f to %4.2f +/- %4.2f and trial to trial variance on stimulus: %4.2f +/- %4.2f,   n neurons %d'%(\
              stats_matrix[cx,mask,2].mean(), stats_matrix[cx,mask,2].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[cx,mask,3].mean(), stats_matrix[cx,mask,3].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[cx,mask,4].mean(),stats_matrix[cx,mask,4].std()*1/np.sqrt(n_all_neurons),\
              stats_matrix[:,mask,:].shape[1]) )
    
        
    if 0:  # display to check
    
        fig, ax = plt.subplots(2,4,figsize=(4*6,2*6))
        for cx,comparison in enumerate(taskaspects):
            for cidx in range(2):       # two classes
                axs = ax[0,cx]
                axs.plot(times,trajectory_matrix[4,cx,:,cidx,0].T)
                axs.set_title(comparison)
                if cx==0: axs.set_ylabel('mean'); axs.legend(['class 1','class2'])
                axs = ax[1,cx]
                axs.plot(times,trajectory_matrix[4,cx,:,cidx,1].T)
                if cx==0: axs.set_ylabel('variance')
        fig.suptitle('across neurons variance')
    
    
        fig, ax = plt.subplots(3,4,figsize=(4*6,3*6))
        for cx,comparison in enumerate(taskaspects):
            axs = ax[0,cx]
    #        axs.bar(np.arange(stats_matrix.shape[1]),stats_matrix[cx,:,0])
            axs.hist(stats_matrix[cx,:,0],bins=20)
            axs.set_title(comparison)
            axs.set_xlabel('mean')
            axs = ax[1,cx]
    #        axs.bar(np.arange(stats_matrix.shape[1]),stats_matrix[cx,:,1])
            axs.hist(stats_matrix[cx,:,1],bins=20)
            axs.set_xlabel('variance')
            axs = ax[2,cx]
            axs.hist(stats_matrix[cx,:,0]/stats_matrix[cx,:,1],bins=20)
            axs.set_xlabel('variance/mean')
    
    
        
        fig.suptitle('across trial trajectory variance')
    
        
    
    return










 # FIGURES

def figure1():
    
    # this will be only 1D holding the behavioural sessions and
    #                   1E holding behavioural rate


    datanamesall = ['ME103', 'ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']
    datanames = ['ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']
    # datanames = ['DT009','DT017']

    dprimes_all = []
    L_max_sessiontypes = [0,0,0,0]
    for n,dn in enumerate(datanames):
        print(pathdatamouse + 'trainingbehaviour/' + dn + '.mat')
        # hf = h5py.File(pathdatamouse + 'trainingbehaviour/' + dn + '_hdf5.mat','r')       # we'have resaved in matlab for hdf5 format -v7.3
        # data = hf['BehavStruct']      # data contains all the sessions, visual, audio, then mixed combined...;  if h5py loading: use ...'][0,0] at the end
        data = sp.io.loadmat(pathdatamouse + 'trainingbehaviour/' + dn + '.mat')['BehavStruct'] # this method preservs the order unfortunately unlike hdf5
        # print(data.size)
        # print(data.shape)
        # print(data.dtype.names)
        
        sessionids = np.array(data.dtype.names)
        

        # collect session index list
        vidx = [ s for sx,s in enumerate(sessionids) if s[0]=='v' and s[:2]!='va']
        aidx = [ s for sx,s in enumerate(sessionids) if s[0]=='a' and s[:2]!='av']
        vaidx = [ s for sx,s in enumerate(sessionids) if s[:2]=='va' ]
        avidx = [ s for sx,s in enumerate(sessionids) if s[:2]=='av' ]

        
        sessiontypelabellist = [vidx,aidx,vaidx,avidx]
        
        dprime_sessiontypes = []
        for sllx,sll in enumerate(sessiontypelabellist): # go through each session type (single and multimodal)
            dprime_sessions = []
            L_max_sessiontypes[sllx] = np.max((L_max_sessiontypes[sllx],len(sll)))
            for sx,sl in enumerate(sll):
                # print('sessions: ',sl)
    
                # exclude all "NaN" trials (use valididx)
                if sl[0]=='v':
                    valididx = np.logical_not( np.isnan(data[sl][0][0].squeeze()[0][:,0]) )
                    go = (data[sl][0][0].squeeze()[0][valididx,0] == 45).astype('int')
                elif sl[0]=='a':
                    valididx = np.logical_not( np.isnan(data[sl][0][0].squeeze()[1][:,0]) )
                    go = (data[sl][0][0].squeeze()[1][valididx,0] == 5).astype('int')
                else:
                    print('error',sll)
        
                L = go.shape[0]
                # print('this is go: ', go)
        
                lick = np.zeros( L, dtype='int16' )       # lick data contains each trial empty array or an array with lick timings; "0" is no lick
                aux = data[sl][0][0].squeeze()[2][valididx]
                for l in range(L):
                    # print(aux)
                    # print(aux.shape)
                    # print(type(aux[l]), aux[l])
                    # if np.isnan(aux[l]): print(l, np.nan)
                    # print(l, aux[l][0])
                    
                    
                    # aux = sessiongroups[sx][2][0][l][0].ravel()
                    if len(aux[l][0])>0:
                        # print(l, aux[l][0][0,:])
                        if len(aux[l][0][0,:])>0:
                            lick[l] = aux[l][0][0,-1]>=2.            # if the animals licks at least once after 2 seconds, it is considered a licked trial, and recorded as "1"
                        # else: print(l,' no lick')
                    # else: print(l,' NaN')
                if len(lick)!=L: print('missing data'); continue     # if there are some issues with the data, then skip
                

                # d′ = Z(hit rate) − Z(false alarm rate),
                # where function Z(p), p ∈ [0,1], is the inverse of the cumulative distribution function of the Gaussian distribution.
                n_go = np.sum(go)
                n_nogo = np.sum(1-go)
                hit,miss,corrrej,fal = np.array([ np.sum(go & lick) / n_go, np.sum( go & (1-lick) ) / n_go,\
                                         np.sum( (1-go) & (1-lick) ) / n_nogo, np.sum( (1-go) & lick ) / n_nogo ])
                
                # also one needs to be careful, because if rate is exactly 1 or 0, the cdf is infinite
                h = sp.stats.norm.ppf( hit )
                if hit==1: h = sp.stats.norm.ppf( 0.99 )
                elif hit==0: h = sp.stats.norm.ppf( 0.01 )
                f = sp.stats.norm.ppf( fal )
                if fal==0: f = sp.stats.norm.ppf( 0.01 )
                elif fal==1: f = sp.stats.norm.ppf( 0.99 )
                
                dprime = h - f
                
                # if dprime==-np.inf or dprime==np.inf:
                #     print(dn,sl,h,f,'  <>   ',hit,fal,dprime)
                    
                
                dprime_sessions.append(dprime)
                
                # print(hit,miss,corrrej,fal,' <>  ', h,f,' <>  ', dprime)
                # print(np.sum(go),np.sum(np.isnan(go)),np.sum(lick),np.sum(np.isnan(go)))
                # return
            # print(len(dprime_sessions))
            dprime_sessiontypes.append(dprime_sessions)
        # print(len(dprime_sessiontypes))
        dprimes_all.append(dprime_sessiontypes)






    fig = plt.figure(figsize=(3*8,1*8))

    # fig = plt.figure(constrained_layout=False,figsize=(4*8,6.4*8))
    # grid spec will be:
    #   bulk:     top    left raster, right firing rates
    #          bottom    TCA panels
    #   bulk dimensions:
    #                        4              6
    #                        10
    
    gs = gridspec.GridSpec(1, 3, figure=fig)
    
    ax = [ fig.add_subplot(gs[0:2]), fig.add_subplot(gs[2]) ]



    # fig,ax = plt.subplots(2,1,figsize=(3*6,2*6))
    panel = 'D'
    colors = ['dodgerblue','olivedrab','navy','darkgreen']
    taskcontextlabels = ['visual only','audio only','combined, attend visual','combined, attend audio']
    xc = 0
    xt = []
    xl = []
    for stx in [0,1,2,3]:
        lm = L_max_sessiontypes[stx]
        container = np.nan*np.ones((lm,len(datanames)))
        # axs = ax[min(stx,2)]
        axs = ax[0]
        
        for n,dn in enumerate(datanames):
            ls = len(dprimes_all[n][stx])
            container[-ls:,n] = dprimes_all[n][stx]
        means = np.nanmean(container,axis=1)
        if stx>1:
            container=container[8:,:]
            means=means[8:]
        
        x = xc + np.arange(len(means)) + 1
        xt.append(x)
        if stx<2: xc = x[-1]
        xl.extend( np.arange(len(means)) + 1 )
        
        axs.plot( x, container,'o', color=colors[stx], alpha=0.6)
        axs.plot( x, means, lw=4, color=colors[stx],label=taskcontextlabels[stx])

    axs.legend(frameon=False)
    axs.set_ylim(-0.5,4)
    axs.set_yticks([0,1,2,3])
    axs.set_xticks(np.arange(1,len(xl)))
    axs.set_xticklabels([])
    axs.set_xlim(0,36)
    axs.set_xlabel('training sessions')
    axs.set_ylabel('D` of behaving mice')       # if stx==0: 
    figs.plottoaxis_chancelevel(axs,ch=1.75)
        

    figs.invisibleaxes(axs,which=['top'])
    figs.invisibleaxes(axs,which=['right'])
    
    figs.labelaxis(axs,panel,x=-0.05,y=1.02)





    # 1E beahviour correct/error fraction
    panel = 'E'
    n_mice = len(datanamesall)

    B = np.zeros((n_mice,2))
    for n,dn in enumerate(datanamesall):
        (hit,miss,correctrejection,falsealarm) = preprocess.assessbehaviouralperformance(dn,modality='all',multimodalonly=True)
        B[n,0] = len(miss)/(len(hit)+len(miss))
        B[n,1] = len(falsealarm)/(len(correctrejection)+len(falsealarm))

    axs = ax[1]

    axs.bar(x=np.arange(n_mice)-0.15,height=B[:,0],width=0.3,color='darkorange',label='miss')
    axs.bar(x=np.arange(n_mice)+0.15,height=B[:,1],width=0.3,color='red',label='false alarm')
    axs.legend(frameon=False)
    axs.set_ylabel('error fraction')
    axs.set_xlabel('mice')
    axs.set_xticks([])
    # axs.set_xticks(np.arange(n_mice))
    # axs.set_xticklabels(datanamesall,rotation=45)

    figs.invisibleaxes(axs,which=['top'])
    figs.invisibleaxes(axs,which=['right'])

    
    figs.labelaxis(axs,panel,x=-0.15,y=1.02)

  



    fig.tight_layout()



    save = 0 or globalsave
    if save:
        # ext = '.pdf'
        # ext = '.png'
        # fig.savefig(resultpath+'Fig1D,E_training,behaviordprime'+ext)
        fig.savefig(resultpath+'Fig1D,E_training,behaviordprime-horizontal'+ext)
        # fig.savefig(resultpath+'Fig1D,E_training,behaviordprime,extended.png')

















def figure2():


    # raster    
    dn = 'DT019'
    block = preprocess.loaddatamouse(dn,T,continuous_method=continuous_method,normalize=False)        # use raw firing rates: normalize=False
    n_neuron = block.segments[0].analogsignals[0].shape[1]
    times = block.segments[0].analogsignals[0].times
    
    taskaspects = ['visual','audio','context','choice']
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [ \
                            [ [ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                            [ [ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                            [ [  blv,  [],[] ],   [ bla,   [], [] ] ],\
                            [ [] ]         ]



        




    
    # raster    
    dn = 'DT019'
    block = preprocess.loaddatamouse(dn,T,continuous_method=continuous_method,normalize=False)        # use raw firing rates: normalzie=False
    n_neuron = block.segments[0].analogsignals[0].shape[1]
    times = block.segments[0].analogsignals[0].times
    
    taskaspects = ['visual','audio','context','choice']
    blv,bla = preprocess.getorderattended(dn)
    comparisongroups  = [ \
                            [ [ [2,4],[45],    [] ], [ [2,4],[135],     [] ] ],\
                            [ [ [2,4],  [],[5000] ], [ [2,4],   [],[10000] ]],\
                            [ [  blv,  [],[] ],   [ bla,   [], [] ] ],\
                            [ [] ]         ]


        
    variablecolors = ['navy','darkgreen','mediumvioletred','orange']
    classlabels = [['45°','135°'],['5kHz','10kHz'],['attend visual','attend audio'],['lick','withhold lick']]




    
    # cherry picks
    # first find the neurons whose activity best describes the first and the second class
    nct = len(taskaspects)
    powers = np.empty( (n_neuron,nct) )                           # cells x taskaspects
    bestcells = np.empty( (n_neuron,nct), dtype='int16' )       # cells x taskaspects
    c_db = []
    for cx,comparison in enumerate(taskaspects):
        acrossdecoder = pickle.load(open('../cache/subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
        wx = int((len(acrossdecoder)-7)/n_neuron)
        c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
        
    c_db = np.array(c_db) # [comparisongroup,neurons,trajectory,stats]
    c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=2)  # average over stimulus timeinterval    # this is a comparison group by neuron by   stats matrix
#    norms_c_db = np.linalg.norm(c_db_means,axis=0)
    powers = c_db_means[:,:,0].T
    for cx in range(nct):
        bestcells[:,cx] = np.argsort(powers[:,cx])[::-1]      # use first and last neuron.  do this for all task aspects
    selected_cell_ids = np.zeros((nct,2))        # best for each class        taskaspects   x    
    selected_cell_ids = bestcells[[0, -1],:].T          # find the cells that respond best for each variable (i.e. taskaspect)
    # selected_cell_ids[3,:] = bestcells[[2,-3],3]       # choice should be separate from visual and audio, so find the third best
    selected_cell_ids[1,:] = bestcells[[1,-2],1]       # audio separate from visual and choice, so find the nextest bestest
    selected_cell_ids[2,:] = bestcells[[1,-2],2]       # context separate from visual and choice, so find the nextest bestest
    print('cell ids selected:',selected_cell_ids)
    trajectory_matrix = np.zeros( (nct,nct+1,len(times),2,2,2) )           # taskaspects   x     selected neurons+all      x        trajectories        x     classes            x mean,s.e.m.
    
    # now find the trajectories for the conditions for each variable for the best neurons
    for cx,comparison in enumerate(taskaspects):
    # collect neural responses
        if not comparison=='choice': # visual, audio, context:
            acrossresponses = preprocess.collect_stimulusspecificresponses(block,comparisongroups[cx])
        else:  # choice:
            acrossresponses = preprocess.collect_stimulusspecificresponses_choice(block,dn)

        for bx in range(2):          # two selected neurons
            for nx,n_id in enumerate(selected_cell_ids[:,bx]):
                for cidx in range(2):       # two classes
                    trajectory_matrix[cx,nx,:,cidx,bx,0] =  np.array(acrossresponses[cidx])[:,:,n_id].mean(axis=0)
                    trajectory_matrix[cx,nx,:,cidx,bx,1] =  2*np.array(acrossresponses[cidx])[:,:,n_id].std(axis=0)/\
                                    np.sqrt(len(acrossresponses[cidx]))
        for cidx in range(2):       # two classes
            trajectory_matrix[cx,4,:,cidx,0,0] =  np.array(acrossresponses[cidx])[:,:,:].mean(axis=(0,2))
            trajectory_matrix[cx,4,:,cidx,0,1] =  2*np.array(acrossresponses[cidx])[:,:,:].std(axis=(0,2))/\
                                    np.sqrt(len(acrossresponses[cidx])*acrossresponses[cidx][0].shape[1])

    final_selection = np.array([0,0,0,0],dtype='int16')
    scs = np.array([ selected_cell_ids[k,final_selection[k]] for k in range(4) ],dtype='int16')
    # selected cells (cell class 1-2 by taskaspects): 2, 1, 1, 2
    # ordered cell ids:
    scs_labels = ['a','b','c','d']
    oscs = n_neuron-1-scs     # since kilosort2 we have the units ordered by depth



    
    # cherry picked neural trial averaged firing rates (or PSTHs?) responding to task-variables one or other values

    variablecolors = ['navy','darkgreen','mediumvioletred','orange']
    classlabels = [['45°','135°'],['5kHz','10kHz'],['attend visual','attend audio'],['lick','withhold lick']]
    
    if 0:     # explore all possible cells for cherry picks
        fig,ax = plt.subplots(8,4,figsize=(32,32))
        for cx,comparison in enumerate(taskaspects):
            for nx,cellsensitivity in enumerate(taskaspects):
                for bx in range(2):
                    axs = ax[nx*2+bx,cx]            # variables times cells matrix
                    for  cidx in [1,0]:
                        axs.plot(times, trajectory_matrix[cx,nx,:,cidx,bx,0],linewidth=2,color=variablecolors[cx],alpha=0.5+0.5*cidx)
                    axs.set_xlim(T['starttime']+200*pq.ms,T['endtime']-200*pq.ms)
                    axs.set_xticks([0,3000])
                    axs.set_ylim(0,np.max(trajectory_matrix[:,:,:,:,:,0].ravel())*1.1)         #axs.get_ylim()[1]
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    figs.plottoaxis_chancelevel(axs,0)
                    if nx==0 and bx==0: axs.set_title(comparison+' variable')
                    if cx==0:
                        if bx==0: axs.set_ylabel(cellsensitivity+' sensitive cells\ncell %d sesnsitive to\n%s'%(bx+1,classlabels[cx][bx]))
                        else: axs.set_ylabel('cell %d response to\n%s'%(bx+1,classlabels[cx][bx]))
                    if nx==3 and bx==1: axs.set_xlabel('time from stimulus onset [ms]')
                    axs.legend(classlabels[cx])
                
        fig.suptitle(dn+'\nEstimated instantenous firing rates [Hz] of selected neurons;\n neurons sensitive to task variables, active in the two values of these variables, selection criteria by decoder coefficients')

        return





    # load PCA signal variance:
    C_pc_signal_variance, n_neurons_pca = pickle.load(open('../cache/subspaces/noise+signal-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'rb'))













    # FIGURE PLOT





    fig = plt.figure(constrained_layout=False,figsize=(4*8,6.4*8))
    # grid spec will be:
    #   bulk:     top    left raster, right firing rates
    #          bottom    TCA panels
    #   bulk dimensions:
    #                        4              6
    #                        10
    
    gs = gridspec.GridSpec((3+4+4)*3+4+2*3, 1, figure=fig)
    # [3*3,1,4*3,1,1*3,1,4*3]
    gsras = gs[0:8].subgridspec(1, 1)
    gsfr = gs[10:22].subgridspec(4, 4)
    # gsfra = gs[22:25].subgridspec(1, 4)
    gspca = gs[24:33].subgridspec(3, 4)
    gssignal = gs[35:40].subgridspec(1, 4)      # we will need one lowest panel row for signal variance early pc saturation
    
    axras = fig.add_subplot(gsras[0])
    axfr = np.empty((5,4),dtype=object)
    axpca = np.empty((3,4),dtype=object)
    axsignal = np.empty((4),dtype=object)
    
    for i in range(4):
        for j in range(4):
            axfr[i,j] = fig.add_subplot(gsfr[i,j])
    # for j in range(4):
    #     axfr[4,j] = fig.add_subplot(gsfra[j])
    for i in range(3):
        for j in range(4):
            axpca[i,j] = fig.add_subplot(gspca[i,j])
    for j in range(4):
        axsignal[j] = fig.add_subplot(gssignal[0,j])





    # raster for multiple subsequent trials
    panel = 'A'

    stimlabels = [' 45°  5 kHz','135° 10 kHz',' 45°  5 kHz',' 45° 10 kHz','135°  5 kHz']
    # raster plot
    axs = axras
    trxstart = 47      # csv 49th line: 47


    # print('segments:',len(block.segments),'neurons in spiketrain',len(block.segments[0].spiketrains))
    # print('trials',trxstart+1)
    # print([  len(block.segments[trxstart+t].spiketrains[0])  for t in [0,1,2,3,4] ])
    # j = 0
    # for st0,st1,st2,st3,st4 in zip(block.segments[trxstart].spiketrains,block.segments[trxstart+1].spiketrains,\
    #                                          block.segments[trxstart+2].spiketrains,block.segments[trxstart+3].spiketrains,\
    #                                          block.segments[trxstart+4].spiketrains):
    #     print(np.r_[st0, st1+6*pq.s, st2+12*pq.s, st3+18*pq.s, st4+24*pq.s]*pq.ms)


    sp =       [ neo.SpikeTrain(   times=np.r_[st0, st1+6*pq.s, st2+12*pq.s, st3+18*pq.s, st4+24*pq.s]*pq.ms,\
                                t_start=-1500*pq.ms, t_stop = 5*6000*pq.ms   ) \
                  for st0,st1,st2,st3,st4 in zip(block.segments[trxstart].spiketrains,block.segments[trxstart+1].spiketrains,\
                                             block.segments[trxstart+2].spiketrains,block.segments[trxstart+3].spiketrains,\
                                             block.segments[trxstart+4].spiketrains) ]
    # ordered display
    # osp = [ sp[-i] for s in sp ]
    osp = sp         # we already have units sorted ascending as depth increases in kilosort2
    osp.reverse()
    
    # axs.eventplot(osp,colors=celltypecolorlist[::-1],linelengths=celltypelinelengths[::-1])#,linewidths=celltypelinewidths[::-1])
    axs.eventplot(osp,colors='black')
    # axs.invert_yaxis()

#    axs.invert_yaxis()   # instead this, which does not work, we reorder with [::-1] and len(sp)-n
    # axs.plot([-1500,-500],[101.5,101.5],'k')
    # axs.text(-1480,102.1,'1 s')
#    axs.plot([-1500,-500],[-0.5,-0.5],'k')
#    axs.text(-1480,-1.1,'1 s')
    
    # for nx,n in enumerate(len(sp)-oscs):
    for nx,n in enumerate(oscs):
        axs.plot([-1800,-1550],[n,n],'k')
#        if nx==2: axs.text(-2500,n,'%s'%(n+1))
        axs.text(-2500,n-0.5,'%s'%(scs_labels[nx]))
    
    # print('sp len',len(sp))
    # ue = len(sp)-1
    # axs.plot([-1800,-1550],[ue,ue],'dodgerblue')
    # axs.text(-2500,ue-0.5,'%s'%('!'),color='dodgerblue')
    
    axs.set_ylim(-1,n_neuron)
    axs.set_xlim(-1900,5*6000-1500)
    for trtx,trt in enumerate([0,6,12,18,24]):
        TR = {'stimstarttime':0+trt*1000, 'stimendtime':3000+trt*1000}
        figs.plottoaxis_stimulusoverlay(axs,TR)
        axs.text(trt*1000+50,len(osp)+2.1,'trial %d: '%(trxstart+trtx)+stimlabels[trtx])
    axs.set_ylim(-1,n_neuron+2)
    axs.axis('off')

    # axs.text(-2300,105,panel,fontsize=24,fontweight='bold')
    figs.labelaxis(axras,panel,-0.025,1.0125)





    
    # firing rates of selected cells
    panel = 'B'
    
# firing rates of selected cells
    
    for cx,comparison in enumerate(taskaspects):
        for nx in range(len(taskaspects)):
            bx = final_selection[nx]
            # cellsensitivity = taskaspects[nx] + ' sensitive cell\n[Hz]'#%(scs_labels[nx])
            cellsensitivity = 'cell %s\nfiring rate [Hz]'%(scs_labels[nx])
            axs = axfr[nx,cx]            # variables times cells matrix
            for  cidx in [0,1]:          # response to given class
                axs.plot(times, trajectory_matrix[cx,nx,:,cidx,bx,0],linewidth=2,color=variablecolors[cx],alpha=1.-2./3.*cidx)
                axs.fill_between(times,trajectory_matrix[cx,nx,:,cidx,bx,0]-trajectory_matrix[cx,nx,:,cidx,bx,1],\
                                        trajectory_matrix[cx,nx,:,cidx,bx,0]+trajectory_matrix[cx,nx,:,cidx,bx,1],\
                                  color=variablecolors[cx],alpha=(1.-2./3.*cidx)/2.)
            axs.set_xlim(T['starttime']+200*pq.ms,T['endtime']-200*pq.ms)
            axs.set_xticks([0,3000])
            axs.set_yticks([0,20])
            axs.set_ylim(axs.get_ylim())
            # axs.set_ylim(0,30)#    np.max(trajectory_matrix[:,:,:,:,:,0].ravel())*1.1)
            # figs.plottoaxis_notickmarks(axs)
            if nx==0 and cx==0:
                # axs.plot([3300, 3800],[16,16],'k')
                # axs.plot([3300, 3300],[16,21],'k')
                # axs.text(3350,13.5,'500 ms',fontsize='small')
                # axs.text(3100,17.5,'5 Hz',rotation=90,fontsize='small')
                axs.text(0.02,0.9,'2 s.e.m.', fontsize='x-small',transform=axs.transAxes)
            figs.plottoaxis_stimulusoverlay(axs,T)
#            figs.plottoaxis_chancelevel(axs,0)
            # axs.spines['left'].set_visible(False)
            axs.spines['right'].set_visible(False)
            axs.spines['top'].set_visible(False)
            # axs.spines['bottom'].set_visible(False)
            if nx==0:
                axs.legend(classlabels[cx],loc='upper right',frameon=False)#    [ 'response to '+classlabels[cx][i] for i in range(2) ],loc=(0.5,0.75) ,frameon=False)
                axs.set_title(comparison,pad=30)
            if nx==len(taskaspects)-1:
                axs.set_xlabel('time [ms]')
#                axs.text(-2000,32,panels[cx],fontsize=24,fontweight='bold')
            # if nx==0 and cx==0:
            #     figs.labelaxis(axs,panels[0])
            # if nx==4 and cx==0:
            #     figs.labelaxis(axs,panels[1])
            #     axs.plot([-800, -800],[6,7],'k')
            #     axs.text(-1000,6.2,'1 Hz',rotation=90,fontsize='small')
            # if nx==4:
            #     axs.set_ylim(3.5,7.5)#
            #     axs.set_yticks([4,5,6,7]); axs.set_yticklabels(['4','','','7'])
            #     axs.spines['left'].set_visible(True)
                   
            if cx==0:
                axs.set_ylabel(cellsensitivity,labelpad=20)
                
            if cx==nx:
                # pos = [pos.x0 - 0.4, pos.y0 - 0.3 ,  pos.width + 0.8, pos.height + 0.6]
                # axf = fig.add_axes(pos)
                axs.add_patch(plt.Rectangle((-0.2,-0.2),1.4,1.4, ec=variablecolors[cx], fc=variablecolors[cx], alpha=0.1, transform=axs.transAxes))
                
                # outergs = gridspec.GridSpec(1, 1)
                # outergs.update(bottom=(i//2)*.47+0.01,left=(i%2)*.5+0.02, 
                #    top=(1+i//2)*.47-0.01,  right=(1+i%2)*.5-0.02)
                # outerax = fig.add_subplot(outergs[0])
                 
                
                # axf.set_facecolor('lightgrey')
                # axf.set_axis('off')
                
    figs.labelaxis(axfr[0,0],panel)
            






    # trajectories projected onto pca axis
    panel = 'C'
    taskcolors = [ ['navy','darkgreen','purple','saddlebrown'],\
                   ['steelblue','c','mediumvioletred','orange'] ]
    classnames = [['45°','135°'],['5kHz','10kHz'],['attend visual','attend audio'],['lick','withhold lick']]

    n_pc = 3
    projected_dynamics,distances,times = pickle.load(open('../cache/subspaces/subspacedynamics,3D-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'rb'))
        
        

    pcmaxs = 0*np.ones(n_pc)
    pcmins = 0*np.ones(n_pc)

    alpha = 0.8
    # time parametrization
    skip_idx = 20
    # t_all = times[skip_idx:-skip_idx]
    t_all = times
    t_0 = times[skip_idx:T['stimstart_idx']+1]
    t_1 = times[T['stimstart_idx']:T['stimend_idx']+1]
    t_oo = times[T['stimend_idx']:-skip_idx]

    
    # fig,ax = plt.subplots(n_pc+3,len(taskaspects),figsize=(len(taskaspects)*8,(n_pc+3)*8))
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
            
            
            
            # 1D 1 PC
            for px in range(n_pc):

                axs = axpca[px,cx]
                

                # gather for identical axes limits
                pcmaxs[px] = np.max(( pcmaxs[px], np.max(trial[skip_idx:-skip_idx,px])   ))
                pcmins[px] = np.min(( pcmins[px], np.min(trial[skip_idx:-skip_idx,px])   ))


                axs.plot( times, trial[:,px], lw=2,color=variablecolors[cx],alpha=1.-2./3.*aix,label=classnames[cx][aix] )

                axs.fill_between(times, trial[:,px]-trial_e[:,px], trial[:,px]+trial_e[:,px], \
                                 color=variablecolors[cx],alpha=(1.-2./3.*aix)/2)



                figs.setxt(axs)
                axs.set_xlim(skip_idx*T['dt']+T['starttime'], -skip_idx*T['dt']+T['endtime'])
                

                if aix==1:
                    if px==0: axs.legend(frameon=False)
                if cx==0: axs.set_ylabel('PC %d'%(px+1))
                
                    
                

            axs.set_xlabel('time [ms]')


    pcmins *=1.1
    pcmaxs *=1.1
    for cx,comparison in enumerate(taskaspects):
        for px in range(n_pc):
            axs = axpca[px,cx]
            axs.set_ylim(pcmins[px],pcmaxs[px])
            figs.plottoaxis_chancelevel(axs)
            figs.plottoaxis_stimulusoverlay(axs,T)
            

    

    figs.labelaxis(axpca[0,0],panel)


    



    # Panels for PCA signal variance
    panel = 'D'
    n_pc_display = n_neurons_pca
    max_v = 0
    for cx,comparison in enumerate(taskaspects):
        axs = axsignal[cx]
        # C_signal = C_pc_signal_variance[cx,:n_pc_display,1]/(C_pc_signal_variance[cx,:n_pc_display,1].sum())
        # C_signal = C_pc_signal_variance[cx,:n_pc_display,1]/(C_pc_signal_variance[cx,0,2])
        C_signal = C_pc_signal_variance[cx,:n_pc_display,1].cumsum()/(C_pc_signal_variance[cx,:n_pc_display,:2].sum(axis=1).cumsum()) # to date
        
        max_v = np.max([np.max(C_signal),max_v])

        bottoms = 0 # will hold previous bar tops for stacked bars

        axs.plot(np.arange(n_pc_display)+1,C_signal,'o-',lw=2,color=variablecolors[cx],label=comparison+' signal\ncumulative variance proportion')
        figs.plottoaxis_chancelevel(axs,ch=1)
        
        axs.set_xlabel('# PC projections')

        axs.legend(frameon=False,loc='lower right')
        
        axs.set_xticks([1,10,20,31])
        if cx==0:
            axs.set_ylabel('relative cumulative\nsignal strength')
        
    for cx,comparison in enumerate(taskaspects):
        axs = axsignal[cx]
        figs.invisibleaxes(axs)

        axs.set_ylim(0,max_v*1.05)
        # axs.set_ylim(0,0.2)


    # show total
    axsignal_inset = axsignal[0].inset_axes([0.3,0.7,0.6,0.3],transform=axsignal[0].transAxes)
    axsignal_inset.plot(np.arange(n_pc_display)+1, C_pc_signal_variance[cx,:n_pc_display,1].cumsum(), 'o-', color=variablecolors[0])
    axsignal_inset_total = axsignal_inset.twinx()
    axsignal_inset_total.plot(np.arange(n_pc_display)+1, C_pc_signal_variance[cx,:n_pc_display,:2].sum(axis=1).cumsum(),'o-',color='grey')
    
    axsignal_inset.tick_params(axis='y', labelcolor=variablecolors[0])
    axsignal_inset_total.tick_params(axis='y', labelcolor='grey')

    for i in [0,1]:
        axsignal_inset.text(10,0.5-i*0.2,['cum. signal var.','cum. total var.'][i],color=[variablecolors[0],'grey'][i],fontsize='small')

    axsignal_inset.set_xticks([1,10,20,31])

    figs.invisibleaxes(axsignal_inset,which=['top'])
    figs.invisibleaxes(axsignal_inset_total,which=['top'])
    
    
    
    figs.labelaxis(axsignal[0],panel)




   









    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Fig2_raster,singlecells,pcaimages'+ext)




























def figure3():

    variablecolors = ['navy','mediumvioletred']    
    
    
    # load data
    
    
    
    # A  decoders
    dn = 'DT019'
    examine = 'allexpcond'
    comparison = 'visual'
    acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
    visualtimecourse = acrossdecoder[:2]
    comparison = 'context'
    acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
    contexttimecourse = acrossdecoder[:2]
    

    # get baseline shuffle:
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    # get the coefficients
    c_db = []
    for cx,comparison in enumerate(['visual','context']):
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
        n_neuron = 31    # DT019
        wx = int((len(acrossdecoder)-7)/n_neuron)
        c_db.append(  np.reshape(np.array(acrossdecoder[7:]), (wx,n_neuron,acrossdecoder[7].shape[0],acrossdecoder[7].shape[1]) ).mean(axis=0)    )
    c_db = np.array(c_db) # [comparisongroup,neurons,trajectory,stats]
    c_db_means = np.array(c_db)[:,:,T['stimstart_idx']:T['stimend_idx'],:].mean(axis=2)  # average over stimulus timeinterval    # this is a comparison group by neuron by   stats matrix
    c_db_order = np.argsort(c_db_means[0,:,0])




    
    # B,C, Fx, Gx +behav x axis
    examine = 'allexpcond'
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032'] # with JRC
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032'] # with ks2
    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting 
    
    timegroups = []
    collect_stats_v = []
    collect_stats_c = []
    for n,dn in enumerate(datanames):
        visualtimegroups = []; contexttimegroups = []
#        timestarts_idx = int((np.arange(0,6001,1500)*pq.ms / np.array(T['dt'])).magnitude)     # cutpoints
        timestarts_idx = np.arange(0,601,150,dtype='int16')
        comparison = 'visual'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        collect_stats_v.append([ acrossdecoder[1][:,0].mean(), acrossdecoder[1][:,2].mean() ])
        for ti in range(4):
            visualtimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )
        comparison = 'context'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        collect_stats_c.append([ acrossdecoder[1][:,0].mean(), acrossdecoder[1][:,2].mean() ])
        
        for ti in range(4):
            contexttimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )
        timegroups.append([visualtimegroups,contexttimegroups])
    timegroups = np.array(timegroups)


    # STATS section
    # print('visual and context means and 2 sem')    
    # print(timegroups.mean(axis=0))
    # print(timegroups.std(axis=0)*2/np.sqrt(len(datanames)))

    # collect_stats_m = np.array(collect_stats_v).mean(axis=0)
    # collect_stats_sem = np.array(collect_stats_v).std(axis=0)*2/np.sqrt(len(datanames))
    # print('visual mean total accuracy: %4.2f (%4.2f) +/- %4.2f (%4.2f)'%(collect_stats_m[0],collect_stats_sem[0],collect_stats_m[1],collect_stats_sem[1]))
    # collect_stats_m = np.array(collect_stats_c).mean(axis=0)
    # collect_stats_sem = np.array(collect_stats_c).std(axis=0)*2/np.sqrt(len(datanames))
    # print('context mean total accuracy: %4.2f (%4.2f) +/- %4.2f (%4.2f)'%(collect_stats_m[0],collect_stats_sem[0],collect_stats_m[1],collect_stats_sem[1]))
    # print('n animals %d'%len(datanames))

    # return



    


    # L       # behavioural performance
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']       # with JRC spike sorting
    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting
    
    n_mice = len(datanames)
    
    Xa = []
    Ya = []
    perf = np.zeros((n_mice,6))
    expertise = np.zeros(n_mice)

    for n,dn in enumerate(datanames):
        print(dn)
        a,b,b_s,_,_,_ = preprocess.loadtrainingbehaviouraldata(dn)       # a: all data, b: behav, b_s behav sem.
        perf[n,:] = [ b[1],b[3], b[5], 2*b_s[1], 2*b_s[3], 2*b_s[5] ]      # behaviour and its s.e.m.


        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        expertise[n] = acrossdecoder[1][ T['start_idx']:T['stimstart_idx'], 0 ].mean()

        Ya.append( [a[1],a[3], np.concatenate( [a[1],a[3]] ) ] )
        Xa.append( [ expertise[n]*np.ones(len(Ya[-1][0])), expertise[n]*np.ones(len(Ya[-1][1])), expertise[n]*np.ones(len(Ya[-1][2])) ] )


    Yl = [ np.concatenate( [  Ya[n][d] for n in range(n_mice)  ] ) for d in range(3) ]
    Xl = [ np.concatenate( [  Xa[n][d] for n in range(n_mice)  ] ) for d in range(3) ]







    # I,J        # subspace projections
    # datanamesefg = ['ME108','ME110','ME112','ME113','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']             # with JRC
    datanamesefg = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032'] # with ks2

    
    # old projections
#     orthooderlabel = 'cxvich'
#     dbnbasiscomps_all = []
#     projections_all = []
# #    projections_vi_all = []
# #    projections_cx_all = []
#     perfidx_all = []
#     for n,dn in enumerate(datanamesefg):
#         dbnbasiscomps,projections,activitycolors,basiscolors,basisvectors,perfidx =\
#             pickle.load(open('../cache/subspaces/dbprojections+projections,%s-%s_%s_%s-%s-%dms%dms_%s.pck'%('visual',orthooderlabel,examine,'all',continuous_method,T['dt'].magnitude,T['bin'].magnitude,dn),'rb'))
#         (ixgohit,ixnogocorrrej,ixgomiss,ixnogofal) = perfidx
#         dbnbasiscomps_all.append(dbnbasiscomps)
#         projections_all.append(projections)
# #        projections_vi_all.append(projections_vi)
# #        projections_cx_all.append(projections_cx)
#         perfidx_all.append(perfidx)
        
    activitycolors = ['darkgreen','darkred','lime','orangered']



    
    # new projections: visual is parallel context is semi: show on orthonormalized closest to it
    # dbnvcoords is visual,context (and choice which is not interesting for this plot)
    depth_idx = 150  # averaging onto 1500 ms into stimulus onset
    # projected_dynamics_all = [] # this will be (mice)(taskaspects,classes,[trials,timecourse,dbnvcoords])
    projections_all = []   # this will be (mice)(taskaspects)(classes)(trials,dbnvcoords)
    basis_all = []
    for n,dn in enumerate(datanamesefg):
        projected_dynamics, basis = pickle.load(open('../cache/subspaces/subspacedynamics,projected+dbnv,vicxch-%s_%s-%dms.pck'%(dn,continuous_method,T['dt'].magnitude),'rb'))
        projected_dynamics = np.array(projected_dynamics)
        projection = [ np.array(projected_dynamics[k])[:,T['stimstart_idx']:T['stimstart_idx']+depth_idx,:].mean(1) for k in [0,1,2,3]]
        projections_all.append(projection)
        
        basis_all.append(basis)








    # dbnv angles
    angles_all = []
    angles_highres_all = []
    for n,dn in enumerate(datanames):
        angles = pickle.load(open('../cache/subspaces/angles,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_all.append(angles)
        angles_highres = pickle.load(open('../cache/subspaces/angles,highres,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_highres_all.append(angles_highres)

    angles_all = np.array(angles_all)

    times_angle = np.arange(0,angles_all.shape[3])*T['dt'] + T['offsettime']





















    #        FIGURE

#    fig = plt.figure(num=2)
#    fig.clf()
#    fig, ax = plt.subplots(5,4,figsize=(36,45))

    # the substructure is compplicated due to the joint distribution axis
    ratio = 7        # ratio for marginal histogram subplots
    f1n=3; f2n=4     # full span of panel grid  vertical and horizontal
    res = 18 # resolution multiplier;    all panels can wiggle 4 directions
    marg = 2
    sizemul = (10/2)/(res/2/marg)
    
    fig = plt.figure(constrained_layout=False,figsize=(f2n*8*sizemul,f1n*8*sizemul))
    gs = fig.add_gridspec(f1n*res, f2n*res)
    ax = []
    ax_sides = []
    for j in np.arange(0,f1n*res,res):
        axa=[]
        for k in np.arange(0,f2n*res,res):
            # default grid:
            jr = slice(j+marg,j+res-marg)
            kr = slice(k+marg,k+res-marg)
            
            # handle big spaces between groups:
            if k==0*res:
                kr = slice(k,k+res-2*marg)      # 1st column
            # if j==1*res and k>0*res:
            #     jr = slice(j+marg-2,j+res-marg-2)      # visual and context rows closer
            if j==3*res:                   # last column
                jr = slice(j+marg+3,j+res-marg+3)
                
            # handle special axis issues
            if j==2*res and k==1*res:      # this is for a joint and marginal histograms of context vs. visual
                gs_marginals = gs[jr,kr].subgridspec(ratio+1, ratio+1)
                ax_joint = fig.add_subplot(gs_marginals[1:, :-1])
                ax_marginals = [ fig.add_subplot(gs_marginals[0 , :-1], sharex=ax_joint), \
                                 fig.add_subplot(gs_marginals[1:,  -1], sharey=ax_joint) ]     #, sharex,sharey=ax_joint)   ]
                axa.append( ax_joint )
            elif j==2*res and k==2*res:     # this is for the angle distribution
                axa.append( fig.add_subplot(gs[jr,kr], projection='polar') )
            elif j in [0*res, 2*res] and k==0*res:   # this is for the cartoons
                axa.append( fig.add_subplot(gs[jr,kr], projection='3d') )
            else:
                axa.append( fig.add_subplot(gs[jr, kr]) )
        ax.append(axa)
    ax = np.array(ax)

    print('figure 3, axis.shape',ax.shape)














    # A,B schematics


    panel = 'A'
    axs = ax[0,0]
    drawsubspaces(axs,0)
    figs.labelaxis(axs,panel,D2=True)
    

    panel = 'H'
    axs = ax[2,0]
    drawsubspaces(axs,1)
    figs.labelaxis(axs,panel,D2=True)
            

    # panel = 'K'   # this was the previous PCA schematics
    # axs = ax[3,0]
    # drawsubspaces(axs,3)
    # figs.labelaxis(axs,panel,D2=True)




    # B,E    # decoder trajectories
    panels = ['B','E']

    for bx,vartimecourse in enumerate([visualtimecourse,contexttimecourse]):
    
        axs = ax[bx,1]
        figs.plottoaxis_decoderrocauc(vartimecourse,axs,colorlist=['',variablecolors[bx]],plottrain=False)       # plot the test performance
        
        # axs.plot(vartimecourse[1].times,shuffle_m,color='red')
        # axs.fill_between(vartimecourse[1].times,shuffle_m-shuffle_e,shuffle_m+shuffle_e,color='red',alpha=0.3)
        
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,chances[dn])
        figs.setxt(axs)
        axs.set_yticks([0.5,1.0])# axs.set_yticklabels([0.5,1.0])
    #    axs.set_xlabel('[ms]')
        axs.set_ylabel(['visual','context'][bx]+' accuracy                ')
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)
    
        # if bx==0:    
        #     axins = axs.inset_axes([-0.42,0.55, 0.4,0.4],transform=axs.transAxes)
        #     drawschematics(axins,1); axins.axis('off')
        
        figs.labelaxis(axs,panels[bx])


#     axs = ax[1,1]
#     figs.plottoaxis_decoderrocauc(contexttimecourse,axs,colorlist=['',variablecolors[1]],plottrain=False)       # plot the test performance
#     figs.plottoaxis_stimulusoverlay(axs,T)
#     figs.plottoaxis_chancelevel(axs,0.5)
#     figs.setxt(axs)
#     axs.set_yticks([0.5,1.0])# axs.set_yticklabels([0.5,1.0])
# #    axs.set_xlabel('[ms]')
#     axs.set_ylabel('context'+' accuracy               ')
#     axs.spines['right'].set_visible(False)
#     axs.spines['top'].set_visible(False)
    
#     figs.labelaxis(axs,panels[2])







    # C,F
    timecourselabels = ['PRE','early ON','late ON','POST']
    panels = ['C','F']
    for bx in range(2):
        axs = ax[0+bx,1+1]        # first visual, below context
        ci = timegroups[:,bx,:].std(axis=0)*2/np.sqrt(len(datanames))
        m = timegroups[:,bx,:].mean(axis=0)
        ci = np.c_[m-ci,m+ci]
#        print(timegroups.shape,m.shape,ci.shape)
        artist = axs.boxplot(x=timegroups[:,bx,:], positions=-750+timestarts_idx[:4]*10, notch=True,usermedians=m,conf_intervals=ci,\
                             whis=[5,95],labels=timecourselabels,widths=350)
    
#        print(artist.keys())
        for element in artist.keys():
            plt.setp(artist[element], color=variablecolors[bx],linewidth=2)
    
        axs.set_yticks([0.5,1.0])
        axs.set_ylim(0.45,1.01)
        axs.set_xlim(-1500,4500)
        figs.plottoaxis_stimulusoverlay(axs,T)
        m_c = np.array(list(chances.values())).mean()
        e_c = np.array(list(chances.values())).std()/np.sqrt(len(datanamesefg))
        figs.plottoaxis_chancelevel(axs,m_c+e_c)
#        axs.set_ylabel('accuracy')
        figs.labelaxis(axs,panels[bx])
#            axs.set_title('n=%d'%14)
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)







    # coefficients
    panels = ['D','G']

    classnames = [['45°','135°'],['attend\nvisual','attend\naudio']]

    for cx,comparison in enumerate(['visual','context']):
        axs = ax[cx,3]
        axs.barh(y=(np.arange(n_neuron)+1),width=c_db_means[cx,c_db_order,0],color=variablecolors[cx] )
        axs.set_yticks([1,10,20,31])
        # axs.invert_yaxis()
        axs.set_xlim(-0.5,0.5)
        axs.set_xticks([-0.5,0.5])
        # axs.set_xticklabels([classnames[cx][0],'',classnames[cx][1]])
        axs.text(-0.4,25,classnames[cx][0])
        axs.text( 0.3,25,classnames[cx][1])

        axs.set_xlabel('neuron weights [SU]',labelpad=-10)
        axs.set_ylabel('neuron id, visual-ordered')
     
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)

        figs.labelaxis(axs,panels[cx])














    # SUBSPACES
    panel = 'I'
    
    n_showmice = datanamesefg.index('DT019')            # DT014 
    basisaspects = ['visual','context']
    basisaspectcolors = ['navy','mediumvioletred']   # np.array(taskcolors[0])[np.array([0,2,3],dtype=np.int16)]
    # basisaspectcolors = ['dodgerblue','fuchsia']   # np.array(taskcolors[0])[np.array([0,2,3],dtype=np.int16)]

    projections = projections_all[n_showmice]
    basis = basis_all[n_showmice]
    # ixgohit,ixnogocorrrej,ixgomiss,ixnogofal = perfidx_all[n_showmice]


    axs = ax[2,1]

    # plot the trial averages as points (av45,av135,iv45,iv135)
    for k in [0,1,2,3]:
        axs.plot(projections[k][:,0], projections[k][:,1], 'o',color=activitycolors[k],alpha=0.8)

    # basis vectors
    x0 = -1.6
    y0 = -2.0
    for bx,basisaspect in enumerate(basisaspects):
        # axs.plot([0,B[bx,0]],[0,B[bx,1]],lw=3,color=basisaspectcolors[bx])
        axs.plot([x0+0,x0+basis[0,bx]],[y0+0,y0+basis[1,bx]],lw=3,color=basisaspectcolors[bx])
    
    axins = axs.inset_axes([-0.2,0.4, 1,1],transform=axs.transAxes)
    axins.axis('off')
    xl0 = -0.2; yl0 = 0.78
    m = 0.08
    for j in range(2):
        for k in range(2):
            axins.add_patch(plt.Rectangle((xl0+m*j,yl0+m*k),m,m,\
              ec=activitycolors[j+2*k],fc=activitycolors[j+2*k],alpha=0.15,transform=axs.transAxes))
            axins.add_patch(plt.Circle((xl0+m/2+m*j,yl0+m/2+m*k),0.015,\
              ec=activitycolors[j+2*k],fc=activitycolors[j+2*k],transform=axs.transAxes))
    axins.text(xl0+m,  yl0+m*2,'45°',ha='right',va='bottom',transform=axs.transAxes)
    axins.text(xl0+m,  yl0+m*2,'135°',ha='left',va='bottom',transform=axs.transAxes)
#    axins.text(xl0+m*2+0.025, yl0+m*3/2,'ignore',ha='left',va='center',transform=axs.transAxes)
#    axins.text(xl0+m*2+0.025, yl0+m/2,'attend',ha='left',va='center',transform=axs.transAxes)
    axins.text(xl0-0.025, yl0+m*3/2,'ignore',ha='right',va='center',transform=axs.transAxes)
    axins.text(xl0-0.025, yl0+m/2,'attend',ha='right',va='center',transform=axs.transAxes)

    
    axs.set_xlim(axs.get_xlim())
    axs.set_ylim(axs.get_ylim())    

    figs.plottoaxis_crosshair(axs)

    
    axs.text(0.04, -0.04,'visual DV [SU]', ha='left',   va='center', transform=axs.transAxes)
    axs.text(-0.04, 0.04,'context DV [SU]', ha='center', va='bottom',rotation=90, transform=axs.transAxes)




    # plot the cov ellipsoids over the per trial activities
    for k in [0,1,2,3]:
        figs.confidence_ellipse(projections[k][:,0], projections[k][:,1], axs, n_std=2.0, facecolor=activitycolors[k],alpha=0.15)



    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    
    axs.axis('off')


    for k in range(2):
        # ax_marginals[0].hist(np.r_[projections[k][:,0],projections[2+k][:,0]],bins=15,color=['green','red'][k],alpha=0.7)
        # ax_marginals[1].hist(np.r_[projections[2*k][:,0],projections[2*k+1][:,0]],bins=15,color=['black','lightgrey'][k],alpha=0.7, orientation='horizontal')


        kde = sp.stats.gaussian_kde(np.r_[projections[k][:,0],projections[2+k][:,0]])
        x = np.arange(-2.4,+3.2+0.01,0.02)
        ax_marginals[0].plot(x,kde(x),color=['green','red'][k],lw=2,alpha=0.7)
        if k==1: ax_marginals[0].plot(x,np.zeros(len(x)),'k--',alpha=0.5)
        
        kde = sp.stats.gaussian_kde(np.r_[projections[2*k][:,1],projections[2*k+1][:,1]])
        x = np.arange(-2.4,+3.2+0.01,0.02)
        ax_marginals[1].plot(kde(x),x,color=['black','lightgrey'][k],lw=2,alpha=0.7)
        if k==1: ax_marginals[1].plot(np.zeros(len(x)),x,'k--',alpha=0.5)


    for margaxx in range(2):
        ax_marginals[margaxx].spines['right'].set_visible(False)
        ax_marginals[margaxx].spines['top'].set_visible(False)
        ax_marginals[margaxx].spines['left'].set_visible(False)
        ax_marginals[margaxx].spines['bottom'].set_visible(False)
        ax_marginals[margaxx].get_xaxis().set_visible(False)
        ax_marginals[margaxx].get_yaxis().set_visible(False)
   
        figs.plottoaxis_crosshair(ax_marginals[margaxx])

    figs.labelaxis(ax_marginals[0],panel,x=-0.45)


















    panel = 'J'

    # polar histogram of bases visual and context in all anumals
    
    angles = np.zeros((len(datanamesefg)))    #  {bv vi,ch}  x  context x mice
    # edges = np.linspace(-np.pi/2,np.pi/2,12+1)
    edges = np.linspace(0,np.pi,24+1)
    
    width=(edges[1]-edges[0])/2
    n_bins = len(edges)-1
    anglecounts = np.zeros((n_bins,2))        # {basisvectors vi,ch}   x context
    
    for n,dn in enumerate(datanamesefg):
        # basis = basis_all[n]
        # x = basis[0,1]
        # y = basis[1,1]
        # angles[n] = np.arctan(y/x)
        angles[n] = angles_all[n,0,2,T['stimstart_idx']:T['stimstart_idx']+150].mean()/180*np.pi
    aux = np.histogram(angles,bins=edges,density=False)
    anglecounts = aux[0]
    print(angles)
    print(anglecounts)
    # anglestats = [  angles.mean(axis=2), angles.std(axis=2)*2/np.sqrt(len(datanamesefg))   ]
    # anglestats_t_p = [ sp.stats.ttest_ind( angles[chvix,0,:], angles[chvix,1,:] )[1] for chvix in range(2) ]

    

    color = 'rebeccapurple'
    axs = ax[2,2]
    axs.bar(edges[:-1]+width, anglecounts, width=width*2,color=color, alpha=0.7)
#                axs.plot(edges[:-1]+width, anglecounts[:,chvix,cx],'o-',\
#                        color=basiscolors[basisindices[chvix][cx]],alpha=0.7)
    # axs.errorbar(anglestats[0],9,xerr=anglestats[1][chvix,cx],color=colors[chvix][cx])
    # axs.plot(anglestats[0][chvix,cx],9,'o',color=colors[chvix][cx])
    # axs.text(anglestats[0][chvix,:].mean(),9,'p=%5.3f'%anglestats_t_p[chvix])
    axs.legend(['angle between visual\nand context DVs'])
    # anglestats_t_p
    # axs.set_xlim(-np.pi/2,np.pi/2)
    axs.set_xlim(0,np.pi)
    axs.set_ylim(0,10)
    axs.set_xticks(edges[::6])
    axs.set_yticklabels([])
    # axs.set_xlabel(['visual basis','choice basis'])
    
    figs.labelaxis(axs,panel)









    # dynamics of dbnv angles between context and visual along the trial, all mice
    panel = 'K'
    
    pair = [0,2]
    axs = ax[2,3]
    singlecolor = 'rebeccapurple'
    
    for n,dn in enumerate(datanames):
        x = neph.smooth(angles_all[n,pair[0],pair[1]],kernelwidth=6,mode='same')
        x[:T['stimstart_idx']] = np.nan

        if n==0: axs.plot(times_angle,x,lw=0.8,color=singlecolor,alpha=0.2,label='single mice')
        else: axs.plot(times_angle,x,lw=0.8,color=singlecolor,alpha=0.2)

    m = angles_all[:,pair[0],pair[1]].mean(axis=0)
    e = angles_all[:,pair[0],pair[1]].std(axis=0)/np.sqrt(len(datanames))
    
    m[:T['stimstart_idx']] = np.nan
    
    axs.plot(times_angle,m,color=singlecolor,lw=2,label='mean of %d mice'%len(datanames))
    axs.fill_between(times_angle,m-2*e,m+2*e,color=singlecolor,alpha=0.2)#,label='2 s.e.m.')


    
    axs.legend(frameon=False)
    axs.set_xlim(T['offsettime']+200*pq.ms,T['endtime']-200*pq.ms)
    figs.setxt(axs)
    axs.set_yticks([0,45,90,135,180])
    axs.set_ylim(0,180)
    figs.plottoaxis_stimulusoverlay(axs,T)
    figs.plottoaxis_chancelevel(axs,90)
    
    axs.set_title('visual to context DV',)
    axs.set_ylabel('angle [deg]')
    # axs.set_xlabel('time from stimulus onset')

    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)    








# delete empty axes:
    for axs in [ax[1,0]]:    axs.axis('off')



    
#    fig.suptitle('V1 context dependent activity with identical stimuli,'+\
#                 'A decoder timecourse, B dec. all mice\nC expertise vs. representation, C orientation selectivity vs. decoder decision coefficients,\n'+\
#                 'D crosstime decoding from trained in pre, 1 mouse, E all mice, F expertise vs. performance\n'+\
#                 'G cx vs. vi decision normals DT018, H same w/all mice, G distrib of cx and vi basis angles w/all mice\n'+\
#                 'M decoder performance on activity projected onto principal components during spontaneous and evoked trial region, '+\
#                 'N explained variance in firing rate by task variables, glm')


    # fig.tight_layout()

    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Fig3_context,orthogonal'+ext)




























#       FIG 4                 GRADIENTS

def figure4():

    

  # B,C, Fx, Gx +behav x axis
    examine = 'allexpcond'
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032'] # with JRC
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032'] # with ks2
    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting 
    n_mice = len(datanames)


    # get baseline shuffle:
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    timegroups = []
    timegroups_se = []
    for n,dn in enumerate(datanames):
        visualtimegroups = []; contexttimegroups = []
        visualtimegroups_se = []; contexttimegroups_se = []
        # timestarts_idx = int((np.arange(0,6001,1500)*pq.ms / np.array(T['dt'])).magnitude)     # cutpoints
        timestarts_idx = np.arange(0,601,150,dtype='int16')
        comparison = 'visual'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        for ti in range(4):
            visualtimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )
            visualtimegroups_se.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 2 ].mean()   )
        comparison = 'context'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        for ti in range(4):
            contexttimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )
            contexttimegroups_se.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 2 ].mean()   )
        timegroups.append([visualtimegroups,contexttimegroups])
        timegroups_se.append([visualtimegroups_se,contexttimegroups_se])
    timegroups = np.array(timegroups)
    timegroups_se = np.array(timegroups_se)



    crossdecoder_matrix_allaspects_all = []   #   (mice)(taskaspects)(trainingpoints)(testtimecourse,stats)
    for n,dn in enumerate(datanames):
        crossdecoder_matrix_allaspects = pickle.load(open('../cache/continuous/crossdecodematrix-%s-%dms_%s.pck'%(continuous_method,T['dt'],dn),'rb'))
        crossdecoder_matrix_allaspects_all.append(crossdecoder_matrix_allaspects)
    crossdecoder_matrix_allaspects_all = np.array(crossdecoder_matrix_allaspects_all)
    sh = crossdecoder_matrix_allaspects_all.shape



    # this collects crossdecoder tests:
    if 1:
        fdr,bdr = pickle.load(open('../cache/continuous/crossdecodematrix,decayconstant-%s-%dms.pck'%(continuous_method,T['dt']),'rb'))
        # will have dimensions: (mice,taskaspects,trainingpoints,{tau,multiplier,additive})
        # fdr[:,:,:,0] = np.log(2) / fdr[:,:,:,0]
        fdr[:,:,:,0] = 100 * fdr[:,:,:,0]    # to make per 100 ms
        bdr[:,:,:,0] = 100 * bdr[:,:,:,0]    # to make per 100 ms


        if 0:           # use this to just calculate stats, print, and leave and close the door
            print('stats, crossdecoding decay rate:')
            print('steady state stimulus context, visual comparison')
            spix = ((np.array([1000,2000])*pq.ms-T['starttime'])/T['dt']).astype(np.int32) # stimulusperiod indices to consider for steady state
            print(spix)
            v = fdr[:,0,spix[0]:spix[1],0].mean(axis=1)
            c = fdr[:,2,spix[0]:spix[1],0].mean(axis=1)
            _,p = sp.stats.ttest_rel(v,c)
            print('n=%d, p=%4.4f'%(len(v),p))
            
            print('pre stimulus context')
            spix = ((np.array([-1500,-500])*pq.ms-T['starttime'])/T['dt']).astype(np.int32) # stimulusperiod indices to consider for steady state
            print(spix)
            c = fdr[:,2,spix[0]:spix[1],0].mean(axis=1)
            print('n=%d, m=%4.4f+/-%4.4f'%(len(c),c.mean(),c.std()/np.sqrt(len(c))))

            print('stimulus onset context')
            spix = ((np.array([-200,200])*pq.ms-T['starttime'])/T['dt']).astype(np.int32) # stimulusperiod indices to consider for steady state
            print(spix)
            c = fdr[:,2,spix[0]:spix[1],0].mean(axis=1)
            print('n=%d, m=%4.4f+/-%4.4f'%(len(c),c.mean(),c.std()/np.sqrt(len(c))))


            return
    else:
            
    
    
        # forward decay rate;  modeled as     exp(-t/tau - c) + B    =   A * exp(-t/tau) + B,      A = exp(-c)
        # log(a[index:index+window])
        # fig,axs = plt.subplots(sh[1],sh[0],figsize=(n_mice*8,sh[1]*8))
        fdr_window = 50
        fdr_skip_window=0
        fdr = np.zeros( (sh[0],sh[1],sh[2],3) )     # (mice,taskaspects,trainingpoints,{tau,multiplier,additive})
        bdr = np.zeros( (sh[0],sh[1],sh[2],3) )     # (mice,taskaspects,trainingpoints,{tau,multiplier,additive})
        for n in range(sh[0]):
            for cx in range(sh[1]):
                for ix in range(sh[2]-fdr_window-fdr_skip_window):
                    #fdr[n,cx,ix] = neph.decay(crossdecoder_matrix_allaspects_all[n,cx,ix,ix:ix+fdr_window,0])
                    m,c,b = neph.decay(crossdecoder_matrix_allaspects_all[n,cx,ix,ix+fdr_skip_window:ix+fdr_skip_window+fdr_window,0])
                    tau = -m
                    # A = np.exp(-c)
                    A = c
                    fdr[n,cx,ix,:] = (tau,A,b)
                    # do backward for reverse approaching stimulus entry
                    m,c,b = neph.decay(crossdecoder_matrix_allaspects_all[n,cx,-ix+fdr_window,-(ix)+fdr_window+fdr_skip_window:-1:(-ix)+fdr_skip_window,0])
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
        # fdr = fdr*100
    
        # fdr[:,0,:T['stimstart_idx']] = np.nan    # make visual pre not available for plot
        # norm with width:
        # fdr = fdr / T['dt'] *100  #  where tau is 1/ 1ms    and here we take decay per 100 ms
        # return


    
    # dbnv angles
    angles_all = []
    angles_highres_all = []
    for n,dn in enumerate(datanames):
        angles = pickle.load(open('../cache/subspaces/angles,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_all.append(angles)
        angles_highres = pickle.load(open('../cache/subspaces/angles,highres,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_highres_all.append(angles_highres)

    angles_all = np.array(angles_all)

    times_angle = np.arange(0,angles_all.shape[3])*T['dt'] + T['offsettime']



    if 0:   #  stats for DV angles, time shifted cross test accuracy stats
        H = np.array(angles_highres_all)   # (mice,taskvariable,shift,reference)
        print('stats, DV angles')

        cx = 2   # we deal with context only here
        band = (100*pq.ms/T['dt']).astype(np.int32)            # leave out this much around the diagonal
        i1 = T['stimstart_idx']       # indices
        i2 = T['stimend_idx']
        
        n_mice = H.shape[0]
        
        r = np.triu(H[:,cx,:i1,:i1],k=band).reshape(n_mice,-1)            # rest
        s = np.triu(H[:,cx,i1:i2,i1:i2],k=band).reshape(n_mice,-1)        # stimulus
        
        c = H[:,cx,:i1,i1:i2]
        c[:,:band,:band] = 0
        c = c.reshape(n_mice,-1)
        
        R = []; S = []; C = []
        
        for n in range(n_mice):
            R.append(r[n,:][r[n,:]>0])
            S.append(s[n,:][s[n,:]>0])
            C.append(c[n,:][c[n,:]>0])
        R = np.array(R); S = np.array(S); C = np.array(C); 
        # plt.hist(R[7,:],bins=100,alpha=0.3)
        # plt.hist(S[7,:],bins=100,alpha=0.3)
        # plt.hist(C[7,:],bins=100,alpha=0.8)
        
        
        mR = R.mean(axis=1)
        mS = S.mean(axis=1)
        mC = C.mean(axis=1)
        
        print('  pre: %4.4f+/-%4.4f'%(mR.mean(),mR.std()/np.sqrt(n_mice)))
        print('   on: %4.4f+/-%4.4f'%(mS.mean(),mS.std()/np.sqrt(n_mice)))
        print('cross: %4.4f+/-%4.4f'%(mC.mean(),mC.std()/np.sqrt(n_mice)))
        print('%4.1f+/-%4.1f°, σ=%4.1f°'%(mR.mean(),mR.std()/np.sqrt(n_mice),mR.std()))
        print('%4.1f+/-%4.1f°, σ=%4.1f°'%(mS.mean(),mS.std()/np.sqrt(n_mice),mS.std()))
        print('%4.1f+/-%4.1f°, σ=%4.1f°'%(mC.mean(),mC.std()/np.sqrt(n_mice),mC.std()))
        return









    # Figure
    
    fig,ax = plt.subplots(4,3,figsize=(3*8,4*8))
    
    
    
    
   # A
    panel = 'A'
    colors = ['navy','mediumvioletred']
    xlabels = ['visual PRE','mean context accuracy\npre stimulus']
    ylabels = ['visual early ON','mean context accuracy\n on stimulus']
    bx = 1
    axs = ax[0,0]           # top context, bottom visual
    x = timegroups[:,bx,0]
    y = timegroups[:,bx,1]
    x_se = timegroups_se[:,bx,0]
    y_se = timegroups_se[:,bx,1]
    axs.scatter(x[0],y[0],s=150,marker='o',color=colors[bx],alpha=0.8,label='indiv. mice')
    axs.scatter(x[1:],y[1:],s=150,marker='o',color=colors[bx],alpha=0.8)
    for n in range(len(x)):
        axs.errorbar(x[n],y[n],x_se[n],y_se[n],color='grey',alpha=0.8)
    axs.set_xlim(0.45,1.01)
    axs.set_ylim(0.45,1.01)
    axs.set_xticks([0.5,1.0])
    axs.set_yticks([0.5,1.0])

    m_c = np.array(list(chances.values())).mean()
    e_c = np.array(list(chances.values())).std()/np.sqrt(len(datanames))
    figs.plottoaxis_crosshair(axs,m_c+e_c,m_c+e_c)
    axs.set_xlabel(xlabels[bx])
    axs.set_ylabel(ylabels[bx])

    l = sp.stats.linregress(x,y)
    if l[3]<0.05:
        line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
        axs.plot(line,l[0]*line+l[1],color=colors[bx],linewidth=2,label='linear fit')
    else:
        line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
        axs.plot(line,l[0]*line+l[1],'--',color=colors[bx],linewidth=2)
        
    xoffs = 0.2
#        if sx in [3,4,5,7]: xoffs = 25
    if l[3]<0.001:
        axs.text(line[1]-xoffs,0.47,'p<%5.3f, $R^2$=%4.2f'%((0.001,l[2]**2)),color=colors[bx])
    else:
        axs.text(line[1]-xoffs,0.47,'p=%5.3f, $R^2$=%4.2f'%((l[3],l[2]**2)),color=colors[bx])
    
    axs.legend(frameon=False,loc=2)
    
    figs.labelaxis(axs,panel)
#            axs.set_title('n=%d'%14)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)    
    


    # B, time-decay presence, visual and context, from accuracy crosstest decays
    panel = 'B'


    dn = 'DT019'    
    n_cherry = datanames.index(dn)
    
    trainingpoints = np.array([-750,750])*pq.ms
    trainingpoints_idx = ((trainingpoints-T['offsettime'])/T['dt']).astype(np.int16)
    taskaspects = ['visual','context']
    times = np.arange(-1500,4400,10)*pq.ms

    for tx,trainingpoint_idx in enumerate(trainingpoints_idx):
        axs = ax[0,1+tx]
        for cix,taskix in enumerate([0,2]):
            if tx==0 and cix==0: continue
            x = neph.smooth(crossdecoder_matrix_allaspects_all[n_cherry,taskix,trainingpoint_idx,:,0],kernelwidth=5,mode='same')
            e = neph.smooth(crossdecoder_matrix_allaspects_all[n_cherry,taskix,trainingpoint_idx,:,2],kernelwidth=5,mode='same')
            if taskix==0:
                x[:T['stimstart_idx']] = np.nan   # visual pre is nonexistent
                e[:T['stimstart_idx']] = np.nan   # visual pre is nonexistent
            axs.plot(times,x,lw=2,color=colors[cix],label='%s accuracy'%taskaspects[cix])
            axs.fill_between(times,x-e,x+e,color=colors[cix],alpha=0.3)
    
        axs.plot([trainingpoints[tx],trainingpoints[tx]],[0.5,1.0],lw=1,color='black',label='trained')
    
        axs.set_xlim(-1400,4400)
        axs.set_ylim(0.5,1.0)

        axs.legend(frameon=False,loc=1)
        axs.set_title('trained at %s stimulus'%['pre','on'][tx])
        axs.set_ylabel('test accuracy')
        if tx==0: axs.set_xlabel('time from stimulus onset')

        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,chances[dn])
        
        axs.set_yticks([0.5,1.0])


        if tx==0: figs.labelaxis(axs,panel)

        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)    







    
    # 2D accuracy test matrices    
    panels = ['C','F']

    dn = 'DT019'    
    n_cherry = datanames.index(dn)
    
    taskaspects = ['visual','context']

    # cmaps = ['PuBu','PuRd']
    mapres = np.arange(0.5,1.,0.01)
    cmap1 = plt.cm.PuBu(mapres)
    cmap1[mapres<chances[dn],:] = np.array([0.94,0.94,0.94,1.0])
    cmap1 = clrs.ListedColormap(cmap1, name='PuBuE', N=cmap1.shape[0])
    cmap2 = plt.cm.PuRd(mapres)
    cmap2[mapres<chances[dn],:] = np.array([0.94,0.94,0.94,1.0])
    cmap2 = clrs.ListedColormap(cmap2, name='PuRdE', N=cmap2.shape[0])
    cmaps = [cmap1, cmap2]

    for cix,taskix in enumerate([0,2]):

        axs = ax[1+cix,0]
        
        # cmap = figs.getcorrcolormap('correlation')
        x = crossdecoder_matrix_allaspects_all[n_cherry,taskix,:,:,0]
        if cix==0:
            x[:T['stimstart_idx'],:] = np.nan
            x[:,:T['stimstart_idx']] = np.nan
        cf = axs.pcolormesh(x,vmin=0.5,vmax=1,cmap=cmaps[cix])
        
        axs.set_aspect('equal')
        
        ticks=[150,450]
        ticklabels=['0 ms','3000 ms']
        axs.set_xticks(ticks)
        axs.set_xticklabels(ticklabels)                
        axs.set_yticks(ticks)
        axs.set_yticklabels(ticklabels,rotation=90)                
        
        axs.set_title('%s accuracy'%taskaspects[cix])
        
        axs.set_xlabel('test timecourse')
        axs.set_ylabel('train timecourse')
    
        # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        fig.colorbar(cf,ax=axs,ticks=np.arange(0.5,1.1,0.25))
    
    
        figs.labelaxis(axs,panels[cix])








    # explanation of forward decay rate, inset into the visual square
    taskix = 0    # choose the visual line
    axs = ax[1,0]         # pin the visual accuracy matrix
    axins = axs.inset_axes([1.0,-0.5, 0.3,0.3],transform=axs.transAxes)
    # axins = axs.inset_axes([1.6,0.6, 0.4,0.4],transform=axs.transAxes)
    # axins = axs.inset_axes([75,15,5,5 ])
    
    # draw the highlight
    # axs.add_patch(plt.Rectangle((15,14),50,2, ec='black', fc=[0,0,0,0], lw=1))
    
    # draw the inset; first the data, then the fit, then the text to show the decay constant
    data = crossdecoder_matrix_allaspects_all[n_cherry,taskix,T['stimstart_idx']+5,T['stimstart_idx']+5:T['stimstart_idx']+55,0]
    axins.plot(data,color=colors[taskix])
    # fit
    m,c,b = neph.decay(data)
    # draw the fit
    xlattice = np.linspace(0,51,500)
    # axins.plot(xlattice,np.exp(m*xlattice)*c+b,color='black',lw=3)
    axins.plot(xlattice,m*xlattice+c+0.08,color='black',lw=3)
    axins.text(-2.5,0.1,'decay rate:\naccuracy loss\nover time',transform=axins.transAxes)
    
    axins.set_xticks([0,50])
    axins.set_xticklabels(['1000','1500'])
    # axins.set_yticks([0.5,1])
    # combine the two images with clear lines

    axs.add_patch(patches.ConnectionPatch(xyA=(250,250), xyB=(0,1),coordsA="data", coordsB="axes fraction", axesA=axs, axesB=axins))
    axs.add_patch(patches.ConnectionPatch(xyA=(300,250), xyB=(1,1),coordsA="data", coordsB="axes fraction", axesA=axs, axesB=axins))

    # axs.plot([65,75],[14,15])
    # axs.plot([65,75],[14+2,15+5])
    # axins.axis('off')
    









    

    # forward decay rate    
    # panels = [['E','F'],['H','I']]
    panels = ['D','G']
    # EF single mice,   # H,I all mice
    times = np.arange(sh[2])*T['dt']+T['offsettime']

    for cix,taskix in enumerate([0,2]):
        
        # axs = ax[1+cix,1]
        
        # # x = neph.smooth(fdr[n_cherry,taskix,:],5) * 1/pq.ms
        # x = fdr[n_cherry,taskix,:,0].copy()
        # if taskix==0: x[:T['stimstart_idx']] = np.nan / pq.ms           # visual pre is nonexistent
        # # x = x*fdr[n_cherry,taskix,:,1].copy()
        # axs.plot(times,x,lw=2,color=colors[cix],label='single mouse')

        # axs.legend(frameon=False,loc=2)

        # # axs.set_ylim(0,500)

        # axs.set_xlim(-1200,3800)
        # figs.setxt(axs)
        # figs.plottoaxis_stimulusoverlay(axs,T)

        # # axs.set_yticks([0,0.1])

        # axs.set_ylabel('accuracy proportion\ndecay constant [1/100ms]')


        # figs.labelaxis(axs,panels[cix][0])

        # axs.spines['right'].set_visible(False)
        # axs.spines['top'].set_visible(False)    


        axs = ax[1+cix,1]
        
        # plot all mice
        for n in range(sh[0]):
            x = neph.smooth(fdr[n,taskix,:,0],10,mode='same')
            if taskix==0: x[:T['stimstart_idx']] = np.nan   # visual pre is nonexistent
            if n==0: axs.plot(times,x,lw=0.5,color=colors[cix],alpha=0.4,label='single mice, smoothed')
            else: axs.plot(times,x,lw=0.5,color=colors[cix],alpha=0.4)
        
        # plot mean
        dr=fdr[:,taskix,:,0]   #,bdr[:,taskix,:,0]]:
        m = neph.smooth(dr.mean(axis=0),kernelwidth=10,mode='same')
        e = neph.smooth(dr.std(axis=0)/np.sqrt(sh[0]),kernelwidth=10,mode='same')*2
        if taskix==0:
            m[:T['stimstart_idx']] = np.nan / pq.ms   # visual pre is nonexistent
            e[:T['stimstart_idx']] = np.nan / pq.ms  # visual pre is nonexistent
        axs.plot(times,m,lw=3,color=colors[cix],label='mean of %d mice'%sh[0])
        axs.fill_between(times,m-e,m+e,color=colors[cix],alpha=0.3,label='2 s.e.m.')
    
        axs.set_yticks([0,0.5,1,1.5])

        axs.legend(frameon=False,loc=2)
        
        axs.set_ylim(0,1.5)
    
        axs.set_xlim(-1200,3800)
        figs.setxt(axs)

        figs.plottoaxis_stimulusoverlay(axs,T)
    
    
        axs.set_title('%s accuracy decay'%taskaspects[cix])
        axs.set_ylabel('rate [1 accuracy/100 ms]')

    
        figs.labelaxis(axs,panels[cix])
    
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)










    # DV dynamics: rotation angles each to its own
    panels = ['E','H']

    #  n_cherry use the same mouse for angles

    for cix,cx in enumerate([0,2]):     # go through visual and context
        axs = ax[1+cix,2]
        
        cmap = figs.getcorrcolormap('correlation')
        x = angles_highres_all[n_cherry][cx,:,:]
        if cix==0:
            x[:T['stimstart_idx'],:] = np.nan
            x[:,:T['stimstart_idx']] = np.nan
        cf = axs.pcolormesh(x,vmin=0,vmax=180,cmap=cmap)
        
        axs.set_aspect('equal')
        
        ticks=[150,450]
        ticklabels=['0 ms','3000 ms']
        axs.set_xticks(ticks)
        axs.set_xticklabels(ticklabels)                
        axs.set_yticks(ticks)
        axs.set_yticklabels(ticklabels,rotation=90)                
        
        # axs.set_title(taskaspects[cix])
        axs.set_title('%s\nshifted DV angles'%taskaspects[cix])
        axs.set_xlabel('shift')
        axs.set_ylabel('reference')

        # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        fig.colorbar(cf,ax=axs,ticks=np.arange(0,181,30))

        figs.labelaxis(axs,panels[cix])








    panels = ['I','J','K']


    n_mice = len(datanames)

    dispersions = []
    datanames_dispersion = []
    n_th = 25         # min number of neurons per mice for aggregate displays
    
    fs = np.ceil(np.sqrt(n_mice)).astype(np.int16)
    cmap = figs.getcorrcolormap('correlation')
    
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


        if n == datanames.index('DT019'):
    
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
    
        
            axs = ax[3,0]
            
            
            
            # axs.plot(times,c_db_matrix)
            # cf = axs.pcolormesh(times,np.arange(n_neurons+1),c_db_matrix[:,order].T,vmin=-0.5,vmax=+0.5,cmap=cmap)
            cf = axs.pcolormesh(times,np.arange(n_neurons+1),cellmatrix.T,vmin=-0.5,vmax=+0.5,cmap=cmap)
            fig.colorbar(cf,ax=axs,ticks=np.arange(-1,1.1,0.5))
            # axs.hlines(changepoint,T['starttime'],T['endtime'],color='white',lw=2)
            axs.vlines( (T['stimstarttime'],T['stimendtime']),0,n_neurons,color='white',lw=2)
            # axs.set_ylim(0,n_neurons)
            axs.set_yticks([])
            figs.setxt(axs)
            axs.set_xlabel('time from stimulus onset')
            axs.set_ylabel('units')
            figs.labelaxis(axs,panels[0])
            
            # add the dispersions for all neurons:
        m_o = c_db_matrix[T['stimstart_idx']:T['stimend_idx'],:].mean(axis=0)
        m_p = np.r_[ c_db_matrix[0:T['stimstart_idx'],:],  c_db_matrix[T['stimend_idx']:,:]  ]
        m_p = m_p.mean(axis=0)
            
            
        if n_neurons>n_th:
            dispersions.append(np.c_[m_p,m_o])
            datanames_dispersion.append(dn)
        


    H = np.concatenate([ np.abs(dispersions[n].reshape((-1,1))) for n in range(len(datanames_dispersion)) ])
    x_res = np.arange(0,0.501,0.01)
    hg,_ = np.histogram(H, bins=x_res, weights=H)
    hg /= np.sum(hg)
    chg = np.cumsum(hg)
    w = np.where(chg>0.1)
    print(w[0][0],x_res[w[0][0]])
    

    # th = 0.075         # minimum value of decoder coefficient for aggregate displays
    th = x_res[w[0][0]]


        
    axs = ax[3,1]
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
    axs.vlines([-th,th],-0.5,0.5,color='white',)
    axs.hlines([-th,th],-0.5,0.5,color='white')

    # axs.text(0.22,0.08,'th=%4.2f'%th,color='white',verticalalignment='center')

    s_label = ['constant','sign change','on-off']
    axs.text(-0.45,0.07,s_label[1],color='white',verticalalignment='center')
    axs.text(-0.45,0.0,s_label[2],color='white',verticalalignment='center')
    axs.text(-0.45,-0.07,s_label[0],color='white',verticalalignment='center')


    figs.labelaxis(axs,panels[1])


    axs = ax[3,2]
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
    

    figs.labelaxis(axs,panels[2])

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)








    fig.tight_layout()


    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Fig4_context,subspace-dynamics'+ext)
























#       FIG 5                 initial conditions

def figure5():


    #prepare data

    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting 
    n_mice = len(datanames)


    dn_Single = 'DT014'
    n_single = datanames.index('DT014')
 
    projected_dynamics_all = []       #       (mice)(tasks,classes,[trials,2])
    B_all = []
    differences = [] # (mice)(class)(mean/s.e.m)(trajectory)



    C_all = []
    for n,dn in enumerate(datanames):
        projected_dynamics,V,B = pickle.load( open('../cache/subspaces/projections,dbnv-visual,context_%s.pck'%(dn),'rb'))
        projected_dynamics_all.append(projected_dynamics)
        B_all.append(B)



        dx = 0     # we are intersted only in visual dbnv projections in the two contexts and behavioural choices
        cx_av = 0  # attend visual
        cx_iv = 1  # ignore visual
        aix_h = 0  # used indices for correct trials only: hit,
        aix_c = 2  #                                       correct rejection

        differences.append([])
        C = []
        for aix in [aix_h,aix_c]:
    
            # mean
            trial_av = np.array(projected_dynamics[cx_av][aix]).mean(axis=0)[:,dx]
            trial_iv = np.array(projected_dynamics[cx_iv][aix]).mean(axis=0)[:,dx]
            
            nf = np.max(np.hstack((np.abs(trial_av),np.abs(trial_iv)))) # normalization factor to compare to maximum signal
            # nf = 1
            trial_av /= nf
            trial_iv /= nf
        
            # s.e.m.
            trial_av_s = np.array(projected_dynamics[cx_av][aix]).std(axis=0)[:,dx]/nf
            trial_iv_s = np.array(projected_dynamics[cx_iv][aix]).std(axis=0)[:,dx]/nf
            trial_av_e = trial_av_s / np.sqrt(len(trial_av_s))
            trial_iv_e = trial_iv_s / np.sqrt(len(trial_iv_s))
            
            
            differences[-1].append(    [trial_av - trial_iv, (trial_av_s + trial_iv_s)*2, (trial_av_e + trial_iv_e)*2 ] )
    

            C.append([trial_av,trial_iv])

        C_all.append(C)      # this will be for (n_mice, stimulus, context, trajectory)
        
    C_all = np.array(C_all)
    
    if 0:     # use only when printing stats
        print('stats for lack of contextual visual displacement in activity space: ')
        print(C_all.shape)
        n_trajectory = trial_av.shape[0]
        p = np.zeros(C_all.shape[3])
        
        for aix in [0,1]:       # 45 135
            for t in range(n_trajectory):
                _, p[t] = sp.stats.ttest_rel(C_all[:,aix,0,t],C_all[:,aix,1,t])
            print(p[T['stimstart_idx']:].mean(),p[T['stimstart_idx']:].std()/np.sqrt(n_trajectory-T['stimstart_idx']))

        return
    
    
    
    # draw figure
    
    
    

    taskaspects = ['attend visual','ignore visual']

    classnames = [['45°','135°'],['attend','ignore']]
    
    perfnames = [ ['45° & hit','135° & correct rejection'],\
                  ['45° & visual hit','135° & visual correct rejection'], \
                  ['45° & audio hit','135° & audio correct rejection'] ]     # (context,visual)

    behavlabels = ['hits','correct rejections']

    # taskcolors = [ ['darkgreen','darkred'],['darkcyan','darkorange'] ]
    taskcolors = [ ['darkgreen','darkred'],['mediumseagreen','orangered'] ]
    

    basisaspects = ['visual','context']
    n_ba = len(basisaspects)
    basisaspectcolors = ['dodgerblue','fuchsia']   # np.array(taskcolors[0])[np.array([0,2,3],dtype=np.int16)]

    tracolor = 'black'


    alpha = 0.8
    alfac = 5
    # time parametrization
    skip_idx = 20
    times = np.arange(-1500,4510,10)*pq.ms
    
    t_all = times[skip_idx:-skip_idx]
    t_pre = times[skip_idx:T['stimstart_idx']+1]
    t_on = times[T['stimstart_idx']:T['stimend_idx']+1]
    t_post = times[T['stimend_idx']:-skip_idx]




    
    fig,ax = plt.subplots(3,2,figsize=(2*8.5, 3*8) )
    
    
    
    
    for cx,comparison in enumerate(taskaspects): # attend visual, ignore visual
        # visual is 0, context is 2, so cx index*2 will be the good one
        # similarly we don't need miss and false alarm trials for this figure, so skip index 1 and 3
        classresponses = [ cr for crx,cr in enumerate(projected_dynamics_all[n_single][cx]) if crx in [0,2] ]
        
        for aix,classresponse in enumerate(classresponses): # gp through classes  45h,m and 135c,f
            # we only need visual and context
            
            # correct sign for visual for this mouse:
            d = -1
            # if cx==0: d = -1
            # else: d = 1
        
            # trial average:
            
            trials = d*np.array(classresponse)            # all single trials
            trial_av = trials.mean(axis=0)        # trial averages
            trial_e = trials.std(axis=0)*2/np.sqrt(len(trials))
            
            # trials_pre = trials[:,skip_idx:T['stimstart_idx']+1,:]
            # trials_on = trials[:,T['stimstart_idx']:T['stimend_idx']+1,:]
            # trials_post = trials[:,T['stimend_idx']:-skip_idx,:]

            trial_av_pre = trial_av[skip_idx:T['stimstart_idx']+1,:]
            trial_av_on = trial_av[T['stimstart_idx']:T['stimend_idx']+1,:]
            trial_av_post = trial_av[T['stimend_idx']:-skip_idx,:]


            panels = ['A','B']

            # 2D 2 dbnv
            axs = ax[0,aix]
            
            # axs.set_title(perfnames[0][aix])
            axs.set_title(classnames[0][aix])
            
   
            # activities (each single trial, criss crossing all over the place)
            # axs.plot( trials_pre[:,:,0].T, trials_pre[:,:,1].T, lw=0.5,color=taskcolors[cx][aix],alpha=alpha/alfac)
            # axs.plot( [trials_pre[:,0,0]], [trials_pre[:,0,1]], 'd',color=taskcolors[cx][aix],alpha=alpha/alfac )
            # axs.plot( trials_on[:,:,0].T, trials_on[:,:,1].T, lw=0.5,color=tracolor,alpha=alpha/alfac )
            # axs.plot( [trials_on[:,0,0]], [trials_on[:,0,1]], 'o',color=tracolor,alpha=alpha/alfac )
            # axs.plot( trials_post[:,:,0].T, trials_post[:,:,1].T, '--',lw=0.25,color=taskcolors[cx][aix],alpha=alpha/alfac )

            # trial averaged activities
            axs.plot( [trial_av_pre[0,0]], [trial_av_pre[0,1]], 'd',markersize=15,color=taskcolors[cx][aix],alpha=alpha)#,label='-1500 ms' )
            axs.plot( trial_av_pre[:,0], trial_av_pre[:,1], lw=2,color=taskcolors[cx][aix],alpha=alpha)#,label='trial av. pre' )
            axs.plot( [trial_av_on[0,0]], [trial_av_on[0,1]], 'o',markersize=15,color=taskcolors[cx][aix],alpha=alpha)#, label='0 ms')
            axs.plot( trial_av_on[:,0], trial_av_on[:,1], lw=5,color=taskcolors[cx][aix],alpha=alpha,label='%s'%comparison )
            # axs.plot( trial_av_post[:,0], trial_av_post[:,1], '--',lw=2,color=taskcolors[cx][aix],alpha=alpha )
            
            
            
            # basis vectors
            if cx==0:
                axs.legend(frameon=False)
                for bx,basisaspect in enumerate(basisaspects):
                    # xoff = -np.sign(B[0,bx]) * 1.9
                    # yoff = -np.sign(B[1,bx]) * 1.9
                    xoff = -1.9
                    yoff = -1.9
                    axs.plot([xoff,xoff+d*B[0,bx]],[yoff,yoff+B[1,bx]],lw=3,color=basisaspectcolors[bx])

            
            axs.set_xlabel('visual DV')
            if aix==0: axs.set_ylabel('orthogonalized context DV')
            
            figs.labelaxis(axs,panels[aix])



            for bx,basisaspect in enumerate([basisaspects[0]]):
            

                # plot 1D projections only
                panels = ['C','D']
                axs = ax[1,aix]

                # single trial activities
                # axs.plot( t_pre, trials_pre[:,:,bx].T, lw=0.5,color=tracolor,alpha=alpha/alfac)
                # for l in range(len(trials_pre)):
                #     axs.plot( t_pre[0], trials_pre[l,0,bx], 'd',color=tracolor,alpha=alpha/alfac )
                #     axs.plot( t_on[0],trials_on[l,0,bx], 'o',color=tracolor,alpha=alpha/alfac )
                # axs.plot( t_on, trials_on[:,:,bx].T, lw=0.5,color=tracolor,alpha=alpha/alfac )
                # axs.plot( t_post, trials_post[:,:,bx].T, '--',lw=0.5,color=tracolor,alpha=alpha/alfac )

                # trial averaged activities
                axs.plot( t_pre[0], trial_av_pre[0,bx],  'd',markersize=15,color=taskcolors[cx][aix],alpha=alpha)
                axs.plot( t_pre, trial_av_pre[:,bx], lw=2,color=taskcolors[cx][aix],alpha=alpha)
                axs.plot( [t_on[0]], [trial_av_on[0,bx]],  'o',markersize=15,color=taskcolors[cx][aix],alpha=alpha)
                axs.plot( t_on, trial_av_on[:,bx], lw=5,color=taskcolors[cx][aix],alpha=alpha,label='%s'%comparison)
                axs.plot( t_post, trial_av_post[:,bx],'--', color=taskcolors[cx][aix],alpha=alpha)

                axs.fill_between( times,trial_av[:,bx]-trial_e[:,bx],trial_av[:,bx]+trial_e[:,bx],\
                         color=taskcolors[cx][aix], alpha=0.2)
                

                axs.set_ylim(-2,2)
                axs.set_yticks([-2,0,2])
                figs.setxt(axs)
                if cx==1:
                    figs.plottoaxis_stimulusoverlay(axs,T)
                    figs.plottoaxis_chancelevel(axs)

                axs.legend(frameon=False)

                if aix==0: axs.set_ylabel(basisaspect+' DV')
                # axs.set_xlabel('trial time from stimulus onset [ms]')

                # axs.set_title(behavlabels[aix])

                axs.spines['right'].set_visible(False)
                axs.spines['top'].set_visible(False)    

                figs.labelaxis(axs,panels[aix])



    
    signs = np.zeros((2,len(datanames),1))
    for aix in [0,1]:
        
        panels = ['E','F']
        axs = ax[2,aix]
        

        for n,dn in enumerate(datanames):
        
            if n==0: axs.plot(times,differences[n][aix][0],color=tracolor,lw=2,alpha=0.15,label='single mice, smoothed')
            else: axs.plot(times,differences[n][aix][0],color=tracolor,lw=2,alpha=0.15)
            
            signs[aix,n] = np.sign(differences[n][aix][0].sum())
            
            # axs.fill_between(times,differences[n][aix][0]-differences[n][aix][2],\
            #                        differences[n][aix][0]+differences[n][aix][2],\
            #                        color=np.array([0,0,(n+1.)/n_mice]),alpha=0.02)



        m = np.array( differences )[:,aix,0,:].mean(axis=0)
        e = np.array( differences )[:,aix,0,:].std(axis=0)/np.sqrt(n_mice)*2
        # below is the flip corrected to check
        # m = (signs[aix,:,:] * np.array( differences )[:,aix,0,:]).mean(axis=0)
        # e = (signs[aix,:,:] * np.array( differences )[:,aix,0,:]).std(axis=0)/np.sqrt(n_mice)*2
        
        
            
            
        
        
        axs.plot(times,m,color=taskcolors[0][aix],lw=3,alpha=0.9,label='mean of %d mice'%n_mice)
        axs.fill_between(times,m-e,m+e,color=taskcolors[0][aix],alpha=0.2)
        
        
        
        axs.legend(frameon=False,loc=2)
        
        # axs.set_title(behavlabels[aix])
        # if aix==0: axs.set_ylabel('contextual difference index of\nvisual dbnv projection')
        if aix==0: axs.set_ylabel('context modulation of visual response')
        axs.set_xlabel('time from stimulus onset')
        
        axs.set_ylim(-1,1)
        axs.set_xlim(times[0],times[-1])
        axs.set_yticks(np.arange(-1,1.1,1))
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs)
        
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)    
        
        
        figs.labelaxis(axs,panels[aix])

    print(signs[:,:,0])

    for aix in [0,1]:
        s = 2
        axs = ax[0,aix]
        axs.legend(frameon=False)
        axs.set_xlim(-s,s)
        axs.set_ylim(-s,s)
        axs.set_xticks([-s,0,s])
        axs.set_yticks([-s,0,s])

        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)






    fig.tight_layout()


    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Fig5_visual,contextinvariant'+ext)




















#         FIG 6        CHOICE




def figure6():
    
    
    variablecolors = ['navy','darkorange']
    
    # A 
    dn = 'DT019'
    n_sel = 8
    examine = 'allexpcond'
    comparison = 'choice'
    acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
    choicetimecourse = acrossdecoder[:2]
    
    
    # B, C, choice tribe             D attend - ignore
    examine = 'allexpcond'
    # datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']                # JRC
    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting

    # get baseline shuffle:
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    # here limit by number of neurons
    # nneuronsdict = preprocess.getnneurons(datanames)
    # print(nneuronsdict)
    # datanames = [ dn for dn in nneuronsdict if nneuronsdict[dn]>20 ]
    # print(datanames)
    
    timegroups = []
    for n,dn in enumerate(datanames):
        
        contexttimegroups = []; choicetimegroups = []; diffbothtimegroups = []; diffignoretimegroups = []
#        timestarts_idx = int((np.arange(0,6001,1500)*pq.ms / np.array(T['dt'])).magnitude)     # cutpoints
        timestarts_idx = np.arange(0,601,150,dtype='int16')

        comparison = 'context'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        for ti in range(4):
            contexttimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )

        comparison = 'choice'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        for ti in range(4):
            choicetimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )

        comparison = 'visual'
        acrossdecoderattend = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('attendignore','all',dn,continuous_method,'attend '+comparison,'all'),'rb'))
        acrossdecoderignore = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('attendignore','all',dn,continuous_method,'ignore '+comparison,'all'),'rb'))
        for ti in range(4):
            attendsignal = acrossdecoderattend[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()
            ignoresignal = acrossdecoderignore[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()
            diffbothtimegroups.append(      attendsignal - ignoresignal    )

        crossdecoder = pickle.load(open('../cache/continuous/responsedecodes,crosscontext_%s-%s.pck'%(dn,continuous_method),'rb'))
        for ti in range(4):
            ignoresignal = crossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()
            attendsignal = crossdecoder[2][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()
            diffignoretimegroups.append(      attendsignal - ignoresignal    )

        timegroups.append([contexttimegroups, choicetimegroups, diffbothtimegroups, diffignoretimegroups])
    timegroups = np.array(timegroups)


    if 0:     # stats
        m = timegroups[:,1,1:3].mean(axis=0)
        e = timegroups[:,1,1:3].std(axis=0)/np.sqrt(timegroups.shape[0])
        _,p = sp.stats.ttest_rel(timegroups[:,1,1],timegroups[:,1,2])
        print('stats, choice')
        print('comparison between early and late accuracy, ttest p=%4.4f'%(p/2))
        print(m,e)



        return
    




    # new projections: choice is parallel, context is semi: show on orthonormalized closest to it
    datanamesefg = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting
    # choice is split into taken in attend visual and attend audio only
    depth_idx = 150  # averaging onto 1500 - 3000 ms into stimulus onset !!!!! unlike visual context, we need the late part!
    projections_all = []   # this will be (mice)(attends)(taskconditionlist)(trials,dbnvcoords)       will have 8 conditionlist
    basis_all = []
    whichbasis = [1,0]          #   0 = choice,  1 = context
    for n,dn in enumerate(datanamesefg):
        projections_attends = []
        basis_attends = []
        for cx,comparison in enumerate(['choiceattendvisual','choiceattendaudio']):
            projected_dynamics, basis = pickle.load(open('../cache/subspaces/subspacedynamics,projected+dbnv,cxchvi,%s-%s_%s-%dms.pck'%(['av','aa'][cx],dn,continuous_method,T['dt'].magnitude),'rb'))
            projected_dynamics = np.array(projected_dynamics)
            projection = []
            for k in np.arange(4)+cx*4:     # choose the projection only that are appropriate for the context
                if len(projected_dynamics[k])>0:
                    projection.append(  np.array(projected_dynamics[k])[:,T['stimstart_idx']+depth_idx:T['stimstart_idx']+2*depth_idx, whichbasis ].mean(1) )
                else: projection.append( np.zeros((0)) )
            projections_attends.append(projection)
            basis_attends.append(basis[whichbasis,:])
        projections_all.append(projections_attends)
        basis_all.append(basis_attends)
    # print(len(projections_all),len(projections_all[0]),len(projections_all[0][0]),len(projections_all[0][0][0]))
    # print(len(basis_all),len(basis_all[0]),len(basis_all[0][0]),len(basis_all[0][0][0]))


    # colors should corerspond to the 8 conditoin list:   av hmcf  and aa hmcf

    # activitycolors = [['darkgreen','darkorange','navy','darkred'],\
    #                   ['lime','gold','dodgerblue','orangered']]
    activitycolors = ['darkorange','firebrick','firebrick','darkorange']












    # dbnv angles
    angles_all = []
    angles_highres_all = []
    for n,dn in enumerate(datanames):
        angles = pickle.load(open('../cache/subspaces/angles,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_all.append(angles)
        angles_highres = pickle.load(open('../cache/subspaces/angles,highres,alongDBNVs-VACC3_%s-%dms_%s'%(continuous_method,T['bin'].magnitude,dn),'rb'))
        angles_highres_all.append(angles_highres)

    angles_all = np.array(angles_all)

    times_angle = np.arange(0,angles_all.shape[3])*T['dt'] + T['offsettime']


    if 0:     # stats
        print('stats, choice angle to context')
        for cx in [0,1]:
            angles = np.zeros((len(datanamesefg)))
            for n,dn in enumerate(datanamesefg):
                #choose choice (4 and 5 for av and aa) and context (2); we show here the last 1500 ms during stimulus
                angles[n] = angles_all[n,2,4+cx,T['stimstart_idx']+150:T['stimstart_idx']+300].mean()
    
    
            m = angles.mean()
            e = angles.std()/np.sqrt(angles.shape[0])
            print('%s angle = %4.2f+/-%4.2f°'%(['attend','ignore'][cx],m,e))
        return









    # behaviour conditioned visual discrimination in attend visual, for all mice
    # actvities_all = [n_mice][visualattend,context][performance][trials] = 15 x 2 x 4 x trials
    activities_all,activities_aggregate = pickle.load(open('../cache/subspaces/dbnvprojections-visual,behaviourdifference-%s-%dms.pck'%(continuous_method,T['dt']),'rb'))











    #        FIGURE

#    fig = plt.figure(num=5)
#    fig.clf()
#    fig, ax = plt.subplots(2,3,figsize=(24,16))  # 36 36
    f1n = 3; f2n = 2
    fig = plt.figure(constrained_layout=False,figsize=(f2n*8,f1n*8))
    gs = fig.add_gridspec(f1n, f2n)
    ax = []
    ax_joint = []
    ax_marginals = []
    for j in np.arange(f1n):
        axa=[]
        for k in np.arange(f2n):


    # handle special axis issues
    # this is for a joint and marginal histograms of context vs. choice
            if j==1 and k in [0,1]:
                ratio = 7
                gs_marginals = gs[j,k].subgridspec(ratio+1, ratio+1)
                ax_joint.append(fig.add_subplot(gs_marginals[1:, :-1]))
                ax_marginals.append([ fig.add_subplot(gs_marginals[0 , :-1], sharex=ax_joint[k]), \
                                 fig.add_subplot(gs_marginals[1:,  -1], sharey=ax_joint[k]) ])
                axa.append( ax_joint[k] )
            elif j==2 and k in [0,1]:
                axa.append( fig.add_subplot(gs[j, k], projection='polar') )
            else:
                axa.append( fig.add_subplot(gs[j, k]) )
        ax.append(axa)
    ax = np.array(ax)











    # A
    panel = 'A'

    axs = ax[0,0]
    figs.plottoaxis_decoderrocauc(choicetimecourse,axs,colorlist=['',variablecolors[1]],plottrain=False)       # plot the test performance
    figs.plottoaxis_stimulusoverlay(axs,T)
    figs.plottoaxis_chancelevel(axs,chances[dn])
    figs.setxt(axs)
    axs.set_yticks([0.5,1.0])# axs.set_yticklabels([0.5,1.0])
    axs.set_ylabel('choice accuracy            ')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    # axins = axs.inset_axes([-0.42,0.55, 0.4,0.4],transform=axs.transAxes)
    # drawschematics(axins,1); axins.axis('off')
    
    figs.labelaxis(axs,panel,x=-0.3)








    # B
    panel = 'B'
    timecourselabels = ['PRE','early ON','late ON','POST']
    axs = ax[0,1]
    ch = 1 # show choice
    ci = timegroups[:,ch,:].std(axis=0)*2/np.sqrt(len(datanames))
    m = timegroups[:,ch,:].mean(axis=0)
    ci = np.c_[m-ci,m+ci]
#    print(timegroups.shape,m.shape,ci.shape)
    artist = axs.boxplot(x=timegroups[:,ch,:], positions=-750+timestarts_idx[:4]*10, notch=True,usermedians=m,conf_intervals=ci,\
                whis=[5,95],labels=timecourselabels,widths=350)

    for element in artist.keys():
        plt.setp(artist[element], color=variablecolors[1],linewidth=2)

    axs.set_yticks([0.5,1.0])
    axs.set_ylim(0.45,1.01)
    axs.set_xlim(-1500,4500)
    figs.plottoaxis_stimulusoverlay(axs,T)
    m_c = np.array(list(chances.values())).mean()
    e_c = np.array(list(chances.values())).std()/np.sqrt(len(datanames))
    figs.plottoaxis_chancelevel(axs,m_c+e_c)
#    axs.set_ylabel('accuracy')
#    axs.set_title('n=%d'%14)
    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)










   
    
    



    # search-hyperjump keyword: asdf
    # SUBSPACES for Context vs. Choice
    panel = 'C'
    
    n_showmice = datanamesefg.index('DT019')
    basisaspects = ['context','choice']
    basisaspectcolors = ['mediumvioletred','gold']   # np.array(taskcolors[0])[np.array([0,2,3],dtype=np.int16)]

    projections = projections_all[n_showmice]
    basis = basis_all[n_showmice]

    # flip: good only for DT019
    for cx in [0,1]:
        for k in [0,1,2,3]:
            if len(projections[0][k])>0:
                for i in [0,1]:
                    projections[cx][k][:,i] = -projections[cx][k][:,i]
        for bx in [0,1]:
            for i in [0,1]:
                basis[cx][i,bx] = -basis[cx][i,bx]


    for cx,comparison in enumerate(['attend visual','attend audio']):
        axs = ax[1,cx]
        print(cx)
        
        # plot the trial averages as points, separate for av,aa
        for k in range(4):
            print('trials:',len(projections[cx][k]))
            if len(projections[cx][k])>0:
                facecolors = [activitycolors[k],'none',activitycolors[k],'none']
                # ['o','x','o','+'][k]
                markersize = [8,10,8,10][k]

                axs.plot(projections[cx][k][:,0], projections[cx][k][:,1], 'o',markersize=markersize,color=activitycolors[k],mfc=facecolors[k],alpha=0.8)
    
    
                # axs.plot(projections[bp[b,0],preon,0][ixgohit],      projections[bp[b,1],preon,0][ixgohit],      'o',markersize=8,color=activitycolors[0],alpha=0.8)
                # axs.plot(projections[bp[b,0],preon,0][ixgomiss],     projections[bp[b,1],preon,0][ixgomiss],     'X',markersize=10,color=activitycolors[4],alpha=0.9)
                # axs.plot(projections[bp[b,0],preon,1][ixnogocorrrej], projections[bp[b,1],preon,1][ixnogocorrrej], 'o',markersize=8,color=activitycolors[1],alpha=0.8)
                # axs.plot(projections[bp[b,0],preon,1][ixnogofal],     projections[bp[b,1],preon,1][ixnogofal],     'P',markersize=10,color=activitycolors[5],alpha=0.9)
    
                # plot the cov ellipsoids over the per trial activities
                figs.confidence_ellipse(projections[cx][k][:,0], projections[cx][k][:,1], axs, n_std=2.0, facecolor=activitycolors[k],alpha=0.15)
    
    
    
        # basis vectors
        x0 = -3.2
        y0 = -2
        for bx,basisaspect in enumerate(basisaspects):
            # axs.plot([0,B[bx,0]],[0,B[bx,1]],lw=3,color=basisaspectcolors[bx])
            axs.plot([x0+0,x0+basis[cx][0,bx]],[y0+0,y0+basis[cx][1,bx]],lw=3,color=basisaspectcolors[bx])
    
    
    
    
    
        # create the legend
    
    
        if cx==0:
            axins = axs.inset_axes([-0.2,0.4, 1,1],transform=axs.transAxes)
            axins.axis('off')
            xl0 = -0.2; yl0 = 0.75
            m = 0.08
            for j in range(2):
                for k in range(2):
                    axins.add_patch(plt.Rectangle((xl0+m*j,yl0+m*k),m,m,\
                      ec=activitycolors[j],fc=activitycolors[j],alpha=0.15,transform=axs.transAxes))
                    axins.add_patch(plt.Circle((xl0+m/2+m*j,yl0+m/2+m*k),0.015,\
                      ec=activitycolors[j],fc=['none',activitycolors[j]][k],transform=axs.transAxes))
            axins.text(xl0+m,  yl0+m*2,'lick ',ha='right',va='bottom',transform=axs.transAxes)
            axins.text(xl0+m,  yl0+m*2,' no lick',ha='left',va='bottom',transform=axs.transAxes)
            # axins.text(xl0+m*2+0.025, yl0+m*3/2,'ignore',ha='left',va='center',transform=axs.transAxes)
            # axins.text(xl0+m*2+0.025, yl0+m/2,'attend',ha='left',va='center',transform=axs.transAxes)
            axins.text(xl0-0.025, yl0+m*3/2,'correct',ha='right',va='center',transform=axs.transAxes)
            axins.text(xl0-0.025, yl0+m/2,'error',ha='right',va='center',transform=axs.transAxes)
    
    

    
        
        axs.text(0.04, -0.04,'choice DV [SU]', ha='left',   va='center', transform=axs.transAxes)
        if cx==0: axs.text(-0.04, 0.04,'context DV [SU]', ha='center', va='bottom',rotation=90, transform=axs.transAxes)

    
    
    

        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)

        axs.axis('off')
        axs.set_xlim(-3.3,3.3)
        axs.set_ylim(-2.1,2.1)
        figs.plottoaxis_crosshair(axs)


        # cumulative histograms:    
        # marginalized context (choice horizontal)
        for kx,k in enumerate([[0,3],[1,2]]):
            data = []
            # for cx in [0,1]:
            for i in k:
                if len(projections[cx][i])>0:
                    data.append(projections[cx][i][:,0])
            data = np.concatenate(data)
            # ax_marginals[cx][0].hist(data,bins=np.arange(-2.4,+3.2+0.01,0.2),color=['orange','red'][kx],alpha=0.7)
    
            kde = sp.stats.gaussian_kde(data)
            x = np.arange(-3.3,+3.3+0.01,0.02)
            ax_marginals[cx][0].plot(x,kde(x),color=['darkorange','firebrick'][kx],lw=2,alpha=0.7)
        # ax_marginals[cx][0].plot(x,np.zeros(len(x)),'k--',alpha=0.5)
        figs.plottoaxis_crosshair(ax_marginals[cx][0])
    
    
            
        # marginalized choice, (context vertical)
        # for cx in [0,1]:
        data = []
        for i in [0,1,2,3]:
            if len(projections[cx][i])>0:
                data.append(projections[cx][i][:,1])
        data = np.concatenate(data)
        # ax_marginals[cx][1].hist(data,bins=np.arange(-2.4,+2.4+0.01,0.2),color=['purple','fuchsia'][cx],alpha=0.7, orientation='horizontal')

        kde = sp.stats.gaussian_kde(data)
        x = np.arange(-2.4,+2.4+0.01,0.02)
        ax_marginals[cx][1].plot(kde(x),x,color=['indigo','deeppink'][cx],lw=2,alpha=0.7)
        ax_marginals[cx][1].plot(np.zeros(len(x)),x,'k--',alpha=0.5)
    
        data = [] # this is for the shadow histogram from the other panel   =1-cx
        for i in [0,1,2,3]:
            if len(projections[1-cx][i])>0:
                data.append(projections[1-cx][i][:,1])
        data = np.concatenate(data)
        # ax_marginals[cx][1].hist(data,bins=np.arange(-2.4,+2.4+0.01,0.2),color=['purple','fuchsia'][cx],alpha=0.7, orientation='horizontal')

        kde = sp.stats.gaussian_kde(data)
        x = np.arange(-2.4,+2.4+0.01,0.02)
        ax_marginals[cx][1].plot(kde(x),x,color=['darkgrey','darkgrey'][cx],lw=2,alpha=0.7)
        
        figs.plottoaxis_crosshair(ax_marginals[cx][1])
        # ax_marginals[cx][1].plot(np.zeros(len(x)),x,'k--',alpha=0.5)
    
    
    
        for margaxx in [0,1]:
            ax_marginals[cx][margaxx].spines['right'].set_visible(False)
            ax_marginals[cx][margaxx].spines['top'].set_visible(False)
            ax_marginals[cx][margaxx].spines['left'].set_visible(False)
            ax_marginals[cx][margaxx].spines['bottom'].set_visible(False)
            ax_marginals[cx][margaxx].get_xaxis().set_visible(False)
            ax_marginals[cx][margaxx].get_yaxis().set_visible(False)
       
    # fontsize=24,
        ax_marginals[cx][0].text(-0.3,1.05,['attend','ignore'][cx],color=['indigo','deeppink'][cx],transform=axs.transAxes)

        if cx==0: figs.labelaxis(ax_marginals[cx][0],panel,x=-0.35)



















    panel = 'D'

    # polar histogram of bases visual and context in all animals
    
    angles = np.zeros((len(datanamesefg)))    #  {bv vi,ch}  x  context x mice
    # edges = np.linspace(-np.pi/2,np.pi/2,12+1)
    edges = np.linspace(0,np.pi,24+1)
    
    width=(edges[1]-edges[0])/2
    n_bins = len(edges)-1

    for cx,comparison in enumerate(['attend','ignore']):

        anglecounts = np.zeros((n_bins,2))        # {basisvectors vi,ch}   x context
        
        for n,dn in enumerate(datanamesefg):
            #choose choice (4 and 5 for av and aa) and context (2); we show here the last 1500 ms during stimulus
            angles[n] = angles_all[n,2,4+cx,T['stimstart_idx']+150:T['stimstart_idx']+300].mean()/180*np.pi
        aux = np.histogram(angles,bins=edges,density=False)
        anglecounts = aux[0]
        # anglestats = [  angles.mean(axis=2), angles.std(axis=2)*2/np.sqrt(len(datanamesefg))   ]
        # anglestats_t_p = [ sp.stats.ttest_ind( angles[chvix,0,:], angles[chvix,1,:] )[1] for chvix in range(2) ]
    
        
    
        color = ['sienna','salmon'][cx]
        axs = ax[2,cx]
        axs.bar(edges[:-1]+width, anglecounts, width=width*2,color=color, alpha=0.7)
        # kde = sp.stats.gaussian_kde(angles)
        # x = np.arange(0,edges[-1],np.pi/180)  #edges[:-1]+width
        # axs.plot(x,kde(x),color=color,lw=2,alpha=0.7)
        
        
        
        #                axs.plot(edges[:-1]+width, anglecounts[:,chvix,cx],'o-',\
        #                        color=basiscolors[basisindices[chvix][cx]],alpha=0.7)
            # axs.errorbar(anglestats[0],9,xerr=anglestats[1][chvix,cx],color=colors[chvix][cx])
            # axs.plot(anglestats[0][chvix,cx],9,'o',color=colors[chvix][cx])
            # axs.text(anglestats[0][chvix,:].mean(),9,'p=%5.3f'%anglestats_t_p[chvix])
        # axs.legend([['attend'],['ignore']][cx],frameon=False)
        if cx==0: axs.text(-0.3,0.7,'choice-context\nDV angle',transform=axs.transAxes)
        axs.text(-0.3,0.9,['attend','ignore'][cx],color=['indigo','deeppink'][cx],transform=axs.transAxes)
        
        # anglestats_t_p
        # axs.set_xlim(-np.pi/2,np.pi/2)
        axs.set_xlim(0,np.pi)
        # axs.set_ylim(0,10)
        axs.set_xticks(edges[::6])
        axs.set_yticklabels([])
        # axs.set_xlabel(['visual basis','choice basis'])
    
        if cx==0: figs.labelaxis(axs,panel,x=-0.35,y=1.)


















# delete empty axes:
    # for axs in [ax[2,0]]:    axs.axis('off')




    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Fig6_choice'+ext,bbox_inches='tight')


















def figure7():

    datanames = ['ME110','ME113','DT009','DT030','DT031','DT032']
    dn_ch = 'DT009'
    n_cherry = datanames.index(dn_ch)
    n_mice = len(datanames)

    # get baseline shuffle:
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    decoderwidth = 5
    patchwidth = 3
    
    max_runspeed = 100.*pq.cm/pq.s
    delta_runspeed = 5.*pq.cm/pq.s
    runspeedlevels = np.linspace(0*pq.cm/pq.s, max_runspeed, int(max_runspeed/delta_runspeed)+1)
    n_runspeedlevels = len(runspeedlevels)-1

    runshifts_times = np.arange(-200,201,100)*pq.ms
    
    cx = 0
    comparison = 'context'




    # load data
    # decoders
    ad_r = []
    ad_n = []
    r_d = []
    for n,dn in enumerate(datanames):
        acrossdecoders_runshifts = []
        acrossdecoders_noeq_runshifts = []
        for rshx,rsh in enumerate(runshifts_times):
            
            acrossdecoders = []
            acrossdecoders_noeq = []
            
            
    
        
            acrossdecoder_m = pickle.load(open('../cache/locomotion/responsedecodes,locomotionequalized,shift%+dms_angles-%s-%s-%s.pck'%(rsh.magnitude,dn,continuous_method,comparison),'rb'))
            acrossdecoder_noeq_m = pickle.load(open('../cache/locomotion/responsedecodes,locomotionequalized-control,shift%+dms_angles-%s-%s-%s.pck'%(rsh.magnitude,dn,continuous_method,comparison),'rb'))
                
            acrossdecoders.append(acrossdecoder_m)
            acrossdecoders_noeq.append(acrossdecoder_noeq_m)
    
            acrossdecoders_runshifts.append(acrossdecoders)
            acrossdecoders_noeq_runshifts.append(acrossdecoders_noeq)
        
        ad_r.append(acrossdecoders_runshifts)
        ad_n.append(acrossdecoders_noeq_runshifts)
    
    
        # locomotion distribution, at 0 ms shift, only context
        r_dists = pickle.load(open('../cache/locomotion/runspeeddistributions_%s-%s-%s.pck'%(dn,continuous_method,comparison),'rb'))
        r_d.append(r_dists)








    # FIGURE
    
    
    fig,ax = plt.subplots(2,3,figsize=(3*8,2*8))
    
    
    
    # patching together 3  50 ms width runspeed distributions
    panel = 'A'
    axs = ax[0,0]
    timepoints = [150,155,160]
    timecourseindices = np.array(timepoints,dtype=np.int16)//decoderwidth

    for tx,t in enumerate(timecourseindices):
        color = np.array(mcs.to_rgb('mediumvioletred'))*((tx+1)/(patchwidth+2))
        axs.plot(runspeedlevels[:-1],r_d[n_cherry][t,:,0],lw=2,color=color,label='at %d ms'%(t*decoderwidth*T['dt']-1500*pq.ms))
        # axs.plot(runspeedlevels[:-1],r_d[n_cherry][t,:,1],'--',lw=2,color=np.array(mcs.to_rgb('mediumvioletred'))*((tx+1)/(patchwidth+2)),label='av %dms'%(t*decoderwidth*T['dt']))

    t = timecourseindices[0]    
    axs.plot(runspeedlevels[:-1],r_d[n_cherry][t:t+patchwidth,:,0].sum(axis=0),lw=4,color=np.array(mcs.to_rgb('mediumvioletred')),label='sum')

    axs.set_ylabel('# trials')
    axs.set_xlabel('runspeed [cm/s]')

    axs.legend(frameon=False, loc='upper left')

    axs.set_ylim(0,90)

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    figs.labelaxis(axs,panel)






    # selected number of trials
    panel = 'B'
    axs = ax[0,1]
    timepoints = [150,155,160]
    timecourseindices = np.array(timepoints,dtype=np.int16)//decoderwidth
    
    colors_classes = ['rebeccapurple','darkcyan']
    t = timecourseindices[0]
    for clx in [0,1]:
        axs.plot(runspeedlevels[:-1],r_d[n_cherry][t:t+patchwidth,:,clx].sum(axis=0),lw=4,color=colors_classes[clx],label=['attend visual','ignore visual'][clx])

    r_d_min = np.min(r_d[n_cherry][t:t+patchwidth,:,:].sum(axis=0),axis=1)
    axs.fill_between(runspeedlevels[:-1],np.zeros(len(r_d_min)),r_d_min,ec=None,fc='mediumvioletred',alpha=0.15,label='# selected trials')

    axs.set_ylabel('# trials')
    axs.set_xlabel('runspeed [cm/s]')
    
    axs.text(0,95,'sum',fontsize='small')

    axs.legend(frameon=False)

    axs.set_ylim(0,90)

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    figs.labelaxis(axs,panel)







    # schematics: runspeed equalization selection activity projections
    panel = 'C'
    
    axs = ax[0,2]

    np.random.seed(152)
    x1 = np.random.randn(15)/7+0.5
    y1 = np.random.randn(15)/12+0.8
    x2 = np.random.randn(20)/7+0.5
    y2 = np.random.randn(20)/12+0.2
    np.random.seed()

    axs.plot(x1,y1,'o',          markersize=12,markeredgewidth=2,markeredgecolor=colors_classes[0],markerfacecolor='none')
    axs.plot(x1[10:],y1[10:],'o',markersize=12,markeredgewidth=2,markerfacecolor='mediumvioletred',markeredgecolor='none',alpha=0.2,label='selected trial')
    axs.plot(x2,y2,'o',          markersize=12,markeredgewidth=2,markeredgecolor=colors_classes[1],markerfacecolor='none')
    axs.plot(x2[10:],y2[10:],'o',markersize=12,markeredgewidth=2,markerfacecolor='mediumvioletred',markeredgecolor='none',alpha=0.2)

    axs.legend(frameon=False,loc='upper left')

    axs.set_xlim(-0.1,1.1)
    axs.set_ylim(-0.1,1.1)
    axs.set_xticklabels([])
    axs.set_yticklabels([])
    axs.set_ylabel('context DV')
    axs.set_xlabel('context nullspace')
    
    # axs.text(0,1.08,'filled: selected trial',fontsize='small')

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    figs.labelaxis(axs,panel)






    
    # decoder acc. single mouse
    panel = 'D'
    axs = ax[1,0]


    s = len(runshifts_times)//2+1-1  # find the zero shift position

    figs.plottoaxis_decoderrocauc(ad_n[n_cherry][s][cx][:2],axs,colorlist=['white','black'],plottrain=False,onlysem=True,label='original')       # plot the performance
    figs.plottoaxis_decoderrocauc(ad_r[n_cherry][s][cx][:2],axs,colorlist=['white','mediumvioletred'],plottrain=False,onlysem=True,label='locomotion dist. equalized')       # plot the performance
    figs.setxt(axs)
    figs.plottoaxis_stimulusoverlay(axs,T)
    figs.plottoaxis_chancelevel(axs,chances[dn_ch])
    axs.legend(frameon=False)
    axs.set_ylabel('context decoder accuracy')



    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    figs.labelaxis(axs,panel)









    # dec. acc. all mice
    panel = 'E'
    axs = ax[1,1]


    s = len(runshifts_times)//2+1-1  # find the zero shift position
    
    times = ad_n[n][s][cx][1].times
    differences = np.zeros((n_mice,len(ad_r[0][s][cx][1][:,0])))
    for n,dn in enumerate(datanames):
        differences[n,:] = (ad_r[n][s][cx][1][:,0]-ad_n[n][s][cx][1][:,0]).squeeze()
    
    M = differences.mean(axis=1)
    S = differences.std(axis=1)
    print('stats: difference   %4.4f+/-%4.4f,  %4.4f+/-%4.4f'%(M.mean(),M.std()/np.sqrt(len(M)),\
                                                               S.mean(),S.std()/np.sqrt(len(S))))


    if 0:
        for n,dn in enumerate(datanames):
        
            if n==0: axs.plot(times,differences[n,:],color='purple',lw=2,alpha=0.2,label='single mice')
            else: axs.plot(times,differences[n,:],color='purple',lw=2,alpha=0.2)
             
            # axs.fill_between(times,differences[n][aix][0]-differences[n][aix][2],\
            #                        differences[n][aix][0]+differences[n][aix][2],\
            #                        color=np.array([0,0,(n+1.)/n_mice]),alpha=0.02)
    
    
        m = differences.mean(axis=0)
        e = differences.std(axis=0)/np.sqrt(n_mice)*2
        
        axs.plot(times,m,color='purple',lw=3,alpha=0.9,label='mean of %d mice'%n_mice)
        axs.fill_between(times,m-e,m+e,color='purple',alpha=0.2)
    
        axs.set_ylim(-0.2,0.2)        
        axs.set_yticks([-0.1,0,0.1])
        axs.set_xlim([times[0],times[-1]])
        figs.setxt(axs)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs)
    
        # axs.set_ylabel('             loc.d.eq. - orig.',labelpad=-20)
        axs.set_ylabel('difference',labelpad=-20)
        
        axs.legend(frameon=False)
    
    if 1:
        axs.boxplot(x=differences.T, notch=True,whis=[5,95],labels=None)#,labels=['mouse %d'%d for d in range(n_mice)])
        axs.set_ylim(-0.15,0.15)
        axs.set_yticks([-0.1,0,0.1])
        figs.plottoaxis_chancelevel(axs)
            
        axs.set_ylabel('difference',labelpad=-20)
        axs.set_xticks([])
        axs.set_xlabel('mice')


    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    figs.labelaxis(axs,panel)







    # runspeed shift has no effect
    panel = 'F'
    axs = ax[1,2]

    
    timepoints = np.array([-500,0,500,2000])*pq.ms
    timecourseindices = np.array((timepoints-T['offsettime'])/T['dt'],dtype=np.int16)//decoderwidth//patchwidth
    timepointlabels = ['%d ms'%tp for tp in timepoints ]
    
    
    for tx,t in enumerate(timecourseindices):
        tc = np.array(mcs.to_rgb('mediumvioletred'))
        tc /= max(tc)
        colors = np.array(tc) * (tx+1) /  (len(timecourseindices)+2)
        shifteddecoder_acc = [ ad_n[n_cherry][s][cx][1][t,0] for s in range(len(runshifts_times)) ]
        axs.plot(runshifts_times,shifteddecoder_acc,lw=3,color=colors,label=timepointlabels[tx] )

    axs.legend(frameon=False)
    
    axs.set_ylim(0.45,1.1)
    axs.set_yticks([0.5,1])
    axs.set_xlim(runshifts_times[0],runshifts_times[-1])
    figs.plottoaxis_chancelevel(axs)
    axs.set_xlabel('shift [ms]')
    axs.set_ylabel('context decoder accuracy')

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    figs.labelaxis(axs,panel)



    # fig.tight_layout()






    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Fig7_locomotioninvariance'+ext,bbox_inches='tight')






    return










def supplementaryfigure1():    # drift control


    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting
    n_mice = len(datanames)
    
    
    example_mouse = 'DT014'
    
    show_single_example = False
    
    # datanames = ['DT008','DT009']
    datanames = ['DT014','DT009']


    spike_times_all = []
    amplitudes_all = []
    pc_features_all = []
    pc_noises_all = []
    firing_rates_all = []

    for n,dn in enumerate(datanames): # go through the animals


        if show_single_example:
            if not n==datanames.index(example_mouse): continue

    
        # get block boundaries
        # block = preprocess.loaddatamouse(dn,T,continuous_method,normalize=True,recalculate=False)
        # ixs1,ixm2,ixs3,ixm4 = getblockidents(dn,block)
        trialsdata = preprocess.loadexperimentdata(dn)
        blocktimes = []
        for b in [1,3]: # start and end of multimodal blocks
            trialtimes = trialsdata['start'][trialsdata['block']==b]
            blocktimes.append( [trialtimes.values[0],trialtimes.values[-1]+4500]   )
        blocktimes = (np.array(blocktimes)*pq.ms).rescale('min')
        
    
    
        pathdatamouse_ks2sorted = '../../../data/mice/raw/ks/default/'
        pathdatamouse_brainarea = 'V1/'
        



        pc_feature_values = []
        amplitude_values = []
        spikes_timing = []
        spikes_incluster = []
        cluster_id_list_good = []
        cluster_channel_good = []
        cluster_id_list_drift = []
        cluster_channel_drift = []
        
        
        
        if show_single_example: shanklist = [1]
        else: shanklist = [0,1]
        
        for sx in shanklist:
            pathshank = 'shank%d/'%(sx+1)
        
        
            datafilepath = pathdatamouse_ks2sorted + pathdatamouse_brainarea + dn + pathshank
        
        
            pc_feature_values.append(np.load(datafilepath+'pc_features.npy'))
            # [nSpikes, nFeaturesPerChannel, nPCFeatures] single
            # matrix giving the PC values for each spike.
            # The channels that those features came from are specified in pc_features_ind.npy.
            # E.g. the value at pc_features[123, 1, 5] is the projection of the 123rd spike onto the 1st
            # PC on the channel given by pc_feature_ind[5].
            amplitude_values.append(np.load(datafilepath+'amplitudes.npy'))
            spikes_timing.append(np.load(datafilepath+'spike_times.npy'))
            spikes_incluster.append(np.load(datafilepath+'spike_clusters.npy'))
        
        
            cluster_info = pd.read_csv(datafilepath+'cluster_info.tsv',sep='\t',usecols=['id','KSLabel','channel','depth','group'])    # KSLabel or the richer cluster_info.npy contains quality measures as well
        
        
            
            # compare good and drifting units
            sua_good_ix = np.logical_and( cluster_info['KSLabel']=='good',   \
                         np.logical_and( cluster_info['group']!='noise', cluster_info['group']!='drift')    )
            sua_drift_ix = np.logical_and( cluster_info['KSLabel']=='good',   \
                             np.logical_and( cluster_info['group']!='noise', cluster_info['group']=='drift')    )
                        
            print('#clusters (all,good,drift)',len(cluster_info),sum(sua_good_ix),sum(sua_drift_ix))
            
            cluster_id_list_good.append( cluster_info['id'].loc[ sua_good_ix ].values )
            cluster_channel_good.append( cluster_info['channel'].loc[ sua_good_ix ].values )
        
            cluster_id_list_drift.append( cluster_info['id'].loc[ sua_drift_ix ].values )
            cluster_channel_drift.append( cluster_info['channel'].loc[ sua_drift_ix ].values )
    
    
        
    
        # assess drift criteria features
    
        spike_times = []
        amplitudes = []
        pc_features = []
        pc_noises = []
        firing_rates = [];   fr_bin = 10*pq.s
        for cluster_id_lists in (cluster_id_list_good,cluster_id_list_drift):
            spike_time = []        
            amplitude = []
            pc_feature = []
            pc_noise = []
            firing_rate = []
            for sx,cluster_id_list in enumerate(cluster_id_lists):   # go through the two shanks
                for u,cluster_id in enumerate(cluster_id_list):
                    print('%s shank %d, cluster %d'%(dn, sx+1,u+1))
                    # for a given sua, extract spike times, amplitudes and features
                    sua_indices = np.where(spikes_incluster[sx]==cluster_id)
                    spike_time.append( ( spikes_timing[sx][sua_indices] / T['ks2samplingrate'] ).rescale('min') )
                    amplitude.append( ( amplitude_values[sx][sua_indices] )  )
                    pc_feature.append( ( pc_feature_values[sx][sua_indices,:,0] ).squeeze()  )   # 3 PCs as a row vector for the best channel
                    aux = ( pc_feature_values[sx][sua_indices,:,1:] ).squeeze()
                    pc_noise.append( aux  )   # 3 PCs as for all background, i.e. 3 by 31 matrix
                    # calculate firing rate
                    n_bins = (spike_time[u][-1].rescale('s')[0]//fr_bin).magnitude.astype(np.int32)+1
                    fr = np.zeros(n_bins)
                    for t in range(n_bins):
                        ind = np.logical_and( spike_time[u].rescale('s')>=t*fr_bin, spike_time[u].rescale('s')<(t+1)*fr_bin )
                        fr[t] = len( spike_time[u][ind]  )/fr_bin.rescale('s')
                    fr_times = np.arange(0,n_bins*fr_bin.rescale('s').magnitude.astype(np.int32),fr_bin.rescale('s').magnitude.astype(np.int32))
                    firing_rate.append(  [fr_times,fr]  )
            
            

            spike_times.append(spike_time)
            amplitudes.append(amplitude)
            pc_features.append(pc_feature)
            pc_noises.append(pc_noise)
            firing_rates.append(firing_rate)
    
        spike_times_all.append(spike_times)
        amplitudes.append(amplitudes)
        pc_features_all.append(pc_features)
        pc_noises_all.append(pc_noises)
        firing_rates_all.append(firing_rates)
    
    # these will be [mice][{good,drift}][clusters][.....]    

    

    # print stats
    # number of clusters
    num_units = np.array( [ [ len(st) for st in mouse ]   for mouse in spike_times_all ] )
    print(num_units)
    
    print('drift stats')
    print('number of units: mean of per animal (%4.2f+/-%4.2f,%4.2f+/-%4.2f), total  (%d,%d) sum %d   (good,drift)'%\
          ( num_units.mean(axis=0)[0], (num_units.std(axis=0)*2/np.sqrt(n_mice))[0],     \
            num_units.mean(axis=0)[1], (num_units.std(axis=0)*2/np.sqrt(n_mice))[1],     \
            num_units.sum(axis=0)[0],num_units.sum(axis=0)[1], num_units.sum()   )
          )
        
        




    if 0:     # display all, to choose from

        for t in [0,1]:
            n_clusters = len(spike_times[t])
            fig,ax = plt.subplots(3,n_clusters,figsize=(6*n_clusters,12*3))
            for n in range(n_clusters):
                axs = ax[0,n]
                axs.plot(spike_times[t][n],amplitudes[t][n],'.',markersize=0.5)         #alpha=5/firing_rates[t][n][1].mean()
                # axs.set_xlim(500000,510000)
                axs.set_ylim(0,40)
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[0])
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[1])
                if n==0: axs.set_ylabel('amplitude')
                axs.set_title(pathshank[:-1]+' sua #%d'%(n+1))
    
                axs = ax[1,n]
                axs.plot(spike_times[t][n],pc_features[t][n],'.',markersize=0.5)
                axs.set_ylim(-20,20)
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[0])
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[1])
                if n==0: axs.set_ylabel('pc features')
    
                axs = ax[2,n]
                axs.plot(firing_rates[t][n][0]*pq.s.rescale('min'),firing_rates[t][n][1])
                axs.set_ylim(0,80)
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[0])
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[1])
                if n==0: axs.set_ylabel('firing rate')
    
    
    
    
            fig.suptitle(['good units','drifting units'][t])


    


    if 0:        # actual figure display
    
            panels = ['A','B','C','D']
            ind = [[0,0,1,1],[0,3,3,1]]
            titles = ['good','good','drifting','drifting']
            
            
            fig,ax = plt.subplots(5,4,figsize=(9*4,6*5))
            for i in range(4):
                t = ind[0][i]
                n = ind[1][i]
                
                # print(i,len(amplitudes[t][n]), len(pc_noises[t][n][:,0,0])*31)

                s = np.random.permutation(len(amplitudes[t][n]))[:15000]

                
                axs = ax[0,i]
                axs.plot(spike_times[t][n],amplitudes[t][n],'.',markersize=0.2)         #alpha=5/firing_rates[t][n][1].mean()
                # axs.set_xlim(500000,510000)
                axs.set_ylim(0,40)
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[0])
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[1])
                if n==0: axs.set_ylabel('amplitude')
                # axs.set_title(titles[i])
                figs.labelaxis(axs,panels[i],x=-0.2,y=1.15)

    
                for px in [0,1,2]:
                    axs = ax[1+px,i]
                    axs.plot(spike_times[t][n][s],pc_noises[t][n][s,px,:],'.',markersize=0.2,color='black',alpha=0.6)
                    axs.plot(spike_times[t][n],pc_features[t][n][:,px],'.',markersize=0.2)
                    axs.set_ylim(-20,20)
                    figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[0])
                    figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[1])
                    if n==0: axs.set_ylabel('pc features %d'%(px+1))
                    
                    axins = axs.inset_axes([1.1,0,0.3,1],transform=axs.transAxes)
                    axins.hist(pc_noises[t][n][:,px,:].ravel(),bins=np.arange(-20,20.1,0.1),orientation='horizontal',color='black',alpha=0.7)
                    axins.hist(pc_features[t][n][:,px],bins=np.arange(-20,20.1,0.1),orientation='horizontal')
                    axins.set_ylim(-20,20)
                    axins.set_xlim(0,    [5000,2000,2000,2000][i]     )
                    
                    
    
                axs = ax[4,i]
                axs.plot(firing_rates[t][n][0]*pq.s.rescale('min'),firing_rates[t][n][1])
                axs.set_ylim(0,80)
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[0])
                figs.plottoaxis_stimulusoverlay(axs,offphase=blocktimes[1])
                if n==0: axs.set_ylabel('firing rate [Hz]')
                
                axs.set_xlabel('session time [min]')
                

            
            fig.tight_layout()


            save = 0 or globalsave
            if save:
                fig.savefig(resultpath+'Supp1_driftcontrol'+ext,bbox_inches='tight')










def supplementaryfigure2():
    # show context decodable from broad spiking putative excitatory neurons only
    datanames = ['DT017','DT018','DT019','ME110','ME113']
    dn_ch = 'DT019'
    n_cherry = datanames.index(dn_ch)

    # get baseline shuffle:
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    comparison = 'context'
    celltype = 'broad'
    
    acrossdecoders_full = []
    acrossdecoders_broad = []
    for n,dn in enumerate(datanames):
        acrossdecoder_full = pickle.load(open(cacheprefix+'continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
        acrossdecoders_full.append(acrossdecoder_full)
        acrossdecoder_broad = pickle.load(open(cacheprefix+'continuous/celltyperestricted-%s,%s-%s-%s-.pck'%(comparison,celltype,dn,continuous_method),'rb'))
        acrossdecoders_broad.append(acrossdecoder_broad)


    print(acrossdecoders_broad[0][1].mean(axis=0))


    fig,ax = plt.subplots(1,2,figsize=(2*8,1*8))

    axs = ax[0]
    panel = 'A'

    figs.plottoaxis_decoderrocauc(acrossdecoders_full[n_cherry][:2],axs,colorlist=['','black'],plottrain=False,onlysem=True,smooth=[1])       # plot the performance
    figs.plottoaxis_decoderrocauc(acrossdecoders_broad[n_cherry][:2],axs,colorlist=['','mediumvioletred'],plottrain=False,onlysem=True,smooth=[1])       # plot the performance
    # axs.legend(['all units','broad spiking units'],frameon=False)
    figs.plottoaxis_stimulusoverlay(axs,T)
    figs.plottoaxis_chancelevel(axs,chances[dn_ch])
    figs.setxt(axs)
    axs.set_xlim(-1400*pq.ms,4400*pq.ms)
    axs.set_ylim(0.45,1.01)
    axs.set_yticks([0.5,1.0])
    axs.set_ylabel('context accuracy')
    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)



    axs = ax[1]
    panel = 'B'

    tpidx = np.array([0,150,300,450,596],dtype=np.int16)
    A = np.array([ acrossdecoders_broad[n][1][:,0].squeeze() for n in range(len(datanames)) ])
    B = np.array([      A[:,tpidx[tpx]:tpidx[tpx+1]].mean(axis=1)     for tpx in range(4)  ]).T
    ci = B.std(axis=0)*2/np.sqrt(len(datanames))
    m = B.mean(axis=0)
    ci = np.c_[m-ci,m+ci]
    artist = axs.boxplot(x=B, positions=-750+tpidx[:4]*10, notch=True,usermedians=m,conf_intervals=ci,\
                         whis=[5,95],widths=350)
    for element in artist.keys():
        plt.setp(artist[element], color='mediumvioletred',linewidth=2)

    figs.setxt(axs)
    axs.set_xlim(-1400*pq.ms,4400*pq.ms)
    axs.set_yticks([0.5,1.0])
    axs.set_ylim(0.45,1.01)
    figs.plottoaxis_stimulusoverlay(axs,T)
    m_c = np.array(list(chances.values())).mean()
    e_c = np.array(list(chances.values())).std()/np.sqrt(len(datanames))
    figs.plottoaxis_chancelevel(axs,m_c+e_c)

    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)


    




    fig.tight_layout()
    
    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Supp2_context,broadspikingonly'+ext)
























def supplementaryfigure3():
    # nullspace projections


    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']
    n_mice = len(datanames)
    # n_cherry = datanames.index('DT022')
    dn_ch = 'DT017'
    n_cherry = datanames.index(dn_ch)

    # get baseline shuffle:
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    fig,ax = plt.subplots(1,2,figsize=(2*8,1*8))
    
    comparison = 'context'
    
    
    
    
    panel='A'
    axs = ax[0]
    
    dn = datanames[n_cherry]
        
    acrossdecoder = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
    acrossdecoder_nullspace = pickle.load(open(cacheprefix+'subspaces/nullspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
    acrossdecoder_visualspace = pickle.load(open(cacheprefix+'subspaces/visualspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
    acrossdecoder_visualspace[1][:T['stimstart_idx'],:] = np.nan
    times = acrossdecoder_nullspace[1].times

    colors = ['mediumvioletred','darkgoldenrod','rebeccapurple']
    # labels = ['full neural space, rank N','context DV nullspace projection, rank N-1','visual DV projection, rank 1']
    labels = ['from full neural space','from projection onto context nullspace','from projection onto visual DV']
    for px,responsedecode in enumerate([acrossdecoder,acrossdecoder_nullspace,acrossdecoder_visualspace]):
        if px==1: continue
        axs.plot(times,responsedecode[1][:,0],linewidth=2,color=colors[px],alpha=0.9,label=labels[px])
        axs.fill_between(times, (responsedecode[1][:,0]-responsedecode[1][:,2]).squeeze(),\
                            (responsedecode[1][:,0]+responsedecode[1][:,2]).squeeze(),\
                            alpha=0.33,color=colors[px])
    axs.legend(frameon=False,loc='upper left')
    axs.set_xlim([times[0],times[-5]])
    axs.set_ylim([0.45,1.01])
    figs.setxt(axs)
    figs.plottoaxis_stimulusoverlay(axs,T)
    figs.plottoaxis_chancelevel(axs,chances[dn_ch])
    
    axs.set_ylabel('context dec. accuracy')
    axs.set_xlabel('time from stimulus onset')


    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)





    panel='B'
    axs = ax[1]
    
    bars_full = [] # will be (mice,projections,mean-sem)
    for n,dn in enumerate(datanames):
        acrossdecoder_full = pickle.load(open(cacheprefix+'subspaces/responsedecodes,subspaces-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
        acrossdecoder_nullspace = pickle.load(open(cacheprefix+'subspaces/nullspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))
        acrossdecoder_visualspace = pickle.load(open(cacheprefix+'subspaces/visualspacedecodes-%s_%s_%s.pck'%(dn,comparison,continuous_method),'rb'))

        
        bars = []
        for bx,acrossdecoder in enumerate([ acrossdecoder_full, acrossdecoder_nullspace, acrossdecoder_visualspace ]):
            bar = [ acrossdecoder[1][T['stimstart_idx']:T['stimend_idx'],0].mean(axis=0), acrossdecoder[1][T['stimstart_idx']:T['stimend_idx'],0].std(axis=0) ]
            bars.append(bar)
        bars_full.append(bars)
        
    bars = np.array(bars_full).squeeze()
    
    
    for bx in [0,1,2]:
        if bx==1: continue
        axs.bar(x=np.arange(n_mice)-0.15+bx*0.15,height=bars[:,bx,0]-0.5,yerr=bars[:,bx,1], width=0.3, color=colors[bx], alpha=0.75, label=labels[bx])
    axs.legend(frameon=False)
    axs.set_yticks(np.arange(0,0.6,0.1))
    axs.set_yticklabels(np.arange(0,0.6,0.1)+0.5)
    axs.set_ylim([-0.05,0.55])


    xranges = np.arange(n_mice)
    for n in range(n_mice):
        chs = chances[datanames[n]]
        figs.plottoaxis_chancelevel( axs, chs-0.5, xs=[ xranges[n]-0.5, xranges[n]+0.5] )


    axs.set_xticklabels([])
    axs.set_ylabel('on stim. context dec. acc')
    axs.set_xlabel('mice')    

    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    fig.tight_layout()


    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Supp3_context,nullspaceprojections'+ext)















def supplementaryfigure4():
    # number of neurons control
    
    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting 
    dn = 'DT014'
    singleanimalindex = datanames.index(dn)
    examine = 'allexpcond'
    comparison = 'context'
    acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
    contexttimecourse = acrossdecoder[:2]

    n_neurons_mice_dict = preprocess.getnneurons(datanames)
    n_neurons_mice = [n_neurons_mice_dict[k] for k in n_neurons_mice_dict]

    # get baseline shuffle:
    n_resample = 10
    chances = pickle.load(open(cacheprefix+'subspaces/chances,allmice,resampled-full,r%d-%s.pck'%(n_resample,continuous_method),'rb'))


    timegroups = []
    for n,dn in enumerate(datanames):
        visualtimegroups = []; contexttimegroups = []
#        timestarts_idx = int((np.arange(0,6001,1500)*pq.ms / np.array(T['dt'])).magnitude)     # cutpoints
        timestarts_idx = np.arange(0,601,150,dtype='int16')
        comparison = 'visual'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        for ti in range(4):
            visualtimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )
        comparison = 'context'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%(examine,'all',dn,continuous_method,comparison,'all'),'rb'))
        for ti in range(4):
            contexttimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )
        timegroups.append([visualtimegroups,contexttimegroups])
    timegroups = np.array(timegroups)
    
    











    variablecolors = ['navy','mediumvioletred']    
    
    panels = ['A','B','C']

    # determine partial correlation between pre on and n_neurons
    # variables will be pre context, on context and number of neurons
    variables_forpartial = np.concatenate( [timegroups[:,1,:2], np.array(n_neurons_mice)[:,np.newaxis] ], axis=1)
    r_partial,res_lsq = nedi.partialcorrelation(variables_forpartial)
    print(variables_forpartial.shape, r_partial.shape,res_lsq.shape)

    
    # plot correlations with each other
    fig,ax = plt.subplots(1,3,figsize=(3*8,1*8))
    colors = ['mediumvioletred','mediumvioletred','darkred']
    ylabels = ['decoder accuracy\ncontext PRE','decoder accuracy\ncontext early ON','residuals\ncontext early ON' ]
    xlabels = ['number of units','number of units','residuals\ncontext PRE']
    for bx in [0,1,2]:

        axs = ax[bx]           # top context
        # x = timegroups[:,1,0]           # context pre
        # y = timegroups[:,bx,1]          #   early
        if bx<2:
            x = variables_forpartial[:,2]
            y = variables_forpartial[:,bx]
        if bx==2:
            # x = variables_forpartial[:,0]
            # y = variables_forpartial[:,1]
            x = res_lsq[0,1,0,:]
            y = res_lsq[0,1,1,:]
            
        
        axs.scatter(x,y,s=150,marker='o',color=colors[bx],alpha=0.8)
        if bx==2: axs.set_xlim(-0.7,0.7); axs.set_xticks([-0.5,0,0.5])
        else: axs.set_xlim(1,np.max(variables_forpartial[:,2])+3);
        if bx==2: axs.set_ylim(-0.7,0.7); axs.set_yticks([-0.5,0,0.5]) 
        else: axs.set_ylim(0.45,1.01); axs.set_yticks([0.5,1.0])

        m_c = np.array(list(chances.values())).mean()
        e_c = np.array(list(chances.values())).std()/np.sqrt(len(datanames))
        if bx==2: figs.plottoaxis_crosshair(axs,0,0,color='black')
        else: figs.plottoaxis_crosshair(axs,y=m_c+e_c,color='black')

        axs.set_xlabel(xlabels[bx])
        axs.set_ylabel(ylabels[bx])

        l = sp.stats.linregress(x,y)
        if l[3]<0.05:
            line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],color=colors[bx],linewidth=2)
        else:
            line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],'--',color=colors[bx],linewidth=2)
            
        xoffs = [40,40,0.7][bx]
        yoffs = [0.52,0.52,-0.4][bx]
#        if sx in [3,4,5,7]: xoffs = 25
        if l[3]<0.001:
            axs.text(line[1]-xoffs,yoffs,'p<%5.3f, $\\rho$=%4.2f'%((0.001,l[2])),color=colors[bx])
        else:
            axs.text(line[1]-xoffs,yoffs,'p=%5.3f, $\\rho$=%4.2f'%((l[3],l[2])),color=colors[bx])
        
        
        
#            axs.set_title('n=%d'%14)
        # axs.spines['right'].set_visible(False)
        # axs.spines['top'].set_visible(False)
        # axs.set_title('context continuity\n%d animals'%len(datanames),pad=40)

        figs.labelaxis(axs,panels[bx])
    
    fig.tight_layout()






    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Supp4_context,controlnumberofneurons'+ext,bbox_inches='tight')
















def supplementaryfigure5():
    # show that discrimination of visual stimuli decreases: project to dbnv

    panels = ['A','B']

    datanames = ['ME103','ME110','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT020','DT021','DT022','DT030','DT031','DT032']                 # with ks2 spike sorting 
    n_mice = len(datanames)

    # activities_all will be [n_mice][comparison][performance] = 15 x 4 x 4
    activities_all,activities_aggregate = pickle.load(open('../cache/subspaces/dbnvprojections-visual,behaviourdifference-%s-%dms.pck'%(continuous_method,T['dt']),'rb'))

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


    # plot each mouse separately
    n_panels = int(np.ceil(np.sqrt(n_mice)))

    cx = 0    # choose the context conditioned block (0), or context (1) for displaying projected activities
    lx = 1     # labels for the conditions  0 visual        1 context


    
    dprimes = np.zeros((n_mice,2))   #    (n_mice x {45,135})
    for n,dn in enumerate(datanames):
        dprimes[n,:] = [ nedi.dprime_normal(activities_all[n][cx][0], activities_all[n][cx][1]),\
                         nedi.dprime_normal(activities_all[n][cx][2], activities_all[n][cx][3]) ]


    if 0:
        fig,ax = plt.subplots(n_panels,n_panels,figsize=(12*n_panels,8*n_panels))
        for n,dn in enumerate(datanames):
            print(dn,'activity len',[len(activities_all[n][cx][i]) for i in [0,1,2,3]])
    
            # print(dn,'activity vector@ hmcf',activities_all[n][cx])
    
            # find significantly different mice with behaviour
            s45,p45 = sp.stats.ttest_ind(activities_all[n][cx][0], activities_all[n][cx][1], equal_var=False)
            s135,p135 = sp.stats.ttest_ind(activities_all[n][cx][2], activities_all[n][cx][3], equal_var=False)
    
            
            axs = ax[int(np.floor(n/n_panels)),n%n_panels]
            
        
            for k in [2,3]:            # range(len(activities_all[n][cx]))   # hit, miss, corr rej, false al
                axs.hist(activities_all[n][cx][k],bins=10,color=colors[k],label=labels[lx][2+k],alpha=0.5)
            # axs.set_xlim(-2,2)
            axs.plot([0,0],[0,axs.get_ylim()[1]],'--',color='grey')
            axs.legend(frameon=False)
            # axs.set_title('%s; 45: p=%5.4f, 135: p=%5.4f'%(dn,p45,p135,))
            axs.set_title('%s; p=%5.4f'%(dn,p135))
    
    
        save = 0
        if save:
            fig.savefig(resultpath+'Supp3_visualdiscrimination,behaviourconditioned-singlemice_%s_%dms%dms'%(continuous_method,T['dt'].magnitude,T['bin'].magnitude)+ext)
            
            
            
    if 1:            
        fig,ax = plt.subplots(1,2,figsize=(2*8,1*8))

        # find out statistically if performance makes a difference in discrimination peaks
        cx = 0    # only one projection
        s45,p45 = sp.stats.ttest_ind(activities_aggregate[cx][0], activities_aggregate[cx][1], equal_var=False)
        s135,p135 = sp.stats.ttest_ind(activities_aggregate[cx][2], activities_aggregate[cx][3], equal_var=False)
        # s45,p45,s135,p135 = 0,0,0,0
        print(p45,p135)
    
        # selection = np.abs(dprimes)<1
        for d in [1]:       # separate 45 and 135 degrees
            axs = ax[0]
            for k in d*2+np.array([0,1]):
                # print(k,'hmcf'[k]+':',len(activities_aggregate[cx][k]))
                axs.hist(activities_aggregate[cx][k],bins=20,color=colors[k],label=labels[lx][2+k],alpha=0.5)
                m = np.mean(activities_aggregate[cx][k])
                e = 2*np.std(activities_aggregate[cx][k])/np.sqrt(len(activities_aggregate[cx][k]))
                axs.plot([m,m],[0,48],lw=3,color=colors[k])
                axs.plot([m-e,m+e],[10,10],lw=7,color=colors[k])
            axs.set_xlim(-3,3)
            axs.set_ylim(0,60)
            axs.plot([0,0],[0,axs.get_ylim()[1]],'--',color='grey')
            axs.legend(frameon=False)
            # axs.set_title('%d mice; %d: p=%5.4f'%(n_mice,[45,135][d],(p45,p135)[d]))
            axs.text(1,10,'p=%.2e'%((p45,p135)[d]))
        
            # axs = ax[1]
            # axs.hist( dprimes[:,d], color=colors[2*d], bins=5, label=['45 hit-miss','135 corr.rej.-fal.al.'][d], alpha=0.5 )
        axs.legend(frameon=False,fontsize='small')
        axs.set_xlabel('activity projected\nonto visual DV')
        axs.set_ylabel('# of trials')
        
        figs.labelaxis(axs,panels[0])

        
        
        axs = ax[1]
        
        for n,dn in enumerate(datanames):
            s135,p135 = sp.stats.ttest_ind(activities_all[n][cx][2], activities_all[n][cx][3], equal_var=False)
            m = [activities_all[n][cx][2].mean(),activities_all[n][cx][3].mean()]
            e = [activities_all[n][cx][2].std()*2/np.sqrt(len(activities_all[n][cx][2])),\
                 activities_all[n][cx][3].std()*2/np.sqrt(len(activities_all[n][cx][3]))]
            
            if p135>0.05:
                if m[0]<m[1]: color = 'green'
                elif m[0]>m[1]: color = 'red'
            else: color = 'grey'

            print(m,e)    
            axs.errorbar( np.array([0,1])+n*0.001, m, e, color=color, lw=2 )

        axs.set_xlim(-0.1,1.1)

        figs.plottoaxis_crosshair(axs,0)

        axs.set_xticks([0,1]); axs.set_xticklabels(['correct\nrejection','false\nalarm'])
        axs.set_yticks(np.arange(-3,3.1))
        axs.set_ylabel('# of trials')
        axs.set_ylabel('activity projected\nonto visual DV')

        figs.labelaxis(axs,panels[1])
        
        
        # for k in [0,1]:
        #     x = dprimes[:,k]
        #     y = behavfrac[:,k]

        #     axs.scatter(x,y,marker='o',s=150,color=colors[2*k],alpha=0.8,label=['hit/45','corr.rej./135'][k])
        #     # for n,dn in enumerate(datanames):
        #     #     axs.text(x,y,dn,fontsize='x-small',color=colors[2*k],alpha=0.3)
    
    
        #     l = sp.stats.linregress(x,y)
        #     print(k,l)
        #     if l[3]<0.05:
        #         line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
        #         axs.plot(line,l[0]*line+l[1],color=colors[2*k],linewidth=2)
        #     else:
        #         line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
        #         axs.plot(line,l[0]*line+l[1],'--',color=colors[2*k],linewidth=2)
                
        #     xoffs = 1
        # #        if sx in [3,4,5,7]: xoffs = 25
        #     if l[3]<0.05:
        #         axs.text(line[1]-xoffs,0.47+k*0.1,'p<%5.3f, $\\rho$=%4.2f'%((0.001,l[2])),color=colors[2*k])
        #     else:
        #         axs.text(line[1]-xoffs,0.47+k*0.1,'p=%5.3f, $\\rho$=%4.2f'%((l[3],l[2])),color=colors[2*k])
        
        # axs.legend(frameon=False)
        # axs.set_xlabel('signed d-prime')
        # axs.set_ylabel('behaviour success')
        
        # fig.suptitle(taskaspects[dbprojections[0]])
        
        fig.tight_layout()
        
        save = 0 or globalsave
        if save:
            fig.savefig(resultpath+'Supp5_visualdiscrimination,behaviourconditioned-singlemice'+ext)













# UNUSED SUPPLEMENTARY
def supplementaryfigure4_unused():



    datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    
    # load rate statistics
    rss = pickle.load(open('../cache/phys/rate,trajectories,celltypes,allmice-%s.pck'%(continuous_method),'rb'))
    # rss dimensions: [n,cx,lidx,cidx][trials][trajectory,neurons]
    
    dn = 'DT017'
    n = 7
    
    rs = rss[n,:,:,:]
    layerlist = [1,2,3,4]       # choose which layer to focus on: 0 all, 1 superf, 2 input, 3 deep5, 4 deep6
    celltypelist = [0,1]        # choose which spike width to use: 0 narrow, 1 broad
    
    
    # data to assess:
    # pyr: {all trials} x {trajectories} x {pyrcells}
    # inh: {all trials} x {trajectories} x {inhcells}
    
    wx = int((1000*pq.ms/T['dt']).magnitude)
    
    fr_p = []
    # for trx in range(len(rs[0,0,0]))
    # pyr:
    aux1 = np.concatenate( [ np.array(rs[0,k,1]) for k in range(rs.shape[2])  ], axis=2 )
    aux2 = np.concatenate( [ np.array(rs[1,k,1]) for k in range(rs.shape[2])  ], axis=2 )
#    fr_p = pq.Quantity( np.concatenate( [aux1,aux2]), units = pq.ms )
    fr_p = np.concatenate( [aux1,aux2] )

    # inh:
    aux1 = np.concatenate( [ np.array(rs[0,k,0]) for k in range(rs.shape[2])  ], axis=2 )
    aux2 = np.concatenate( [ np.array(rs[1,k,0]) for k in range(rs.shape[2])  ], axis=2 )
#    fr_i = pq.Quantity( np.concatenate( [aux1,aux2]), units = pq.ms )
    fr_i = np.concatenate( [aux1,aux2] )
    
    # fr statistics    
    m_p = fr_p.mean(axis=(0,2))
    sem_p = fr_p.std(axis=(0,2))/np.sqrt(fr_p.shape[0]*fr_i.shape[2])
    m_i = fr_i.mean(axis=(0,2))
    sem_i = fr_i.std(axis=(0,2))/np.sqrt(fr_i.shape[0]*fr_i.shape[2])
    
    # histogram data
    h_p_sa = fr_p[:,T['stimstart_idx']-wx:T['stimstart_idx'],:].mean(axis=1).reshape(-1)
    h_p_ea = fr_p[:,T['stimstart_idx']:T['stimstart_idx']+wx,:].mean(axis=1).reshape(-1)
    h_i_sa = fr_i[:,T['stimstart_idx']-wx:T['stimstart_idx'],:].mean(axis=1).reshape(-1)
    h_i_ea = fr_i[:,T['stimstart_idx']:T['stimstart_idx']+wx,:].mean(axis=1).reshape(-1)
    
    







    # D  cross decoders: across trajectory and blocks
    dn = 'DT017'
    comparison = 'context'
    trainingoffsets = np.array([0,375,750,1225])*pq.ms
    trainingpoints = T['starttime'] + trainingoffsets
    trainingoffsetidx = [[0,27],[27,27+27]]
    trx = 3
    offstimdecoders = []
    for trx,toff in enumerate(trainingpoints):
        offstimdecoders.append(pickle.load(open('../cache/continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,comparison,trx),'rb')))
    offstimdecoders_sa = []
    offstimdecoders_ea = []
    for n,dn in enumerate(datanames):
        trx = 2
        aux = pickle.load(open('../cache/continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,comparison,trx),'rb'))
        aux = aux[1][0:27,0].mean()
        offstimdecoders_sa.append(aux)
        trx = 3
        aux = pickle.load(open('../cache/continuous/offstimdecodes-%s-%s-%s_t%d.pck'%(dn,continuous_method,comparison,trx),'rb'))
        aux = aux[1][27:27+27,0].mean()
        offstimdecoders_ea.append(aux)



    # E behaviour
    datanames = ['ME108','ME110','ME112','ME113','DT008','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    n_mice = len(datanames)
    
    Xa = []
    Ya = []
    perf = np.zeros((n_mice,6))
    expertise = np.zeros(n_mice)

    for n,dn in enumerate(datanames):
        print(dn)
        a,b,b_s,_,_,_ = preprocess.loadtrainingbehaviouraldata(dn)       # a: all data, b: behav, b_s behav sem.
        perf[n,:] = [ b[1],b[3], b[5], 2*b_s[1], 2*b_s[3], 2*b_s[5] ]      # behaviour and its s.e.m.


        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        expertise[n] = acrossdecoder[1][ T['start_idx']:T['stimstart_idx'], 0 ].mean()

        Ya.append( [a[1],a[3], np.concatenate( [a[1],a[3]] ) ] )
        Xa.append( [ expertise[n]*np.ones(len(Ya[-1][0])), expertise[n]*np.ones(len(Ya[-1][1])), expertise[n]*np.ones(len(Ya[-1][2])) ] )


    Yl = [ np.concatenate( [  Ya[n][d] for n in range(n_mice)  ] ) for d in range(3) ]
    Xl = [ np.concatenate( [  Xa[n][d] for n in range(n_mice)  ] ) for d in range(3) ]









    # FIGURE




    f1n=2; f2n=4     # full span of panel grid  vertical and horizontal

    fig = plt.figure(constrained_layout=False,figsize=(f2n*8,f1n*8))
    gs = fig.add_gridspec(f1n, f2n)
    ax = []
    ax_sides = []
    for j in np.arange(f1n):
        axa=[]
        for k in np.arange(f2n):
            if j==1 and k in [1,2]:
                gs_double = gs[j,k].subgridspec(1, 4)
                ax_main = fig.add_subplot(gs_double[0,:-1])
                ax_sides.append(fig.add_subplot(gs_double[0,-1]))
                axa.append( ax_main )
            else:
                axa.append( fig.add_subplot(gs[j, k]) )
        ax.append(axa)
    ax = np.array(ax)











    




    panels = ['A','B']
    times = np.arange(-1500*pq.ms,4501*pq.ms,T['dt'])
    
    axs = ax[0,0]
    figs.plottoaxis_plottrajectory(axs,times,m_p,2*sem_p,'k',linewidth=2)
    axs.legend(['broad spiking cells'])
    figs.labelaxis(axs,panels[0])

    axs = ax[0,1]
    figs.plottoaxis_plottrajectory(axs,times,m_i,2*sem_i,'r',linewidth=2)
    axs.legend(['narrow spiking cells'])

    axs = ax[0,2]
    axs.hist( (h_p_ea-h_p_sa).squeeze(),bins=50,density=True,color='k')
    axs.legend(['broad spiking\ncells EA-SA'])
    figs.labelaxis(axs,panels[1])

    axs = ax[0,3]
    axs.hist( (h_i_ea-h_i_sa).squeeze(),bins=50,density=True,color='r')
    axs.legend(['narrow spiking\ncells EA-SA'])

    for ix in [0,1]:
        axs = ax[0,ix]
        axs.set_xlim(-1200,4300)
        axs.set_ylim(0,11)
        figs.setxt(axs)
        axs.set_yticks([0,5,10])
        figs.plottoaxis_stimulusoverlay(axs,T)
        axs.set_ylabel('[Hz]')
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)

        axs = ax[0,ix+2]
        axs.set_xlim(-20,20)
        axs.set_ylim(0,0.20)
        axs.set_xticks([-10,0,10])
        axs.set_yticklabels([])
        axs.set_xlabel('[Hz]')
        axs.set_ylabel('hist. freq. [AU]')
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)










    panel = 'C'
    title = 'decoding schematics'   #,'on time\nparadigm','across time, within trial\nparadigm','across time across block\nparadigm']

    axs = ax[1,0]
    drawschematics(axs,0)
    
    axs.set_title(title)
    axs.spines['left'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    axs.spines['bottom'].set_visible(False)
    axs.set_xticks([])
    axs.set_yticks([])

    figs.labelaxis(axs,panel)

    










    panel = 'D'
    for trx in [2,3]:
        axs = ax[1,1+trx-2]
        
        colors = ['indigo','mediumvioletred']
#        ft = slice()
#        figs.plottoaxis_offdecoderrocauc(offstimdecoders[trx],axs,T,colors=colors)

        axs.set_ylim(0.45,1.001)
        axs.set_xlim(-1500,2400)
        ylim = axs.get_ylim()

        # plot and indicate training:

        axs.plot(offstimdecoders[trx][0].times,offstimdecoders[trx][0][:,0],linewidth=5,color=colors[0],alpha=1)

        axs.fill_between([ offstimdecoders[trx][0].times[0], offstimdecoders[trx][0].times[-1],],\
                        [ylim[0],ylim[0]],[ylim[1],ylim[1]],color=colors[0],alpha=0.33)

        labels=['1','2','3','4','5']
        # highlight test area
        ft = slice(trainingoffsetidx[3-2][0],trainingoffsetidx[3-2][1])
        tt = 1
        t = offstimdecoders[trx][tt].times
        axs.plot(t[ft],offstimdecoders[trx][tt][ft,0],linewidth=2,color=colors[tt],label=labels[tt])
        axs.fill_between(t[ft], (offstimdecoders[trx][tt][ft,0]-offstimdecoders[trx][tt][ft,2]).squeeze(),\
                            (offstimdecoders[trx][tt][ft,0]+offstimdecoders[trx][tt][ft,2]).squeeze(),\
                            alpha=0.4,color=colors[tt])
        axs.fill_between(t[ft], (offstimdecoders[trx][tt][ft,0]-offstimdecoders[trx][tt][ft,1]).squeeze(),\
                            (offstimdecoders[trx][tt][ft,0]+offstimdecoders[trx][tt][ft,1]).squeeze(),\
                            alpha=0.2,color=colors[tt])
        axs.plot(t,offstimdecoders[trx][tt][:,0],linewidth=1,color=colors[tt],label=labels[tt],alpha=0.2)
        axs.fill_between(t, (offstimdecoders[trx][tt][:,0]-offstimdecoders[trx][tt][:,2]).squeeze(),\
                            (offstimdecoders[trx][tt][:,0]+offstimdecoders[trx][tt][:,2]).squeeze(),\
                            alpha=0.1,color=colors[tt])
        axs.fill_between(t, (offstimdecoders[trx][tt][:,0]-offstimdecoders[trx][tt][:,1]).squeeze(),\
                            (offstimdecoders[trx][tt][:,0]+offstimdecoders[trx][tt][:,1]).squeeze(),\
                            alpha=0.03,color=colors[tt])
        
        
        
        
        
        
    #    axs.legend(['train','test'],loc=4)
        figs.plottoaxis_stimulusoverlay(axs,T)
        figs.plottoaxis_chancelevel(axs,0.5,'chance level')
        axs.set_xticks([0,2000])
        axs.set_yticks([0.5,1])
        
        axs.text(trainingpoints[trx]-550*pq.ms,0.98,'train',color=colors[0])
        axs.text(trainingpoints[trx]+500*pq.ms,0.65,'test',color=colors[1])

        axs.set_xlabel('[ms]')
        
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)
        
        
#    axins = inset_axes(axs, width="40%", height="40%", loc=6)

        axins = axs.inset_axes([-0.4,0.55, 0.4/3*4,0.4],transform=axs.transAxes)
        drawschematics(axins,5+trx-2)
        axins.axis('off')

        if trx==2: figs.labelaxis(axs,panel,x=-0.37)
  

    # context cross decoder: tasks

    for cx in range(2):
        axs = ax_sides[cx]
        artist = axs.boxplot([offstimdecoders_sa,offstimdecoders_ea][cx],notch=True,widths=0.3,\
                             labels=[['PRE->PRE','PRE->ON'][cx]])
    
        for element in artist.keys():
            plt.setp(artist[element], color=colors[1],linewidth=2)
    
        axs.set_xticks([])
        axs.set_yticks([0.5,1])
        axs.set_yticklabels([])
    
        axs.set_ylim(0.45,1.001)
        figs.plottoaxis_chancelevel(axs,0.5)
#        axs.set_title('n=%d'%len(datanames))
    
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)







    # J, behaviour

#    norminv(fraction trials	with	correct	licking)	– norminv(fraction	trials	with	incorrect	licking
    panel = 'E'
    colors = ['darkorange','darkred','firebrick']
    ylabels = 'error fraction'

    for dx in [0,1]:     # for misses+fals do both separately
        axs = ax[1,3]
        
        x = np.array([ np.mean(Xa[n][dx]) for n in range(len(Xa)) ])
        y = np.array([ np.mean(Ya[n][dx]) for n in range(len(Xa)) ])
        
        axs.scatter(x,y,s=150,marker='o',color=colors[dx],alpha=0.8,label=['miss','false alarm','miss & false alarm'][dx])
        for n in range(n_mice): axs.errorbar(x[n],y[n],yerr=perf[n,dx+3],color=colors[dx],alpha=0.8)
        axs.set_xlim(0.45,1.01)
#        axs.set_ylim(-0.04,1.04)
        axs.set_ylim(-0.02,0.5)
        if dx==0: figs.plottoaxis_crosshair(axs,0.5,0)
        axs.set_xlabel('context PRE', labelpad=0.1, verticalalignment='bottom')

        axs.set_xticks([0.5,1])
        axs.set_yticks([0,0.5])
        axs.set_yticklabels(['0','0.5'])
        

        l = sp.stats.linregress(Xl[dx],Yl[dx])
        if l[3]<0.05:
            line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],color=colors[dx],linewidth=2)
        else:
            line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],'--',color=colors[dx],linewidth=2)
            
        xoffs = 0.05  #0.35
#        yoffs = y.mean()*1.4 - 0.06     #/2.5 + 0.35
        yoffs = [0.05, 0.36, 0.27][dx]
        if l[3]<1e-5:
            axs.text(line[0]+xoffs,yoffs,'p<%s, slope=%5.3f, R$^2$=%6.4f'%(('10e$^{-5}$',l[2],l[2]**2)),color=colors[dx])
        elif l[3]<0.001:
            axs.text(line[0]+xoffs,yoffs,'p<%5.3f, slope=%5.3f, R$^2$=%6.4f'%((0.001,l[2],l[2]**2)),color=colors[dx])
        else:
            axs.text(line[0]+xoffs,yoffs,'p=%5.3f, slope=%5.3f, R$^2$=%6.4f'%((l[3],l[2],l[2]**2)),color=colors[dx])
        
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)

    axs.legend(frameon=False)
#    axs.set_title('%d mice, %d trials'%(n_mice,np.sum(Xl[2])))
    axs.set_ylabel(ylabels, labelpad=0, verticalalignment='top')
    figs.labelaxis(axs,panel)





    panel = 'J' # multiple mouse visual-context subspaces
    
    vib = basisvectors.index(0)   # look for visual coordinate to separete into two figures
    cxb = basisvectors.index(1)   # look for context coordinate to separete into two figures

    labels = ['visual comp. orth. to context','context']
    bp = np.array([[0,1],[0,2],[1,2]]) # basis pairs for the 2D unfolding
    b=0   # choisen coordinate axes: visual and context
    # use the same pre/on from 'G' panel
    for n,dn in enumerate(datanamesefg):
        # if not n in [9,10,11,12,13]: continue

        dbnbasiscomps = dbnbasiscomps_all[n]
        projections = projections_all[n]

        # decide on + or - corr between vi and cx:
        sx = dbnbasiscomps[bp[b,1],vib]>0
        axs = ax[2,1+1+sx]
        
        

        # In the following plots we have to flip the vectors to the upper-right quadrant (np.sign mirroring of coordinates).
        # Uncomment below to display the flipper coordinates (x component of visual, y of context):
        #        print(dn,'        vi:',np.sign(dbnbasiscomps[bp[b,0],0]),np.sign(dbnbasiscomps[bp[b,1],0]),\
        #                 '        cx:', np.sign(dbnbasiscomps[bp[b,0],1]),np.sign(dbnbasiscomps[bp[b,1],1]))
        
        L = np.max(  np.hstack(  [ np.abs(projections[j,preon,k]) for j in range(3) for k in range(4)  ]   ) ) * 1.
        

    
        # plot the attend and the ignore
        for k in [0,1,2,3]:
            axs.plot(projections[bp[b,0],preon,k],\
                      projections[bp[b,1],preon,k],\
                      'o',color=activitycolors[k],alpha=0.1)
    
        # plot the cov ellipsoids
        for k in range(4):
            figs.confidence_ellipse(projections[bp[b,0],preon,k],\
                                    projections[bp[b,1],preon,k],\
                                    axs, n_std=2.0, facecolor=activitycolors[k],alpha=0.05)

        axs.set_xlabel(labels[bp[b,0]])
        axs.set_ylabel(labels[bp[b,1]])
    

# create the legend
#    axins = axs.inset_axes([-0.35,0.55, 0.4/3*4,0.4],transform=axs.transAxes)
#    axins.axis('off')
#
#    xl0 = 0.2; yl0 = 0.0
#    for j in range(2):
#        for k in range(2):
#            axins.add_patch(plt.Rectangle((xl0+0.05*j,yl0+0.05*k),0.05,0.05,\
#              ec=activitycolors[j+2*k],fc=activitycolors[j+2*k],alpha=0.2,transform=axs.transAxes))
#            axins.add_patch(plt.Circle((xl0+0.025+0.05*j,yl0+0.025+0.05*k),0.01,\
#              ec=activitycolors[j+2*k],fc=activitycolors[j+2*k],alpha=0.4,transform=axs.transAxes))
#    axins.text(xl0+0.05,  yl0+0.1,'45°',ha='right',va='bottom',transform=axs.transAxes)
#    axins.text(xl0+0.05,  yl0+0.1,'135°',ha='left',va='bottom',transform=axs.transAxes)
#    axins.text(xl0+0.125, yl0+0.075,'ignore',ha='left',va='center',transform=axs.transAxes)
#    axins.text(xl0+0.125, yl0+0.025,'attend',ha='left',va='center',transform=axs.transAxes)


    axs = ax[2,2]

    axs.text(0.04, -0.04,'visual DV [SU]', ha='left',   va='center', transform=axs.transAxes)
    axs.text(-0.04, 0.04,'context DV [SU]', ha='center', va='bottom',rotation=90, transform=axs.transAxes)


    axs.set_xlim(axs.get_xlim())
    axs.set_ylim(axs.get_ylim())


    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    axs.axis('off')



    for n,dn in enumerate(datanamesefg):
        dbnbasiscomps = dbnbasiscomps_all[n]
        # sx = dbnbasiscomps[bp[b,1],vib]>0
        # plot skewed basis:
        for r,color in zip(basisvectors,basiscolors):
            if not r in [0,1]: continue
            x0 = axs.get_xlim()[0]+0.03
            y0 = axs.get_ylim()[0]+0.03
            axs.plot( [x0,x0+(dbnbasiscomps[bp[b,0],r])],\
                      [y0,y0+(dbnbasiscomps[bp[b,1],r])],\
                      color=color,linewidth=4,alpha=0.9 )

    axs = ax[2,2]
    axs.set_xlim(axs.get_xlim())
    axs.set_ylim(axs.get_ylim())

    figs.labelaxis(axs,panel)








    ax[1,1].axis('off')


    fig.suptitle('Supplementary Figure 3')








    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Supp3_fr,evoked+crosstimedecoders+glm'+ext)





























def supplementaryfigure5_unused():
    
    
    datanames = ['ME108','ME110','ME112','ME113','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']



    # FIGURE


    fig,ax = plt.subplots(3,4,figsize=(36,24))



    # behaviour single:

#    dn = 'DT014'
#    
#    blv,bla = preprocess.getorderattended(dn)
#    windowwidth = 100*pq.ms
#    timeaveragewidth = 100*pq.ms
#    timeaverage_range = ((np.array([windowwidth+T['dt'], windowwidth+timeaveragewidth+T['dt']])/T['dt']).magnitude).astype('int16')
#    
#    n_max_trial = 30
#    triallist = np.arange(n_max_trial)+1
#    titles = ['pre stimulus','on stimulus']
#    colorlists = [ ['darkorange','chocolate'], ['tomato','firebrick'] ]
#
#
#    _,_,_,_,cuetrainig,cue_est = preprocess.loadtrainingbehaviouraldata(dn)
#
#
#    
#
#    for bl in [0,1]:
#        blockorder = int(([blv[0],bla[0]][bl]-1)/2)
#        ptrx = 0 # pre/on             trainingpoint in enumerate(trainingpoints):
#
#        # behavioural
#        axs = ax[0, blockorder]
#        axs.plot(triallist,cue_est[bl,0,:],lw=2,color='darkblue')
##        axs.fill_between(np.arange(30)+1,cue_est[bl,0,:]-cue_est[bl,2,:],cue_est[bl,0,:]+cue_est[bl,2,:],\
##                         color='darkblue',alpha=0.2)     # this was with symmetric variance estimate
#        axs.fill_between(triallist,cue_est[bl,3,:],cue_est[bl,5,:],\
#                         color='darkblue',alpha=0.4)      # this is with proper beta estimate
#        axs.fill_between(triallist,cue_est[bl,2,:],cue_est[bl,4,:],\
#                         color='darkblue',alpha=0.2)      # this is with proper beta estimate
#        axs.set_xticks([1,15,30])
#        axs.set_ylim(-0.04,1.04)
#        axs.set_yticks([0,0.5,1])
#        if bl==0: axs.set_ylabel('across sessions performance')
#        axs.set_title('%s'%( ['initial','transition'][blockorder]  ))
#        if bl==0:
#            figs.labelaxis(axs,panels[0])
#            axs.legend(['MLE/MAP','67% h.dens.cred.int.','95% h.dens.cred.int.'],frameon=False)
#
#
#
#
#        # this is to show if there is a trend
#        axs = ax[0,blockorder]
#        color = 'green'
#        n_samples = 1000
#        x = np.concatenate([ np.ones(n_samples)*(t+1) for t in range(n_max_trial) ])
#        y = []
#        for t in range(n_max_trial):
#            y.append(sp.stats.beta.rvs(cue_est[bl,6,t],cue_est[bl,7,t],size=n_samples))
#        y = np.concatenate(y)
#
##        axs.plot(x,y,'o',alpha=0.15)        # uncomment to show the sampling
#
#        l = sp.stats.linregress(x,y)
#        if l[3]<0.05:
#            line = np.linspace(start=x.min(),stop=x.max(),num=2)
#            axs.plot(line,l[0]*line+l[1],color=color,linewidth=2)
#        else:
#            line = np.linspace(start=x.min(),stop=x.max(),num=2)
#            axs.plot(line,l[0]*line+l[1],'--',color=color,linewidth=2)
#            
#        xoffs = 22
#    #        if sx in [3,4,5,7]: xoffs = 25
#        if l[3]<0.001:
#            axs.text(line[1]-xoffs,0.17,'intercept %5.3f, slope %5.3f\np<%5.3f, $R^2$=%4.2f'%((l[1],l[0],0.001,l[2]**2)),color=color)
#        else:
#            axs.text(line[1]-xoffs,0.17,'intercept %5.3f, slope 5.3f\np=%5.3f, $R^2$=%4.2f'%((l[1],l[0],l[3],l[2]**2)),color=color)
#
#















# behaviour multiple mouse


    panel = 'A'
    datanames = ['ME108','ME110','ME112','ME113','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    n_mice = len(datanames)

    n_max_trials = 30
    
    aucdeltat = 750*pq.ms     # for the order length

    crossgradients = 0.5*np.ones( (2,2,n_mice,n_max_trials) )      # starting values for empty trials, this is raw trial numbers
#    crossgradients2 = 0.5*np.ones( (2,2,n_mice,trialcompress) )      # starting values for empty trials, this is timecompressed
    blockorder = np.zeros((n_mice,2),dtype='int16')      # this will mean 0: attend visual first, 1: attend audio first
    cx_decauc_pre = np.zeros((n_mice))       # this is context PRE: timeaveraged decoder accuracy
    cuetrainings = []
    cue_ests = []

    
    

    for n,dn in enumerate(datanames):
    
        export_crossgradient = pickle.load(open('../cache/continuous/offstimgradient,cue,decodeprobs,savgol_%s-%s.pck'%(dn,continuous_method),'rb'))
#        export_crossgradients[n][bl][ptrx,:,0]
        
        blv,bla = preprocess.getorderattended(dn)
        
        for bl in [0,1]:          # block visual, audio
            n_cuetrials = export_crossgradient[bl].shape[1]
            blockorder[n,bl] = int(([blv[0],bla[0]][bl]-1)/2)
            for ptrx in [0,1]:    # pre, on
                crossgradients[blockorder[n,bl],ptrx,n,:n_cuetrials] = export_crossgradient[bl][ptrx,:n_max_trials,2]


        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # lower confidence bound
        cx_decauc_pre[n] = neph.get_auctimenormalized(decsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]        



        _,_,_,_,cuetraining,cue_est = preprocess.loadtrainingbehaviouraldata(dn)

        cuetrainings.append(cuetraining)
        # [mice][init/trans,sessions,h+c/m+f,trials]           =(14)(2,8,2,30)

        cue_ests.append(cue_est)       # ix 6,7 (from 0) will contain the beta parameters



    # display aggregate behaviour and neural context decoder probability for all mice
    resolution = 100
    n_max_trial = 30
    window = 5             # first few and last few trials
    triallist = np.arange(n_max_trial)+1
    probmap = np.zeros((resolution+1,n_max_trial))
    mle = np.zeros(n_max_trial)
    x = np.linspace(0,1,resolution+1)    # beta probability mesh
    probmapwindowaverage = np.zeros( (2,2,resolution+1) )      # blockorder, first/last, probability density
    probmapwindowaverage_s = np.zeros( (2,2,resolution+1) )      # same, just s.e.m.
    R = np.zeros((2,n_mice,2))      # blockorder, mouse, tangent/intercept
    labelcount = [0,0,0]
    for bl in [0,1]:

        axs = ax[0,bl]
        
        for n,dn in enumerate(datanames):
#            if n!=11: continue
#            print(datanames[n])

#            print(cue_ests[n])
            
            for t in range(n_max_trial):
#                print('beta a b', cue_ests[n][bl,6,t], cue_ests[n][bl,7,t])
                
#                if cue_ests[n][bl,6,t]==np.nan: print(dn,bl,t,'nan a')
#                if cue_ests[n][bl,7,t]==np.nan: print(dn,bl,t,'nan b')
                
                # collect the posterior values at the resolution:
#                posterior = sp.stats.beta.cdf(x[1:], cue_ests[n][bl,6,t], cue_ests[n][bl,7,t] ) - \
#                            sp.stats.beta.cdf(x[:-1], cue_ests[n][bl,6,t], cue_ests[n][bl,7,t] )
                            
                posterior = sp.stats.beta.pdf( x, cue_ests[n][bl,6,t], cue_ests[n][bl,7,t] )
#                print(posterior)
#                posterior[posterior==np.nan]=0
                
                probmap[:,t] += posterior
            
            
            
            
            # regression to test:
#            n_samples = 100
#            x_list = np.concatenate([ np.ones(n_samples)*(t+1) for t in range(n_max_trials) ])
#            y = []
#            for t in range(n_max_trials):
#                samples = sp.stats.beta.rvs(cue_ests[n][bl,6,t],cue_ests[n][bl,7,t],size=n_samples)
#                y.append(samples)
#                
#            y = np.concatenate(y)
#    
##                axs.plot(x_list,y,'o',alpha=0.15,color=color)        # uncomment to show the sampling
#    
#            l = sp.stats.linregress(x_list,y)
#            line = np.array([triallist[0],triallist[-1]])
#            negativeslope = (l[0]<0).astype('int16')
#            color = ['chartreuse','salmon'][ negativeslope ]
#            label = ['+ slope',' - slope'][ negativeslope ]
#            if labelcount[negativeslope]>0: label=None
#            if l[3]<0.05:
#                axs.plot(line,l[0]*line+l[1],color=color,linewidth=1.5,label=label)
#                labelcount[negativeslope] += 1
#            else:
#                color = 'gold'
#                label = ' 0 slope'
#                if labelcount[2]>0: label=None
#                labelcount[2] += 1
#                axs.plot(line,l[0]*line+l[1],color=color,linewidth=1.5,label=label)
#            R[bl,n,:] = l[:2]




    
        probmap = probmap / n_mice / n_max_trial                 # make an actual average/mixture

        # save the first and last density means
        probmapwindowaverage[bl,0,:] = probmap[:,:window].mean(axis=1)
        probmapwindowaverage[bl,1,:] = probmap[:,-window:].mean(axis=1)
        probmapwindowaverage_s[bl,0,:] = 2*probmap[:,:window].std(axis=1)/np.sqrt(window)
        probmapwindowaverage_s[bl,1,:] = 2*probmap[:,-window:].std(axis=1)/np.sqrt(window)

        cumprobmap = np.cumsum(probmap,axis=0) / np.sum(probmap,axis=0)#          * (x[1]-x[0])
        
        
        # get contours for the mle and highest density credible intervals 67 and 95%
        c_l = np.zeros((2,n_max_trial))
        c_u = np.zeros((2,n_max_trial))
        mle = np.zeros(n_max_trial)
        for t in range(n_max_trial):
            mle[t] = x[np.argmax(probmap[:,t])]
            c_l[0,t], c_u[0,t] = neba.hdci(ppf=lambda q: neba.empiricalppf(q,cumprobmap[:,t],x), alpha=19/20,resolution=resolution)  # 95% hdci
            c_l[1,t], c_u[1,t] = neba.hdci(ppf=lambda q: neba.empiricalppf(q,cumprobmap[:,t],x), alpha=2/3,resolution=resolution)  # 67% hdci


#        print(c_l,c_u)
#        s = axs.imshow(probmap,aspect='auto',origin='lower',cmap='viridis',vmax=0.2,vmin=0)
        mX,mY = np.meshgrid( triallist, x )
        s = axs.contourf(mX,mY,probmap,levels=100,vmin=0,vmax=0.15)
        axs.plot(triallist,mle,lw=2,color='white',label='MLE/MAP')
        axs.plot(triallist,c_l[1,:],'--',color='white',label='67% h.dens.cred.int.')
        axs.plot(triallist,c_u[1,:],'--',color='white')
        axs.plot(triallist,c_l[0,:],'--',color='grey',label='95% h.dens.cred.int.')
        axs.plot(triallist,c_u[0,:],'--',color='grey')
        
        if bl==0:
            leg = axs.legend(frameon=False)
            for text in leg.get_texts(): text.set_color('white')
        
        axs.set_title('%s'%( ['initial\nbehavioural performance','transition\nbehavioural performance'][bl]  ))
        axs.set_xticks([1,15,30])
        axs.set_yticks([0,0.5,1])
        axs.set_yticklabels(['0','0.5','1'])
#        axcb = fig.colorbar(s,ax=axs)#,ticks=[0,0.1,0.2])

    
    
        if bl==0:
            axs.set_ylabel('posterior across sessions and mice')
            figs.labelaxis(axs,panel)
    
    
    
    





















    panel = 'B'

    # load context PRE:
    contextpre = []
    for n,dn in enumerate(datanames):
        contexttimegroups = []
#        timestarts_idx = int((np.arange(0,6001,1500)*pq.ms / np.array(T['dt'])).magnitude)     # cutpoints
        timestarts_idx = np.arange(0,601,150,dtype='int16')
        comparison = 'context'
        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,comparison,'all'),'rb'))
        for ti in range(4):
            contexttimegroups.append(      acrossdecoder[1][ timestarts_idx[ti]:timestarts_idx[ti+1], 0 ].mean()   )
        contextpre.append(contexttimegroups[0])

    contextpre = np.array(contextpre)







    # collect regression behavioural data as long lists for each animal (i.e. at each context PRE value)
    initials = []
    transitions = []
    for n,dn in enumerate(datanames):
        initials.append( [ np.reshape(cuetrainings[n][0][:,0,:window],-1),\
                           np.reshape(cuetrainings[n][0][:,0,-window:],-1) ]   )
        transitions.append( [ np.reshape(cuetrainings[n][1][:,0,:window],-1),\
                              np.reshape(cuetrainings[n][1][:,0,-window:],-1) ] )

        
    X_list = [[],[]]       # as below: context pre values
    Y_list = [[],[]]       # these will store the per session first and last data for regression
    Y_m_list = []
    for bl in [0,1]:          # initial and transition block
        period = [initials,transitions][bl]         # shortcut call
        
        Y = []
        Ys = []
        for i in [0,1]:          # first and last, on the same figure
            

            
            # behaviour for all mice
            X_list[i].append(   np.concatenate( [ [  contextpre[n] for d in range( len(period[n][i]) ) ] \
                               for n in range(n_mice)  ] )  )
            Y_list[i].append( np.concatenate(  [ [  period[n][i][d] for d in range( len(period[n][i]) ) ] \
                               for n in range(n_mice)  ]  ) )
                    # Y_list will be Y_list[block][first/last][sessions*window]

            
            Y.append(  [ (period[n][i]).mean()  \
                               for n in range(n_mice)  ]  )
            Ys.append(  [ (period[n][i]).std()/np.sqrt( len(period[n][i]) )  \
                               for n in range(n_mice)  ]  )
        Y_m_list.append(Y)         # this will be [bl][i]

    # uses Y_list from above
    # Y_list is Y_list[block][first/last][sessions*window]
    colors = ['seagreen','mediumblue']

    for bl in [0,1]:    # initial and transition block
        for i in [0,1]: # first and last
            axs = ax[0,bl+2]
            x = X_list[i][bl]
            y = Y_list[i][bl]
            y_m = Y_m_list[bl][i]
                        # or using: #   R[bl,:,1]        # 0 slope 1 intercept
            axs.scatter(contextpre,y_m,s=150,marker='o',color=colors[i],alpha=0.8)
            axs.set_xlim(0.49,1.01)
            axs.set_ylim(0.49,1.01)
            axs.set_xticks([0.5,1.0])
            axs.set_yticks([0.5,1.0])
            axs.set_title(['initial\n','transition\n'][bl])
            axs.set_xlabel('context PRE', labelpad=0.1, verticalalignment='bottom')
    #        axs.set_ylabel(ylabels[bx], labelpad=0, verticalalignment='top')
    
            l = sp.stats.linregress(x,y)
            if l[3]<0.05:
                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                axs.plot(line,l[0]*line+l[1],color=colors[i],linewidth=2)
            else:
                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                axs.plot(line,l[0]*line+l[1],'--',color=colors[i],linewidth=2)
                
            xoffs = 0.35
    #        if sx in [3,4,5,7]: xoffs = 25
            if l[3]<0.001:
                axs.text(line[1]-xoffs,line[1]*0.49-i/30,'p<%5.3f, $R^2$=%4.2f'%((0.001,l[2]**2)),color=colors[i])
            else:
                axs.text(line[1]-xoffs,line[1]*0.49-i/30,'p=%5.3f, $R^2$=%4.2f'%((l[3],l[2]**2)),color=colors[i])
            
        if bl==0:
            axs.legend(['first %d trials'%window,'last %d trials'%window], frameon=False,loc=2)
            axs.set_ylabel('across sessions MAP per mice',verticalalignment='top',labelpad=0)
            figs.labelaxis(axs,panel)
#            axs.set_title('n=%d'%14)
        figs.plottoaxis_crosshair(axs,0.5,0.5)
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)

















    # neural data
    # multiple mouse


    panel = 'D'
    n_mice = len(datanames)

    trialcompress = 30
    
    aucdeltat = 750*pq.ms     # for the order length

    crossgradients = 0.5*np.ones( (2,2,n_mice,100) )      # starting values for empty trials
    crossgradients2 = 0.5*np.ones( (2,2,n_mice,trialcompress) )      # starting values for empty trials
    blockorder = np.zeros((n_mice,2),dtype='int16')      # this will mean 0: attend visual first, 1: attend audio first
    decauc_pre = np.zeros((n_mice))
    
    ls = np.zeros((2,2,n_mice,2))
    for n,dn in enumerate(datanames):
    
        export_crossgradient = pickle.load(open('../cache/continuous/offstimgradient,cue,decodeprobs,savgol_%s-%s.pck'%(dn,continuous_method),'rb'))
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


        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # lower confidence bound
        decauc_pre[n] = neph.get_auctimenormalized(decsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]        
        
                
    
    orders = np.argsort(ls[1,1,:,1])[::-1]     # cue 3, onstimulus,   all mice,    slope/intercept
    orders = np.argsort(ls[1,1,:,1]+ls[1,1,:,0]*15+crossgradients2[blockorder[n,1],1,n,-1]+\
                        crossgradients2[blockorder[n,0],1,n,-1]\
                        )[::-1]     # cue 3, prestimulus,   all mice,    slope/intercept
    orders = np.argsort(decauc_pre)[::-1]
    
    
    blocklabels = ['initial\nneural decoder','transition\nneural decoder']
    timeperiodlabels = [' PRE',' ON']
        
    for blo in [0,1]:          # block cue 1, cue 3
        for ptrx in [0,1]:    # pre, on
            
            G = crossgradients[blo,ptrx,orders,:]
            
            axs = ax[2,blo+2*ptrx]

            G = crossgradients2[blo,ptrx,orders,:]
            
            im = axs.imshow(G,vmin=0.0,vmax=1.0,cmap='RdBu',aspect='auto')
            figs.plottoaxis_notickmarks(axs)


            axs.set_title(blocklabels[blo]+timeperiodlabels[ptrx])
            

    axs = ax[2,0]
    figs.labelaxis(axs,panel)
    axs.set_ylabel('mice, ordered by context PRE')
    axs.set_xlabel('unimodal cueing session trials')














# neural all mice, beta plot
    crossgradients = 0.5*np.ones( (2,2,n_mice,n_max_trials) )      # starting values for empty trials, this is raw trial numbers
#    crossgradients2 = 0.5*np.ones( (2,2,n_mice,trialcompress) )      # starting values for empty trials, this is timecompressed
    blockorder = np.zeros((n_mice,2),dtype='int16')      # this will mean 0: attend visual first, 1: attend audio first
    max_trials_ofallmice = [0,0]
    export_crossgradient = []
    for n,dn in enumerate(datanames):
    
        export_crossgradient.append(pickle.load(open('../cache/continuous/offstimgradient,cue,decodeprobs,savgol_%s-%s.pck'%(dn,continuous_method),'rb')))
#        export_crossgradients[n][bl][ptrx,:,0]
        
        blv,bla = preprocess.getorderattended(dn)
        
        for bl in [0,1]:          # block visual, audio
            n_cuetrials = export_crossgradient[n][bl].shape[1]
            blockorder[n,bl] = int(([blv[0],bla[0]][bl]-1)/2)
            max_trials_ofallmice[blockorder[n,bl]] = max(n_cuetrials, max_trials_ofallmice[blockorder[n,bl]])
            for ptrx in [0,1]:    # pre, on
                crossgradients[blockorder[n,bl],ptrx,n,:n_cuetrials] = export_crossgradient[n][bl][ptrx,:n_max_trials,2]


        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        decsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # lower confidence bound
        cx_decauc_pre[n] = neph.get_auctimenormalized(decsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]        



        _,_,_,_,cuetraining,cue_est = preprocess.loadtrainingbehaviouraldata(dn)

        cuetrainings.append(cuetraining)
        # [mice][init/trans,sessions,h+c/m+f,trials]           =(14)(2,8,2,30)

        cue_ests.append(cue_est)       # ix 6,7 (from 0) will contain the beta parameters


    panel = 'C'
    ptrx = 0  # only PRE
    
    probability_resolution = 100

    for bl in [0,1]:
        
        beta_a = np.ones(max_trials_ofallmice[bl])
        beta_b = np.ones(max_trials_ofallmice[bl])
        x = np.linspace(0,1,probability_resolution+1)
        triallist = np.arange(max_trials_ofallmice[bl])+1
        sea = np.meshgrid( x, np.linspace(1,max_trials_ofallmice[bl]) )
        pr = np.zeros( (probability_resolution+1,max_trials_ofallmice[bl]) )
        
        for n in range(n_mice):
            p_m = export_crossgradient[n][blockorder[n,bl]][ptrx,:,0]
            p_s = export_crossgradient[n][blockorder[n,bl]][ptrx,:,1]
            for trx in range(len(p_m)):
                beta_a[trx] += p_m[trx]*1
                beta_b[trx] += (1-p_m[trx])*1
        for trx in range(len(beta_a)):
            pr[:,trx] = sp.stats.beta.pdf( x, beta_a[trx], beta_b[trx] )
        pr = pr / np.sum(pr,axis=0)
        cumprobmap = np.cumsum(pr,axis=0)               #          * (x[1]-x[0])
        
        
        # get contours for the mle and highest density credible intervals 67 and 95%
        c_l = np.zeros((2,max_trials_ofallmice[bl]))
        c_u = np.zeros((2,max_trials_ofallmice[bl]))
        mle = np.zeros(max_trials_ofallmice[bl])
        for t in range(max_trials_ofallmice[bl]):
            mle[t] = x[np.argmax(pr[:,t])]
            c_l[0,t], c_u[0,t] = neba.hdci(ppf=lambda q: neba.empiricalppf(q,cumprobmap[:,t],x), alpha=19/20,resolution=resolution)  # 95% hdci
            c_l[1,t], c_u[1,t] = neba.hdci(ppf=lambda q: neba.empiricalppf(q,cumprobmap[:,t],x), alpha=2/3,resolution=resolution)  # 67% hdci




        axs = ax[1,bl]
        handle = axs.contourf(triallist,x,pr,levels=100)
#        axs.plot(pr)       # debug
#        fig.colorbar(handle,ax=axs)

        axs.plot(triallist,mle,lw=2,color='white',label='MLE/MAP')
        axs.plot(triallist,c_l[1,:],'--',color='white',label='67% h.dens.cred.int.')
        axs.plot(triallist,c_u[1,:],'--',color='white')
        axs.plot(triallist,c_l[0,:],'--',color='grey',label='95% h.dens.cred.int.')
        axs.plot(triallist,c_u[0,:],'--',color='grey')
        if bl==0:
            leg = axs.legend(frameon=False)
            for text in leg.get_texts(): text.set_color('white')
    

#        axs.set_xticks( [ [1,30,max_trials_ofallmice[bl]], [1,30,60,max_trials_ofallmice[bl]] ][bl] )
        axs.set_xticks( [1,15,30] )
        axs.set_xlim([1,30])
        axs.set_yticks([0,0.5,1])
#        axs.set_xlim([0.1,max_trials_ofallmice[bl]+0.9])

        axs.set_title(['initial\nneural decoder PRE','transition\nneural decoder PRE'][bl])
        if bl==0:
            axs.set_ylabel('correct context probability')
            figs.labelaxis(axs,panel)





















# These are pairwise statistics for each mouse

#
#    panel = 'E'
#        
#    X_list = [[],[]]       # as below: context pre values
#    Y_list = [[],[]]       # these will store the per session first and last data for regression
#    Y_m_list = []
#    for bl in [0,1]:          # initial and transition block
#        period = [initials,transitions][bl]         # shortcut call
#        
#        axs = ax[2,bl]
#        
#        Y = []
#        Ys = []
#        for i in [0,1]:          # first and last, on the same figure
#            
#
#            
#            # behaviour for all mice
#            X_list[i].append(   np.concatenate( [ [  contextpre[n] for d in range( len(period[n][i]) ) ] \
#                               for n in range(n_mice)  ] )  )
#            Y_list[i].append( np.concatenate(  [ [  period[n][i][d] for d in range( len(period[n][i]) ) ] \
#                               for n in range(n_mice)  ]  ) )
#                    # Y_list will be Y_list[block][first/last][sessions*window]
#
#            
#            Y.append(  [ (period[n][i]).mean()  \
#                               for n in range(n_mice)  ]  )
#            Ys.append(  [ (period[n][i]).std()/np.sqrt( len(period[n][i]) )  \
#                               for n in range(n_mice)  ]  )
#        Y_m_list.append(Y)         # this will be [bl][i]
#            
#            
#        m = np.array(Y)
#        s = np.array(Ys)
#        print(m.shape,s.shape)
#        d = m[1]-m[0]
#        ci = 2*(s[1]+s[0])
#        crit = (d>ci).astype('int16')
##        axs.plot(i*np.ones(n_mice),Y,'o-')
#        for n,dn in enumerate(datanames):
#            if crit[n]:
#                axs.errorbar(x=[0,1],y=m[:,n],yerr=s[:,n],color='darkgreen',capsize=3)
#            else:
#                axs.errorbar(x=[0,1],y=m[:,n],yerr=s[:,n],color='darkred',capsize=3)
#            axs.text(-0.1,m[0,n],dn, fontsize='x-small',horizontalalignment='right',verticalalignment='center')
#        axs.set_yticks([0.5,1])
#        axs.set_ylim(0.44,1.01)
#        axs.set_xticks([0,1])
#        axs.set_xticklabels(['first %d'%window,'last %d'%window])
#        axs.set_xlim(-0.5,1.5)
#        figs.plottoaxis_chancelevel(axs,0.5)
#        if bl==0:
#            axs.set_ylabel('behaviour\nfraction correct')
#        
#            figs.labelaxis(axs,panel)
#
#
#
#
#
#    if 0:
#    # visualize HDCI change between first window and last window difference
#        panel = 'Z+'
#        initials = []
#        transitions = []
#        for n,dn in enumerate(datanames):
#            initials.append( [ np.reshape(cuetrainings[n][0][:,0,:window],-1),\
#                               np.reshape(cuetrainings[n][0][:,0,-window:],-1) ]   )
#            transitions.append( [ np.reshape(cuetrainings[n][1][:,0,:window],-1),\
#                                  np.reshape(cuetrainings[n][1][:,0,-window:],-1) ] )
#    
#    
#    
#        colors = ['darkturquoise','goldenrod']
#        labels = ['first %d'%window,'last %d'%window]
#        
#        for bl in [0,1]:          # initial and transition block
#            period = [initials,transitions][bl]         # shortcut call
#            
#            axs = ax[3,bl+2]
#            
#            Y = []
#            Ys = []
#            for i in [0,1]:          # first and last, on the same figure
#                
#                axs.plot(x,probmapwindowaverage[bl,i,:],lw=2,color=colors[i],label=labels[i])
#                axs.fill_between(x,probmapwindowaverage[bl,i,:]-probmapwindowaverage_s[bl,i,:],\
#                                   probmapwindowaverage[bl,i,:]+probmapwindowaverage_s[bl,i,:],\
#                                   color=colors[i],alpha=0.2)
#    
#                
#                # behaviour for all mice
#            
#            axs.set_ylim(0,0.2)
#            axs.set_yticks([0,0.1,0.2])
#            axs.set_xticks([0,0.5,1])
#            
#            if bl==0: axs.legend(frameon=False)
#    #        if bl==0: figs.labelaxis(axs,panel)
#        
#
#
#







    fig.suptitle('Supplementary Figure 4')






    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Supp4_context,cueblockgradients'+ext)














#       FIG   Supp5     Confounders

def supplementaryfigure6_unused():

    
    datanamesefg = ['ME108','ME110','ME112','ME113','DT009','DT014','DT017','DT018','DT019','DT021','DT030','DT031','DT032']
    variablecolors = ['navy','mediumvioletred']    

#   N                  # explained variance from predictive GLM
    modeldisplayP = np.array([17,0+4,2+4, 18,8+0+4,8+2+4],dtype='int16')              # find the individuals subtracted models, and compare to the full model;    with runspeed: 18, 0+4+8, 2+4+8
    modeldisplayQ = np.array([4,5,6,7,12,13,14,15,16,17],dtype='int16')
    contextcomp_idx = [ [2,8], [6,9] ]       # first two is without run, last two is with run
    
#    M = np.zeros(n_mice)
    n_mice = len(datanamesefg)
    cxpre = np.zeros(n_mice)

    for n,dn in enumerate(datanamesefg):
        mc,r,a,d = pickle.load(open('../cache/glm/glm_predictive,modelselection_%s'%(dn),'rb'))
        if n==0:
            explainedvariance = np.zeros((n_mice,len(modeldisplayP)))
            aic = np.zeros((n_mice,len(modeldisplayP)))
            R = np.zeros((n_mice,len(modeldisplayQ)))
            A = np.zeros((n_mice,len(modeldisplayQ)))

        explainedvariance[n,:] = r[modeldisplayP]
        aic[n,:] = a[modeldisplayP]

        R[n,:] = r[modeldisplayQ]
        A[n,:] = a[modeldisplayQ]

        acrossdecoder = pickle.load(open('../cache/continuous/responsedecodes,angles-%s_%s-%s-%s-%s,%s.pck'%('allexpcond','all',dn,continuous_method,'context','all'),'rb'))
        cxpre[n] = neph.get_auctimenormalized(acrossdecoder[1][:,0],T['starttime'],T['stimstarttime'])[0]
        # cxpre holds context pre variable for each mice










    # load decoders for task variables from runspeed
    taskaspects = ['visual','audio','context','choice']
    aucdeltat = 750*pq.ms     # for the order length
    aucdeltaix = int((aucdeltat / T['dt']).magnitude)
    features = []         # this will be (n_mice,taskaspect,trajectoryposition) = (14,4,4)
    for n,dn in enumerate(datanamesefg):
        blv,bla = preprocess.getorderattended(dn)
        featurevector = []
        for cx,comparison in enumerate(taskaspects):
            acrossdecoder = pickle.load(open('../cache/phys/rundecodes-%s-%s-%s.pck'%(dn,continuous_method,comparison),'rb'))

            rundecsignal = acrossdecoder[1][:,0]     #-acrossdecoder[1][:,2]       # change to 2 for lower confidence bound or mean to the left
            rundecsignal_var = acrossdecoder[1][:,1].mean()**2       # this is the decoder cross validation variance
            
            rundecauc_pre = neph.get_auctimenormalized(rundecsignal,T['stimstarttime']-aucdeltat,T['stimstarttime'])[0]
            rundecauc_early = neph.get_auctimenormalized(rundecsignal,T['stimstarttime'],T['stimendtime']/2+T['stimstarttime']/2)[0]
            rundecauc_late = neph.get_auctimenormalized(rundecsignal,T['stimendtime']/2+T['stimstarttime']/2,T['stimendtime'])[0]
            rundecauc_post = neph.get_auctimenormalized(rundecsignal,T['stimendtime'],T['stimendtime']+aucdeltat)[0]

            rundecauc_pre_ci = np.sqrt(rundecsignal[T['stimstart_idx']-aucdeltaix:T['stimstart_idx']].var()+rundecsignal_var)/np.sqrt(aucdeltaix)
            rundecauc_early_ci = np.sqrt(rundecsignal[T['stimstart_idx']:int(T['stimstart_idx']/2+T['stimstart_idx']/2)].var()+rundecsignal_var)/np.sqrt(aucdeltaix)
            rundecauc_late_ci = np.sqrt(rundecsignal[int(T['stimstart_idx']/2+T['stimend_idx']/2):T['stimstart_idx']].var()+rundecsignal_var)/np.sqrt(aucdeltaix)
            rundecauc_post_ci = np.sqrt(rundecsignal[T['stimend_idx']:T['stimend_idx']+aucdeltaix].var()+rundecsignal_var)/np.sqrt(aucdeltaix)

            featurevector.append( [  rundecauc_pre, rundecauc_early, rundecauc_late, rundecauc_post,\
                                     rundecauc_pre_ci, rundecauc_early_ci, rundecauc_late_ci, rundecauc_post_ci ])

        features.append(np.array(featurevector))

#    print(features)

    features = np.array(features)    # (mice,taskaspect,preearlylatepost+same w/ci)




    fig,ax = plt.subplots(1,3,figsize=(3*8,1*8))










    # predictive
    panel = 'A'
#    axs = ax[3,2]
#    axs.set_ylim(0,0.2)
#    axs.set_xticks(np.arange(n_mice)); axs.set_xticklabels(datanamesefg,rotation=45)
#    axs.set_yticks([0,0.1,0.2])
#    axs.set_ylabel('explained variance by visual predictor')

    ordercxpre = np.argsort(cxpre)#[::-1]

    axs = ax[0]
#    axs.bar(np.arange(n_mice)-0.2,explainedvariance[ordercxpre,0],width=0.4,color=variablecolors[0],label='visual')
#    axs.bar(np.arange(n_mice)+0.2,explainedvariance[ordercxpre,1],width=0.4,color=variablecolors[1],label='context')
    
    axs.boxplot(x=explainedvariance, positions=[0,1,2,4,5,6], notch=True)
    
    
#    axs.boxplot(explainedvariance,labels=['visual','context'], notch=True)
#    xlabels = [ datanamesefg[n] for n in ordercxpre ]
#    axs.set_ylim(0,0.5)
    axs.set_xlim(-1,7)
    axs.set_xticks([0,1,2,4,5,6])
    axs.set_xticklabels(['full','- visual','- context','full w/run','- visual','- context'],rotation=30)
#    axs.text(1,-0.02,'without locomotion',horizontalalignment='center')
#    axs.text(5,-0.02,'with locomotion',horizontalalignment='center')
    
#    axs.set_xticklabels([])
    axs.set_yticks([0,0.2,0.4])
    axs.set_ylabel('explained variance')
    figs.plottoaxis_chancelevel(axs,0)

    axs.legend()

    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)










    if 0:
        # Predictive
        panel = 'B'
        # loglikelihood
        x = cxpre
        y = np.c_[ A[ :, contextcomp_idx[0][1] ] - A[ :, contextcomp_idx[0][0] ],\
                   A[ :, contextcomp_idx[1][1] ] - A[ :, contextcomp_idx[1][0] ] ]
    
    
    #    labels = ['LL cx+ - cx-   r- ','LL cx+ - cx-   r+ '] 
        labels = ['no locomotion','with locomotion']
        colors = ['darkred','orange']
    
        axs = ax[1]
    
        for m in range(2):
            for n,dn in enumerate(datanamesefg):
                
                axs.scatter(x[n],y[n,m],s=150,marker='o',color=colors[m],alpha=0.8)
    #            axs.text(x[n]*1.02,y[n,m]*0.995,datanames[n],color='k',alpha=0.3,fontsize=10)
    
            l = sp.stats.linregress(x,y[:,m])
            if l[3]<0.13:
                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                axs.plot(line,l[0]*line+l[1],color=colors[m],linewidth=2,label=labels[m])
            else:
                line = np.linspace(start=x.min()*0.8,stop=x.max()*1.2,num=2)
                axs.plot(line,l[0]*line+l[1],'--',color=colors[m],linewidth=2,label=labels[m])
                
            xoffs = 0.65
    #            if sx in [3,4,5,7]: xoffs = 25
            if l[3]<0.001:
                axs.text(line[1]-xoffs,11-2*m,'p<%5.3f, $R^2$=%4.2f'%((0.001,l[2]**2)),color=colors[m])
            else:
                axs.text(line[1]-xoffs,11-2*m,'p=%5.3f, $R^2$=%4.2f'%((l[3],l[2]**2)),color=colors[m])
        
        
        axs.set_xlim(0.45,1.04)
        axs.set_ylim(-22,15)
        axs.set_xticks([0.5,1])
        axs.set_yticks([-20,-10,0,10])
        
        figs.plottoaxis_crosshair(axs,0.5,0)
        axs.legend()    
    
        axs.set_xlabel('context PRE', labelpad=0.1, verticalalignment='bottom')
        axs.set_ylabel('AIC difference between\nwith and without context')
    
        figs.labelaxis(axs,panel)
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)













    # runspeed context decoders are near chance
    panel = 'B'
    axs = ax[1]
    axs.boxplot(features[:,2,np.array([0,2],dtype='int16')], notch=True)

    axs.set_yticks([0.5,1.0])
    axs.set_ylim(0.45,1.01)
    figs.plottoaxis_chancelevel(axs,0.5)
    axs.set_xticklabels(['PRE','ON'])
    
    axs.set_ylabel('locomotion decoder acc.', labelpad=0.1, verticalalignment='bottom')
    
    figs.labelaxis(axs,panel)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    
    







    # show runspeed and context pre commonalities
    panel = 'C'


        
    alpha = 1-1*0.15
    
    axs = ax[2]
    for k,timing in enumerate([0,2]):
        color = ['fuchsia','purple'][k]
    # run cx pre,       and run cx post
        X = cxpre                      # context pre _neural_
        Y = features[:,2,timing]            # context pre [0] and post [4] _runnning speed_
    
        for n,mousename in enumerate(datanamesefg):
            # features dimensions: (mice,taskaspect,preearlylatepost)
            
            # plot context pre neural vs. context pre running speed
            axs.scatter(X[n],Y[n],s=150,marker='o',color=color,alpha=alpha)
    #                        axs.text(features.iloc[n,xx]*1.05,features.iloc[n,yx]*0.98,mousename,color=colors[ynx],alpha=0.1)
    #        if ynx==0 or sx==8: axs.text(features.iloc[n,xx]*1.02,features.iloc[n,yx]*0.995,mousename,color='k',alpha=0.3,fontsize=10)
    #                        if ynx==0: axs.text(features.iloc[n,xx],tribe_stds[:,0,:,pox].min()*0.95,mousename,color='k',alpha=0.3,rotation=90,fontsize=10)
            
            axs.errorbar(X[n],Y[n],yerr=features[n,2,timing+4],color=color,alpha=alpha)
            
    
        l = sp.stats.linregress(X,Y)
        if l[3]<0.13:
            line = np.linspace(start=X.min()*0.8,stop=X.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],color=color,alpha=alpha,linewidth=2,label='context %s'%['PRE','ON'][k])
        else:
            line = np.linspace(start=X.min()*0.8,stop=X.max()*1.2,num=2)
            axs.plot(line,l[0]*line+l[1],'--',color=color,alpha=alpha,linewidth=2,label='context %s'%['PRE','ON'][k])
            
        xoffs = 0.4; yoffs = 0
        xoffs = 0.6; yoffs = 0.3
        if l[3]<0.001:
            axs.text(line[1]-xoffs,l[0]*line[1]+l[1]*1.02+yoffs,'p<%5.3f, $R^2$=%4.2f'%((0.001,l[2]**2)),color=color)
        else:
            axs.text(line[1]-xoffs,l[0]*line[1]+l[1]*1.02+yoffs,'p=%5.3f, $R^2$=%4.2f'%((l[3],l[2]**2)),color=color)
        
        axs.set_xlim(0.45,1.01); axs.set_ylim(0.45,1.01)
        axs.set_xticks([0.5,1.0]); axs.set_yticks([0.5,1.0])
        figs.plottoaxis_crosshair(axs,0.5,0.5)
        axs.set_xlabel('neural decoder acc.\ncontext PRE',labelpad=0.1, verticalalignment='center')
        axs.set_ylabel('locomotion decoder acc.', labelpad=0.1, verticalalignment='bottom')
#        axs.legend(['locomotion context %s'%['PRE','POST'][k]])
#        axs.set_title(trajectorygrouplabels[ynx])

        if k==1:
            axs.legend(frameon=False)
            figs.labelaxis(axs,panel)
            axs.spines['right'].set_visible(False)
            axs.spines['top'].set_visible(False)


















#    fig.suptitle('Figure 7',verticalalignment='bottom')


    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Fig7_predictglm,locomotion'+ext)



















def supplementaryfigure7_unused():









    panel = 'E'
    n_cherry = datanames.index('DT019')
    cx = 0
    axs = ax[3,0]
    
    axs.hist(activities_all[n_cherry][cx][3],bins=np.arange(-2.6666666,0.1,0.3333333),density=True,\
             edgecolor=activitycolors[3],alpha=0.25,lw=2,facecolor='white',orientation='horizontal')
    axs.hist(activities_all[n_cherry][cx][2],bins=np.arange(-2.6666666,0.1,0.3333333),density=True,\
             edgecolor='darkred',facecolor=activitycolors[2],alpha=0.25,orientation='horizontal')

    for k in [1,0]:
        data = activities_all[n_cherry][cx][k+2]
        kde = sp.stats.gaussian_kde(data)
        x = np.arange(-3.3,1.1,0.02)
        axs.plot(kde(x),x,color=activitycolors[k+2],lw=3,alpha=0.9,label=['correct rejection','false alarm'][k])

    axs.legend(frameon=False,loc='upper right')
    axs.set_ylim(-3.3,1.1)

    axs.set_ylabel('projected activity\nto visual DV')
    axs.set_xticklabels([])
    axs.set_xlabel('hist. freq. [AU]')

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    figs.labelaxis(axs,panel,x=-0.3)





            
        
    panel = 'F'
    
    cx = 0    # choose the context conditioned block (0), or context (1) for displaying projected activities
    lx = 1     # labels for the conditions  0 visual        1 context
    dprimes = np.zeros((len(datanames),2))   #    (n_mice x {45,135})
    colors = ['black','darkgreen','darkred']
    axs = ax[3,1]
    
    for n,dn in enumerate(datanames):
        print(dn,'activity len',[len(activities_all[n][cx][i]) for i in [0,1,2,3]])

        # print(dn,'activity vector@ hmcf',activities_all[n][cx])

        # find significantly different mice with behaviour
        s45,p45 = sp.stats.ttest_ind(activities_all[n][cx][0], activities_all[n][cx][1], equal_var=False)
        s135,p135 = sp.stats.ttest_ind(activities_all[n][cx][2], activities_all[n][cx][3], equal_var=False)
        

        dprimes[n,:] = [ nedi.dprime_normal(activities_all[n][cx][0], activities_all[n][cx][1]),\
                         nedi.dprime_normal(activities_all[n][cx][2], activities_all[n][cx][3]) ]
    
        m = np.array([ np.mean(activities_all[n][cx][2]), np.mean(activities_all[n][cx][3]) ])
        e = np.array([ np.std(activities_all[n][cx][2])/np.sqrt(len(activities_all[n][cx][2]))*2, \
              np.std(activities_all[n][cx][3])/np.sqrt(len(activities_all[n][cx][3]))*2 ])
        

        if n==n_cherry: lwm = 1.5
        else: lwm = 1

        if p135>0.05: color = colors[0]; alpha = 0.3; lw = 1
        else:
            if m[0]<m[1]: color = colors[1]; alpha = 0.7; lw = 3*lwm
            else: color = colors[2]; alpha = 0.5; lw = 2*lwm
    
        # for k in range(len(activities_all[n])):
            # axs.hist(activities_all[n][cx][k],bins=10,color=colors[k],label=taskcontext+labels[lx][2+k],alpha=0.5)
        axs.errorbar([0,1],m,e,lw=lw,color=color,alpha=alpha)
        
           
    axs.set_ylim(-3.3,1.1)

    # axs.set_ylabel('projection to visual DV')
    axs.set_xticks([0,1])
    axs.set_xticklabels(['correct rejection','false alarm'])
        

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)


    figs.labelaxis(axs,panel)




    fig.suptitle('Supplementary Figure 6')






    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'Supp6_choice'+ext)



















def unusedcollections():
    # fig 2 G:  multiple mice activity projection histograms

    panel = 'G'
    axs = ax[2,2]
    varidx = 0   # visual
    varidy = 1   # y coords of visual
    P_cx = [[],[]]
    for k in range(2):
        for n in range(len(projections_all)):
            P_cx[k].extend(projections_all[n][varidx,preon,k%2] * np.sign(dbnbasiscomps_all[n][varidx,0]) )
            P_cx[k].extend(projections_all[n][varidx,preon,k%2+2] * np.sign(dbnbasiscomps_all[n][varidx,0]) )
    for k in range(2):
        axs.hist(P_cx[k],bins=50,color=['darkgreen','red'][k],alpha=0.7)
    axs.legend(['  45°','135°'],loc=2)
    axs.plot( [0,0+3/2* np.sqrt( (dbnbasiscomps[0,varidx])**2 + (dbnbasiscomps[1,varidx])**2 ) ], [70,70],\
              color=basiscolors[0],linewidth=4,alpha=0.9 )
    axs.text(-0.05,72,'decision normal',color=basiscolors[0])

    axs.set_xlabel('visual decision boundary normal')
    figs.labelaxis(axs,panel)
    axs.set_title('n=%d'%13)


    axs = ax[2,3]
    varidx = 1   # context
    P_cx = [[],[]]
    for k in range(2):
        for n in range(len(projections_all)):
            P_cx[k].extend(projections_all[n][varidx,preon,2*k] * np.sign(dbnbasiscomps_all[n][varidx,1]) )
            P_cx[k].extend(projections_all[n][varidx,preon,2*k+1] * np.sign(dbnbasiscomps_all[n][varidx,1]) )
    for k in range(2):
        axs.hist(P_cx[k],bins=50,color=['silver','dimgrey'][k],alpha=0.7)
    axs.legend(['attend visual','ignore visual'],loc=2)
    axs.plot( [0,0+3/2* np.sqrt( (dbnbasiscomps[0,varidx])**2 + (dbnbasiscomps[1,varidx])**2 ) ], [55,55],\
              color=basiscolors[3],linewidth=4,alpha=0.9 )
    axs.text(-0.05,56.6,'decision normal',color=basiscolors[3])
    
    axs.set_xlabel('context decision boundary normal')
    axs.set_title('n=%d'%13)
    
    
    


    
    # Figure 5 D, multiple mouse histograms    

    panel = 'D'
    preon = 1
    
    axs = ax[1,0]
    varidx = 0   # visual
    varidy = 1   # y coords of visual
    vi = 0
    ch = 2
    P_vi = [[],[]]
    for k in range(2):
        for n in range(len(projections_all)):
            P_vi[k].extend(projections_all[n][varidx,preon,k%2] * np.sign(dbnbasiscomps_all[n][varidx,0]) )
            P_vi[k].extend(projections_all[n][varidx,preon,k%2+2] * np.sign(dbnbasiscomps_all[n][varidx,0]) )
    for k in range(2):
        axs.hist(P_vi[k],bins=50,color=['darkgreen','red'][k],alpha=0.7)
    axs.legend(['  45°','135°'],loc=2)
    axs.plot( [0,0+3/2* np.sqrt( (dbnbasiscomps[0,varidx])**2 + (dbnbasiscomps[1,varidx])**2 ) ], [60,60],\
              color=basiscolors[0],linewidth=4,alpha=0.9 )
    axs.set_ylim(0,65)
    axs.text(-0.05,65,'decision normal',color=basiscolors[0])

    axs.set_xlabel('visual decision boundary normal')
    axs.text(-3.5,65,panel,fontsize=24,fontweight='bold')
    axs.set_title('n=%d'%13)



    axs = ax[1,1]
    varidx = 2   # choice
    P_ch = [[],[]]
    for n in range(len(projections_all)):
        ixgohit,ixnogocorrrej,ixgomiss,ixnogofal = perfidx_all[n]
        P_ch[0].extend(projections_all[n][varidx,preon,0][ixgohit] * np.sign(dbnbasiscomps_all[n][varidx,0]))
        P_ch[1].extend(projections_all[n][varidx,preon,0][ixgomiss] * np.sign(dbnbasiscomps_all[n][varidx,0]))
        P_ch[1].extend(projections_all[n][varidx,preon,1][ixnogocorrrej] * np.sign(dbnbasiscomps_all[n][varidx,0]))
        P_ch[0].extend(projections_all[n][varidx,preon,1][ixnogofal] * np.sign(dbnbasiscomps_all[n][varidx,0]))
    for k in range(2):
        axs.hist(P_ch[k],bins=50,color=['blue','brown'][k],alpha=0.7)
    axs.legend(['lick','no lick'],loc=2)
    axs.plot( [0,0+3/2* np.sqrt( (dbnbasiscomps[0,varidx])**2 + (dbnbasiscomps[1,varidx])**2 ) ], [60,60],\
              color=basiscolors[4],linewidth=4,alpha=0.9 )
    axs.set_ylim(0,65)
    axs.text(-0.05,56.6,'decision normal',color=basiscolors[4])
    
    axs.set_xlabel('context decision boundary normal')
    axs.set_title('n=%d'%13)









def figrn():

    dn = 'DT017'
#    selected_cell_ids = [2,10]
    selected_cell_ids = [7,26]
    
    block = preprocess.loaddatamouse(dn,T,continuous_method=continuous_method,normalize=False)        # use raw firing rates: normalzie=False
    n_neuron = block.segments[0].analogsignals[0].shape[1]
    times = block.segments[0].analogsignals[0].times
    
    ti = np.arange(0, 36*6,6)*pq.s

    V = [ block.segments[42+trx].annotations['visual'] for trx in range(30)]

    Sp = [ np.concatenate([block.segments[42+trx].spiketrains[n][0].magnitude*pq.ms + ti[trx] for trx in range(30)]) for n in range(len(block.segments[42].spiketrains)) ]

    fig,axs = plt.subplots(1,1,figsize=(54,8))
    axs.eventplot([  Sp[i] for i in selected_cell_ids ], color = 'navy',linelengths=0.618)
    
    axs.set_xlim(-1000,176000)
    axs.set_ylim(axs.get_ylim())
    
    for k in range(30):
        figs.plottoaxis_stimulusoverlay(axs,{'stimstarttime':0*pq.ms+ti[k], 'stimendtime':3000*pq.ms+ti[k]} )
        axs.text(ti[k]*1000+100*pq.ms,-0.5 + (V[k]==45)*2,'%d°'%V[k],verticalalignment='center')

    axs.set_yticks([])
    axs.set_xticks(np.arange(0,150001,30000))
    axs.set_xticklabels(['0 s','30 s','60 s','90 s','120 s', '150 s'])

    axs.spines['left'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)

    save = 0 or globalsave
    if save:
        fig.savefig(resultpath+'raster-%s-c%d,%d'%(dn,selected_cell_ids[0],selected_cell_ids[1])+ext)










def main():
    # drawschematics()

    figure1()
    # figure2()     # raster, firing rates, PCA
    # figure3()     # visual and context orthogonal
    # figure4()     # context dynamics
    # figure5()     # visual discrimination independent of context
    # figure6()     # choice
    # figure7()     # locomotion-invariance


    # statshelper()


    # supplementaryfigure1()    # drift control
    # supplementaryfigure2()    # show that without interneurons, context is still decodable in excitatory only
    # supplementaryfigure3()    # show that decoding from visual DV, information is limited.
    # supplementaryfigure4()    # control for number of neurons
    # supplementaryfigure5()    # visual discrimination is dependent on animal performance





    # test schematics
    # fig,axs = plt.subplots(1,1,figsize=(8,8))
    # axs = fig.add_subplot(1,1,1, projection='3d') 
    # drawsubspaces(axs,0)



    
    # supplementaryfigure4()
    # supplementaryfigure5()
    # supplementaryfigure6()
    # supplementaryfigure7()     # confounds, GLMs

    # figrn()

    # unusedcollections()
    
    
    
    return


if __name__ == '__main__':
    main()
    plt.show()