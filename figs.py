# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 20:08:30 2019

@author: mahajnal
"""


import numpy as np
import scipy as sp


import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


import quantities as pq
import neo

import neurodiscover as nedi
import neurophysiology as neph



plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'legend.fontsize': 14})
plt.rcParams.update({'lines.linewidth': 1})





colors1=['dodgerblue','darkorange','darkred']
colors2=['dodgerblue','chocolate','darkorange']







def getcorrcolormap(which='correlation'):

    if which=='firingrate':
        cd    = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.0),
                           (0.666, 0.5, 0.5),
                           (1.0, 0.9, 0.9)),
            
                 'blue':  ((0.0, 0.9, 0.9),
                           (0.333, 0.5, 0.5),
                           (0.5, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),
        
                 'green':  (( 0.0, 0.6, 0.6),
                           (0.25, 0.0, 0.0),
                           (0.75, 0.0, 0.0),
                           ( 1.0, 0.6, 0.6))
                 }        
    elif which=='variance':
        cd    = {'red':   ((0.0, 0.9, 0.9),
                           (0.333, 0.5, 0.5),
                           (0.5, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),
        
                 'blue':  (( 0.0, 0.6, 0.6),
                           (0.25, 0.0, 0.0),
                           (0.75, 0.0, 0.0),
                           ( 1.0, 0.6, 0.6)),
        
                 'green':  ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.0),
                           (0.666, 0.5, 0.5),
                           (1.0, 0.9, 0.9))
                 }
    if which=='correlation':
        cd    = {'red':   ((0.0, 0.9, 0.9),
                           (0.333, 0.5, 0.5),
                           (0.5, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),
        
                 'blue':  (( 0.0, 0.6, 0.6),
                           (0.25, 0.0, 0.0),
                           (0.75, 0.0, 0.0),
                           ( 1.0, 0.6, 0.6)),
        
                 'green':  ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.0),
                           (0.666, 0.5, 0.5),
                           (1.0, 0.9, 0.9))
                 }
    return colors.LinearSegmentedColormap('BuBaRe', cd)







# utilities


def labelaxis(ax,panel='A',x=-0.1,y=1.05,fontsize=24,D2=False):
#    return     # use this to avoid printing any panel labels
    if not D2:
        ax.text(x,y,panel,fontsize=fontsize,fontweight='bold',transform=ax.transAxes)
    else:
        ax.text2D(x,y,panel,fontsize=fontsize,fontweight='bold',transform=ax.transAxes)
    
    
    
def plottoaxis_notickmarks(ax):
    ax.set_xticks([])
    ax.set_xticklabels('')
    ax.set_yticks([])
    ax.set_yticklabels('')
    

def setxt(axs,ticks=[0,3000],labels=['0 ms','3000 ms']):
    axs.set_xticks(ticks)
    if labels!=None: axs.set_xticklabels(labels)



def plottoaxis_stimulusoverlay(ax,T=None,dmts=False,offphase=None):
    if (np.array(offphase)==None).any():
        start = T['stimstarttime']
        stop = T['stimendtime']
    else:
        start = offphase[0]
        stop = offphase[1]
    ylim=ax.get_ylim()
    ax.fill_between([ start,stop ],[ylim[0],ylim[0]],[ylim[1],ylim[1]],color='gray',alpha=0.33)
    if dmts==True:
        ax.fill_between([T['stim2starttime'],T['stim2endtime'] ],[ylim[0],ylim[0]],[ylim[1],ylim[1]],color='gray',alpha=0.3)
    elif dmts>1:
        ax.fill_between([T['stim2starttime'],T['stim2starttime']+dmts ],[ylim[0],ylim[0]],[ylim[1],ylim[1]],color='gray',alpha=0.2)
        ax.fill_between([T['stim2starttime']+dmts,T['stim2starttime']+dmts+500*pq.ms ],[ylim[0],ylim[0]],[ylim[1],ylim[1]],color='gray',alpha=0.3)


def plottoaxis_chancelevel(ax,ch=0.0,label=None,lw=2,xs=None):
    if xs==None:
        xlim=ax.get_xlim()
        ax.plot([xlim[0],xlim[1]],[ch,ch],'k--',linewidth=lw,alpha=0.2,label=label)
    else:
        ax.plot([xs[0],xs[1]],[ch,ch],'k--',linewidth=lw,alpha=0.2,label=label)
        

def plottoaxis_crosshair(ax,x=0.0,y=0.0,color='black'):
    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    ax.plot([xlim[0],xlim[1]],[y,y],'k--',color=color,linewidth=2,alpha=0.2)
    ax.plot([x,x],[ylim[0],ylim[1]],'k--',color=color,linewidth=2,alpha=0.2)




def invisibleaxes(axs,which=['right','top']):
    if 'right' in which: axs.spines['right'].set_visible(False)
    if 'top' in which: axs.spines['top'].set_visible(False)
    if 'left' in which: axs.spines['left'].set_visible(False); axs.set_yticks([])
    if 'bottom' in which: axs.spines['bottom'].set_visible(False); axs.set_xticks([])









# activity plots


def rasterplot(block,trial_ids=[13,16]):
    
    for trial in block.segments[  trial_ids[0]:(trial_ids[1]+1)  ]:
        sp = [ st[0].rescale('s').magnitude for st in trial.spiketrains ]
        
        t = trial.analogsignals[0].times.rescale('s').magnitude
    #            r = np.asarray([fr.ravel() for fr in trial.analogsignals]).T
        r = trial.analogsignals[0]
        plt.figure(figsize=(32,8))
        plt.subplot(1,2,1)
    #            for n in sp: print(len(sp),len(n),n);plt.eventplot(n,color='r')
        plt.eventplot(sp,color='r')
        plt.title(trial.description +'  '+ trial.name )
        plt.subplot(1,2,2)
        plt.plot(t,r)
        plt.title(block.description+' block %d'%trial.annotations['block'])



def plottoaxis_rate(ax,t,x,colorlist='darkred'):
    ax.plot(t,x,color=colorlist,alpha=0.5)


def plottoaxis_plottrajectory(ax,t,m,s=None,colorlist='mediumvioletred',alpha=0.5,linewidth=1,label=None,fill=True):
    ax.plot(t,m,linewidth=linewidth,color=colorlist,label=label,alpha=alpha)
    if fill:
        ax.fill_between(t, m-s*2, m+s*2, color=colorlist,alpha=1/3*alpha**2)#,label=label)


def plottoaxis_plottrajectories(ax,t,ms,ss=None,colormapname='Reds',alpha=0.33,linewidth=1,label=None):
    n_colors = ms.shape[1]
    if colormapname=='Reds':
        colorlist = plt.cm.Reds( np.linspace(1.0, 0.33, n_colors) )
    elif colormapname=='Blues':
        colorlist = plt.cm.Blues( np.linspace(1.0, 0.33, n_colors) )
    for ch in range(ms.shape[1]):
        plottoaxis_plottrajectory(ax,t,ms[:,ch],ss[:,ch],colorlist=colorlist[ch],alpha=alpha,linewidth=linewidth,label=label)




def plottoaxis_timecoursecolorsurface(surface,ax,order=1):
    # order 1: means and unnormalized; order 2=
#    t = surface.times
    
    if order==1:
        cmap = getcorrcolormap('firingrate')
        ax.pcolormesh(surface.T,vmin=-2,vmax=2,cmap=cmap)
    elif order==2:
        cmap = getcorrcolormap('firingrate')
        ax.pcolormesh(surface.T,vmin=-2,vmax=2,cmap=cmap)
    elif order==3:
        cmap = getcorrcolormap('correlation')
        ax.pcolormesh(surface.T,vmin=-0.66666,vmax=0.66666,cmap=cmap)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_yticks([])




# other statistics

def plottoaxis_halftrialvshalftrial(R1,R2,ax,order=1,alpha=1):
    colors = ['orange','sienna','seagreen']
    ax.plot(R1,R2,'o',color=colors[order-1],alpha=alpha)
















# decoders


def plottoaxis_difference(responsedifference,ax,colorlist='darkcyan',topsem=False):
    dm,de = nedi.get_maxnormdifference(responsedifference)
    ax.plot(dm.times,dm,linewidth=3,color=colorlist,label='mean difference')
    ax.fill_between(dm.times, de.squeeze(), dm.squeeze(), alpha=0.5,color=colorlist,label='95% s.e.m.')
    if topsem: ax.fill_between(dm.times, dm.squeeze(), 2*dm.squeeze()-de.squeeze(), alpha=0.5,color=colorlist,label='95% s.e.m.')
    ax.set_xlim([dm.times[0],dm.times[-1]])
    ax.set_ylim([-0.1,None])

#    ax.set_xlabel('time from trial onset [ms]')
#    ax.set_ylabel('z-score difference')






    
def plottoaxis_decoderrocauc(responsedecode,ax,colorlist=colors1,label=None,plottrain=True,plotcrosstest=False,onlysem=False,sem=True,smooth=[],lw=2):
    # this only plots the decoder performances
    t = responsedecode[0].times
    labels = ['train','test','crosstest']
    if plottrain:
        whatplots = [0,1]
    else:
        whatplots = [1]
    if plotcrosstest:
        whatplots.extend([2])

    for sx in smooth:
        s = neph.smooth(responsedecode[sx][:,0],mode='same')
        ax.plot(t,s,linewidth=lw,color=colorlist[sx],alpha=0.5)

    for tt in whatplots:
        if tt in smooth: alphabase=0.05
        else: alphabase=0.2
        ax.plot(t,responsedecode[tt][:,0],linewidth=lw,color=colorlist[tt],alpha=alphabase*5,label=label)
        if sem:
            ax.fill_between(t, (responsedecode[tt][:,0]-responsedecode[tt][:,2]).squeeze(),\
                                (responsedecode[tt][:,0]+responsedecode[tt][:,2]).squeeze(),\
                                alpha=alphabase*2,color=colorlist[tt])
        if not onlysem:   # if we don't need the clutter for full variance, only sem, ignore this:
            ax.fill_between(t, (responsedecode[tt][:,0]-responsedecode[tt][:,1]).squeeze(),\
                            (responsedecode[tt][:,0]+responsedecode[tt][:,1]).squeeze(),\
                            alpha=alphabase,color=colorlist[tt])
    ax.set_xlim([t[0],t[-1]])
    ax.set_ylim([0.45,1.01])




def plottoaxis_decodercoeffs(responsedecode,ax):
    # this only plots the decoder coefficients in time
    n_colors = len(responsedecode)
    colorlist = plt.cm.viridis( np.linspace(0, 0.8, n_colors) )
    
    t = responsedecode[0].times

    for cx,coeff in enumerate(responsedecode):
        ax.plot(t,coeff[:,0],linewidth=2,color=colorlist[cx],label='coeff %d'%(cx+1),alpha=0.6)
        ax.fill_between(t, (coeff[:,0]-coeff[:,2]).squeeze(),\
                            (coeff[:,0]+coeff[:,2]).squeeze(),\
                            alpha=0.1,color=colorlist[cx])
        ax.fill_between(t, (coeff[:,0]-coeff[:,1]).squeeze(),\
                            (coeff[:,0]+coeff[:,1]).squeeze(),\
                            alpha=0.05,color=colorlist[cx])

    ax.set_xlim([t[0],t[-1]])






def plottoaxis_decoderangles(responsedecode,ax):
    # this only plots the angles
    n_colors = len(responsedecode)
    colorlist = plt.cm.viridis( np.linspace(0, 0.8, n_colors) )
    
    t = responsedecode[0].times

    for angx,angles in enumerate(responsedecode):
        ax.plot(t,angles[:,0],linewidth=2,color=colorlist[angx],label='DB to PC%d'%(angx+1),alpha=0.6)
        ax.fill_between(t, (angles[:,0]-angles[:,2]).squeeze(),\
                            (angles[:,0]+angles[:,2]).squeeze(),\
                            alpha=0.1,color=colorlist[angx])
        ax.fill_between(t, (angles[:,0]-angles[:,1]).squeeze(),\
                            (angles[:,0]+angles[:,1]).squeeze(),\
                            alpha=0.05,color=colorlist[angx])

    ax.set_xlim([t[0],t[-1]])


    






def plottoaxis_offdecoderrocauc(offresponsedecode,ax,T,colors=colors2,labels=['train','test']):
    for tt in range(2):
        t = offresponsedecode[tt].times
        print(tt,len(offresponsedecode),offresponsedecode[0].shape)
        ax.plot(t,offresponsedecode[tt][:,0],linewidth=3,color=colors[tt],label=labels[tt])
        ax.fill_between(t, (offresponsedecode[tt][:,0]-offresponsedecode[tt][:,2]).squeeze(),\
                            (offresponsedecode[tt][:,0]+offresponsedecode[tt][:,2]).squeeze(),\
                            alpha=0.4,color=colors[tt])
        ax.fill_between(t, (offresponsedecode[tt][:,0]-offresponsedecode[tt][:,1]).squeeze(),\
                            (offresponsedecode[tt][:,0]+offresponsedecode[tt][:,1]).squeeze(),\
                            alpha=0.2,color=colors[tt])
    ax.set_xlim([T['starttime'],T['endtime']])
    ax.set_ylim([0.25,1.01])



def plottoaxis_offgradientdecoder(n_cuetrials,dtm,dtsem,dtm_c_s,ax,colorlist=['darkorange','chocolate'],labels=['raw','savgol']):

    t = np.arange(n_cuetrials)+1
    
    ax.plot(t,dtm[:,0],color=colorlist[0],linewidth=2,alpha=0.7,label=labels[0])
    ax.fill_between(t, dtm[:,0]-dtsem[:,0], dtm[:,0]+dtsem[:,0],color=colorlist[0],alpha=0.3)
#        ax.fill_between(np.arange(n_cuetrials)+1, dtm[:,0]-dtstd[:,0], dtm[:,0]+dtstd[:,0],color=colorlist[0],alpha=0.1)
    ax.plot(t,dtm_c_s,color=colorlist[1],linewidth=2,label=labels[1])
    
    ax.set_ylim([-0.05,1.05])
    


def plottoaxis_performancebarplot(n_cuetrials,bp_vals,bpx,ax,colorlist):
    t = np.arange(n_cuetrials)+1
    ax.bar( t, bp_vals[bpx], color=[colorlist[t] for t in bpx ], alpha=1  )
    ax.set_ylim([-2,2])
    ax.set_yticks([])





def plottoaxis_normalpdf(axs,m=0,s=1,y0=0,yh=1,color='grey',alpha=1.0):
    x = np.linspace(m-4*s,m+4*s,101)
    y = sp.stats.norm.pdf(x,m,s)
    y = y/np.max(y)
#    axs.plot(x,y0+y*yh,color=color,lw=2)
    axs.fill_between(x,y0,y0+y*yh,color=color,alpha=alpha,zorder=2)










def plottoaxis_jointtimecourserunspeedneural(surface,ax):
#    cmap = getcorrcolormap('firingrate')
    ax.pcolormesh(surface)#,cmap=cmap)







# subspaces


def plottoaxis_subspacestrajectory(X,coeffs,coords,ax=None):
    
    Xa = np.array(X[0])
    Xb = np.array(X[1])
        
    c1s = sp.signal.savgol_filter(coeffs[coords[0]][:,0], 15, 5, axis=0)
    c2s = sp.signal.savgol_filter(coeffs[coords[1]][:,0], 15, 5, axis=0)
    
    xa1 = np.array(Xa).mean(axis=0)[:,coords[0]]           # average over trials
    xa2 = np.array(Xa).mean(axis=0)[:,coords[1]]
    
    xb1 = np.array(Xb).mean(axis=0)[:,coords[0]]
    xb2 = np.array(Xb).mean(axis=0)[:,coords[1]]

    ax.plot(  c1s, c2s, linewidth=2,alpha=0.8,color='purple',label='coeff' )
    
    ax.plot(  xa1, xa2, linewidth=2,alpha=0.8,color='darkgoldenrod', label='class 1 pc' )
    
    ax.plot(  xb1, xb2, linewidth=2,alpha=0.8,color='gold', label='class 2 pc' )
    
    return




def plottoaxis_neuralandboundary(X,coeffs,tpc_coords,cpc_coords=[0,1],ax=None):
    
    tpc_coords = np.array(tpc_coords)
    cpc_coords = np.array(cpc_coords)

#    print('X shape, decision boundary shape',X[0][0].shape,np.array(coeffs).shape)

    PsXa, PsXb, PXa, PXb, PVa, PVb, Pcoeffs = nedi.get_orthogonalprojections_ontopcaxes(X,coeffs,tpc_coords)


    ax.plot(  PXa[:,0], PXa[:,1],'o', linewidth=2,alpha=0.4,color='fuchsia', label='class 1' )
    ax.plot(  PXb[:,0], PXb[:,1],'o', linewidth=2,alpha=0.4,color='darkturquoise', label='class 2' )

    ax.plot(  PsXa[:,0], PsXa[:,1],'o', linewidth=2,alpha=0.4,color='darkgoldenrod', label='coords class 1' )
    ax.plot(  PsXb[:,0], PsXb[:,1],'o', linewidth=2,alpha=0.4,color='gold', label='coords class 2' )

    
    ax.plot(  [0, PVa[0,cpc_coords[0]]],  [0, PVa[1,cpc_coords[0]]], linewidth=3,alpha=0.6,color='darkmagenta',label='class 1 pc 1' )
    ax.plot(  [0, PVa[0,cpc_coords[1]]],  [0, PVa[1,cpc_coords[1]]], linewidth=3,alpha=0.6,color='orchid',label='class 1 pc 2' )

    ax.plot(  [0, PVb[0,cpc_coords[0]]],  [0, PVb[1,cpc_coords[0]]], linewidth=3,alpha=0.6,color='darkblue',label='class 2 pc 1' )
    ax.plot(  [0, PVb[0,cpc_coords[1]]],  [0, PVb[1,cpc_coords[1]]], linewidth=3,alpha=0.6,color='royalblue',label='class 2 pc 2' )



    ax.plot(  [0, Pcoeffs[0]],  [0, Pcoeffs[1]], linewidth=3,alpha=1.0,color='red',label='decision' )

#    ax.arrow(  0,0, PVa[0,cpc_coords[0]],  PVa[1,cpc_coords[0]], width=0.05, alpha=0.6,color='darkmagenta',label='class 1 pc 1' )
#    ax.arrow(  0,0, PVa[0,cpc_coords[1]],  PVa[1,cpc_coords[1]], width=0.05, alpha=0.6,color='orchid',label='class 1 pc 2' )
#
#    ax.arrow(  0,0, PVb[0,cpc_coords[0]],  PVb[1,cpc_coords[0]], width=0.05, alpha=0.6,color='navy',label='class 2 pc 1' )
#    ax.arrow(  0,0, PVb[0,cpc_coords[1]],  PVb[1,cpc_coords[1]], width=0.05, alpha=0.6,color='royalblue',label='class 2 pc 2' )
#
#
#
#    ax.arrow(  0,0, Pcoeffs[0], Pcoeffs[1], width=0.05, alpha=1.0,color='lime',label='decision' )

    
    return














def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of `x` and `y`
    
    See how and why this works: https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
    
    This function has made it into the matplotlib examples collection:
    https://matplotlib.org/devdocs/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py
    
    Or, once matplotlib 3.1 has been released:
    https://matplotlib.org/gallery/index.html#statistics
    
    I update this gist according to the version there, because thanks to the matplotlib community
    the code has improved quite a bit.
    Parameters
    ----------
    x, y : array_like, shape (n, )
        Input data.
    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.
    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.
    Returns
    -------
    matplotlib.patches.Ellipse
    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)




def plottoaxis_gmm(axs, X, Y, means, covariances, colors=['m','b','r','y','g','c']):
    n_components = len(means)
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, colors[:n_components])):
        v, w = np.linalg.eigh(covar)
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / np.linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y == i):
            continue
        axs.scatter(X[Y == i, 0], X[Y == i, 1], .8, color=color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ell = Ellipse(mean, v[0], v[1], 180. + angle, color=color)
        ell.set_clip_box(axs.bbox)
        ell.set_alpha(0.5)
        axs.add_artist(ell)






        
if __name__ == '__main__':
    print('cameleon disappears in colors...')