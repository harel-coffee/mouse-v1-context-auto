# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:01:55 2019

@author: mahajnal
"""







import numpy as np
import scipy as sp
import quantities as pq
import pandas as pd

import os
import pickle
# import dask.bag as db
import multiprocessing
from os import cpu_count
n_cpu = cpu_count()

#import matplotlib.image as mimg
import matplotlib.pyplot as plt
import matplotlib.lines
import matplotlib.colors as clrs
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

import neo

import neurophysiology as neph
import neurodiscover as nedi

import figs


# import preprocess










globalrecalculate = 0
globaldoplot = 0
globalsave = 0


# kilosort3 
# continuous_method = 'ks3ifr'; binsize = 100*pq.ms     # dt 10 ms, bin 100 ms

# kilosort2.5,   note, that for any analysis for 2.5 or 3 version is recommended to symlink to ks2ifr neo files
# continuous_method = 'ks25ifr'; binsize = 100*pq.ms     # dt 10 ms, bin 100 ms


# kilosort2
continuous_method = 'ks2ifr'; binsize = 100*pq.ms     # dt 10 ms, bin 100 ms
# continuous_method = 'ks2count'; binsize = 10*pq.ms     # dt 1 ms, bin 1 ms      for DT014
# continuous_method = 'ks2count'; binsize = 1*pq.ms     # dt 1 ms, bin 1 ms      for DT014



# janelia JRC
# continuous_method = 'instfr'; binsize = 100*pq.ms     # dt 10 ms, bin 100 ms
# continuous_method = 'count'; binsize = 200*pq.ms       # only bin is important, dt := bin





# define time windows for each trials:
# starttime: recording start for trial; similar stimulus: stimulus start time
# endtime: timepoint where recording ends in a trial
# epoch: interesting time window to calculate with and display; outside this will be the baseline, between start+trim:end-trim
# stim: stimulus presentation start and end times
T = {'dt':10*pq.ms,'bin':binsize,'ks2samplingrate':25*pq.kHz,\
     'starttime':-1500*pq.ms,'endtime':4500*pq.ms,'trimtime':0*pq.ms,\
     'epochstarttime':-500*pq.ms,'epochendtime':3500*pq.ms,\
     'stimstarttime':0*pq.ms,'stimendtime':3000*pq.ms,'rewardtime':2000*pq.ms}
if continuous_method=='count': T['dt'] = T['bin']
if continuous_method=='ks2count': T['dt'] = T['bin']
if continuous_method=='ks25count': T['dt'] = T['bin']
if continuous_method=='ks3count': T['dt'] = T['bin']
T['offsettime'] = T['starttime']
T['offset_idx'] = int((T['starttime']/T['dt']).magnitude)
T['start_idx'] = 0
T['end_idx'] = int((T['endtime']/T['dt']).magnitude) - T['offset_idx']
T['stimstart_idx'] = int((T['stimstarttime']/T['dt']).magnitude) - T['offset_idx']
T['stimend_idx'] = int((T['stimendtime']/T['dt']).magnitude) - T['offset_idx']
T['reward_idx'] = int((T['rewardtime']/T['dt']).magnitude) - T['offset_idx']
T['epochstart_idx'] = int((T['epochstarttime']/T['dt']).magnitude) - T['offset_idx']
T['epochend_idx'] = int((T['epochendtime']/T['dt']).magnitude) - T['offset_idx']
T['trim_idx'] = int((T['trimtime']/T['dt']).magnitude)
dt = T['dt'] # shortcut
# decoding time-course feature T['width_idx']
T['decodingwindow-widthtime'] = 50*pq.ms
T['width_idx'] = int(( T['decodingwindow-widthtime'] / T['dt'] ).magnitude)
T['videofps'] = 20.




# data storages

# # pre 2022
# pathdatamouse = '../../../data/ucla/2018-2020,jrc+behav+training-mat,events-csv/'
pathdatamouse = '../cache/events/'
# trialsfolder = 'trials-generated/'
trialsfolder = ''
pathdatamouse_ks2sorted = '../../../data/ucla/2018-2020,wigner,ks+phy/kilosort2/'
pathdatamouse_ks3sorted = '../../../data/ucla/2018-2020,wigner,ks+phy/kilosort3/'

# from 2021
# pathdatamouse = '../../../data/ucla/2021-2022,gonogo/'
# trialsfolder = ''
# pathdatamouse_ks2sorted = pathdatamouse
# pathdatamouse_ks25sorted = pathdatamouse
# pathdatamouse_ks3sorted = pathdatamouse





# cache storages

cacheprefix = '../cache/'



# results savepoints
# __________
# discovery:

resultpathprefix = '../results/'
# resultpathseries = 'phys/'
# resultpathseries = 'behaviour/'
# resultpathseries = 'symmetric/'
resultpathseries = 'motion/pca/'
# resultpathseries = 'motion/subspaces/'
# resultpathseries = 'firingrate/'
# resultpathseries = 'differences,temporal/'
# resultpathseries = 'differences,count/'
# resultpathseries = 'runspeed,temporal/'
# resultpathseries = 'subspaces/'
# resultpathseries = 'subspaces-acc/'
# resultpathseries = 'decoding/'
# resultpathseries = 'predictive/'
# resultpathseries = 'bayesian/'
# resultpathseries = 'vlgp/'
# resultpathseries = 'tca/'    # 'tca,instfr/'
# resultpathseries = 'pca/'
# resultpathseries = 'actionreward/'
# resultpathseries = 'tribe/'
# resultpathseries = 'anime/'


resultpath = resultpathprefix + resultpathseries
ext = '.png'
publish = False




# _____________
# publications:

#resultpath = '../../../publish/technicalreports/hfsp-report2019/'
#resultpath = '../../../publish/conferences/ccn2019/'
#resultpath = '../../../publish/journals/journal2019summer/figures/'
# resultpath = '../../../publish/conferences/mitt2022/poster/'
# ext = '.pdf'
#publish = True




print('Saving %s will be into '%ext,resultpath)












# # rate below which neurons are considered non-participating, and exlusion from trials; None to have no threshold
# activerate_threshold = 1.0 #None #1.0

# # sampling rate in the original data; do not change!
# f_s = 1000.

# # sampling bin in ms
# binsize = 100.



# only needed for runoff

# # runspeed filtering lowpass cutoff frequency [Hz] (1Hz  = 1/1000 ms)
# f_rs_lp = 4.9           #  (0.15 Hz =>  6.7s,  0.5 Hz -> 2s, 2 Hz -> 500 ms)

# # neural data filters
# f_n_bp_gamma = np.asarray([45,75]) # bandpass
# f_n_lp = 0.49 # cutoff for lowpass (this is only for long display in selectbyneurons 2/b variant figure)


layernames = ['superf','input','deep5','deep6']
spikewidthnames = ['narrow','broad']










