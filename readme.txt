Contains code for figures and statistics for the publication:

Márton Albert Hajnal, Karen Safaryan, Duy Tran, Michael Einstein, Mauricio Vallejo Martelo, Pierre-Olivier Polack, Peyman Golshani, Gergő Orbán
"Continuous multiplexed population representations of task context in the mouse primary visual cortex"


makefigure.py: will produce figure files from precalculated and cached results data

caches can be recalculated using the following files:

config.py: setup file for modules, sampling constants, experiment trial structure, output files, etc.
run.py: command center, choose which mice to include in a given analysis and choose which calculation routines to run
aggregateallmice.py: performs calculations and exports caches for data with multiple mice

preprocess.py: handles mouse experiment specific loaders for raw and cached data
neurophysiology.py: math and helper routines to preprocess spike train data into smoothed spike counts
neurodiscover.py: statistics and simple ai routines for decoders

All other files contain above referenced calculation routines in thematic groups

The files contain more calculations and routines than that were used in the publication.
The makefigure.py file should guide the reader as to which calculation routines are needed by searching
for the cached results files loaded within makefigure.py as pickle dump files in the calculation files



The following folder structure is necessary to be created before being able to cache calculation results and generate output:
../cache/*
../results/*
where * should be inferred from the recalculate = 1 branches of the calculation functions inside each calculation file
or the pickle.load(open()) functions from makefigure.py file



The raw input data structures are from files in .mat and phy2 via kilosort,
which are huge, and not shared, so Supp Fig 1 cannot be recreated.

We provide an internal version of the raw data in neo format (.nio files):
https://neo.readthedocs.io/en/stable/install.html



Download data from the accompanying repository:
doi://10.5281/zenodo.7065334
https://zenodo.org/record/7065334

neural-events,cache.zip, unpack in the parent folder of the code:
- cache/events/*.csv: behavioural event files in csv format
- cache/events/trainingbehaviour/*.mat: training behavioural event files in mat cell format
- cache/neo/*.nio: trial segmented spike counted neural activity, run speed
- cache/phys/*.pck: spike times and info of single activity units in pickle list format
                    provided as a convenience (also available within .nio files):



Questions about the code and data should be addressed to: hajnal.marton@wigner.mta.hu
