Contains code for figures and statistics for the publication:

Márton Albert Hajnal, Duy Tran, Michael Einstein, Mauricio Vallejo Martelo, Karen Safaryan, Pierre-Olivier Polack, Peyman Golshani, Gergő Orbán
"Continuous multiplexed population representations of task context in the mouse primary visual cortex"


makefigure.py: will produce figure files from precalculated and cached results data
config.py: setup file for modules, sampling constants, experiment trial structure, output files, etc.
run.py: command center, choose which mice to include in a given analysis and choose which calculation routines to run
aggregateallmice.py: performs calculations and exports caches for data with multiple mice

preprocess.py: handles mouse experiment specific loaders for raw and cached data
neurophysiology.py: math and helper routines to preprocess spike train data into smoothed spike counts
neurodiscover.py: statistics and simple ai routines for decoders

All other files contain above referenced calculation routines in thematic groups



The following folder structure is necessary to be created before being able to cache calculation results and generate output:
../cache/*
../results/*
where * should be inferred from the recalculate = 1 branches of the calculation functions inside each calculation file


The raw input data structures are from .mat and phy2 via various JRClust and kilosort versions.

The internal data structure is neo:
https://neo.readthedocs.io/en/stable/install.html
../cache/neo/*.nio cached smoothed data files are available upon reasonable request

The following data repository contains:
- spike times of single activity units in pickle list format
- event files in csv format
https://zenodo.org/record/5044249
phy2 format spike times are available upon reasonable request


The files contain more calculations and routines than that were used in the publication.
The makefigure.py file should guide the reader as to which calculation routines are needed by searching
for the cached results files loaded within makefigure.py as pickle dump files in the calculation files.


Questions about the code should be addressed to: hajnal.marton@wigner.mta.hu
