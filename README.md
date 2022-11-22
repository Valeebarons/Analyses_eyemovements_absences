# Analyses_eyemovements_absences
Here you will find the analyses of EEG and eye movements of patients with absence seizures used for our manuscript. 



The analyses are developed in:
- MATLAB (version R2018b) and the freely available EEGlab (version 2021.0) and Fieldtrip (version 20220714)
- Python 3.8

# EEG 
- absence_quantification.m
quantitative EEG features (i.e., amplitude, peak frequency, slope of peak frequency, centre of gravity) of absence seizures for 10 specific pediatric patients with absence seizure. 

- connectivityanalysis.m
analyses of bivariate connectivity using phase lag index at 3 Hz during absence seizures for 10 specific patients. 
unweighted graph analyses using PLI: hubness, cluster coefficient, path length. Use functions pathlength.m, plot_clustercoef.m

- dipfitFT_allpat.m 
source reconstruction of each time point of every absence seizures of 10 specific patients using dipole fitting via Fieldtrip.

- FEF.m 
from extracted positions of dipoles (from dipfitFT_allpat.m), derive dipoles lying within frontal eye field (FEF). 

(credits V. Barone, MC. Piastra)

# ET  

- eyemovements.py
preprocessing of ET from Tobii Nano Pro. 
extrapolation of eye movements variables (e.g., fixation, saccadic latency, saccadic duration, processing speed).
visualization of eye movements of one trial for one specific subject

(Credits: V. Barone, J. Boons)
