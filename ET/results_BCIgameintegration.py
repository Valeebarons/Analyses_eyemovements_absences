# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:21:32 2023

@author: VB


Data analysis for pilot data from KH with BCI integration on monkey game
Measures eye movements and reaction times during interictal EEG activity (monkey trials) and during absences (red dots)


"""

import sys
#print(sys.path)
sys.path.append(r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\functions_ET\\') #location to scripts to import
import functions_BCIintegration
from IPython import get_ipython
# Make plots interactive
ipython = get_ipython() 
ipython.magic("matplotlib qt") # or %matplotlib qt 
# when error 'NoneType' object has no attribute 'to_csv' modify first column in tsv file_trial (take recording timestamp column, got to Data-->texttocolumns, save)


# 1) You will use two OpenSesame (OS) files: 
    #-1 monkey trials file: DELETE MANUALLY TRIALS IN OPENSESAME FILES WITH NO REAL RESPONSE (RT <1 (0. something), these are the trials in which press of trial before is recorded for the next trial too) and with dots (absences)
#    -2 absence file: create OS file called patnum_absences, in which you keep only the trials during absences (dot) to calculate RT 


#2) PREPROCESSING 
filename = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\pat5_pressall.tsv'
TobiiFolder = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\\'
subject_nr = 5
functions_BCIintegration.getRightStructure_paper5(filename, TobiiFolder, subject_nr)
#if there are some null trials in the OS file (response time 0. something) then you need to go through the function getRightstructure
#run it until the eye smoothing and delete the indeces end of the monkey trials to equal the number of indeces start. in fact, you'll see some of the
#indeces end are lower than indeces start, but idneces end need to be always higher. 


#3) VISUALIZE TRIALS MONKEY

#plot
subject_nr=5
filename = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\trials\\subj5_monkey_trial.tsv'
OSFile = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\OS_files\pat5_pressall.csv'

for i in range(41,52):
   functions_BCIintegration.plotYoverX_monkey(filename, OSFile,subject_nr,i, False)
    

#4) VISUALIZE TRIALS RED DOTS
subject_nr=5
filename = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\trials\\subj5_absences_trial.tsv'
## ALL ABSENCES
functions_BCIintegration.plotYoverX_absences(filename,subject_nr)

# 5) percentage eye movemements correct vs incorrect + total RT during absences + response rate during absences
OSFile = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\OS_files\pat5_absences.csv'
filename =  r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\trials\\subj5_absences_trial.tsv'
subject_nr=5
resrate_abs, average_RT, correct_eyeposition= functions_BCIintegration.absencetrials_output(filename,OSFile, subject_nr)



#7) RT subcomponents for monkey trial + response rate save as excel for each trial


# first open the tsv file in filename, select first column and then go to Data --> Text to column --> Finish --> save 
filename = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\trials\\subj5_monkey_trial.tsv'
TobiiFolder =  r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\trials\\'
OSFile = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\OS_files\pat5_pressall.csv'
subject_nr = 5
#start = '19' #why string before?
#stop = '22'
#AOI=[960-200,540+200 , 960+200, 540-200] 
functions_BCIintegration.ValuesToExcel_paper5(filename, OSFile,TobiiFolder, subject_nr)

