# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:26:58 2023

@author: VB

functions for results_BCIgameintegration.py
"""
from PyTrack.formatBridge import generateCompatibleFormat
from PyTrack.Stimulus import Stimulus
import pandas as pd
import numpy as np
import math
from math import atan2, degrees, sqrt
import os
import time
import ast   
import random
from matplotlib import pyplot
from numpy import arange
from matplotlib.patches import Rectangle, Circle
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D   
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from statistics import mean, pvariance, variance, pstdev, stdev, median

#################################################### PREPROCESSING ##############
def getRightStructure_paper5(filename, TobiiFolder, subject_nr):
    '''


    Parameters
    ----------
    filename : folder + name of tsv files from tobii without header
    TobiiFolder: folder where tsv files without header are stored
    subject_nr: number of subject, to save new preprocessed data file

    Returns
    -------
    raw_data : Returns the raw data in the correct format as .tsv file
    
    saves preprocessed data in the same directory as tsv

    '''
##remember to always remove the 'extra text' from the excel file (keep headers!), plus decimal separator = . and thousands separator = ,

    raw_data = pd.read_csv(filename, delimiter='\t', header=0)
    raw_data.columns = ['Recording timestamp', 'Event Value', 'Gaze2d_Left.x', 'Gaze2d_Left.y', 'ValL', 'Gaze2d_Right.x', 'Gaze2d_Right.y', 'ValR', 'GazeX', 'GazeY', 'PupilDiam_Left', 'PupValL', 'PupilDiam_Right', 'PupValR']
#Pytrack expects integers of time in milliseconds. timestamps from tobii are in microseconds
    raw_data['Recording timestamp'] = raw_data['Recording timestamp'].astype('float64')*1000
#raw_data['Recording timestamp'] = raw_data['Recording timestamp'].astype('float64')
    print('Converting data to right format...')
#Pygaze expects that the markers are numbers, so we convert strings to numbers first.
       
##Check if the Stimulus is done       
    a = raw_data['Event Value'].str.contains('StimulusOUT')
    for index, row in enumerate(a):
       if row == True:
          raw_data.at[index, 'Event Value'] = 22
    a = raw_data['Event Value'].str.contains('Stimulus OUT')
    for index, row in enumerate(a):
       if row == True:
          raw_data.at[index, 'Event Value'] = 22
#Check if a Stimulus appeared        
    a = raw_data['Event Value'].str.contains('StimulusIN')
    for index, row in enumerate(a):
        if row == True:
            #remove first 300 ms coz of trigger delay --> 0.3*fs. 
            #not do it for pat 1, acute TBI- 
          # raw_data.at[index+18, 'Event Value'] = 19
          raw_data.at[index, 'Event Value'] = 19
           
    a = raw_data['Event Value'].str.contains('Stimulusdownright Absence IN')
    for index, row in enumerate(a):
        if row == True:
           raw_data.at[index, 'Event Value'] = 50
           
    a = raw_data['Event Value'].str.contains('Stimulusdownright Absence OUT')
    for index, row in enumerate(a):
        if row == True:
           raw_data.at[index, 'Event Value'] = 51            

    a = raw_data['Event Value'].str.contains('Stimuluscenter Absence IN')
    for index, row in enumerate(a):
        if row == True:
           raw_data.at[index, 'Event Value'] = 60
           
    a = raw_data['Event Value'].str.contains('Stimuluscenter Absence OUT')
    for index, row in enumerate(a):
        if row == True:
           raw_data.at[index, 'Event Value'] = 61  
    a = raw_data['Event Value'].str.contains('Stimulusupleft Absence IN')
    for index, row in enumerate(a):
        if row == True:
           raw_data.at[index, 'Event Value'] = 70
           
    a = raw_data['Event Value'].str.contains('Stimulusupleft Absence OUT')
    for index, row in enumerate(a):
        if row == True:
           raw_data.at[index, 'Event Value'] = 71               
              

###Delete all the nans
    raw_data['Event Value'] = raw_data['Event Value'].replace(np.nan, 0)
    ### Find all values that are negative or outside of the screen and set them to non-valid ###
    df = raw_data
    jj = 0    
    print('Finding invalid datapoints...')
    while jj < len(df):
        if df.loc[jj, 'Gaze2d_Left.x'] <0.0 or df.loc[jj, 'Gaze2d_Left.x'] > 1920.0:
            df.loc[jj, 'Gaze2d_Left.x'] = -1
        if df.loc[jj, 'Gaze2d_Left.y'] <0.0 or df.loc[jj, 'Gaze2d_Left.y'] > 1080.0:
            df.loc[jj, 'Gaze2d_Left.y'] = -1
        if df.loc[jj, 'Gaze2d_Right.x'] <0.0 or df.loc[jj, 'Gaze2d_Right.x'] > 1920.0:
            df.loc[jj, 'Gaze2d_Right.x'] = -1
        if df.loc[jj, 'Gaze2d_Right.y'] <0.0 or df.loc[jj, 'Gaze2d_Right.y'] > 1080:
            df.loc[jj, 'Gaze2d_Right.y'] = -1
        if df.loc[jj, 'PupilDiam_Left'] <0.0:
            df.loc[jj, 'PupilDiam_Left'] = -1
        if df.loc[jj, 'PupilDiam_Right'] <0.0: 
            df.loc[jj, 'PupilDiam_Right'] = -1
        jj +=1
    
    df = df.fillna(-1)  
    ###Delete all the nans
    raw_data['Event Value'] = raw_data['Event Value'].replace(np.nan, 0)
    ### Find all values that are negative or outside of the screen and set them to non-valid ###
    df = raw_data
    jj = 0
    print('Finding invalid datapoints...')
    while jj < len(df):
        if df.loc[jj, 'Gaze2d_Left.x'] <0 or df.loc[jj, 'Gaze2d_Left.x'] > 1920.0:
            df.loc[jj, 'Gaze2d_Left.x'] = -1
        if df.loc[jj, 'Gaze2d_Left.y'] <0 or df.loc[jj, 'Gaze2d_Left.y'] > 1080.0:
            df.loc[jj, 'Gaze2d_Left.y'] = -1
        if df.loc[jj, 'Gaze2d_Right.x'] <0 or df.loc[jj, 'Gaze2d_Right.x'] > 1920.0:
            df.loc[jj, 'Gaze2d_Right.x'] = -1
        if df.loc[jj, 'Gaze2d_Right.y'] <0 or df.loc[jj, 'Gaze2d_Right.y'] > 1080:
            df.loc[jj, 'Gaze2d_Right.y'] = -1
        if df.loc[jj, 'PupilDiam_Left'] <1.0:
            df.loc[jj, 'PupilDiam_Left'] = -1
        if df.loc[jj, 'PupilDiam_Right'] <1.0:
            df.loc[jj, 'PupilDiam_Right'] = -1
        jj +=1

    df = df.fillna(-1)
    ###Interpolate small gaps ###
    jj = 1
    print('Interpolating left eye-x...')
    while jj < len(df)-3:
        if df.loc[jj, 'Gaze2d_Left.x'] == -1:
            if sum(df.loc[jj:jj+4, 'Gaze2d_Left.x']) ==-4:
                print('Too large to interpolate')
                while df.loc[jj, 'Gaze2d_Left.x'] == -1 and jj < len(df)-3:
                    jj +=1

            else:
                if df.loc[jj, 'Gaze2d_Left.x'] ==-1 and df.loc[jj+1, 'Gaze2d_Left.x'] == -1 and df.loc[jj+2, 'Gaze2d_Left.x'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Left.x']
                    value_after_x = df.loc[jj+3, 'Gaze2d_Left.x']
                    steps = 4

                    df.loc[jj, 'Gaze2d_Left.x'] = (((1/steps)*value_before_x)+((3/steps)*value_after_x))
                    df.loc[jj+1,'Gaze2d_Left.x'] = (((2/steps)*value_before_x)+((2/steps)*value_after_x))
                    df.loc[jj+2, 'Gaze2d_Left.x'] = (((3/steps)*value_before_x)+((1/steps)*value_after_x))
                elif df.loc[jj, 'Gaze2d_Left.x'] ==-1 and df.loc[jj+1, 'Gaze2d_Left.x'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Left.x']
                    value_after_x = df.loc[jj+2, 'Gaze2d_Left.x']
                    steps = 3
                    df.loc[jj,'Gaze2d_Left.x'] = ((2/steps)*value_before_x+((1/steps)*value_after_x))

                    df.loc[jj+1,'Gaze2d_Left.x'] = (((1/steps)*value_before_x)+((2/steps)*value_after_x))

                elif df.loc[jj, 'Gaze2d_Left.x'] ==-1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Left.x']
                    value_after_x = df.loc[jj+1, 'Gaze2d_Left.x']

                    steps = 2
                    df.loc[jj,'Gaze2d_Left.x'] = (value_before_x+value_after_x)/steps

        jj +=1


    jj = 1
    print('Interpolating left eye-y...')
    while jj < len(df)-3:
        if df.loc[jj, 'Gaze2d_Left.y'] == -1:
            if sum(df.loc[jj:jj+4, 'Gaze2d_Left.y']) ==-4:
                print('Too large to interpolate')
                while df.loc[jj, 'Gaze2d_Left.y'] == -1 and jj < len(df)-3:
                    jj +=1

            else:
                if df.loc[jj, 'Gaze2d_Left.y'] ==-1 and df.loc[jj+1, 'Gaze2d_Left.y'] == -1 and df.loc[jj+2, 'Gaze2d_Left.y'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Left.y']
                    value_after_x = df.loc[jj+3, 'Gaze2d_Left.y']
                    steps = 4

                    df.loc[jj, 'Gaze2d_Left.y'] = (((1/steps)*value_before_x)+((3/steps)*value_after_x))
                    df.loc[jj+1,'Gaze2d_Left.y'] = (((2/steps)*value_before_x)+((2/steps)*value_after_x))
                    df.loc[jj+2, 'Gaze2d_Left.y'] = (((3/steps)*value_before_x)+((1/steps)*value_after_x))
                elif df.loc[jj, 'Gaze2d_Left.y'] ==-1 and df.loc[jj+1, 'Gaze2d_Left.y'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Left.y']
                    value_after_x = df.loc[jj+2, 'Gaze2d_Left.y']
                    steps = 3
                    df.loc[jj,'Gaze2d_Left.y'] = ((2/steps)*value_before_x+((1/steps)*value_after_x))

                    df.loc[jj+1,'Gaze2d_Left.y'] = (((1/steps)*value_before_x)+((2/steps)*value_after_x))

                elif df.loc[jj, 'Gaze2d_Left.y'] ==-1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Left.y']
                    value_after_x = df.loc[jj+1, 'Gaze2d_Left.y']

                    steps = 2
                    df.loc[jj,'Gaze2d_Left.y'] = (value_before_x+value_after_x)/steps

        jj +=1

    jj = 1
    print('Interpolating right eye-y...')
    while jj < len(df)-3:
        if df.loc[jj, 'Gaze2d_Right.y'] == -1:
            if sum(df.loc[jj:jj+4, 'Gaze2d_Right.y']) ==-4:
                print('Too large to interpolate')
                while df.loc[jj, 'Gaze2d_Right.y'] == -1 and jj < len(df)-3:
                    jj +=1

            else:
                if df.loc[jj, 'Gaze2d_Right.y'] ==-1 and df.loc[jj+1, 'Gaze2d_Right.y'] == -1 and df.loc[jj+2, 'Gaze2d_Right.y'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Right.y']
                    value_after_x = df.loc[jj+3, 'Gaze2d_Right.y']
                    steps = 4

                    df.loc[jj, 'Gaze2d_Right.y'] = (((1/steps)*value_before_x)+((3/steps)*value_after_x))
                    df.loc[jj+1,'Gaze2d_Right.y'] = (((2/steps)*value_before_x)+((2/steps)*value_after_x))
                    df.loc[jj+2, 'Gaze2d_Right.y'] = (((3/steps)*value_before_x)+((1/steps)*value_after_x))
                elif df.loc[jj, 'Gaze2d_Right.y'] ==-1 and df.loc[jj+1, 'Gaze2d_Right.y'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Right.y']
                    value_after_x = df.loc[jj+2, 'Gaze2d_Right.y']
                    steps = 3
                    df.loc[jj,'Gaze2d_Right.y'] = ((2/steps)*value_before_x+((1/steps)*value_after_x))

                    df.loc[jj+1,'Gaze2d_Right.y'] = (((1/steps)*value_before_x)+((2/steps)*value_after_x))

                elif df.loc[jj, 'Gaze2d_Right.y'] ==-1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Right.y']
                    value_after_x = df.loc[jj+1, 'Gaze2d_Right.y']

                    steps = 2
                    df.loc[jj,'Gaze2d_Right.y'] = (value_before_x+value_after_x)/steps

        jj +=1

    jj = 1
    print('Interpolating right eye-x...')
    while jj < len(df)-3:
        if df.loc[jj, 'Gaze2d_Right.x'] == -1:
            if sum(df.loc[jj:jj+4, 'Gaze2d_Right.x']) ==-4:
                print('Too large to interpolate')
                while df.loc[jj, 'Gaze2d_Right.x'] == -1 and jj < len(df)-3:
                    jj +=1

            else:
                if df.loc[jj, 'Gaze2d_Right.x'] ==-1 and df.loc[jj+1, 'Gaze2d_Right.x'] == -1 and df.loc[jj+2, 'Gaze2d_Right.x'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Right.x']
                    value_after_x = df.loc[jj+3, 'Gaze2d_Right.x']
                    steps = 4

                    df.loc[jj, 'Gaze2d_Right.x'] = (((1/steps)*value_before_x)+((3/steps)*value_after_x))
                    df.loc[jj+1,'Gaze2d_Right.x'] = (((2/steps)*value_before_x)+((2/steps)*value_after_x))
                    df.loc[jj+2, 'Gaze2d_Right.x'] = (((3/steps)*value_before_x)+((1/steps)*value_after_x))
                elif df.loc[jj, 'Gaze2d_Right.x'] ==-1 and df.loc[jj+1, 'Gaze2d_Right.x'] == -1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Right.x']
                    value_after_x = df.loc[jj+2, 'Gaze2d_Right.x']
                    steps = 3
                    df.loc[jj,'Gaze2d_Right.x'] = ((2/steps)*value_before_x+((1/steps)*value_after_x))

                    df.loc[jj+1,'Gaze2d_Right.x'] = (((1/steps)*value_before_x)+((2/steps)*value_after_x))

                elif df.loc[jj, 'Gaze2d_Right.x'] ==-1:
                    value_before_x = df.loc[jj-1, 'Gaze2d_Right.x']
                    value_after_x = df.loc[jj+1, 'Gaze2d_Right.x']

                    steps = 2
                    df.loc[jj,'Gaze2d_Right.x'] = (value_before_x+value_after_x)/steps

        jj +=1

    ### Smoothen the signal ###
    dfSmooth = df
    jj = 1
    print('Started smoothening')
    while jj < len(df)-1:
        if not any(x == -1 for x in df.loc[jj-2:jj+2, 'Gaze2d_Left.x']):
            dfSmooth.loc[jj, 'Gaze2d_Left.x'] = sum(df.loc[jj-4:jj+4, 'Gaze2d_Left.x'])/9

        if not any(x == -1 for x in df.loc[jj-2:jj+2, 'Gaze2d_Left.y']):
            dfSmooth.loc[jj, 'Gaze2d_Left.y'] = sum(df.loc[jj-4:jj+4, 'Gaze2d_Left.y'])/9

        if  not any(x == -1 for x in df.loc[jj-2:jj+2, 'Gaze2d_Right.x']):
            dfSmooth.loc[jj, 'Gaze2d_Right.x'] = sum(df.loc[jj-4:jj+4, 'Gaze2d_Right.x'])/9

        if  not any(x == -1 for x in df.loc[jj-2:jj+2, 'Gaze2d_Right.y']):
            dfSmooth.loc[jj, 'Gaze2d_Right.y'] = sum(df.loc[jj-4:jj+4, 'Gaze2d_Right.y'])/9

        if not any(x == -1 for x in df.loc[jj-2:jj+2, 'PupilDiam_Left']):
            dfSmooth.loc[jj, 'PupilDiam_Left'] = sum(df.loc[jj-4:jj+4, 'PupilDiam_Left'])/9

        if not any(x == -1 for x in df.loc[jj-2:jj+2, 'PupilDiam_Right']):
            dfSmooth.loc[jj, 'PupilDiam_Right'] = sum(df.loc[jj-4:jj+4, 'PupilDiam_Right'])/9

        jj+=1
    
    dfSmooth['GazeX'] = (dfSmooth['Gaze2d_Left.x']+dfSmooth['Gaze2d_Right.x'])/2
    dfSmooth['GazeY'] = (dfSmooth['Gaze2d_Left.y']+dfSmooth['Gaze2d_Right.y'])/2
    
    
    
    # save file between trials
    aa= list(dfSmooth['Event Value'])
# when indices end and start are different i look in opensesame file the trials in which response time is <1 ms and i manually remove all the indices end who are higher than corresponding indeces start

    indeces_start= [i for i in range(len(aa)) if aa[i] == 19]
    indeces_end= [i for i in range(len(aa)) if aa[i] == 22]
    
    
    # remove trials where stimulus does not appear (e.g. they press during the cue)
    indeces_end1 = indeces_end 
    for i, x in enumerate(indeces_start):
        if indeces_start[i]>indeces_end1[i]:
            del(indeces_end1[i])
        
    
    mm = pd.DataFrame()
    mm1 = pd.DataFrame()
    for i in range(len(indeces_start)):
        mm = dfSmooth.iloc[indeces_start[i]:indeces_end[i]+1,:]
        mm1 = mm1.append(mm)
    
    outputName=TobiiFolder+'trials\subj'+str(subject_nr)+'_monkey_trial.tsv'
    mm1.to_csv(outputName, sep= '\t', index=False) 
    
    indeces_start_downright= [i for i in range(len(aa)) if aa[i] == 50]
    indeces_end_downright= [i for i in range(len(aa)) if aa[i] == 51]
    
    indeces_start_center= [i for i in range(len(aa)) if aa[i] == 60]
    indeces_end_center= [i for i in range(len(aa)) if aa[i] == 61]
    
    indeces_start_upleft= [i for i in range(len(aa)) if aa[i] == 70]
    indeces_end_upleft= [i for i in range(len(aa)) if aa[i] == 71]
    
    nn= pd.DataFrame()
    nn1 = pd.DataFrame()
    for i in range(len(indeces_start_downright)):
        nn = dfSmooth.iloc[indeces_start_downright[i]:indeces_end_downright[i]+1,:]
        nn1 = nn1.append(nn)
        nn = dfSmooth.iloc[indeces_start_center[i]:indeces_end_center[i]+1,:]
        nn1 = nn1.append(nn)
        nn = dfSmooth.iloc[indeces_start_upleft[i]:indeces_end_upleft[i]+1,:]
        nn1 = nn1.append(nn)
        
    outputName=TobiiFolder+'trials\subj'+str(subject_nr)+'_absences_trial.tsv'
    nn1.to_csv(outputName, sep= '\t', index=False) 
    
    
################################3 VISUALIZE TRIAL MONKEY
def plotYoverX_monkey(filename, OSFile, subject_nr, trialtoplot, saveFigure):
    

######################################################### IF BAD CALIBRATION #################################
    
    '''
    command in console to plot animations: %matplotlib qt
    
    
    filename: name of tsv file after getRightStructure_paper5()
    
    OSFile: name of OSfile after removing 0.1
    Automatic of one output logfile from Opensesame. (.csv). 
    
    Subject_nr: subject number
 
    Options for 'start':
    2. 19 --> Starts at trigger 'Stimulus appears'
    
    Options for 'stop':
    3. 22 --> Starts at trigger 'Stimulus done'
    
    Trialtoplot:
        choose the trial you want to plot

    Options for 'saveFigure':
    1. True: saves figure
    2. False: do not save figure, only plot in Python.


    '''
    if subject_nr ==3:
       start = 19
       stop = 22      
    else:  
       start = '19'
       stop = '22'    
    df1 = pd.read_csv(filename, sep= '\t')
    OSdata = pd.read_csv(OSFile, delimiter = ',')
    aa= df1['Event Value']
    indeces_start = aa[aa == start]
    indeces_end = aa[aa == stop]
    df = pd.DataFrame()

  
    i = trialtoplot
    df = df1.iloc[indeces_start.index[i]:indeces_end.index[i]+1,:] 
    y1 = df['GazeX']-(1920/2) 
    y2 = -(df['GazeY']-(1080/2))
    fig, ax = plt.subplots()
    line, = ax.plot(y1, y2, color='k', label='eye movements')
    AOI_stim = [(list(OSdata['x_pos1'])[i])-180, (list(OSdata['y_pos1'])[i])+180,(list(OSdata['x_pos1'])[i])+200, (list(OSdata['y_pos1'])[i])-200]
    radius = sqrt(abs(((AOI_stim[2] - AOI_stim[0])*(AOI_stim[3] - AOI_stim[1]))/math.pi))
    patch = Circle((AOI_stim[0]+180, AOI_stim[1]-180),
                          radius,
                          color='r', fill=False, linestyle = '-')
    #patch = Rectangle((AOI_stim[0], AOI_stim[1]),
      #                    (AOI_stim[2] - AOI_stim[0]),
      #                    (AOI_stim[3] - AOI_stim[1]),
      #                    color='r', fill=False, linestyle = '-')

    
    ax.set_title('Eye movements trial'+str(i))
    ax.set_ylabel('Y movements')
    ax.set_xlabel('X movements')
    ax.set_ylim(-600,600)
    ax.set_xlim(-1100,1100)
    ax.add_patch(patch)
    
   

    def update(num, x, y, line):
        line.set_data(y1[:num], y2[:num])
        line.set_color('k')
        line.set_linewidth(1.5)
        return line,

   # line2, = ax.plot(y1, y2)
    ani1 = animation.FuncAnimation(fig, update, len(y1), fargs=[y1, y2, line,], interval=30, blit=True)

    plt.legend() 
    plt.show
    if saveFigure == True:
       ani1.save(r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\trials'+str(subject_nr)+'_monkeytrial'+str(i)+'.gif')
       
       
def plotYoverX_absences(filename, subject_nr):
    
    df1 = pd.read_csv(filename, sep= '\t')
    df1['Event Value'] = df1['Event Value'].replace('Absence started', '0')
    df1['Event Value'] = df1['Event Value'].replace('Absence ended', '0')
    aa= df1['Event Value']
    if subject_nr==3 or subject_nr==5:
       res = [eval(i) for i in list(aa)] #make event value a list of integers
    else:
       res=aa 
    hfont = {'fontname':'Times New Roman', 'size':16}    
    if 50 and 51 in res:
       start_downright= [i for i, n in enumerate(res) if n == 50]
       end_downright = [i for i, n in enumerate(res) if n == 51]
       df_downright = pd.DataFrame()
       for kk in range(len(start_downright)): #plot all the absence down right 
           df_downright = df1.iloc[start_downright[kk]: end_downright[kk],:] 
           y1 = df_downright['GazeX']-(1920/2) 
           y2 = -(df_downright['GazeY']-(1080/2))
           fig, ax = plt.subplots()
           line, = ax.plot(y1, y2, color='k', label='eye movements')
           AOI_stim = [0, 0,800, -475]
           patch = Rectangle((AOI_stim[0], AOI_stim[1]),
                                 (AOI_stim[2] - AOI_stim[0]),
                                 (AOI_stim[3] - AOI_stim[1]),
                                 color='r', fill=False, linestyle = '-')
           ax.set_title('Eye movements downright, absence num '+ str(kk+1), **hfont)
           ax.set_ylabel('Y movements', **hfont)
           ax.set_xlabel('X movements', **hfont)
           ax.set_ylim(-600,600)
           ax.set_xlim(-1100,1100)
           plt.yticks(**hfont)
           plt.xticks(**hfont)
           ax.add_patch(patch)
           def update(num, x, y, line):
               line.set_data(y1[:num], y2[:num])
               line.set_color('k')
               line.set_linewidth(1.5)
               return line,
          # line2, = ax.plot(y1, y2)
           ani1 = animation.FuncAnimation(fig, update, len(y1), fargs=[y1, y2, line,], interval=30, blit=True)
           plt.show
        #   if saveFigure == True:
         #     ani1.save(r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\paper5\tobiinoheader\trials'+str(subject_nr)+'_absence_downright'+str(kk+1)+'.gif')
              
        
    if 60 and 61 in res: 
       start_center= [i for i, n in enumerate(res) if n == 60]
       end_center = [i for i, n in enumerate(res) if n == 61]   
       df_center = pd.DataFrame()
       for kk in range(len(start_center)): #plot all the absence down right 
           df_center = df1.iloc[start_center[kk]: end_center[kk],:] 
           y1 = df_center['GazeX']-(1920/2) 
           y2 = -(df_center['GazeY']-(1080/2))
           fig, ax = plt.subplots()
           line, = ax.plot(y1, y2, color='k', label='eye movements')
           AOI_stim = [-400,235, 400,-235]
           patch = Rectangle((AOI_stim[0], AOI_stim[1]),
                                 (AOI_stim[2] - AOI_stim[0]),
                                 (AOI_stim[3] - AOI_stim[1]),
                                 color='r', fill=False, linestyle = '-')
           ax.set_title('Eye movements center, absence num '+ str(kk+1), **hfont)
           ax.set_ylabel('Y movements', **hfont)
           ax.set_xlabel('X movements', **hfont)
           ax.set_ylim(-600,600)
           ax.set_xlim(-1100,1100)
           plt.yticks(**hfont)
           plt.xticks(**hfont)
           ax.add_patch(patch)
           def update(num, x, y, line):
               line.set_data(y1[:num], y2[:num])
               line.set_color('k')
               line.set_linewidth(1.5)
               return line,
          # line2, = ax.plot(y1, y2)
           ani1 = animation.FuncAnimation(fig, update, len(y1), fargs=[y1, y2, line,], interval=30, blit=True)
           plt.show
        
    if 70 and 71 in res:
       start_upleft= [i for i, n in enumerate(res) if n == 70]
       end_upleft = [i for i, n in enumerate(res) if n == 71]        
       df_upleft = pd.DataFrame()
       for kk in range(len(start_center)): #plot all the absence down right 
           df_upleft = df1.iloc[start_upleft[kk]: end_upleft[kk],:] 
           y1 = df_upleft['GazeX']-(1920/2) 
           y2 = -(df_upleft['GazeY']-(1080/2))
           fig, ax = plt.subplots()
           line, = ax.plot(y1, y2, color='k', label='eye movements')
           AOI_stim = [-800, 475, 0,0]
           patch = Rectangle((AOI_stim[0], AOI_stim[1]),
                                 (AOI_stim[2] - AOI_stim[0]),
                                 (AOI_stim[3] - AOI_stim[1]),
                                 color='r', fill=False, linestyle = '-')
           ax.set_title('Eye movements upleft, absence num '+ str(kk+1),**hfont)
           ax.set_ylabel('Y movements',**hfont)
           ax.set_xlabel('X movements',**hfont)
           ax.set_ylim(-600,600)
           ax.set_xlim(-1100,1100)
           plt.yticks(**hfont)
           plt.xticks(**hfont)
           ax.add_patch(patch)
           def update(num, x, y, line):
               line.set_data(y1[:num], y2[:num])
               line.set_color('k')
               line.set_linewidth(1.5)
               return line,
          # line2, = ax.plot(y1, y2)
           ani1 = animation.FuncAnimation(fig, update, len(y1), fargs=[y1, y2, line,], interval=30, blit=True)
           plt.show 
           
           
def absencetrials_output(filename, OSFile, subject_nr):
    
    
    correct_eyeposition = pd.DataFrame()
    df1 = pd.read_csv(filename, sep= '\t')
    df1['Event Value'] = df1['Event Value'].replace('Absence started', '0')
    df1['Event Value'] = df1['Event Value'].replace('Absence ended', '0')
    aa= df1['Event Value']
    if subject_nr==3 or subject_nr==5:
       res = [eval(i) for i in list(aa)] #make event value a list of integers
    else:
       res=aa 
    if 50 and 51 in res:
       start_downright= [i for i, n in enumerate(res) if n == 50]
       end_downright = [i for i, n in enumerate(res) if n == 51]
       df_downright = pd.DataFrame()
       for kk in range(len(start_downright)): #plot all the absence down right 
           df_downright = df1.iloc[start_downright[kk]: end_downright[kk],:] 
           y1 = df_downright['GazeX']-(1920/2) 
           y2 = -(df_downright['GazeY']-(1080/2))
           points = np.transpose(np.vstack((list(y1), list(y2))))    
           AOI_stim = [0, 0,800, -475]
           patch = Rectangle((AOI_stim[0], AOI_stim[1]),
                                 (AOI_stim[2] - AOI_stim[0]),
                                 (AOI_stim[3] - AOI_stim[1]),
                                 color='r', fill=False, linestyle = '-')
           containsr = patch.contains_points(points)
           gaze_aoi_newr = np.asarray(containsr, dtype=int)
           
           if gaze_aoi_newr.any()==False or len(gaze_aoi_newr[gaze_aoi_newr==1])<4:
              correct_eyeposition.loc[kk,'downright']=0
           else:   
              correct_eyeposition.loc[kk,'downright']=1
    if 60 and 61 in res:
        start_center= [i for i, n in enumerate(res) if n == 60]
        end_center = [i for i, n in enumerate(res) if n == 61]
        df_center = pd.DataFrame()
        for kk in range(len(start_center)): #plot all the absence down right 
            df_center = df1.iloc[start_center[kk]: end_center[kk],:] 
            y1 = df_center['GazeX']-(1920/2) 
            y2 = -(df_center['GazeY']-(1080/2))
            points = np.transpose(np.vstack((list(y1), list(y2))))    
            AOI_stim = [-400,235, 400,-235]
            patch = Rectangle((AOI_stim[0], AOI_stim[1]),
                                  (AOI_stim[2] - AOI_stim[0]),
                                  (AOI_stim[3] - AOI_stim[1]),
                                  color='r', fill=False, linestyle = '-')
            containsr = patch.contains_points(points)
            gaze_aoi_newr = np.asarray(containsr, dtype=int)
            
            if gaze_aoi_newr.any()==False or  len(gaze_aoi_newr[gaze_aoi_newr==1])<4:
               correct_eyeposition.loc[kk,'center']=0
            else:   
               correct_eyeposition.loc[kk,'center']=1   
    if 70 and 71 in res:
         start_upleft= [i for i, n in enumerate(res) if n == 70]
         end_upleft = [i for i, n in enumerate(res) if n == 71]
         df_upleft = pd.DataFrame()
         for kk in range(len(start_upleft)): #plot all the absence down right 
             df_upleft= df1.iloc[start_upleft[kk]: end_upleft[kk],:] 
             y1 = df_upleft['GazeX']-(1920/2) 
             y2 = -(df_upleft['GazeY']-(1080/2))
             points = np.transpose(np.vstack((list(y1), list(y2))))    
             AOI_stim = [-800, 475, 0,0]
             patch = Rectangle((AOI_stim[0], AOI_stim[1]),
                                   (AOI_stim[2] - AOI_stim[0]),
                                   (AOI_stim[3] - AOI_stim[1]),
                                   color='r', fill=False, linestyle = '-')
             containsr = patch.contains_points(points)
             gaze_aoi_newr = np.asarray(containsr, dtype=int)
             
             if gaze_aoi_newr.any()==False or  len(gaze_aoi_newr[gaze_aoi_newr==1])<4:
                correct_eyeposition.loc[kk,'upleft']=0
             else:   
                correct_eyeposition.loc[kk,'upleft']=1   
                
    # to extract numbers from RTs         
    def extract_numbers_from_brackets(s):
        numbers = ast.literal_eval(s)
    # Convert the list of strings to a list of integers
        numbers = [int(num) for num in numbers]
        return numbers
    
    OSdata = pd.read_csv(OSFile, delimiter = ',')  
    B = OSdata['responseTime_abs'].apply(extract_numbers_from_brackets)    
    
    #average without zeros of RT during dot trials
    average_RT=[]
    for mm in range(len(B)):
        average_RT1= np.nanmean(B[mm])
        average_RT.append(average_RT1)
    
    #response rate
    resrate_abs =[]
    C = OSdata['resRate_abs'].apply(extract_numbers_from_brackets)       
    for nn in range(len(C)):  
        if len(C[nn]) ==1:
            resrate_abs.append(33) #percentage of response given
        elif len(C[nn]) ==2:    
            resrate_abs.append(67)
        elif len(C[nn]) ==3:    
            resrate_abs.append(100)
        elif len(C[nn]) ==0:    
            resrate_abs.append(0)    
           
           
    return resrate_abs, average_RT, correct_eyeposition
        

def DataOfInterest(filename, start, stop, AOI): #before transfor first column of trials results
    '''


    Parameters
    ----------
    filename : Give filename with only data between stimuli - data on events will be lost (it is -1s so does not matter)
    start : Start trigger for the analysis
    stop : End trigger for the analysis
    AOI : Area of interest

    Returns
    -------
    ConvData : Outputs csv file, which is needed for PyTrack. 
    sensor_dict : Needed for Pyrack, which has some information about screen, sampling frquency and AOI

    '''
    generateCompatibleFormat(filename,
                             device = 'tobii',
                             start = start,
                             stop = stop
                             )
    ConvData = pd.read_csv(filename[:-4]+".csv")
    sensor_dict = {
        "EyeTracker": {
            "Sampling_Freq": 60,
            "Display_width": 1920,
            "Display_height": 1080,
            "aoi": AOI
            }
        }
    return ConvData, sensor_dict



def find_subcomponentsRT_VT(filename, OSFile,OSdata, start, stop, AOI, Plot):

    '''
    the function finds the saccadic latency of the first saccades with a velocity higher than 35 degress/s
    and  its amplitude, just for the trials in which the target has been visited. 
    Parameters
    ----------

    Options for 'filename':
    Choose one of the generated filenames by the function divideData. (.tsv)
    
    Options for 'OSFile':
    Choose one output logfile from Opensesame. (.csv)

    Options for 'start':
        1. 10 --> Starts at trigger 'Start trial'
        2. 19 --> Starts at trigger 'Stimulus appears'
        3. 22 --> Starts at trigger 'Stimulus done'

    Options for 'stop':
        1. 10 --> Starts at trigger 'Start trial'
        2. 19 --> Starts at trigger 'Stimulus appears'
        3. 22 --> Starts at trigger 'Stimulus done'
        4. None --> Each trial runs from start to start

    Options for 'AOI_fixcross':
        AOI_stim=[960-200,540+200 , 960+200, 540-200]
        
    OPtions for Plot:
        True --> plot
        False --> no plots

    Returns
    -------
    first_Saccade : Returns a dataframe with all parameters obtained by the algorithms.
    '''
    
    visitedAoi = []
    visitedFixcross = []
    percentageonAOI = []
    percentageonFixcross = []
    saccade_end_index = []
    stimul_name = []
    index_visitedAoi = []
    rtime_et =[]
    angular_velocity_all = []
    currentStim_tot = [1]
    index_saclat = []
    amplitude_saccade = []
    saclat_ms1 = []
    visual_rs = []
    index_of_wrong_values = []
    mean_velocity_saccade = []
    processing_speed = []
    sum_vrssaclat = []
    visual_rs_ms = []
    percprocspeed = [] 
    distance_target =[]
    percvisualrspeed = []
    #gain = []
    percentage_saclatency = []
    first_saccade = pd.DataFrame(columns=[ 'saccadic latency', 'visual reaction time', 'processing speed', 'saccade onset_index', 'saccadic amplitude', 'percentage of saccadic latency', 'percentage of visual reaction time', 'percentage of processing speed', 'total RT', 'velocity of saccade'])  # initialize array for saccadic latency

    OSdata = pd.read_csv(OSFile, delimiter = ',')
            # Get the Tobii Data
    df, sensor_dict = DataOfInterest(filename,
                                     start=start, stop=stop,
                                     AOI=AOI)
# Find out what the amount of stimuli are
    if df['StimulusName'].iloc[-1][-3] != '_' and df['StimulusName'].iloc[-1][-2] != '_':
       lastStim = int(df['StimulusName'].iloc[-1][-3:])
    elif df['StimulusName'].iloc[-1][-2] != '_':
       lastStim = int(df['StimulusName'].iloc[-1][-2:])
    else:
       lastStim = int(df['StimulusName'].iloc[-1][-1:])
    
# Calculate the number of degrees that correspond to a single pixel. 
 #usually it is a very small value, arond 0.03.
    h = 19.5  # Monitor height in cm
    d = 60  # Distance between monitor and participant in cm. always 57cm on tobii logfile
    r = 1080  # Vertical resolution of the monitor in pixels
# size_monkey=100 # The stimulus size in pixels. 
    deg_per_px = degrees(atan2(.5*h, d)) / (.5*r)
    dt = 1.0/60  # delta t to calculate velocity. = sampling time or time between two samples
    fs = 60

# divide gaze data per stimulus
    stims_y = np.zeros((1000, lastStim+1))
    k = 0
    i = 0
    jj = 1
# Convert the 1D array to 2D array, with a column for each trial. Do the same for stims_x. stim_y and x are the average y and x gaze data per trial, 1000 is just a random number to prelocate the variable that is big enough
    if (OSFile[-21] == '1' and OSFile[-20] == '5' and OSFile[-18] == '7') or (OSFile[-56]=='p' and OSFile[-57]=='m' and OSFile[-12]=='2' and OSFile[-13]=='j') or (OSFile[-21] == '2' and OSFile[-20] == '5' and OSFile[-18] == '7'):
        while jj < len(df):
            if df.loc[jj, 'StimulusName'] != df.loc[jj-1, 'StimulusName']:
                k += 1
                i = 0
            stims_y[i, k] = df.loc[jj, 'GazeRighty'] 
            i += 1
            jj += 1
 
        stims_y[stims_y == 0] = 'nan'

        stims_x = np.zeros((1000, lastStim+1))
        k = 0
        i = 0
        jj = 1
        while jj < len(df):
            if df.loc[jj, 'StimulusName'] != df.loc[jj-1, 'StimulusName']:
               k += 1
               i = 0
            stims_x[i, k] = df.loc[jj, 'GazeRightx'] 
            i += 1
            jj += 1

        stims_x[stims_x == 0] = 'nan'
    elif (OSFile[-21] == '2' and OSFile[-20] == '6' and OSFile[-18] == '7' and OSFile[-24] == 'B' ):
        while jj < len(df):
            if df.loc[jj, 'StimulusName'] != df.loc[jj-1, 'StimulusName']:
                k += 1
                i = 0
            stims_y[i, k] = df.loc[jj, 'GazeLefty'] 
            i += 1
            jj += 1
 
        stims_y[stims_y == 0] = 'nan'

        stims_x = np.zeros((1000, lastStim+1))
        k = 0
        i = 0
        jj = 1
        while jj < len(df):
            if df.loc[jj, 'StimulusName'] != df.loc[jj-1, 'StimulusName']:
               k += 1
               i = 0
            stims_x[i, k] = df.loc[jj, 'GazeLeftx'] 
            i += 1
            jj += 1

        stims_x[stims_x == 0] = 'nan'
    else:   
# Convert the 1D array to 2D array, with a column for each trial. Do the same for stims_x. stim_y and x are the average y and x gaze data per trial, 1000 is just a random number to prelocate the variable that is big enough
        while jj < len(df):
            if df.loc[jj, 'StimulusName'] != df.loc[jj-1, 'StimulusName']:
               k += 1
               i = 0
            stims_y[i, k] = (df.loc[jj, 'GazeLefty'] + df.loc[jj, 'GazeRighty'])/2
            i += 1
            jj += 1
    
        stims_y[stims_y == 0] = 'nan'
        
        stims_x = np.zeros((1000, lastStim+1))
        k = 0
        i = 0
        jj = 1
        while jj < len(df):
            if df.loc[jj, 'StimulusName'] != df.loc[jj-1, 'StimulusName']:
               k += 1
               i = 0
            stims_x[i, k] = (df.loc[jj, 'GazeLeftx'] + df.loc[jj, 'GazeRightx'])/2
            i += 1
            jj += 1

        stims_x[stims_x == 0] = 'nan'
        
   
    
# Loop for every trial and determine if the subject has visited the stimuli
    for index in range(len(stims_y.T)):  # .T = transpose. same as np.transpose()
        nan_array_x = np.isnan(stims_x[:, index])
        nan_array_y = np.isnan(stims_y[:, index])
        non_nan_array_x = ~ nan_array_x
        non_nan_array_y = ~ nan_array_y
   # gaze points per each trial - remove the first 300 ms because waiting time for trigger
        points = np.transpose(np.vstack((stims_x[non_nan_array_x, index], stims_y[non_nan_array_y, index])))
       # remove the first 300 ms because waiting time for trigger
       # points = np.delete(points, [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], axis=0)

        AOI_stim = [960+(list(OSdata['x_pos1'])[index])-200, 540-(list(OSdata['y_pos1'])[index]) +
                    200, 960+(list(OSdata['x_pos1'])[index])+200, 540-(list(OSdata['y_pos1'])[index])-200]
        radius = ((AOI_stim[2] - AOI_stim[0])*sqrt(2))/2
        patch = Circle((960+(list(OSdata['x_pos1'])[index]), 540-(list(OSdata['y_pos1'])[index])),
                          radius,
                          color='r', fill=False, linestyle = '-')
        
        #fix cross to get last sample in fixation cross before saccade begin
        AOI_fixcross=[960-110,540+110 , 960+110, 540-110]
        radius_fix = ((AOI_fixcross[2] - AOI_fixcross[0])*sqrt(2))/2
        patchfix = Circle((960, 540),
                          radius_fix,
                          color='b', fill=False, linestyle = '-')


        containsr = patch.contains_points(points)
        containsfix = patchfix.contains_points(points)
        gaze_aoi_newr = np.asarray(containsr, dtype=int)
        gaze_fixcross_newr = np.asarray(containsfix, dtype=int)
        
        gaze_fixcross_newr = list(gaze_fixcross_newr)
        index_fixcross = [i for i, e in enumerate(gaze_fixcross_newr) if e == 1]
        distance_target1 =(sqrt((((960+(OSdata['x_pos1'][index]))-960)**2)+(((540-(OSdata['y_pos1'][index]))-540)**2)))*deg_per_px
        distance_target.append(distance_target1)
                                                 
        gaze_aoi_newr = list(gaze_aoi_newr)
        index_aoi = [i for i, e in enumerate(gaze_aoi_newr) if e == 1]
        if len(index_aoi) != 0:
           index_fixcross = [i for i, e in enumerate(index_fixcross) if e < index_aoi[0]]
        else:
           index_fixcross = []
      #  msonAOI = 
        if not gaze_aoi_newr:
           percentageonAOI.extend([0.0]) 
        else:   
           percentageonAOI.extend([100*sum(gaze_aoi_newr)/len(gaze_aoi_newr)])
        if not gaze_fixcross_newr: 
           percentageonFixcross.extend([0.0]) 
        else:
           percentageonFixcross.extend([100*sum(gaze_fixcross_newr)/len(gaze_fixcross_newr)])
        
        stimul_name.append('stimulus_'+str(index))
        

        if any(iii == 1 for iii in gaze_aoi_newr) == True:
      # prima facciamo solo questo, se é 1 allora calcolo la velocitá e uso una threshold di 35
            visitedAoi.extend([1])
        else:
            visitedAoi.extend([0])
            
        if any(iii == 1 for iii in gaze_fixcross_newr) == True:
      # prima facciamo solo questo, se é 1 allora calcolo la velocitá e uso una threshold di 35
            visitedFixcross.extend([1])
        else:
            visitedFixcross.extend([0])
            
        trl_cnt = 0     
        for ii in range(len(points)-1):
           # visual angle calculation --> θ = atan (Size/Distance) Size = distanze of points, distance = distance from screen.
           # teorema di pitagora tra due points per trovare distanza tra di loro
           size_in_px = math.sqrt(((points[ii+1,0]-points[ii,0])**2)+((points[ii+1,1]-points[ii,1])**2))
           # 1pix is between 0.01 and 0.03 degrees. SACCADIC AMPLITUDE!
           size_in_deg = size_in_px*deg_per_px
           if visitedAoi[-1] == 1:
              angular_velocity = size_in_deg/dt
           else:
              angular_velocity = 0
              
           angular_velocity_all.append(angular_velocity)

           if angular_velocity> 35 and visitedAoi[-1] == 1 and visitedFixcross[-1] == 1 and len(index_fixcross)>0 and ii> index_fixcross[-1] and len(index_aoi)>1  and trl_cnt ==0:
              #saclat = np.where(points == points[ii,0])
              #saclat = saclat[0][0]
              index_saclat.append(ii)  # per each trial
               # amplitude of saccade = gain (>1 = saccade was too large, <1 = saccade too small)
                   # wrong calculation, it should be when point[0, to when saccade ends]
              saclat_ms = ii*(1000/fs)
              saclat_ms1.append(saclat_ms)
              D = np.where(containsr == True)
              m = 0 
              if D[0].size > 1:
                 m = 0 
                 for j in range(D[0].size-1):  
                     if D[0][j] > ii and D[0][j+1] == D[0][j]+1 and m==0:
                            saccade_end_index.append(D[0][j])                        
                            m +=1     
                     elif j>3 and m==0:
                            saccade_end_index.append(0)
                            m +=1
              else:
                saccade_end_index.append(0)
                m +=1
                     
              amp = sqrt(((points[saccade_end_index[-1],0]-points[index_saclat[-1],0])**2)+((points[saccade_end_index[-1],1]-points[index_saclat[-1],1])**2))
              amp= amp*deg_per_px
              mean_velocity_saccade.append(amp/dt)
              amplitude_saccade.append(amp)
              
              
              trl_cnt += 1    
              
        if trl_cnt ==0:
           amplitude_saccade.append(0)
           saccade_end_index.append(0)
           saclat_ms1.append(0)
           mean_velocity_saccade.append(0)
           index_saclat.append(0)
           
      #RT #its eye response time, the button press is different because we start from after fix cross appears
        rtime_et.append(len(points)*(1000/fs) )  
##################
# calculate velocity and amplitude of saccades only for the trials when they look at the target stimulus

  
    
    rtime_real = OSdata['responseTime1']
    for i in range(len(rtime_et)):
        if rtime_et[i] == 0.0:  #when no button press is required. do it so no 0 is left and SL and VRT can also be evaluated of these trials
           rtime_et[i] = 2166.7
    for i in range(len(rtime_real)):
         if rtime_real[i] == 0.0:  #when no button press is required. do it so no 0 is left and SL and VRT can also be evaluated of these trials
            rtime_real[i] = 2166.7       
    
   # rtime_et =[x for x in rtime_et if x > 300]
   # rtime_real =[x for x in rtime_real if x > 300]

    zip_object = zip(saccade_end_index, index_saclat)
    for s, r in zip_object:
        visual_rs.append(s-r)
        
    for i in range(len(visual_rs)):
        visual_rs_ms.append(visual_rs[i]*(1000/fs))           


  #  for i in range(len(amplitude_saccade)):
   #     gain.append(amplitude_saccade[i]/distance_target[i]) #Gains of <1 indicate the saccade was too small or hypometric; gains of >1 indicate the saccade was too large or hypermetric. 
    
              
#PS         
    for i in range(len(visual_rs_ms)):
        A = visual_rs_ms[i]+saclat_ms1[i]
        sum_vrssaclat.append(A)

#calculate processing speed

    for i in range(len(rtime_et)):
        if rtime_et[i] != 2166.7:
           processing_speed.append(rtime_et[i]-sum_vrssaclat[i])
        else: 
           processing_speed.append(0)
       
    for i in range(len(rtime_et)):
       if visual_rs_ms[i] == 0:
          processing_speed[i] = 0
    processing_speed= [abs(x) for x in processing_speed]
    
    for i in range(len(rtime_et)):
        percprocspeed.append((processing_speed[i]*100)/rtime_et[i])
        percvisualrspeed.append((visual_rs_ms[i]*100)/rtime_et[i])
        percentage_saclatency.append((saclat_ms1[i]*100)/rtime_et[i]) # normal values around 25%
        
 #   gain = list(np.around(np.array(gain),1))
    mean_velocity_saccade= list(np.around(np.array(mean_velocity_saccade),2))
    amplitude_saccade= list(np.around(np.array( amplitude_saccade),2))
    saclat_ms1= list(np.around(np.array(saclat_ms1),1))
    visual_rs_ms= list(np.around(np.array(visual_rs_ms),1))
    processing_speed= list(np.around(np.array(processing_speed),1))
    rtime_et= list(np.around(np.array(rtime_et),1))
    rtime_real= list(np.around(np.array(rtime_real),1))
    percprocspeed= list(np.around(np.array(percprocspeed),1))
    percvisualrspeed= list(np.around(np.array(percvisualrspeed),1))
    percentage_saclatency= list(np.around(np.array(percentage_saclatency),1))
# plot percentage of time spent on AOI
    tot_RT = [100] * len(percentage_saclatency)


    diff_1 = []
    zip_object = zip(tot_RT,percentage_saclatency)
    for s, r in zip_object:
        diff_1.append(s-r)
    
    diff_2 = []
    zip_object = zip(tot_RT,percvisualrspeed)
    for s, r in zip_object:
        diff_2.append(s-r)
    
    diff_3 = []
    zip_object = zip(tot_RT,percprocspeed)
    for s, r in zip_object:
        diff_3.append(s-r)

  #  first_saccade['gain']= gain    
    first_saccade['percentage of saccadic latency'] = percentage_saclatency
    first_saccade['saccadic latency'] = saclat_ms1
    first_saccade['saccade onset_index'] = index_saclat
    first_saccade['saccadic amplitude'] = amplitude_saccade
    first_saccade['visual reaction time'] = visual_rs_ms
    first_saccade['processing speed'] = processing_speed
    first_saccade['percentage of visual reaction time'] = percvisualrspeed
    first_saccade['percentage of processing speed'] =  percprocspeed
    first_saccade['RT ET'] = rtime_et
    first_saccade['Total RT'] = rtime_real
    first_saccade['mean velocity of saccade'] = mean_velocity_saccade
    first_saccade['visited AOI'] = visitedAoi
    
    

    if Plot == False:
   ###plot percentage saccadic latency      
        r = range(len(rtime_et))
        barWidth = 1
# plot bars
        plt.figure(figsize=(10, 7))
        plt.bar(r, percentage_saclatency, color='black',
                edgecolor='white', width=barWidth, label="saccadiclatency")
        plt.bar(r, diff_1, bottom=np.array(percentage_saclatency), color='darkgrey',
                edgecolor='white', width=barWidth, label='totalRT')

        plt.title('SL detection, 80cm')
# Custom X axis
        plt.xticks(r, percentage_saclatency, fontweight='bold')
        plt.ylabel("RT (%)")
        plt.xlabel("saccadic latency  (%)")
        plt.xticks([])
# plt.savefig("stacked1.png")
        
      #  plt.savefig(r'C:\Users\vb\Documents\BCIProject\Eyetracking\EyeTrackingAnalysis\images_paper\SL80cmKemp.png', format='png', dpi=1000)
        plt.show()
        
         ###plot vrs     
        r = range(len(rtime_et))
        barWidth = 1
# plot bars
        plt.figure(figsize=(10, 7))
        plt.bar(r, percvisualrspeed, color='black',
               edgecolor='white', width=barWidth, label="visual reaction speed")
        plt.bar(r, diff_2, bottom=np.array(percvisualrspeed), color='darkgrey',
               edgecolor='white', width=barWidth, label='totalRT')

        #plt.legend()
       # plt.get_legend().remove()
  #   plt.legend(bbox_to_anchor=(0, -0.15), loc=2, prop={'size': 8}, frameon=False)
        plt.title('VRT detection, 80cm')
       # 
# Custom X axis
        plt.xticks(r, percvisualrspeed, fontweight='bold')
        plt.ylabel("RT (%)")
        plt.xlabel("VRT  (%)")
        plt.xticks([])
# plt.savefig("stacked1.png")
        
       # plt.savefig(r'C:\Users\vb\Documents\BCIProject\Eyetracking\EyeTrackingAnalysis\images_paper\VRT80cmKemp.png', format='png', dpi=1000)
        #plt.show()
        
        ###plot  processing speed    
        r = range(len(rtime_et))
        barWidth = 1
# plot bars
        plt.figure(figsize=(10, 7))
        plt.bar(r, percprocspeed, color='red',
               edgecolor='white', width=barWidth, label="processing speed")
        plt.bar(r, diff_3, bottom=np.array(percprocspeed), color='orange',
               edgecolor='white', width=barWidth, label='totalRT')

        plt.legend()
# Custom X axis
        plt.xticks(r, percprocspeed, fontweight='bold')
        plt.ylabel("RT (%)")
        plt.xlabel("processing speed  (%)")
# plt.savefig("stacked1.png")
        plt.show()
        
    if Plot == True:
#stacked bar
        df1 = pd.DataFrame(columns=['SL (%)', 'VRS (%)', 'PS (%)'])
        df1['SL (%)'] = percentage_saclatency
        df1['VRS (%)'] = percvisualrspeed
        df1['PS (%)'] = percprocspeed
        ax = df1.plot.barh(stacked = True, xlabel='Number of valid trials', color={"SL (%)": "grey", "VRS (%)": "brown", "PS (%)": "orange"});
        for p in ax.patches:
             left, bottom, width, height = p.get_bbox().bounds
             ax.annotate(str(int(width)), xy=(left+width/2, bottom+height/2), 
                ha='center', va='center')

        plt.legend(bbox_to_anchor=(0, -0.15), loc=2, prop={'size': 8}, frameon=False)  
   
    
    return first_saccade

def ValuesToExcel_paper5(filename, OSFile,TobiiFolder, subject_nr):
    '''

    Saves all param in one excel file in the same folder as where the tobiifile is located 
    Both Pytrack and Pixels data are saved
    
    Parameters
    ----------
    TobiiFolder : Give the folder in which all the converted .tsv files are for one subject (output of divide data)
    OSFile :  Give the openSesame output file for this subject
    subject_nr : Nr of subject. Thus will be used to indicate the sheet number


    '''
    #In the filearrays, for each trial the values for this particular measure will be extended
    AOI=[960-200,540+200 , 960+200, 540-200]
    start = '19'
    stop = '22'
    correctresponse = []
  #  responsegiven = [] 
    psyes = []
    subj_nr = []
    

                
    OSdata = pd.read_csv(OSFile, delimiter = ',')           

    RTsubcomponents_VT = find_subcomponentsRT_VT(filename, OSFile,OSdata, start, stop, AOI, True)
    
   # U =  list(RTsubcomponents['stimulus num'])
    A= list(RTsubcomponents_VT['Total RT'])
    B= list(RTsubcomponents_VT['RT ET'])
    C = list(RTsubcomponents_VT['saccadic latency'])
    D = list(RTsubcomponents_VT['saccadic amplitude'])
    E =list(RTsubcomponents_VT['visual reaction time'])
    F = list(RTsubcomponents_VT['processing speed'])
    G = list(RTsubcomponents_VT['mean velocity of saccade'])
    H = list(RTsubcomponents_VT['percentage of saccadic latency'])
    I = list(RTsubcomponents_VT['percentage of visual reaction time'])  
    L = list(RTsubcomponents_VT['percentage of processing speed'])
   
    for i in range(len(RTsubcomponents_VT['visited AOI'])):
        if RTsubcomponents_VT['visited AOI'][i] ==0:
           A[i] = 0
           B[i] = 0 
           C[i] = 0 
           D[i] = 0 
           E[i] = 0 
           F[i] = 0 
           G[i] = 0 
           H[i] = 0 
           I[i] = 0 
           L[i] = 0 
          
        #   U.insert(i,0)
          
           
           
    #correct response is given
    for i in range(len(OSdata)): 
        if OSdata['correct_response'][i] == OSdata['responseKey1'][i]: #use 
            correctresponse.append(1)
        elif OSdata['correct_response'][i] != OSdata['responseKey1'][i]:
            correctresponse.append(0)
            
    ## put if response was supposed to be given for PS   
    for i in range(len(OSdata)): 
        if OSdata['correct_response'][i] == 'None':
            psyes.append(0)
        else: 
            psyes.append(1)


    ## add subj number        
    for i in range(len(OSdata)):
        subj_nr.append(subject_nr)
   
    allvar = list(zip(subj_nr, RTsubcomponents_VT['visited AOI'], A,B,C, E,F, D, G, H, I, L, correctresponse, OSdata['accuracy_total'],  OSdata['target_coordinates']))
    outData = pd.DataFrame(allvar, columns = ['Subject Number', 'Visited Target (yes/no)',
                                              'Total Reaction Time (ms)','RT ET','SL (ms) IVT','VRT (ms) IVT','PS (ms) IVT','Saccadic Amplitude IVT','Mean Velocity saccade IVT','SL (%) IVT', 'PS (%) IVT', 'VRT (%) IVT', 'Correct Response given','Total accuracy (%)','Target Coordinates'])        
            
    
    outData.to_excel(TobiiFolder + 'Output_patient'+ str(subject_nr) +".xlsx", sheet_name=str(subject_nr), index = False)
