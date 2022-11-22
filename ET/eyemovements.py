# -*- coding: utf-8 -*-
"""

Preprocessing of Eye tracking data, 
extrapolation of eye tracking variables
Plot of eye movements per trial for each subject


Created on Tue Nov 22 2022

@author: VB
"""

from PyTrack.formatBridge import generateCompatibleFormat
from PyTrack.Stimulus import Stimulus
import pandas as pd
import numpy as np
import math
from math import atan2, degrees, sqrt
import os
import time
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

# PREPROCESSING 
def getRightStructure_pat(TobiiFolder, subject_nr):
    


    '''


    Parameters
    ----------
    TobiiFile : folder of tsv files from tobii without header
    subject_nr: #subject number

    Returns
    -------
    raw_data : Returns the raw data in the correct format as .tsv file
    
    saves preprocessed data in the same directory as tsv

    '''
##remember to always remove the 'extra text' from the excel file (keep headers!), plus decimal separator = . and thousands separator = ,
    for filename1 in os.listdir(TobiiFolder):
        if filename1[4] == str(subject_nr) or filename1[4:6] == str(subject_nr):
           filename = TobiiFolder + filename1

    raw_data = pd.read_csv(filename, delimiter='\t', header=0)
    raw_data.columns = ['Recording timestamp', 'Event Value', 'Gaze2d_Left.x', 'Gaze2d_Left.y', 'ValL', 'Gaze2d_Right.x', 'Gaze2d_Right.y', 'ValR', 'GazeX', 'GazeY', 'PupilDiam_Left', 'PupValL', 'PupilDiam_Right', 'PupValR']
#Pytrack expects integers of time in milliseconds. timestamps from tobii are in microseconds
    raw_data['Recording timestamp'] = raw_data['Recording timestamp'].astype('float64')*1000
#raw_data['Recording timestamp'] = raw_data['Recording timestamp'].astype('float64')
    print('Converting data to right format...')
#Pygaze expects that the markers are numbers, so we convert strings to numbers first.
#Check if the message is Start
    a = raw_data['Event Value'].str.contains('Start experiment2 trial')
    for index, row in enumerate(a):
       if row == True:
          raw_data.at[index, 'Event Value'] = 10
    
    a = raw_data['Event Value'].str.contains('Start trial')
    for index, row in enumerate(a):
       if row == True:
          raw_data.at[index, 'Event Value'] = 10
       
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
           raw_data.at[index, 'Event Value'] = 19
           
    a = raw_data['Event Value'].str.contains('CueIN')
    for index, row in enumerate(a):
        if row == True:
           raw_data.at[index, 'Event Value'] = 17
            

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

    jj = 1
    print('Interpolating left pupil...')
    while jj < len(df)-3:
        if df.loc[jj, 'PupilDiam_Right'] == -1:
            if sum(df.loc[jj:jj+4, 'PupilDiam_Right']) ==-4:
                print('Too large to interpolate')
                while df.loc[jj, 'PupilDiam_Right'] == -1 and jj < len(df)-3:
                    jj +=1

            else:
                if df.loc[jj, 'PupilDiam_Right'] ==-1 and df.loc[jj+1, 'PupilDiam_Right'] == -1 and df.loc[jj+2, 'PupilDiam_Right'] == -1:
                    value_before_x = df.loc[jj-1, 'PupilDiam_Right']
                    value_after_x = df.loc[jj+3, 'PupilDiam_Right']
                    steps = 4

                    df.loc[jj, 'PupilDiam_Right'] = (((1/steps)*value_before_x)+((3/steps)*value_after_x))
                    df.loc[jj+1,'PupilDiam_Right'] = (((2/steps)*value_before_x)+((2/steps)*value_after_x))
                    df.loc[jj+2, 'PupilDiam_Right'] = (((3/steps)*value_before_x)+((1/steps)*value_after_x))
                elif df.loc[jj, 'PupilDiam_Right'] ==-1 and df.loc[jj+1, 'PupilDiam_Right'] == -1:
                    value_before_x = df.loc[jj-1, 'PupilDiam_Right']
                    value_after_x = df.loc[jj+2, 'PupilDiam_Right']
                    steps = 3
                    df.loc[jj,'PupilDiam_Right'] = ((2/steps)*value_before_x+((1/steps)*value_after_x))

                    df.loc[jj+1,'PupilDiam_Right'] = (((1/steps)*value_before_x)+((2/steps)*value_after_x))

                elif df.loc[jj, 'PupilDiam_Right'] ==-1:
                    value_before_x = df.loc[jj-1, 'PupilDiam_Right']
                    value_after_x = df.loc[jj+1, 'PupilDiam_Right']

                    steps = 2
                    df.loc[jj,'PupilDiam_Right'] = (value_before_x+value_after_x)/steps

        jj +=1

    jj = 1
    print('Interpolating right pupil...')
    while jj < len(df)-3:
        if df.loc[jj, 'PupilDiam_Right'] == -1:
            if sum(df.loc[jj:jj+4, 'PupilDiam_Right']) ==-4:
                print('Too large to interpolate')
                while df.loc[jj, 'PupilDiam_Right'] == -1 and jj < len(df)-3:
                    jj +=1

            else:
                if df.loc[jj, 'PupilDiam_Right'] ==-1 and df.loc[jj+1, 'PupilDiam_Right'] == -1 and df.loc[jj+2, 'PupilDiam_Right'] == -1:
                    value_before_x = df.loc[jj-1, 'PupilDiam_Right']
                    value_after_x = df.loc[jj+3, 'PupilDiam_Right']
                    steps = 4

                    df.loc[jj, 'PupilDiam_Right'] = (((1/steps)*value_before_x)+((3/steps)*value_after_x))
                    df.loc[jj+1,'PupilDiam_Right'] = (((2/steps)*value_before_x)+((2/steps)*value_after_x))
                    df.loc[jj+2, 'PupilDiam_Right'] = (((3/steps)*value_before_x)+((1/steps)*value_after_x))
                elif df.loc[jj, 'PupilDiam_Right'] ==-1 and df.loc[jj+1, 'PupilDiam_Right'] == -1:
                    value_before_x = df.loc[jj-1, 'PupilDiam_Right']
                    value_after_x = df.loc[jj+2, 'PupilDiam_Right']
                    steps = 3
                    df.loc[jj,'PupilDiam_Right'] = ((2/steps)*value_before_x+((1/steps)*value_after_x))

                    df.loc[jj+1,'PupilDiam_Right'] = (((1/steps)*value_before_x)+((2/steps)*value_after_x))

                elif df.loc[jj, 'PupilDiam_Right'] ==-1:
                    value_before_x = df.loc[jj-1, 'PupilDiam_Right']
                    value_after_x = df.loc[jj+1, 'PupilDiam_Right']

                    steps = 2
                    df.loc[jj,'PupilDiam_Right'] = (value_before_x+value_after_x)/steps

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
    
    
    
### EXTRACTS EYE PARAMETERS 
def ValuesToExcel_patients(filename, OSFile,TobiiFolder, subject_nr):
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
    #tot_fix_time = total_fixationtime(filename, OSFile, start, stop, AOI, True) ##pensa a come plottare, tipo box plot 
    RTsubcomponents = subcomponentsRTs_pixels(filename, OSFile, start, stop, AOI, False)
   # saccadiclatency_PT = find_saccadiclatency_pytrack(filename, OSFile, start, stop, AOI)
   # RTsubcomponents_PT = RTsubcomponents_pytrack(filename, OSFile, OSdata, start, stop, AOI)
    RTsubcomponents_VT = find_subcomponentsRT_VT(filename, OSFile,OSdata, start, stop, AOI, True)
    
    A = list(RTsubcomponents['saccadic latency'])
    B = list(RTsubcomponents['processing speed'])
    C = list(RTsubcomponents['visual reaction time'])
    D = list(RTsubcomponents['saccadic amplitude'])
    E = list(RTsubcomponents['percentage of saccadic latency']) 
    F = list(RTsubcomponents['percentage of processing speed'])
    G = list(RTsubcomponents['percentage of visual reaction time'])
    H = list(RTsubcomponents['total RT'])
    I = list(RTsubcomponents['velocity of saccade'])
   # L = list(RTsubcomponents_PT['VRT'])
   # N = list(RTsubcomponents_PT['saccade amplitude'])
   # O = list(RTsubcomponents_PT['peak velocity'])
    P = list(RTsubcomponents['time on fixcross(%)'])
    Q = list(RTsubcomponents['time on target(%)'])
    R = list(RTsubcomponents['total fixation duration on AOI target (ms)']) 
    S = list(RTsubcomponents['first fixation duration on AOI target (ms)']) 
    T=  list(RTsubcomponents['fixation duration on AOI fixcross (ms)'])
   # U =  list(RTsubcomponents['stimulus num'])
    M = list(RTsubcomponents_VT['saccadic latency'])
    V = list(RTsubcomponents_VT['saccadic amplitude'])
    Z =list(RTsubcomponents_VT['visual reaction time'])
    W = list(RTsubcomponents_VT['processing speed'])
    X = list(RTsubcomponents_VT['mean velocity of saccade'])
    T1 = list(RTsubcomponents_VT['percentage of saccadic latency'])
    T2 = list(RTsubcomponents_VT['percentage of visual reaction time'])  
    T3 = list(RTsubcomponents_VT['percentage of processing speed'])
    G1 = list(RTsubcomponents_VT['gain'])
    for i in range(len(RTsubcomponents['visited AOI'])):
        if RTsubcomponents['visited AOI'][i] ==0:
           A[i] = 0
           B[i] = 0 
           C[i] = 0 
           D[i] = 0 
           E[i] = 0 
           F[i] = 0 
           G[i] = 0 
           G1[i] = 0 
           H[i] = 0 
           I[i] = 0 
           M[i] = 0 
           P[i] = 0 
           Q[i] = 0 
           R[i] = 0 
           S[i] = 0 
           T[i] = 0 
           V[i] = 0 
           Z[i] = 0 
           W[i] = 0 
           X[i] = 0 
           T1[i] = 0 
           T2[i] = 0 
           T3[i] = 0 
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
   
    allvar = list(zip(subj_nr, G1,P,Q,R,S,T, RTsubcomponents['visited fix'], RTsubcomponents['visited AOI'], H, A, B, C, M,Z,W, D, V, E, F, G,T1,T3,T2,I,X, correctresponse, OSdata['accuracy_total'], psyes,  OSdata['target_coordinates']))
    outData = pd.DataFrame(allvar, columns = ['Subject Number','Gain IVT','Time on fixcross (%) IDT','Time on target (%) IDT', 'Total Fixation Duration (ms) IDT', 'First Fixation Duration (ms) IDT', 'Fixation Duration on fixcross (ms) IDT', 'Visited Fixation Cross (yes/no)','Visited Target (yes/no)',
                                              'Total Reaction Time (ms)','SL (ms) IDT','PS (ms) IDT','VRT (ms) IDT','SL (ms) IVT','VRT (ms) IVT', 'PS (ms) IVT','Saccadic Amplitude (pixels)','Saccadic Amplitude IVT', 'SL (%) IDT','PS (%) IDT','VRT (%) IDT','SL (%) IVT', 'PS (%) IVT', 'VRT (%) IVT','Mean Velocity (pixels/s)','Mean Velocity IVT', 'Correct Response given','Total accuracy (%)','PS can be calculated?','Target Coordinates'])        
            
    
    outData.to_excel(TobiiFolder + 'Output_patient'+ str(subject_nr) +".xlsx", sheet_name=str(subject_nr), index = False)



# PLOT EYE MOVEMENTS OF SINGLE TRIAL FOR SINGLE PATIENTS + EXTRACT MAIN SEQUENCE 
def visualize_eyemovement_trial(trialtoplot, subject_nr, SaveFig, MainSeq):

#plot
   start = 19 #why string before?
   stop = 22
   AOI=[960-200,540+200 , 960+200, 540-200] 
   TobiiFolder = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\OSfiles-Tobiifiles-SAGA\tobiinoheader\monkey\trials\\'
   OSFolder = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\OSfiles-Tobiifiles-SAGA\OSfilesmonkey\\'
#subject_nr = 3
   fs = 60
#trialtoplot =107#num of trial got from EEG. ones during absence + 1 before and 1 after
## fix duration
   OutputFolder = r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\OSfiles-Tobiifiles-SAGA\tobiinoheader\monkey\trials\output\\'
   for filename3 in os.listdir(OutputFolder):
       if filename3[-6] == str(subject_nr) or filename3[-7:-5] == str(subject_nr):
           output1 = OutputFolder + filename3
    
   out = pd.read_excel(output1)


   cross_percentage= list(out['Time on fixcross (%) IDT'])
   target_percentage= list(out['Time on target (%) IDT'])
   target_fixation = list(out['Total Fixation Duration (ms) IDT'])
   cross_fixation = list(out['Fixation Duration on fixcross (ms) IDT'])

   cross_percentage[trialtoplot]
   target_fixation[trialtoplot]
   cross_fixation[trialtoplot]
   target_percentage[trialtoplot]
   print( 'the fixation duration on target is: ' + str(target_fixation[trialtoplot])+ 'ms')

   target_fixation = [i for i in target_fixation if i != 0]  
   cross_fixation = [i for i in cross_fixation if i != 0]  
   print( 'the total average fixation duration is: ' + str(mean(target_fixation))+ 'ms')
   print( 'the std of the total fixation duration is: ' + str(stdev(target_fixation))+ 'ms')
 
#mean(cross_fixation)
#stdev(cross_fixation)

#num of trial from OS is +2 than from ET, coz first row is title, then it starts at 1 while python at 0

   for filename1 in os.listdir(TobiiFolder):
       if filename1[4] == str(subject_nr) or filename1[4:6] == str(subject_nr):
            filename = TobiiFolder + filename1
           
   for filename2 in os.listdir(OSFolder):  
       if filename2[5] == str(subject_nr) or filename2[4:6] == str(subject_nr):
           OSFile = OSFolder + filename2
      
   df1 = pd.read_csv(filename, sep= '\t')
   OSdata = pd.read_csv(OSFile, delimiter = ',')
   aa= df1['Event Value']
   indeces_start = aa[aa == start]
   indeces_end = aa[aa == stop]
   df = pd.DataFrame()


   df = df1.iloc[indeces_start.index[trialtoplot]:indeces_end.index[trialtoplot]+1,:] 
   y1 = df['GazeX']-(1920/2) 
   y2 = -(df['GazeY']-(1080/2))
   
   #######PLOT STATIC ALL
   hfont = {'fontname':'Times New Roman', 'size':14} 
   fig, ax = plt.subplots()
   line, = ax.plot(y1, y2, color='silver', marker='D', markeredgecolor='black', label='eye movement', linewidth = 3)
   AOI_stim = [(list(OSdata['x_pos1'])[trialtoplot])-180, (list(OSdata['y_pos1'])[trialtoplot])+180,(list(OSdata['x_pos1'])[trialtoplot])+200, (list(OSdata['y_pos1'])[trialtoplot])-200]
   radius = sqrt(abs(((AOI_stim[2] - AOI_stim[0])*(AOI_stim[3] - AOI_stim[1]))/math.pi))
   patch = Circle((AOI_stim[0]+180, AOI_stim[1]-180),
                             radius,
                             color='r', fill=False, linestyle = '-',label='target stimulus')

   #AOI_cue = [(list(OSdata['cue_xpos'])[trialtoplot])-80, (list(OSdata['cue_ypos'])[trialtoplot])+80,(list(OSdata['cue_xpos'])[trialtoplot])+80, (list(OSdata['cue_ypos'])[trialtoplot])-80]
   #radius = sqrt(abs(((AOI_cue[2] - AOI_stim[0])*(AOI_cue[3] - AOI_cue[1]))/math.pi))
   #patch1 = Circle((AOI_cue[0]+80, AOI_cue[1]-80),
                           #  radius,
                           #  color='g', fill=False, linestyle = '--', label='cue')

   AOI_fixcross=[960-100,540+100 , 960+100, 540-100]
   radius_fix = ((AOI_fixcross[2] - AOI_fixcross[0])*sqrt(2))/2
   patchfix = Circle((0, 0),
                     radius_fix,
                     color='steelblue', fill=False, linestyle = '--', label='fixation cross')
   


      #ax.set_title('Eye movements ET trial  '+str(trialtoplot+1), **hfont)
      #ax.set_ylabel('Y movements', **hfont)
      #ax.set_xlabel('X movements', **hfont)
   ax.set_ylim(-600,600)
   ax.set_xlim(-1100,1100)
      #ax.add_patch(patch)
      #ax.add_patch(patchfix)
   plt.xticks([], [])
   plt.yticks([], [])
   plt.axis('off')
     # plt.legend(prop={'family':'Times New Roman', 'size':16}, loc= 'lower right') 
     # plt.xticks(**hfont)
      #plt.yticks(**hfo


#determine duration of each trials to determine if trials with absences are above average (cannot check just using OS file because if button is not pressed then RT saved there is 0)
   length_trial = []
   for j in range(0,len(indeces_end)):
       df = df1.iloc[indeces_start.index[j]:indeces_end.index[j]+1,:] 
       y11 = df['GazeX']-(1920/2) 
       length_trial.append((len(y11)*(1/fs))*1000) #in ms 

# delete trials with None response 
   for ii in range(0,len(length_trial)): 
      if OSdata['correct_response'][ii] == 'None':
         length_trial[ii] =0 

# take average duration to see how much far away from trials with absences   
   length_trial[trialtoplot]    
   length_trial1 =[i for i in length_trial if i !=  0]
   
   print( 'the duration of this trial is: ' + str(length_trial[trialtoplot])+ 'ms')
   print( 'the total average duration of the trial is: ' + str(mean(length_trial1))+ 'ms')
   print( 'the total median duration of the trial is: ' + str(median(length_trial1))+ 'ms')
   print( 'the std of the total duration of the trial is: ' + str(stdev(length_trial1))+ 'ms')
   



############################################
# define jitter (i.e.  deviation of the significant instances of a signal from their ideal locations in time)
#plot it and check if it makes sense or not (if it is just normal saccade than no)
#it doesnt make sense so i am not plotting it anymore 
   
   yy1 = list(y1)
   yy2 = list(y2)
#saccades_group = [] 
#saccadex =[]
#saccadey = []
#for jj in range(2,len(y1)-2):
#    if (yy1[jj]>yy1[jj+1]+150 or yy1[jj]<yy1[jj+1]-150 or  yy1[jj]>yy1[jj-1]+150 or yy1[jj]<yy1[jj-1]-150) or (yy2[jj]>yy2[jj+1]+150 or yy2[jj]<yy2[jj+1]-150 or  yy2[jj]>yy2[jj-1]+150 or yy2[jj]<yy2[jj-1]-150) or yy1[jj]>yy1[jj+2]+250 or yy1[jj]<yy1[jj+2]-250 or  yy1[jj]>yy1[jj-2]+250 or yy1[jj]<yy1[jj-2]-250 or yy2[jj]>yy2[jj+2]+250 or yy2[jj]<yy2[jj+2]-250 or  yy2[jj]>yy2[jj-2]+250 or yy2[jj]<yy2[jj-2]-250: 
#       saccadex.append(yy1[jj])
#       saccadey.append(yy2[jj])
#       saccades_group.append(jj)

#line, = ax.plot(saccadex,saccadey, color='yellow', label='jitter', linestyle = 'dashdot')

   target_position = ['(-800, 475)', '(800, -475)','(800, 475)', '(-800, -475)','(800, -250)','(800, 250)','(-800, -250)','(-800, 250)','(600, 250)','(600, -250)','(-600, -250)','(-600, 250)']

   trials_pos0, trials_pos1,trials_pos2, trials_pos3, trials_pos4, trials_pos5, trials_pos6, trials_pos7, trials_pos8 = ([] for i in range(9))
   trials_pos9, trials_pos10,trials_pos11= ([] for i in range(3))
# define averaged ideal line 
   for f,coordinates in enumerate(OSdata['target_coordinates']):
       if coordinates == target_position[0]:
          trials_pos0.append(f)
       elif coordinates == target_position[1]:
          trials_pos1.append(f)
       elif coordinates == target_position[2]:
          trials_pos2.append(f)
       elif coordinates == target_position[3]:
          trials_pos3.append(f)
       elif coordinates == target_position[4]:
          trials_pos4.append(f)
       elif coordinates == target_position[5]:
          trials_pos5.append(f)
       elif coordinates == target_position[6]:
          trials_pos6.append(f)
       elif coordinates == target_position[7]:
          trials_pos7.append(f)
       elif coordinates == target_position[8]:
          trials_pos8.append(f)
       elif coordinates == target_position[9]:
          trials_pos9.append(f)
       elif coordinates == target_position[10]:
          trials_pos10.append(f)
       elif coordinates == target_position[11]:
          trials_pos11.append(f)

#average trials

#function to average lists with different length 
   def column_wise_sum(rows):
       columns = zip_longest(*rows, fillvalue=0)
       return [sum(col)/len(rows) for col in columns]  

#trials to exclude from average determined looking at each individual trial 
   if subject_nr ==3:
       trials_to_exclude = [1,2,3,6,7,8,9,19,21,24,25,27,30,31,36,37,38,39,41,43,44,45,46,48,50,57,58,59,62,64,69,71,73,77,78,79,81,84,85,87,88,93,94,101,103,104,105,111,112,115,118,119] #for absences + trials who do not look like ideal line towards target 
   elif subject_nr ==5:
       trials_to_exclude = [3,4,8,9,12,17,21,22,25,26,37,38,52,53,55,56,57,58,67,74,76,78,79,82,92,96,112,119,130,135]
   elif subject_nr ==7:
        trials_to_exclude = [1,2,3,4,5,6,8,9,11,12,13,14,15,16,18,19,20,21,22,23,24,25,28,29,30,31,32,34,36,37,38,41,42,43,44,45,46,48,49,51,52,53,54,55,57,58,59,60,61,62,64,66,67,68,69,70,71,72,73,74,75,76,77,78,80,82,83,84,85,86,87,88,89,90,91,94,95,97,98,99,100,101,102,103,104,105,106,108,110,111,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,131,132,133,135,137,138, 139,140,141,142 ]
   elif subject_nr ==8:
        trials_to_exclude = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,57,58,60,61,63,66,67,68,71,72,74,75,76,77,79,81,82,85,86,87,89,90,92,93,94,96,98,100,101,103,105,107,108]
   elif subject_nr ==10:
        trials_to_exclude = []
   elif subject_nr ==13:
        trials_to_exclude = []     
        
   
   cnt = 0
   average_x_pos0 = []
   for ff, elem in enumerate(trials_pos0): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos0 = list(df['GazeX']-(1920/2))
                 average_y_pos0 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos0,X_eyes1]
                 rows1 = [average_y_pos0,Y_eyes1]
                 average_x_pos0 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos0 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]
   cnt = 0
   average_x_pos1 = []
   for ff, elem in enumerate(trials_pos1): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos1 = list(df['GazeX']-(1920/2))
                 average_y_pos1 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos1,X_eyes1]
                 rows1 = [average_y_pos1,Y_eyes1]
                 average_x_pos1 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos1 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]
   cnt = 0
   average_x_pos2 = []
   for ff, elem in enumerate(trials_pos2): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos2 = list(df['GazeX']-(1920/2))
                 average_y_pos2 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos2,X_eyes1]
                 rows1 = [average_y_pos2,Y_eyes1]
                 average_x_pos2 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos2 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]
   cnt = 0 
   average_x_pos3 = []
   for ff, elem in enumerate(trials_pos3): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos3 = list(df['GazeX']-(1920/2))
                 average_y_pos3 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos3,X_eyes1]
                 rows1 = [average_y_pos3,Y_eyes1]
                 average_x_pos3 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos3 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)] 
   cnt = 0
   average_x_pos4 = []
   for ff, elem in enumerate(trials_pos4): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos4 = list(df['GazeX']-(1920/2))
                 average_y_pos4 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos4,X_eyes1]
                 rows1 = [average_y_pos4,Y_eyes1]
                 average_x_pos4 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos4 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]    
   cnt = 0
   average_x_pos5 = []
   for ff, elem in enumerate(trials_pos5): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos5 = list(df['GazeX']-(1920/2))
                 average_y_pos5 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos5,X_eyes1]
                 rows1 = [average_y_pos5,Y_eyes1]
                 average_x_pos5 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos5 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]          
   cnt = 0
   average_x_pos6 = []
   for ff, elem in enumerate(trials_pos6): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos6 = list(df['GazeX']-(1920/2))
                 average_y_pos6 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos6,X_eyes1]
                 rows1 = [average_y_pos6,Y_eyes1]
                 average_x_pos6 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos6 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]                 
   cnt = 0 
   average_x_pos7 = []
   for ff, elem in enumerate(trials_pos7): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos7 = list(df['GazeX']-(1920/2))
                 average_y_pos7 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos7,X_eyes1]
                 rows1 = [average_y_pos7,Y_eyes1]
                 average_x_pos7 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos7 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]
   cnt = 0 
   average_x_pos8 = []
   for ff, elem in enumerate(trials_pos8): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos8 = list(df['GazeX']-(1920/2))
                 average_y_pos8 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos8,X_eyes1]
                 rows1 = [average_y_pos8,Y_eyes1]
                 average_x_pos8 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos8 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]          
   cnt = 0 
   average_x_pos9 = []
   for ff, elem in enumerate(trials_pos9): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos9 = list(df['GazeX']-(1920/2))
                 average_y_pos9 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos9,X_eyes1]
                 rows1 = [average_y_pos9,Y_eyes1]
                 average_x_pos9 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos9 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]                 
   cnt = 0 
   average_x_pos10 = []
   for ff, elem in enumerate(trials_pos10): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos10 = list(df['GazeX']-(1920/2))
                 average_y_pos10 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos10,X_eyes1]
                 rows1 = [average_y_pos10,Y_eyes1]
                 average_x_pos10 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos10 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]         
   cnt = 0 
   average_x_pos11 = []
   for ff, elem in enumerate(trials_pos11): 
          if (elem in trials_to_exclude) == False: 
              if cnt == 0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:] 
                 average_x_pos11 = list(df['GazeX']-(1920/2))
                 average_y_pos11 = list(-(df['GazeY']-(1080/2)))
                 cnt = 1 
              elif cnt>0:
                 df = df1.iloc[indeces_start.index[elem]:indeces_end.index[elem]+1,:]  
                 X_eyes1 = list(df['GazeX']-(1920/2))
                 Y_eyes1 = list(-(df['GazeY']-(1080/2)))
                 rows = [average_x_pos11,X_eyes1]
                 rows1 = [average_y_pos11,Y_eyes1]
                 average_x_pos11 = [sum(column)/len(rows) for column in zip_longest(*rows, fillvalue=0)]
                 average_y_pos11 = [sum(column)/len(rows) for column in zip_longest(*rows1, fillvalue=0)]                 
                 
#######PLOT STATIC ALL
   hfont = {'fontname':'Times New Roman', 'size':14} 
   fig, ax = plt.subplots()

   line, = ax.plot(y1, y2, color='black', marker='p', markerfacecolor= 'lightgrey', markeredgecolor='black', label='eye movement', linewidth = 2.5)
   line, =  ax.plot(y1[0:1], y2[0:1], color='white',marker='>', markersize=15,markeredgecolor='black', label='start', linewidth = 1)
   line, =  ax.plot(y1[len(y1)-1:len(y1)], y2[len(y1)-1:len(y1)], color='white',marker='<', markersize=15,markeredgecolor='black', label='end', linewidth = 1)
   
  
   AOI_stim = [(list(OSdata['x_pos1'])[trialtoplot])-180, (list(OSdata['y_pos1'])[trialtoplot])+180,(list(OSdata['x_pos1'])[trialtoplot])+200, (list(OSdata['y_pos1'])[trialtoplot])-200]
   radius = sqrt(abs(((AOI_stim[2] - AOI_stim[0])*(AOI_stim[3] - AOI_stim[1]))/math.pi))
   patch = Circle((AOI_stim[0]+180, AOI_stim[1]-180),
                          radius,
                          color='r', fill=False, linestyle = '-',label='target stimulus', linewidth = 2)

#AOI_cue = [(list(OSdata['cue_xpos'])[trialtoplot])-80, (list(OSdata['cue_ypos'])[trialtoplot])+80,(list(OSdata['cue_xpos'])[trialtoplot])+80, (list(OSdata['cue_ypos'])[trialtoplot])-80]
#radius = sqrt(abs(((AOI_cue[2] - AOI_stim[0])*(AOI_cue[3] - AOI_cue[1]))/math.pi))
#patch1 = Circle((AOI_cue[0]+80, AOI_cue[1]-80),
                        #  radius,
                        #  color='g', fill=False, linestyle = '--', label='cue')

   AOI_fixcross=[960-100,540+100 , 960+100, 540-100]
   radius_fix = ((AOI_fixcross[2] - AOI_fixcross[0])*sqrt(2))/2
   patchfix = Circle((0, 0),
                  radius_fix,
                  color='lightcoral', fill=False, linestyle = '--', label='_None', linewidth = 1.7)

   #ax.set_title('Eye movements ET trial  '+str(trialtoplot+1), **hfont)
   #ax.set_ylabel('Y movements', **hfont)
   #ax.set_xlabel('X movements', **hfont)
   ax.set_ylim(-600,600)
   ax.set_xlim(-1100,1100)
   ax.add_patch(patch)
   ax.add_patch(patchfix)
   plt.xticks([], [])
   plt.yticks([], [])
   plt.scatter(0, 0, s=600, c='lightcoral', marker='+', clip_on=False, linewidth=5,label='fixation cross')
   
  
   plt.legend(prop={'family':'Times New Roman', 'size':14}, loc= 'lower right', frameon=False) 
   plt.show()
   
   #ANIMATION PER 
   def update(num, x, y, line):
       line.set_data(y1[:num], y2[:num])
       line.set_color('black')
       line.set_linewidth(2.5)
       return line,
   line2, = ax.plot(y1, y2)
   ani1 = animation.FuncAnimation(fig, update, len(y1), fargs=[y1, y2, line2], interval=30, blit=True)
  # plt.xticks(**hfont)
   #plt.yticks(**hfont)
   ani1.save(r'C:\Users\vb\Desktop\lala.gif', dpi=300)
   plt.savefig(r'C:\Users\vb\Desktop\V\MOTION\PhD\paper4-absencecharacteristics\plots\pat7_trial18.png', format = 'png', dpi=1000)
    

### PLOT ANYMATED IN 3 TIMES 
   hfont = {'fontname':'Times New Roman', 'size':14} 
   fig, ax = plt.subplots()
   hal = int(len(yy1)/3)
   tt = int((hal*(1/fs))*1000)
   line, = ax.plot(y1[0:hal], y2[0:hal], color='lightcoral', marker='D', markeredgecolor='black', label='t1 = '+ str(tt)+' ms', linewidth = 3)
   line, = ax.plot(y1[hal:hal*2], y2[hal:hal*2], color='blue',marker='D', markeredgecolor='black', label='t2 = '+ str(tt)+' ms', linewidth = 3)
   line, = ax.plot(y1[hal*2:hal*3], y2[hal*2:hal*3], color='#4b0082',marker='D', markeredgecolor='black', label='t3 = ' + str(tt)+' ms', linewidth = 3)

#create ideal trajectory line 
   x_coordinate_fix = 0
   y_coordinate_fix = 0
   x_coordinate_stim = (list(OSdata['x_pos1'])[trialtoplot])
   y_coordinate_stim =(list(OSdata['y_pos1'])[trialtoplot])
#ideal fake line
#line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')

#ideal averaged line 
   if OSdata['target_coordinates'][trialtoplot] == target_position[0]:
       if len(average_x_pos0)>0:
          line1, = ax.plot(average_x_pos0, average_y_pos0, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[1]:
       if len(average_x_pos1)>0:
          line1, = ax.plot(average_x_pos1, average_y_pos1, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[2]:
       if len(average_x_pos2)>0:
          line1, = ax.plot(average_x_pos2, average_y_pos2, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[3]:
       if len(average_x_pos3)>0:
          line1, = ax.plot(average_x_pos3, average_y_pos3, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')
 
   elif OSdata['target_coordinates'][trialtoplot] == target_position[4]:
       if len(average_x_pos4)>0:
          line1, = ax.plot(average_x_pos4, average_y_pos4, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[5]:
       if len(average_x_pos5)>0:
          line1, = ax.plot(average_x_pos5, average_y_pos5, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[6]:
       if len(average_x_pos6)>0:
          line1, = ax.plot(average_x_pos6, average_y_pos6, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[7]:
       if len(average_x_pos7)>0:
          line1, = ax.plot(average_x_pos7, average_y_pos7, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[8]:
       if len(average_x_pos8)>0:
          line1, = ax.plot(average_x_pos8, average_y_pos8, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[9]:
      if len(average_x_pos9)>0:
         line1, = ax.plot(average_x_pos9, average_y_pos9, color='limegreen',linewidth = 4,label= 'ideal line')
      else:
         line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[10]:
       if len(average_x_pos10)>0:
          line1, = ax.plot(average_x_pos10, average_y_pos10, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   
   elif OSdata['target_coordinates'][trialtoplot] == target_position[11]:
       if len(average_x_pos11)>0:
          line1, = ax.plot(average_x_pos11, average_y_pos11, color='limegreen',linewidth = 4,label= 'ideal line')
       else:
          line1, = ax.plot([x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim], color='limegreen',linewidth = 4,label= 'ideal line')   

   polygon = plt.fill(np.append(y1, [x_coordinate_fix, x_coordinate_stim][::-1]), np.append(y2, [y_coordinate_fix, y_coordinate_stim][::-1]),alpha=0.3)

   AOI_stim = [(list(OSdata['x_pos1'])[trialtoplot])-180, (list(OSdata['y_pos1'])[trialtoplot])+180,(list(OSdata['x_pos1'])[trialtoplot])+200, (list(OSdata['y_pos1'])[trialtoplot])-200]
   radius = sqrt(abs(((AOI_stim[2] - AOI_stim[0])*(AOI_stim[3] - AOI_stim[1]))/math.pi))
   patch = Circle((AOI_stim[0]+180, AOI_stim[1]-180),
                          radius,
                          color='r', fill=False, linestyle = '--')

#AOI_cue = [(list(OSdata['cue_xpos'])[trialtoplot])-80, (list(OSdata['cue_ypos'])[trialtoplot])+80,(list(OSdata['cue_xpos'])[trialtoplot])+80, (list(OSdata['cue_ypos'])[trialtoplot])-80]
#radius = sqrt(abs(((AOI_cue[2] - AOI_stim[0])*(AOI_cue[3] - AOI_cue[1]))/math.pi))
#patch1 = Circle((AOI_cue[0]+80, AOI_cue[1]-80),
                        #  radius,
                        #  color='g', fill=False, linestyle = '--', label='cue')

   AOI_fixcross=[960-100,540+100 , 960+100, 540-100]
   radius_fix = ((AOI_fixcross[2] - AOI_fixcross[0])*sqrt(2))/2
   patchfix = Circle((0, 0),
                  radius_fix,
                  color='m', fill=False, linestyle = '--')

   ax.set_title('Eye movements ET trial  '+str(trialtoplot+1), **hfont)
   ax.set_ylabel('Y movements', **hfont)
   ax.set_xlabel('X movements', **hfont)
   ax.set_ylim(-600,600)
   ax.set_xlim(-1100,1100)
   ax.add_patch(patch)
   ax.add_patch(patchfix)
   plt.legend(prop={'family':'Times New Roman', 'size':10}) 
   plt.xticks(**hfont)
   plt.yticks(**hfont)
   def update(num, x, y, line):
       line.set_data(y1[:num], y2[:num])
       line.set_color('w')
       line.set_linewidth(1.5)
       return line,
   line2, = ax.plot(y1, y2)
   ani1 = animation.FuncAnimation(fig, update, len(yy1), fargs=[y1, y2, line2], interval=30, blit=True)

   if SaveFig == True:
      ani1.save(r'C:\Users\vb\Desktop\V\MOTION\PhD\paper4-absencecharacteristics\absences_features\eyemove2.gif', dpi=300)
#plt.savefig(r'C:\Users\vb\Desktop\V\MOTION\PhD\paper4-absencecharacteristics\absences_features\eyemove2.png', format = 'png', dpi=1000)
 
   def computeArea(pos):
       x, y = (zip(*pos))
       return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

# pyplot.fill(a, b) will return a list of matplotlib.patches.Polygon.

# The area of the polygon can be computed as follows:
# (you could also sum the areas of all polygons in the list).
   routeeff = computeArea(polygon[0].xy)
   print('the route efficacy (area) is:'+ str(routeeff))

   plt.show()


#########dot prod
   myline = [[x_coordinate_fix, x_coordinate_stim], [y_coordinate_fix, y_coordinate_stim]]
   realine = [yy1, yy2]

   path_len = 0 
  # dotprod_realine = [] 
   for index in range(len(yy1)-1):
       a = [realine[0][index+1]-realine[0][index],realine[1][index+1]-realine[1][index]] # for each point take the difference between point thereafter and before (so just points)
       b = math.sqrt(a[0]**2+a[1]**2) #lenght of vector 
       path_len += b
      # dotprod_realine.append([a[0]/b, a[1]/b]) #  normalize each distance of vectors for the length of it so we do not have differences related to length/duration of trial because each trial has different numebr of points


   a1 = [myline[0][1]-myline[0][0],myline[1][1]-myline[1][0]]
   b1 = math.sqrt(a1[0]**2+a1[1]**2) 
   #dotprod_myline_norm = [a1[0]/b1, a1[1]/b1]
   path_len1 = path_len/b1  # sum of the lengths of all the vectors divided by the lenght of the ideal vector
   print('the route efficacy (length) is:'+ str(path_len1))

 

  # dot_prod = []
   #for i, j in enumerate(dotprod_realine):
    #    dot_prod.append(dotprod_myline_norm[0]*dotprod_realine[i][0]+dotprod_myline_norm[1]*dotprod_realine[i][1]) #angle between 
     
   #ef = mean(dot_prod) # average of all the angles between vectors #the closer to 1 the better (the two lines overlap more)

 ## cannot account of microsaccades going back and forth on the same directions, because they will just sum up (so u see negative values if look on opposite direction of line)

    

##################################################################################
#MAIN SEQUENCE
   if MainSeq == True:
      target_position = ['(-800, 475)', '(800, -475)','(800, 475)', '(-800, -475)','(600, 250)','(600, -250)','(-600, -250)','(-600, 250)']
            
      listac, listac1, listac2, listac3,listac4, listac5, listac6, listac7, listac8, listac9 = ([] for i in range(10))
      amp_averaged = []
      meanvel_averaged = []
      slIDT_averaged = []
      slIVT_averaged = []
      vrtIDT_averaged = []
      vrtIVT_averaged = []
      ampIVT_averaged = []
                  
      df = out
                
      SL1,SL2,SL3,SL4,SL5,SL6,SL7,SL8 = ([] for i in range(8))
      SLIVT1,SLIVT2,SLIVT3,SLIVT4,SLIVT5,SLIVT6,SLIVT7,SLIVT8 = ([] for i in range(8))
      VRT1, VRT2, VRT3, VRT4, VRT5, VRT6, VRT7, VRT8= ([] for i in range(8))
               
      VRTIVT1, VRTIVT2, VRTIVT3, VRTIVT4, VRTIVT5, VRTIVT6, VRTIVT7, VRTIVT8= ([] for i in range(8))
      amp1, amp2,amp3,amp4, amp5, amp6, amp7, amp8 = ([] for i in range(8))
      ampIVT1, ampIVT2,ampIVT3,ampIVT4, ampIVT5, ampIVT6, ampIVT7, ampIVT8 = ([] for i in range(8))              
      meanvel1, meanvel2, meanvel3, meanvel4,  meanvel5,  meanvel6,  meanvel7,  meanvel8 = ([] for i in range(8))
      A = list(df['Saccadic Amplitude (pixels)'])
      A2 = list(df['Saccadic Amplitude IVT'])
      B= list(df['VRT (ms) IDT'])
      B1= list(df['VRT (ms) IVT'])              
      C = list(df['Mean Velocity IVT'])           
      D = list(df['SL (ms) IDT'])
      D1 = list(df['SL (ms) IVT'])
                   
      listac.append(A)               
      listac2.append(A2)
      listac3.append(B)
      listac4.append(B1)              
      listac6.append(C)               
      listac8.append(D)
      listac9.append(D1)
           
  ### average same target position
           
      for j in range(len(df['Target Coordinates'])):
          if df['Target Coordinates'][j] == target_position[0]:
              ampIVT1.append(df['Saccadic Amplitude IVT'][j])              
              meanvel1.append(df['Mean Velocity IVT'][j])             
              VRTIVT1.append(df['VRT (ms) IVT'][j])
          elif df['Target Coordinates'][j] == target_position[1]:
              ampIVT2.append(df['Saccadic Amplitude IVT'][j])               
              meanvel2.append(df['Mean Velocity IVT'][j])            
              VRTIVT2.append(df['VRT (ms) IVT'][j])
          elif df['Target Coordinates'][j] == target_position[2]:
              ampIVT3.append(df['Saccadic Amplitude IVT'][j])          
              meanvel3.append(df['Mean Velocity IVT'][j])
              VRTIVT3.append(df['VRT (ms) IVT'][j])
          elif df['Target Coordinates'][j] == target_position[3]:
              ampIVT4.append(df['Saccadic Amplitude IVT'][j])      
              meanvel4.append(df['Mean Velocity IVT'][j])
              VRTIVT4.append(df['VRT (ms) IVT'][j])
          elif df['Target Coordinates'][j] == target_position[4]:
              ampIVT5.append(df['Saccadic Amplitude IVT'][j])               
              meanvel5.append(df['Mean Velocity IVT'][j])
              VRTIVT5.append(df['VRT (ms) IVT'][j])
          elif df['Target Coordinates'][j] == target_position[5]:
              ampIVT6.append(df['Saccadic Amplitude IVT'][j])              
              meanvel6.append(df['Mean Velocity IVT'][j])
              VRTIVT6.append(df['VRT (ms) IVT'][j])
          elif df['Target Coordinates'][j] == target_position[6]:
              ampIVT7.append(df['Saccadic Amplitude IVT'][j])             
              meanvel7.append(df['Mean Velocity IVT'][j])              
              VRTIVT7.append(df['VRT (ms) IVT'][j])
          elif df['Target Coordinates'][j] == target_position[7]:              
              ampIVT8.append(df['Saccadic Amplitude IVT'][j])              
              meanvel8.append(df['Mean Velocity IVT'][j])           
              VRTIVT8.append(df['VRT (ms) IVT'][j])
                

      ampIVT1 = mean(ampIVT1)              
      meanvel1 = mean(meanvel1)              
      VRTIVT1 = mean(VRTIVT1)
               
      ampIVT2 = mean(ampIVT2)              
      meanvel2 = mean(meanvel2)                
      VRTIVT2 = mean(VRTIVT2)
               
      ampIVT3 = mean(ampIVT3)               
      meanvel3 = mean(meanvel3)             
      VRTIVT3 = mean(VRTIVT3)
              
      ampIVT4 = mean(ampIVT4)         
      meanvel4 = mean(meanvel4)               
      VRTIVT4 = mean(VRTIVT4)  
              
      ampIVT5 = mean(ampIVT5)                
      meanvel5 = mean(meanvel5)      
      VRTIVT5 = mean(VRTIVT5)
               
      ampIVT6 = mean(ampIVT6)               
      meanvel6 = mean(meanvel6)               
      VRTIVT6 = mean(VRTIVT6)

      ampIVT7 = mean(ampIVT7)        
      meanvel7 = mean(meanvel7)       
      VRTIVT7= mean( VRTIVT7) 
              
      ampIVT8 = mean(ampIVT8)
      meanvel8 = mean(meanvel8)           
      VRTIVT8 = mean(VRTIVT8)
            
              
      vrtIVT_averaged.append([VRTIVT1,VRTIVT2,VRTIVT3,VRTIVT4,VRTIVT5,VRTIVT6,VRTIVT7, VRTIVT8])             
      ampIVT_averaged.append([ampIVT1,ampIVT2, ampIVT3,ampIVT4,ampIVT5,ampIVT6,ampIVT7,ampIVT8])               
      meanvel_averaged.append([meanvel1, meanvel2, meanvel3, meanvel4,  meanvel5,  meanvel6,  meanvel7,  meanvel8])
      
            
      A, A1, A2, B, B1,B2,C,C1,D,D1,F,F1,F2,G,G1,G2,H,I,L = ([] for i in range(19))
       # NUM OF LISTAC changes with the number of subjects included per category 
      for i, e in enumerate(listac):
          F1 += vrtIVT_averaged[i]             
          G2 +=ampIVT_averaged[i]             
          I += meanvel_averaged[i]
   
    
      amplitudeivt = G2         
      meanvelivt = I 
      durationivt = F1            
      amplitudeivtHEAL = G2 
      durationivtHEAL = F1
      meanvelivtHEAL = I 
            ##############save ps to plot all fit lines together
      Ff1 = []   
      for i in range(len(G2)):
          if G2[i] != 0:
             Ff1.append(F1[i])
      Gg2H = [i for i in G2 if i != 0]  
      p, res, _, _, _  =  np.polyfit(Gg2H, Ff1, 1, full=True)
      phealthy_durampivt = p
      mmhealthy_durampivt = [i * p[0] for i in Gg2H]
            #########
      Hh = []   
      for i in range(len(G1)):
          if G1[i] != 0:
             Hh.append(H[i])
                   
           
            #######
      Ggg2H = []   
      for i in range(len(I)):
          if I[i] != 0 and I[i] < 500:
             Ggg2H.append(G2[i]) 
      Ii = [i for i in I if i != 0 and i <500]
      p, res, _, _, _  =  np.polyfit(Ggg2H, Ii, 1, full=True)
      mmhealthy_meanvel = [i * p[0] for i in Ggg2H]
      phealthy_meanvel = p
            ######


###################### IVT  AVER DURATION (per target pos) vs AVER AMP (per target pos)     
      fig, ax = plt.subplots()
      hfont = {'fontname':'Times New Roman', 'size': '17'}
      durationivtHEAL1 = []   
      for i in range(len(amplitudeivtHEAL)):
          if amplitudeivtHEAL[i] != 0:
             durationivtHEAL1.append(durationivtHEAL[i])
      amplitudeivtHEAL1 = [i for i in amplitudeivtHEAL if i != 0]  
     
      durationivt1 = []   
      for i in range(len(amplitudeivt)):
         if amplitudeivt[i] != 0:
            durationivt1.append(durationivt[i])
             
      amplitudeivt1 = [i for i in amplitudeivt if i != 0]  
      plt.plot(amplitudeivtHEAL1, durationivtHEAL1, 'o', markersize=3, color = 'maroon', alpha=0.5, label ='_Hidden label')
     #plt.plot(amplitudeivt1,durationivt1, 'o', markersize=4, color = 'silver', label ='_Hidden label' )   
      plt.plot(Gg2H, mmhealthy_durampivt + phealthy_durampivt[1],  color = 'maroon', linewidth =2, linestyle='-')   
      plt.xlabel("Saccadic Amplitude (°)", **hfont)
     #plt.title("Average per target position of children with AE (N=3) IVT")
      plt.ylabel("Saccadic Duration (ms)", **hfont)
      plt.title("Main Sequence pat 8", **hfont)
    # plt.xticks([0, 2,4,6,8,10,12], ['', '2', '4', '6', '8', '10','12']) to remove zero from x axis
      plt.yticks(**hfont)
      plt.xticks(**hfont)
      plt.xlim(0,12)
 #    plt.ylim(0,300)
  #   plt.legend(prop={'family':'Times New Roman', 'size':7}, frameon=False)
      plt.legend(prop={'family':'Times New Roman', 'size':12}, frameon=False)
    # plt.legend()
      plt.tight_layout()
     #plt.savefig(r'C:\Users\vb\Documents\MainSeq_durvsamp.png', format='png', dpi=1000)
      plt.show()
     # pearson correlation coefficient and p-value
    # pearson_coef, p_value = stats.pearsonr(amplitudeivtTBI1, durationivtTBI1)

######################   PT AMP vs PEAK VEL    

##############################################################
      amplitudeivt2 = []   
      for i in range(len(meanvelivt)):
          if meanvelivt[i] != 0 and meanvelivt[i] < 500:
             amplitudeivt2.append(amplitudeivt[i])
                             
      meanvel_averaged1 = [i for i in meanvelivt if i != 0 and i <500]
           
      amplitudeivtHEAL2 = []   
      for i in range(len(meanvelivtHEAL)):
          if meanvelivtHEAL[i] != 0 and meanvelivtHEAL[i] < 500:
             amplitudeivtHEAL2.append(amplitudeivtHEAL[i])
             
             
      meanvel_averagedHEAL1 = [i for i in meanvelivtHEAL if i != 0 and i <500]
      fig, ax = plt.subplots()
      plt.plot(amplitudeivtHEAL2, meanvel_averagedHEAL1, 'o', markersize=3, color = 'maroon', alpha=0.5, label ='_Hidden label')      
      plt.plot(Ggg2H, mmhealthy_meanvel + phealthy_meanvel[1],  color = 'maroon', linewidth =2, linestyle='-')     
      plt.xlabel("Saccadic Amplitude (°)", **hfont)
     # plt.title('Average per trial of amplitude vs Mean Velocity in Epileptic children with AE (N=3) IVT')
      plt.ylabel('Saccade Mean Velocity (°/s)', **hfont)
      plt.yticks(**hfont)
      plt.xticks(**hfont)
      plt.xlim(0,12)
    # plt.ylim(0, 400)
    # plt.legend(prop={'family':'Times New Roman', 'size':7}, frameon=False)
      plt.tight_layout()
    # plt.savefig(r'C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\images_paper\MeanSeq_velvsamp.eps', format='eps', dpi=1000)
      plt.show()
    
    