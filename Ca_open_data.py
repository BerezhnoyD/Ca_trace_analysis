# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:07:29 2023

@author: Daniil.Berezhnoi


These are the functions to open and normalize the data
from the Inscopix Data Analysis Software

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import signal, ndimage
#import caiman
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.signal import find_peaks
from scipy import interpolate
import imageio
import random


def open_traces(directory, file):
    
    #Open the IDAS DataFrame and save cell traces in numpy format
    
    Ca_df=pd.read_csv(directory+file, header=1, sep=', ', dtype=np.float32, engine='python')
    Ca_cells=pd.read_csv(directory+file[:-4]+'-props.csv', header=0, sep=',', engine='python')
    image = imageio.imread(directory + file[:-4]+'_map.tif')

    Ca_cells = Ca_cells.loc[Ca_cells['Status']=='accepted']
    Ca_cells.drop(columns=['Name','Status','NumComponents','Size','ActiveSegment0'], inplace=True)
    Ca_cells.reset_index(inplace=True)
    
    
    col_dr = []
    for column in Ca_df.columns.unique():
        if 'rejected' in column:
            col_dr.append(column)
        elif 'undecided' in column:
            col_dr.append(column)
            
    Ca_df.drop(columns = col_dr, inplace = True)
    Ca_df.drop(columns='Time(s)/Cell Status', inplace=True)
    Ca = Ca_df.transpose().to_numpy()
    
    return Ca_df, Ca_cells, Ca, image

def split_cellset(cells: np.array, props: pd.DataFrame):
    n = props.loc[props['Red']>1].index.shape[0]
    list1 = props.loc[props['Red']>1].index
    list2 = props.loc[props['Red']==0].index
    list3 = random.choices(list2, k=n)
    
    cell_set1 = cells[list1]
    cell_set2 = cells[list2]
    cell_set3 = cells[list3]
    
    return cell_set1, cell_set2, cell_set3


def make_3darray(cells: np.array, reaches: pd.DataFrame, window=80):
    
    #Make the tensor-3dimensional array with Ncells x Timepoints x Trials 
    # The returned 'tensor' can be adressed by any dimension to get either one trial or one cells
    # Example: cell_tensor[:,:,0] gives you all cell traces in the first trial
    # cell_tensor[0,:,:] gives you traces in all trials for the first cell
    
    for k in enumerate(reaches['insc_peak']):
        if k[0] == 0:
            cell_tensor = np.expand_dims(cells[:,k[1]-window:k[1]+window], axis=2)
        else:
            add_cell_tensor = np.expand_dims(cells[:,k[1]-window:k[1]+window], axis=2)
            cell_tensor = np.append(cell_tensor, add_cell_tensor, axis=2)
    
    return cell_tensor

def make_3darray_category(cells: np.array, reaches: pd.DataFrame, window=80, category='Missed'):
    
    #Make the tensor-3dimensional array with Ncells x Timepoints x Trials 
    # The returned 'tensor' can be adressed by any dimension to get either one trial or one cells
    # Example: cell_tensor[:,:,0] gives you all cell traces in the first trial
    # cell_tensor[0,:,:] gives you traces in all trials for the first cell
    
    for k in enumerate(reaches.loc[reaches['category']==category, 'insc_peak']):
        if k[0] == 0:
            cell_tensor = np.expand_dims(cells[:,k[1]-window:k[1]+window], axis=2)
        else:
            add_cell_tensor = np.expand_dims(cells[:,k[1]-window:k[1]+window], axis=2)
            cell_tensor = np.append(cell_tensor, add_cell_tensor, axis=2)
    
    return cell_tensor


def z_score_data(cells: np.array):
    #make a z-scored copy of data
    
    cells_z = stats.zscore(np.nan_to_num(cells, nan=0.0),axis=1)
    fig = plt.figure()
    im_z=plt.imshow(cells_z, aspect='auto', vmin=0, vmax=cells_z.max())
    plt.title('Ca-signal Z-scored')
    plt.colorbar(im_z)
    
    return cells_z


def max_score_data(cells: np.array):    
    #make a max-scored copy of data
    
    row_max = np.nan_to_num(cells).max(axis=1)
    cells_max = np.nan_to_num(cells) / row_max[:, np.newaxis]
    fig = plt.figure()
    im_m=plt.imshow(cells_max, aspect='auto', vmin=0, vmax=1)
    plt.title('Ca-signal Max-scored')
    plt.colorbar(im_m)
    
    return cells_max


def filter_data(cells: np.array, length1=199, length2=5):
    
    #make a cleared copy of data with removed median trace and smoothing
    
    cells = np.nan_to_num(cells)
    cells_med = signal.medfilt2d(cells, kernel_size=[1,length1])
    cells_subst = cells - cells_med
    cells_filt = ndimage.uniform_filter(cells_subst, size=[1,length2])
    
    a= plt.figure()
    plt.plot(cells[0,:],'b-')
    plt.plot(cells_subst[0,:],'g-.')
    plt.plot(cells_filt[0,:],'r-')
    plt.legend(['raw signal','median substracted','averaged signal'])
    
    return cells_filt

def open_accel(directory, file):
    # Open the data from accelerometer and save changes in orientation and acceleration
    
    table=pd.read_csv(directory+file)
    table['Acc']=table[' Acc x'].abs()+table[' Acc y'].abs()+table[' Acc z'].abs()
    table['Ori']=table[' Ori pitch'].abs()+table[' Ori roll'].abs()+table[' Ori yaw'].abs()
    table['time']=table['IMU Time (s)'].diff()
    table1=table.iloc[::5,:]
    table1=table1.reset_index()
    accel=table1['Acc'].to_numpy()
    ori=table1['Ori'].to_numpy()
    table.Acc.plot()
    print('Accelerometer stops at {} sec'.format(table1['IMU Time (s)'].max()))
    
    return accel, ori

def stretch1d(to_stretch, dimension):
    #stretch one 1d array to match certain dimension
    #iterpolation is used to fill non-present dots
    
    f=interpolate.interp1d(np.arange(0,len(to_stretch)), to_stretch)
    stretch=f(np.linspace(0, len(to_stretch)-1, dimension))
    
    return stretch

def open_gpio(directory, file):
    
    #Open gpio file for Ca recording and get the events for each channel in separate tables
        
    table1 = pd.read_csv(directory + file)
    table1.rename(columns = {'Time (s)':'time', ' Channel Name':'channel',' Value':'value'}, inplace = True)
    
    Flir = table1[table1['channel'] == ' GPIO-1'].reset_index()                   #this is for FLIR synchro signal
    Inscopix = table1[table1['channel'] == ' BNC Sync Output'].reset_index()      #this is for Inscopix synchro signal
    Trigger = table1[table1['channel'] == ' BNC Trigger Input'].reset_index()   #this is for the Trigger from the Arduino

    # if you need to truncate the recording    
    #value = 20
    #Flir = Flir.truncate(before=(Flir['time']-value).abs().argsort()[:1].values[0])                   #this is for FLIR synchro signal
    #Inscopix = Inscopix.truncate(before=(Inscopix['time']-value).abs().argsort()[:1].values[0])     #this is for Inscopix synchro signal
    #Trigger = Trigger.truncate(before=(Trigger['time']-value).abs().argsort()[:1].values[0])     #this is for the Trigger from the Arduino
    
    
    #Making the tables for the Frames in both videos and summary
    Flir.drop(['channel'], axis=1, inplace=True)
    to_map = Flir['value'] <= 10000
    Flir['value'].mask(to_map, other=0, inplace=True)
    Flir['value'].mask(~to_map, other=30000, inplace=True)
    
    Flir = Flir.loc[Flir['value'].diff() < -10000]    #find the big downward leap in the signal as this is the end of the frame capture
    Flir_Frames = pd.DataFrame({'time': Flir.iloc[0][0], 'value': Flir.iloc[0][1]}, index=[1])
    Flir_Frames = pd.concat([Flir_Frames, Flir.loc[Flir['time'].diff() > 0.001]])
    Flir_Frames['value'] = 1
    Flir_Frames.drop(index=Flir_Frames.iloc[0,:].name, inplace=True)
    Flir_Frames.index = np.arange(1, (len(Flir_Frames)+1))
    
    Inscopix.reset_index(drop=True, inplace=True)
    Inscopix.drop(['channel'], axis=1, inplace=True)
    
    Trigger.reset_index(drop=True, inplace=True)
    Trigger.drop(['channel'], axis=1, inplace=True)
    Trigger.drop(Trigger.index[[0, -1]], inplace=True)
    Trigger_Frames = Trigger.query('value == True')
    Trigger_Frames.index = np.arange(1, (len(Trigger_Frames)+1))
    TriggerOff_Frames = Trigger.query('value == False')
    TriggerOff_Frames.index = np.arange(1, (len(TriggerOff_Frames)+1))
    
    Insc_Frames = Inscopix.query('value == True')
    Insc_Frames.index = np.arange(1, (len(Insc_Frames)+1))
    
    Summary = pd.DataFrame({'Inscopix': Insc_Frames.shape[0], 'Flir': Flir_Frames.shape[0], 'Triggers': Trigger_Frames.shape[0]}, index=[0])
    print(Summary)
    
    return Flir_Frames, Insc_Frames,Trigger_Frames, TriggerOff_Frames
    

    