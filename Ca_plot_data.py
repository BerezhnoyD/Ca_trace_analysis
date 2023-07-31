# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 17:36:46 2023

@author: Daniil.Berezhnoi

These are the functions to open and plot the data
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
from scipy.signal import find_peaks, correlate
import scipy.cluster
import seaborn as sns
from scipy.spatial.distance import squareform


def raw_plot(Ca: np.array):
    # Plot the accepted Cells
    
    fig,ax = plt.subplots(figsize=(15,10))
    plt.imshow(Ca, aspect='auto', interpolation='none')
    plt.ylabel('Accepted cells')
    plt.yticks([])
    plt.xlabel('Frames (at 20fps)')
    plt.title('Total Cells = {}'.format(Ca.shape[0]))


def raw_plot_all(Ca: np.array, reaches: pd.DataFrame, mark='insc_peak'):
    # Plot the accepted Cells with marks for reach onsets
    
    fig,ax = plt.subplots(figsize=(15,10))
    plt.imshow(Ca, aspect='auto', interpolation='none')
    #plt.ylabel('Accepted cells')
    plt.yticks(list(range(Ca.shape[0])))
    plt.xlabel('Frames (at 20fps)')
    plt.title('Total Reaches = {}'.format(reaches.shape[0]))
    
    for k in reaches[mark]:
        plt.axvline(x=k, color='k', alpha=0.2)
        

def raw_plot_triggers(Ca: pd.DataFrame, reaches: pd.DataFrame):
    # Plot the accepted Cells with marks for reach onsets and triggers
    
    fig,ax = plt.subplots(figsize=(15,10))
    plt.imshow(Ca, aspect='auto', interpolation='none')
    plt.ylabel('Accepted cells')
    plt.yticks([])
    plt.xlabel('Frames (at 20fps)')
    plt.title('Total Reaches = {}'.format(reaches.shape[0]))
    
    for k in reaches['insc_peak']:
        plt.axvline(x=k, color='k', alpha=0.2)
        
    for k in reaches['insc_trig']:
        plt.axvline(x=k, color='r', alpha=1)

    
def raw_plot_category(cells: np.array, reaches: pd.DataFrame, mark='Grasped'):
    # Plot the accepted Cells with marks for only reinforced reaches
    
    fig,ax = plt.subplots(figsize=(15,10))
    plt.imshow(cells, aspect='auto', interpolation='none', cmap='binary')
    plt.ylabel('Accepted cells')
    plt.yticks([])
    plt.xlabel('Frames (at 20fps)')
    plt.title('Reaches = {}'.format(reaches.loc[reaches['name']=='grasped', 'start'].shape[0]))
    
    for k in reaches.loc[reaches['name']==mark, 'start']:
        plt.axvline(x=k, color='r', alpha=0.2)


def raw_plot_mean(cells: np.array, reaches: pd.DataFrame, mark='Grasped'):
    # Add the Mean Trace and plot it along the others
    
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    im = ax.imshow(cells, aspect='auto', interpolation='none', cmap='binary')
    ax.set_ylabel('Accepted cells')
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title('Reaches = {}'.format(reaches.loc[reaches['category']=='Grasped', 'start'].shape[0]))
    #plt.colorbar(im, orientation='vertical')
    
    for k in reaches.loc[reaches['category']==mark, 'start']:
        ax.axvline(x=k, color='r', alpha=0.2)
    
    ca_mean = cells.mean(axis=0)
    #ax2.imshow(np.expand_dims(ca_mean, axis=0), aspect='auto', interpolation='none')
    ax2.plot(ca_mean)
    
    ax2.set_xlim([0, ca_mean.shape[0]])
    ax2.set_ylim([0, max(ca_mean)])
    ax2.set_xlabel('Frames (at 20fps)')
    ax2.set_ylabel('Mean trace')
    
def raw_plot_accel(cells: np.array, accel: np.array, ori: np.array):
    # Add the Trace for acceleration and plot it along the Ca signal
    
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    ax = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    im = ax.imshow(cells, aspect='auto', interpolation='none', cmap='viridis')
    ax.set_ylabel('Accepted cells')
    ax.set_yticks([])
    ax.set_xticks([])

    ax2.plot(accel)
    
    ax2.set_xlim([0, accel.shape[0]])
    ax2.set_ylim([0, max(accel)])
    ax2.set_ylabel('Accelerometer movement')    
    
    ax3.plot(ori)
    
    ax3.set_xlim([0, ori.shape[0]])
    ax3.set_ylim([0, max(ori)])
    ax3.set_xlabel('Timepoints at 50hz')
    ax3.set_ylabel('Rotation')  


def trials_plot_all(Ca: np.array, reaches: pd.DataFrame, mark='insc_peak'):
    #Plot the Mean population level activity aligned to the event
    
    n=80
    n_list = [0, 0.25*n, 0.5*n, 0.75*n, n, 1.25*n, 1.5*n, 1.75*n, 2*n]
    n_list_2 =['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    
    Ca_mean_array = np.array([])
    all_cells = dict({})
    for i in range(Ca.shape[0]):
        cell_mean_array = np.array([])
        for k in reaches[mark]:
            cell_mean_array = np.append(cell_mean_array, Ca[i, k-n:k+n])   
        cell_mean_array = cell_mean_array.reshape((-1,2*n))
        all_cells.update({i:cell_mean_array})
        Ca_mean_array = np.append(Ca_mean_array, np.nanmean(cell_mean_array, axis=0))
    
    Ca_mean_array = Ca_mean_array.reshape((Ca.shape[0],-1))
    
    
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    im = ax.imshow(Ca_mean_array, aspect='auto', interpolation='none', cmap='viridis')
    ax.set_ylabel('Cells')
    ax.set_xticks([])
    ax.set_title('Mean Activity Aligned to Reach Peak')
    ax.axvline(x=n_list[4], color='r', alpha=0.2)
    #plt.colorbar(im)
    
    
    ax2.plot(Ca_mean_array.mean(axis=0))
    ax2.set_xlim([0, Ca_mean_array.mean(axis=0).shape[0]])
    ax2.set_xticks(n_list)
    ax2.set_xticklabels(labels=n_list_2)
    ax2.axvline(x=n_list[4], color='r', alpha=0.2)
    #ax2.set_xticklabels([])
    ax2.set_ylim([0, max(Ca_mean_array.mean(axis=0))])
    ax2.set_xlabel('Seconds from reach peak')
    ax2.set_ylabel('Mean trace')
    
    return Ca_mean_array

def trials_plot_one(Ca: np.array, reaches: pd.DataFrame, cell=0, mark='insc_peak'):
    #Plot the Single Cell level activity aligned to the event
    
    n=160
    n_list = [0, 0.25*n, 0.5*n, 0.75*n, n, 1.25*n, 1.5*n, 1.75*n, 2*n]
    n_list_2 =['-8', '-6', '-4', '-2', '0', '2', '4', '6', '8']
  
    
    Ca_mean_array = np.array([])
    for k in reaches[mark]:
        Ca_mean_array = np.append(Ca_mean_array, Ca[cell, k-n:k+n])
        
    Ca_mean_array = Ca_mean_array.reshape((-1,2*n))
    
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    im = ax.imshow(ca_mean_array, aspect='auto', interpolation='none', cmap='viridis')
    ax.set_ylabel('Reaches')
    ax.set_xticks([])
    ax.set_title('Mean Activity Cell ' + str(cell))
    
    ax2.plot(Ca_mean_array.mean(axis=0))
    
    ax2.set_xlim([0, ca_mean_array.mean(axis=0).shape[0]])
    ax2.set_xticks(ticks=n_list)
    ax2.set_xticklabels(labels=n_list_2)
    ax2.axvline(x=n_list[4], color='k', alpha=0.2)
    #ax2.set_xticklabels([])
    ax2.set_ylim([0, max(ca_mean_array.mean(axis=0))])
    ax2.set_xlabel('Seconds from reach peak')
    ax2.set_ylabel('Mean trace')
    
    return Ca_mean_array

def trials_spikes_all(spikes: np.array, raster=False):
    # Plot the Spikes in all trials aligned to the event
    # Works with the data in the form of Tensor
    
    n=80
    n_list = [0, 0.25*n, 0.5*n, 0.75*n, n, 1.25*n, 1.5*n, 1.75*n, 2*n]
    n_list_2 =['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    #n_list_3 =np.arange(0,spikes.shape[0],22)
    spikes2 = np.copy(spikes)
    if raster == True:
        i,j=np.where(spikes==1)
        for ii,jj in zip(i,j):
            spikes2[ii,jj]=jj
            
        spike_times = np.where(spikes2==0, np.nan, spikes2)
        frq, edges = np.histogram(spike_times, bins=[0,10,20,30,40,50,60,70,
                                                     80,90,100,110,120,130,140,150,160])
    else:     
        positions=np.where(spikes ==1)[1]
        #positions[(positions <= 40) & (positions >= 20)]
            
    fig = plt.figure(figsize=(6,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    im = ax.imshow(spikes, aspect='auto', interpolation='none', cmap='binary')
    ax.set_ylabel('Cells')
    ax.set_xticks([])
    #ax.set_yticks(n_list_3)
    ax.set_title('Raster Plots Aligned to Reach Start')
    ax.axvline(x=n_list[4], color='r', alpha=0.2)
    
    if raster == True:
        ax2.bar(edges[:-1], frq, width=np.diff(edges), edgecolor="black", align="edge")
        ax2.set_xlim([0, spikes.shape[1]])
        ax2.set_xticks(n_list)
        ax2.set_xticklabels(labels=n_list_2)
        ax2.axvline(x=n_list[4], color='r', alpha=0.2)
        ax2.set_xlabel('Seconds from reach onset')
        ax2.set_ylabel('Spikes sum')
        
    else:
        ax2.boxplot(positions, vert=False)
        ax2.set_xlim([0, spikes.shape[1]])
        ax2.set_xticks(n_list)
        ax2.set_yticks([])
        ax2.set_xticklabels(labels=n_list_2)
        ax2.axvline(x=n_list[4], color='r', alpha=0.2)
        ax2.set_xlabel('Seconds from reach onset')
        ax2.set_ylabel('Spike timing')


def one_trial_all(Ca: np.array, reaches: pd.DataFrame):
    # Plot all cells activity in a single reach
    # This function works with the data in the form of Tensor
    
    n=80
    n_list = [0, 0.25*n, 0.5*n, 0.75*n, n, 1.25*n, 1.5*n, 1.75*n, 2*n]
    n_list_2 =['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    
    for k in range(Ca[2]):
        fig = plt.figure(figsize=(15,10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
        ax = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0])
        im = ax.imshow(Ca[:,:,k], aspect='auto', interpolation='none', cmap='viridis')
        ax.set_ylabel('Cells')
    
        ax.set_xticks([])
        ax.set_title('Mean Activity Aligned to Reach Peak ' + str(k))
        ax.axvline(x=n_list[4], color='r', alpha=0.5)
        #plt.colorbar(im)
        
        
        ax2.plot(Ca[:,:,k].mean(axis=0))
        ax2.set_xlim([0, Ca[:,:,k].mean(axis=0).shape[0]])
        ax2.set_xticks(n_list)
        ax2.set_xticklabels(labels=n_list_2)
        ax2.axvline(x=n_list[4], color='r', alpha=0.2)
        #ax2.set_xticklabels([])
        ax2.set_ylim([0, max(Ca[:,:,k].mean(axis=0))])
        ax2.set_xlabel('Seconds from reach peak')
        ax2.set_ylabel('Mean trace')
        

def plot_mean(arr: np.array):
    
    # Plot all cells in a big dataset aligned from tensor 
    # (Cells in single trial X Length of the window)

    n=arr.shape[1]/2
    n_list = [0, 0.25*n, 0.5*n, 0.75*n, n, 1.25*n, 1.5*n, 1.75*n, 2*n]
    n_list_2 =['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    
    
    fig = plt.figure(figsize=(6,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    #im = ax.imshow(arr, aspect='auto', vmin=0, vmax=1)
    im = ax.imshow(arr, aspect='auto', cmap='viridis', interpolation='none')
    ax.set_ylabel('CellsxTrials')

    ax.set_xticks([])
    ax.set_title('Mean Activity Aligned to Reach Peak')
    ax.axvline(x=n_list[4], color='r', alpha=0.5)
    #plt.colorbar(im)
    
    
    ax2.plot(arr.mean(axis=0))
    ax2.set_xlim([0, arr.mean(axis=0).shape[0]])
    ax2.set_xticks(n_list)
    ax2.set_xticklabels(labels=n_list_2)
    ax2.axvline(x=n_list[4], color='r', alpha=0.2)
    #ax2.set_xticklabels([])
    ax2.set_ylim([0, max(arr.mean(axis=0))])
    ax2.set_xlabel('Seconds from reach peak')
    ax2.set_ylabel('Mean trace')


def plot_scatterbox(spikes1: np.array, spikes2: np.array, plot='sum', yaxis='Events per minute', labels=['PT cells','IT cells']):
    # Plot two sets of parameters as a boxplot + scatterplot
    
    if plot=='sum':
        event_rates1 = spikes1.sum(axis=1)/(spikes1.shape[1]/1200)
        event_rates2 = spikes2.sum(axis=1)/(spikes2.shape[1]/1200)
    elif plot=='mean':
        event_rates1 = np.nanmean(spikes1, axis=1)
        event_rates2 = np.nanmean(spikes2, axis=1)
    
    xs1=[]
    for i in range(event_rates1.shape[0]):
        xs1.append(np.random.normal(loc=1, scale=0.01))
    
    xs2=[]
    for i in range(event_rates2.shape[0]):
        xs2.append(np.random.normal(loc=2, scale=0.01))
        
    plt.boxplot([event_rates1, event_rates2], showfliers=False)
    plt.scatter(xs1, event_rates1)
    plt.scatter(xs2, event_rates2)
    plt.xticks(ticks=[1,2], labels=labels)
    plt.ylabel(yaxis)
    #plt.ylim([0,200])
    print('t-statistic = %6.3f pvalue = %6.4f' %  stats.ttest_ind(event_rates1, event_rates2, equal_var=False, nan_policy='omit'))
    
    return event_rates1, event_rates2
        

def spike_parameters(spikes_data: tuple):
    # Plot all the parameters of the Ca events
    
    n=80
    n_list = [0, 0.25*n, 0.5*n, 0.75*n, n, 1.25*n, 1.5*n, 1.75*n, 2*n]
    n_list_2 =['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    #n_list_3 =np.arange(0,spikes.shape[0],22)
    spikes2 = np.copy(spikes_data[3])      
    i,j=np.where(spikes2==1)
    for ii,jj in zip(i,j):
        spikes2[ii,jj]=jj

    spike_times = np.where(spikes2==0, np.nan, spikes2)
    frq1, edges1 = np.histogram(spikes_data[1], range=(1,spikes_data[1].max()), bins=20)
    frq2, edges2 = np.histogram(spikes_data[2], range=(1,spikes_data[2].max()), bins=20)
    frq3, edges3 = np.histogram(spike_times, bins=[0,10,20,30,40,50,60,70,
                                                 80,90,100,110,120,130,140,150,160])
        
    fig = plt.figure(figsize=(6,10))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1,1,1])
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    
    ax1.set_title('Main parameters for the peaks')
    ax1.bar(edges1[:-1], frq1, width=np.diff(edges1), edgecolor="black", align="edge")
    ax1.set_xlabel('Peak height, dF/F')
    ax1.set_ylabel('Number of spikes')
    ax1.set_xlim([0, 50])
    ax1.set_xticks([0,10,20,30,40,50])
    
    ax2.bar(edges2[:-1], frq2, width=np.diff(edges2), edgecolor="black", align="edge")
    ax2.set_xlabel('Peak half-width (sec)')
    ax2.set_ylabel('Number of spikes')
    ax2.set_xlim([0, 100])
    ax2.set_xticks([0,20,40,60,80,100])
    ax2.set_xticklabels(['0','1s','2s','3s','4s','5s'])
    
    ax3.bar(edges3[:-1], frq3, width=np.diff(edges3), edgecolor="black", align="edge")
    ax3.set_xlim([0, spikes2.shape[1]])
    ax3.set_xticks(n_list)
    ax3.set_xticklabels(labels=n_list_2)
    ax3.axvline(x=n_list[4], color='r', alpha=0.3)
    ax3.set_xlabel('Seconds from reach onset')
    ax3.set_ylabel('Spikes sum')
        


def plot_correlation(corr = np.array, diagonal=True, mask=False, annot=False, Vmin=0, Vmax=0.25, Thresh=0.1):        
    #Plot the correlation as a triange seaborn matrix
    
    #get rid of the diagonal in symmetric matrix
    if diagonal==False:
        np.fill_diagonal(corr,0)
    
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))
    
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    
    # Draw the heatmap with the mask and correct aspect ratio
    # Generate a mask for the upper triangle
    if mask==True:
        mask = np.triu(np.ones_like(corr, dtype=bool))
        sns.heatmap(corr, cmap=cmap, center=0,
                    square=True, annot=annot, mask=mask | (np.abs(corr) <= Thresh), vmin=Vmin, vmax=Vmax, linewidths=.5, cbar_kws={"shrink": .5})
    else:    
        sns.heatmap(corr, cmap=cmap, center=0,
                    square=True, annot=annot, vmin=Vmin, vmax=Vmax, linewidths=.5, cbar_kws={"shrink": .5})
    
    
    fig1 = plt.figure(figsize=(8, 4))
    ax1 = fig1.add_subplot() 
    frq1, edges1 = np.histogram(corr, bins=20)
    ax1.set_title('Mean={0:.2f}    STD={1:.2f}'.format(corr.mean(), corr.std()), {'fontsize': 14,
 'fontweight' : 'bold'})
    ax1.bar(edges1[:-1], frq1, width=np.diff(edges1), edgecolor="black", align="edge")
    ax1.set_xlabel('correlation')
    ax1.set_ylabel('Number of pairs')



def plot_distances(data: np.array):
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    clust = scipy.cluster.hierarchy.linkage(data, method='centroid')
    sns.clustermap(data, cmap=cmap, center=0 ,row_linkage=clust, col_linkage=clust)
    
    #sns.clustermap(dist['ar[init]'],cmap='cubehelix', metric='correlation', method='complete')
    plt.show()
    


def plot_ccs(cell1: np.array, cell2: np.array):
    # Plot the cross-correlation between two traces
    # Input: two cell traces
    
    cell1 = np.nan_to_num(cell1)
    cell2 = np.nan_to_num(cell2)
    corr = correlate(cell1, cell2, mode='same')
    fig, (ax_cell1, ax_cell2, ax_corr) = plt.subplots(3, 1, sharex=True)
    ax_cell1.plot(cell1)  
    ax_cell1.set_title('Cell1')
    ax_cell2.plot(cell2)
    ax_cell2.set_title('Cell2')
    ax_corr.plot(corr/cell1.shape[0])
    ax_corr.axhline(0.5, ls=':')
    ax_corr.set_title('Cross-correlation')
    fig.tight_layout()
    plt.show()

def plot_overlay(ca_data: np.array, behav_data: np.array, alpha=0.5):
    # plot the Ca_data over the behavioral data
    
    fig,ax1=plt.subplots(figsize=(20,5))
    ax2=plt.twinx()
    ax1.plot(ca_data,'k', alpha=0.5)
    ax2.fill_between(np.arange(0, len(behav_data)), behav_data, color='skyblue',alpha=alpha)
    ax1.set_ylabel('Number of events')
    ax2.set_ylabel('Motor activity')
    ax1.set_xlabel('Frames (at 20fps)')
    plt.show()
    