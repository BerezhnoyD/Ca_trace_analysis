# -*- coding: utf-8 -*-
"""
Created on Fri May 26 18:05:05 2023

@author: Daniil.Berezhnoi

These are the functions to detect neural events from Ca traces acquired
with the Inscopix Data Analysis Software


"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import signal, ndimage
#import caiman
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.signal import find_peaks, peak_widths, correlate
from itertools import product


def detect_events_caiman(cells: np.array, p=1, s_min = 0):
    
    # Deconvolution using the OASIS built-in in CaImAn package

    # Extract DF/F values
    # cells.detrend_df_f(quantileMin=8, frames_window=250)
    # cells.select_components(use_object=True)

    # s_min is set automatically if 0, but can be chosen to discard certain small peaks if 0.5 or 1
   
    Ca_Spike_array = np.array([])
    Spike_array = np.array([])
    for k in range(cells.shape[0]):
        cell=np.nan_to_num(cells[k,:])
        cell=caiman.source_extraction.cnmf.deconvolution.constrained_foopsi(fluor=cell, method_deconvolution='oasis', p=p, s_min = s_min)
        Ca_Spike_array = np.append(Ca_Spike_array, cell[0])
        Spike_array = np.append(Spike_array, cell[5])
    
    Ca_Spike_array= Ca_Spike_array.reshape((cells.shape[0],-1))
    Spike_array= Spike_array.reshape((cells.shape[0],-1))
    
    plt.figure()
    plt.imshow(cells,aspect='auto')
    plt.title('Raw Ca data')
    plt.figure()
    plt.imshow(Ca_Spike_array,aspect='auto')
    plt.title('Deconvolved Ca data')
    plt.figure()
    plt.imshow(Spike_array,aspect='auto', vmin=0, vmax=1)
    plt.title('Spikes inferred from Ca data - OASIS')
    
    return Ca_Spike_array, Spike_array


def detect_events(cells: np.array, sd = 2, wid = 10 , dist = 6, rel=0.5):    
    # Finding the Calcium events simply by thresholding
    
    Peak_array = np.array([])
    Peak_height = np.array([])
    Peak_width = np.array([])
    L_HalfPeak = np.array([])
    R_HalfPeak = np.array([])
    Peak_width2 = np.array([])
    L_HalfPeak2 = np.array([])
    R_HalfPeak2 = np.array([])
    for k in range(cells.shape[0]):
        cell=np.nan_to_num(cells[k,:])
        peaks, properties = find_peaks(cell, height=cell.std()*sd, width=wid, distance = dist)
        peaks_full = peak_widths(cell, peaks, rel_height=rel)
        cell = np.zeros_like(cell)
        height = np.zeros_like(cell)
        width = np.zeros_like(cell)          
        left = np.zeros_like(cell) 
        right = np.zeros_like(cell)
        width2 = np.zeros_like(cell)
        left2 = np.zeros_like(cell) 
        right2 = np.zeros_like(cell)
        n=0
        for i in peaks:    
            cell[i] = 1
            height[i] = properties['peak_heights'][n]
            width[i] = properties['widths'][n]
            width2[i] = peaks_full[0][n]
            n += 1
            
        for i in np.rint(properties['left_ips']).astype(np.int32):            
            left[i] = 1
            
        for i in np.rint(properties['right_ips']).astype(np.int32):            
            right[i] = 1
            
        for i in np.rint(peaks_full[2]).astype(np.int32):            
            left2[i] = 1
            
        for i in np.rint(peaks_full[3]).astype(np.int32):            
            right2[i] = 1
            
        Peak_array = np.append(Peak_array, cell)
        Peak_height = np.append(Peak_height, height)
        Peak_width = np.append(Peak_width, width)
        L_HalfPeak = np.append(L_HalfPeak, left)
        R_HalfPeak = np.append(R_HalfPeak, right)
        Peak_width2 = np.append(Peak_width2, width2)
        L_HalfPeak2 = np.append(L_HalfPeak2, left2)
        R_HalfPeak2 = np.append(R_HalfPeak2, right2)
        
    Peak_array= Peak_array.reshape((cells.shape[0],-1)).astype(int)
    Peak_height= Peak_height.reshape((cells.shape[0],-1)).astype(float)
    Peak_width= Peak_width.reshape((cells.shape[0],-1)).astype(float)
    L_HalfPeak= L_HalfPeak.reshape((cells.shape[0],-1)).astype(int)
    R_HalfPeak= R_HalfPeak.reshape((cells.shape[0],-1)).astype(int)
    Peak_width2= Peak_width2.reshape((cells.shape[0],-1)).astype(float)
    L_HalfPeak2= L_HalfPeak2.reshape((cells.shape[0],-1)).astype(int)
    R_HalfPeak2= R_HalfPeak2.reshape((cells.shape[0],-1)).astype(int)
    
    
    L_Peak=np.where(L_HalfPeak2 == 1)
    R_Peak=np.where(R_HalfPeak2 == 1)
    Peak=np.where(Peak_array == 1)
    Ca_Spikes=np.empty_like(cells)
    Ca_Spikes[:]=np.nan
    Areas=np.empty_like(cells)
    Areas[:]=np.nan
    Onsets=np.zeros_like(cells)
    
        
    for a,b,c,d in zip(L_Peak[0],L_Peak[1],R_Peak[1],Peak[1]):
        Ca_Spikes[a,b:c] = cells[a,b:c]
        Areas[a,b] = np.trapz(cells[a,b:c])
        Onsets[a,b:d] = 1
        
    Sum_Onsets=Onsets.sum(axis=0)
    
    plt.figure()
    plt.imshow(Onsets,aspect='auto', vmin=0, vmax=1)
    plt.title('Onsets positions')
    plt.figure()
    plt.imshow(Peak_height,aspect='auto', vmin=0, vmax=Peak_height.max())
    plt.title('Peak heights')
    plt.figure()
    plt.plot(cells[0])
    plt.plot(Peak_array[0])
    plt.plot(L_HalfPeak2[0])
    plt.plot(Ca_Spikes[0])
    
    
    return Peak_array, Peak_height, Peak_width, L_HalfPeak, R_HalfPeak, Peak_width2, L_HalfPeak2, R_HalfPeak2, Ca_Spikes, Areas, Onsets,  Sum_Onsets 


def correlation_cells(data: np.array):
    # Calculating Pearson's correlation
    
    cor_data = np.corrcoef(data)
    
    return cor_data

# From old MESMERIZE package
# https://github.com/kushalkolar/MESmerize/blob/master/mesmerize/analysis/math/cross_correlation.py
# @author: kushal
#
# Chatzigeorgiou Group
# Sars International Centre for Marine Molecular Biology
#
# GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 

def ncc_c(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Must pass 1D array to both x and y
    :param x: Input array [x1, x2, x3, ... xn]
    :param y: Input array [y2, y2, x3, ... yn]
    :return:  Returns the normalized cross correlation function (as an array) of the two input vector arguments "x" and "y"
    :rtype: np.ndarray
    """
    cc = correlate(x.reshape(-1, 1), y.reshape(-1, 1))
    denominator = (np.sqrt(np.dot(x.reshape(1, -1),x.reshape(-1, 1)))*np.sqrt(np.dot(y.reshape(1, -1),y.reshape(-1, 1))))
    cc = cc/denominator
    return cc

def get_omega(x: np.ndarray = None, y: np.ndarray = None, cc: np.ndarray = None) -> int:
    """
    Must pass a 1D array to either both "x" and "y" or a cross-correlation function (as an array) to "cc"
    :param x: Input array [x1, x2, x3, ... xn]
    :param y: Input array [y2, y2, x3, ... yn]
    :param cc: cross-correlation function represented as an array [c1, c2, c3, ... cn]
    :return: index (x-axis position) of the global maxima of the cross-correlation function
    :rtype: np.ndarray
    """
    if cc is None:
        cc = ncc_c(x, y)
    w = np.argmax(cc)
    return w

def get_lag(x: np.ndarray = None, y: np.ndarray = None, cc: np.ndarray = None) -> float:
    """
    Must pass a 1D array to either both "x" and "y" or a cross-correlation function (as an array) to "cc"
    :param x: Input array [x1, x2, x3, ... xn]
    :param y: Input array [y2, y2, x3, ... yn]
    :param cc: cross-correlation function represented as a array [c1, c2, c3, ... cn]
    :return: Position of the maxima of the cross-correlation function with respect to middle point of the cross-correlation function
    :rtype: np.ndarray
    """
    if cc is None:
        w = get_omega(x, y)
        s = w - x.shape[0]
    else:
        w = get_omega(cc=cc)
        s = w - int(cc.shape[0] / 2)
    return float(-s)

def get_epsilon(x: np.ndarray = None, y: np.ndarray = None, cc: np.ndarray = None) -> float:
    """
    Must pass a 1D vector to either both "x" and "y" or a cross-correlation function to "cc"
    :param x: Input array [x1, x2, x3, ... xn]
    :param y: Input array [y2, y2, x3, ... yn]
    :param cc: cross-correlation function represented as an array [c1, c2, c3, ... cn]
    :return: Magnitude of the global maxima of the cross-correlationn function
    :rtype: np.ndarray
    """
    if cc is None:
        cc = ncc_c(x, y)
    e = np.max(cc)
    return e


def get_lag_matrix(curves: np.ndarray = None, ccs: np.ndarray = None) -> np.ndarray:
    """
    Get a 2D matrix of lags. Can pass either a 2D array of 1D curves or cross-correlations
    :param curves: 2D array of 1D curves
    :param ccs:    2D array of 1D cross-correlation functions represented by arrays
    :return: 2D matrix of lag values, shape is [n_curves, n_curves]
    :rtype: np.ndarray
    """
    if ccs is None:
        m = curves.shape[0]
        a = np.zeros((m, m))

        for i, j in product(*(range(m),) * 2):
            a[i, j] = get_lag(curves[i], curves[j])
        return a

    return _compute_from_ccs(ccs, get_lag)


def get_epsilon_matrix(curves: np.ndarray = None, ccs: np.ndarray = None) -> np.ndarray:
    """
    Get a 2D matrix of maximas. Can pass either a 2D array of 1D curves or cross-correlations
    :param curves: 2D array of 1D curves
    :param ccs:    2D array of 1D cross-correlation functions represented by arrays
    :return: 2D matrix of maxima values, shape is [n_curves, n_curves]
    :rtype: np.ndarray
    """
    if ccs is None:
        m = curves.shape[0]
        a = np.zeros((m, m))
        for i, j in product(*(range(m),) * 2):
            a[i, j] = get_epsilon(curves[i], curves[j])
        return a

    return _compute_from_ccs(ccs, get_epsilon)


def _compute_from_ccs(ccs: np.ndarray, func: callable) -> np.ndarray:
    m = ccs.shape[0]
    a = np.zeros(shape=(m, m))
    for i, j in product(*(range(m),) * 2):
        a[i, j] = func(cc=ccs[i, j, :])
    return a

def compute_ccs(a: np.ndarray) -> np.ndarray:
    """
    Compute cross-correlations between all 1D curves in a 2D input array
    :param a: 2D input array of 1D curves, shape is [n_samples, curve_size]
    :rtype: np.ndarray
    """
    n = a.shape[0]
    m = a.shape[1]
    out = np.zeros(shape=(n, n, (2 * m) - 1))

    for i, j in product(*(range(n),) * 2):
        out[i, j, :] = ncc_c(a[i,:], a[j,:]).swapaxes(0,1)
    return out


def compute_cc_data(curves: np.ndarray):
    """
    Compute cross-correlation data (cc functions, lag and maxima matrices)
    :param curves: input curves as a 2D array, shape is [n_samples, curve_size]
    :return:    cross correlation data for the input curves as a CC_Data instance
    :rtype: CC_Data
    
    ccs - Cross correlation functions
    l - Position of the maxima of the cross-correlation function with respect to middle point of the cross-correlation function
    e - Magnitude of the global maxima of the cross-correlationn function
    
    
    """
    ccs = compute_ccs(curves)

    l = get_lag_matrix(ccs=ccs)
    e = get_epsilon_matrix(ccs=ccs)

    return ccs, l, e

#Calculate the Jakkard distance between two cells - only on binary traces
# Jaccard distance https://www.jneurosci.org/content/40/50/9576

def asymmetric_correlation(vec1, vec2):     #formula for asymmetric Jaccard similarity
    
    both=[]
    for i in range(0, len(vec1)):
        both.append(int(all([vec1[i], vec2[i]])))

    j=sum(both)/sum(vec1)

    return j

def symmetric_correlation(vec1, vec2):     #formula for symmetric Jaccard similarity
    
    both=[]
    for i in range(0, len(vec1)):
        both.append(int(all([vec1[i], vec2[i]])))

    j=sum(both)/sum(vec1)+sum(vec2)

    return j



def compute_asymmetric(a: np.ndarray, correlation='asymmetric') -> np.ndarray:
    """
    Compute asymmetric correlation between all 1D curves in a 2D input array
    :param a: 2D input array of 1D curves, shape is [n_samples, curve_size]
    :rtype: np.ndarray
    """
    n = a.shape[0]
    m = a.shape[0]
    out = np.zeros(shape=(n, m))
    
    if correlation=='asymmetric':
        for i, j in product(range(n),range(m)):
            out[i, j] = asymmetric_correlation(a[i,:], a[j,:])
        return out
    
    if correlation=='symmetric':
        for i, j in product(range(n),range(m)):
            out[i, j] = symmetric_correlation(a[i,:], a[j,:])
        return out
    
    if correlation=='normalized':
        for i, j in product(range(n),range(m)):
            a1 = a[i,:]
            a2 = a[j,:]
            out_real = asymmetric_correlation(a1, a2)
            np.random.shuffle(a1)
            np.random.shuffle(a2)            
            out_shuffle = asymmetric_correlation(a1, a2)
            while out_shuffle == 0:
               np.random.shuffle(a1)
               np.random.shuffle(a2)            
               out_shuffle = asymmetric_correlation(a1, a2) 
            out[i, j] = out_real/out_shuffle
        return out



def spike_to_txt(directory, Spike_array):
    #Saving the numpy array for raster plots
    
    spikes=np.where(Spike_array >1)
    spikes=np.stack((spikes[0],spikes[1]/20), axis=0).T.astype(np.float64)
    np.savetxt(directory+'reaching_spikes.txt', spikes, fmt='%.1f',delimiter='\t', newline='\n')


def event_based_sort(data: np.array, events: np.array):
    sort_order = []
    
    for n in range(events.shape[0]):
        found=np.argwhere(events[n]== 1)
        if len(found)==0:
            found=np.array([[0]])
        sort_order.append(found[0][0])
    
    sort_order = pd.DataFrame(sort_order).sort_values(by=0).index
    
    return data[sort_order,:]
        
    
        
    
            
        