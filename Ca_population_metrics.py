# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 17:10:37 2023

@author: Daniil.Berezhnoi

These are functions used for the 
Principal Component Analysis 
and related visualizations

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib import gridspec


def principal_components(data: np.array):

    # Perform principal component analysis - 
    # reduces the dimension of input array to 10 PC representing the population activity
    # Input: Ca traces for all the cells 
    # Output: Principal component traces

          
    x = StandardScaler().fit_transform(data)      
    pca = PCA(n_components=10)
    principalComponents = pca.fit_transform(x)
    
    principalDf = pd.DataFrame(data=principalComponents, columns = ['PC'+str(k) for k in range(1,11)])
    PC_variance = pca.explained_variance_ratio_
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    plt1 = ax.plot(principalComponents[:,0], principalComponents[:,1], principalComponents[:,2], color='grey', label='projection', linewidth=1, alpha=0.5)
    
      
    return principalDf, PC_variance 


def d_prime(pc1: np.array, pc2: np.array):

    #Sensitivity index
    # Calculates the difference between two sets of PCs
    # Input: Matrix with principal components for two classes of trials

        
    pc1_mean = pc1.mean(axis=2)
    pc2_mean = pc2.mean(axis=2)
    pc1_dev = pc1.std(axis=2)
    pc2_dev = pc2.std(axis=2)
    
    d_prime = 0
    for k in range(3):
        mean_difference = pc1_mean[k,:] - pc2_mean[k,:]
        mean_difference = mean_difference **2
        
        deviation_sum = pc1_dev[k,:] **2 + pc2_dev[k,:] **2
        deviation_sum = deviation_sum/2
        
        d_prime = d_prime + (mean_difference/deviation_sum)
    
    fig = plt.figure(figsize=(10,10))
    plt.plot(d_prime, color='black', linestyle='solid', linewidth=3)
    plt.xlabel("Time from reach onset (s)")
    plt.ylabel("Sensitivity index (d')")
    n=40
    n_list = [0, 0.25*n, 0.5*n, 0.75*n, n, 1.25*n, 1.5*n, 1.75*n, 2*n]
    n_list_2 =['-2', '-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2']
    plt.xticks(ticks=n_list, labels=n_list_2)
    plt.vlines(n, 0, d_prime.max(), color='black', linestyle='dashed', linewidth=1 )

    return d_prime


def project_pc(pc = np.array, only_mean = False):    
    
    # Plot all the PC components projections for all the reaches in current category
    # Input: Principal components for the whole recording
    # If only_mean option is used than plot only the average trace for all the trials
 
    pc_mean = pc.mean(axis=2)
    
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])
    ax = fig.add_subplot(gs[0,0], projection='3d')
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    ax2.set_ylabel('PC1')
    ax3.set_ylabel('PC2')
    ax4.set_ylabel('PC3')
    
    if only_mean == False:
        for i in range(pc.shape[2]):
            pc_plot = pc[:,:,i]
            
            x = pc_plot[0,:]
            y = pc_plot[1,:]
            z = pc_plot[2,:]
                    
            plt1 = ax.plot(x,y,z, linewidth=1, color='grey', label='open')
            plt2 = ax2.plot(x, color='grey')
            plt3 = ax3.plot(y, color='grey') 
            plt4 = ax4.plot(z, color='grey')
        
    plt5 = ax.plot(pc_mean[0,:], pc_mean[1,:], pc_mean[2,:], color='green', label='open', linewidth=3)
    plt6 = ax2.plot(pc_mean[0,:], color='green')
    plt7 = ax3.plot(pc_mean[1,:], color='green') 
    plt8 = ax4.plot(pc_mean[2,:], color='green')
