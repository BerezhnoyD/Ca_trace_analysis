# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 15:29:33 2023


Pipeline to run the scripts for spontaneous activity analysis

@author: Daniil.Berezhnoi
"""

from Ca_open_data import *
from Ca_plot_data import *
from Ca_event_detection import *
from Ca_population_metrics import *

#%%

directory = 'C:/Users/daniil.berezhnoi.VAI/Desktop/Current Experiments/Ca Imaging/Ca_IT_analysis/14wo/ncb1907/'
file = 'Oct25_ncb1907.csv'
accel_file = file[:-4]+'-accel.csv'
gpio_file = file[:-4]+'-gpio.csv'



#%% Open the traces

Ca_df, Ca_cells, Ca, image = open_traces(directory, file)
#Flir_Frames, Insc_Frames, Trigger_Frames, TriggerOff_Frames = open_gpio(directory, gpio_file)
#accel,ori= open_accel(directory, accel_file)


#%% align with accelerometer data

Insc_Frames.tail(30)

#%%

accel=stretch1d(accel,24133)
ori=stretch1d(ori,24133)


#%% normalize the traces

Ca_IT, Ca_PTall, Ca_PT = split_cellset(Ca, Ca_cells)
Ca_IT_filter = filter_data(Ca_IT)
Ca_IT_max = max_score_data(Ca_IT_filter)
Ca_IT_z = z_score_data(Ca_IT_filter)

Ca_PT_filter  = filter_data(Ca_PT)
Ca_PT_max = max_score_data(Ca_PT_filter)
Ca_PT_z = z_score_data(Ca_PT_filter)

Ca_PTall_filter  = filter_data(Ca_PTall)
Ca_PTall_max = max_score_data(Ca_PTall_filter)
Ca_PTall_z = z_score_data(Ca_PTall_filter)


#%% plot Ca signal along the accelerometer

raw_plot_accel(Ca_PT_max, accel, ori)
plt.savefig(directory+file[:-4]+'_Ca_accel_PT')
raw_plot_accel(Ca_IT_max, accel, ori)
plt.savefig(directory+file[:-4]+'_Ca_accel_IT')
#raw_plot_accel(Ca_PTall_max, accel, ori)
#plt.savefig(directory+'Ca_accel_PTall')

#%% detect events

PT_events = detect_events(Ca_PT_filter)
IT_events = detect_events(Ca_IT_filter)
#PTall_events = detect_events(Ca_PT_all_filter)

#%% plot the parameters for the spikes as histograms


spike_parameters(IT_events)
plt.savefig(directory+'Ca_param_IT')
spike_parameters(PT_events)
plt.savefig(directory+'Ca_param_PT')
spike_parameters(PTall_events)
plt.savefig(directory+'Ca_param_PTall')

#%% compare two groups of neurons

epm1907_14=plot_scatterbox(PT_events[0],IT_events[0], yaxis='Events per minute')
#plt.savefig(directory+'Events_per_minute')
#plt.close()
plt.figure()
auc1907_14=plot_scatterbox(PT_events[9],IT_events[9], plot='mean',  yaxis='Area under the curve')
#plt.ylim([0,600])
#plt.savefig(directory+'Area_under_curve')
plt.figure()
pkh1907_14=plot_scatterbox(PT_events[1],IT_events[1], plot='mean', yaxis='Peak height')


#%% calculate the correlations

PT_correl=compute_asymmetric(PT_events[10])
IT_correl=compute_asymmetric(IT_events[10])
#PTall_correl=compute_asymmetric(PTall_events[10])
PT_correl_norm=compute_asymmetric(PT_events[10],'normalized')
IT_correl_norm=compute_asymmetric(IT_events[10],'normalized')
np.fill_diagonal(PT_correl,0)
np.fill_diagonal(IT_correl,0)
np.fill_diagonal(PT_correl_norm,0)
np.fill_diagonal(IT_correl_norm,0)


#%% plot the jaccard correlations

plot_correlation(IT_correl, diagonal=True)
#plt.savefig(directory+'PT_correl')
plot_correlation(IT_correl_norm, diagonal=True)
#plt.savefig(directory+'IT_correl')


#%% compare two groups of neurons

acc1907_14=plot_scatterbox(PT_correl, IT_correl, plot='mean',  yaxis='Asymmetric correlation')
#plt.savefig(directory+'Cell_correlations')


#%%

acc_norm1907_14=plot_scatterbox(PT_correl_norm,IT_correl_norm, plot='mean',  yaxis='Norm.Asymmetric correlation')
#plt.savefig(directory+'Cell_correlations_normalized')



#%%

pkh1403_17=plot_scatterbox(PT_events[1],IT_events[1], plot='mean', yaxis='Peak height')
auc1403_17=plot_scatterbox(PT_events[9],IT_events[9], plot='mean', yaxis='Area under curve')
epm1403_17=plot_scatterbox(PT_events[0],IT_events[0], plot='sum')

#plt.savefig(directory+'Peak heights')

#%% Comparing different animals on the same day

plot_scatterbox(epm1279_0331[0],epm1280_0331[0], yaxis='Events per minute')
plt.savefig(directory+'Events_per_minute')
plt.close()
plot_scatterbox(auc1279_0331[0],auc1280_0331[0], plot='mean',  yaxis='Area under the curve')
plt.ylim([0,600])
plt.savefig(directory+'Area_under_curve')
plt.close()
plot_scatterbox(pkh1279_0331[0],pkh1280_0331[0], plot='mean', yaxis='Peak height')
plt.savefig(directory+'Peak heights')
plt.close()
acc=plot_scatterbox(acc1279_0317[0],acc1280_0317[0], plot='mean',  yaxis='Asymmetric correlation')
plt.savefig(directory+'Cell_correlations')
plt.close()
acc_norm=plot_scatterbox(acc_norm1279_0331[0],acc_norm1280_0331[0], plot='mean',  yaxis='Norm.Asymmetric correlation')
plt.savefig(directory+'Cell_correlations_normalized')
plt.close()

#%% Calculating and plotting the Pearsons correlation

pearson_PT=correlation_cells(Ca_PT_filter)
plot_correlation(pearson_PT,diagonal=False)
plt.savefig(directory+'PT_Rcorrel')

pearson_IT=correlation_cells(Ca_IT_filter)
plot_correlation(pearson_IT,diagonal=False)
plt.savefig(directory+'IT_Rcorrel')

pc1280_0512=plot_scatterbox(pearson_PT,pearson_IT, plot='mean',  yaxis='Pearsons correlation') 
plt.savefig(directory+'Pearson_correlations')

#%% Calculating same parameters for each cell reach-related activity

def calculate_per_cell(table_PT, table_IT, trials):
    PT_corr = []
    IT_corr = []
    
    for k in range(0,table_PT.shape[0],trials):
        current_cell_PT=detect_events(table_PT[k:k+trials,:])
        PT_correl=compute_asymmetric(current_cell_PT[10])
        PT_corr.append(np.nanmean(PT_correl))
        
    for l in range(0,table_IT.shape[0],trials):
        current_cell_IT=detect_events(table_IT[l:l+trials,:])
        IT_correl=compute_asymmetric(current_cell_IT[10])
        IT_corr.append(np.nanmean(IT_correl))
        
    return np.array(PT_corr), np.array(IT_corr)

#%% Make a big table for a pair of animals and plot the parameters in dynamics

data1398_PT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])
data1403_PT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])
data1398_IT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])
data1403_IT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])

dates = ['17','20','23']

data1398_days = [-1]
data1403_days = [-1]


acc1398 = np.array([[],[]])        #asymmetric correlation
acc1403 = np.array([[],[]])

acc_norm1398 = np.array([[],[]])   #asymmetric correlation normalized
acc_norm1403 = np.array([[],[]])

auc1398 = np.array([[],[]])        #area under the curve
auc1403 = np.array([[],[]])

epm1398 = np.array([[],[]])        #events per minute
epm1403 = np.array([[],[]])

pkh1398 = np.array([[],[]])        #peak heights
pkh1403 = np.array([[],[]])


for date in dates:
    exec('acc1398 = np.append(acc1398, acc1398_' + date + ',axis=1)')
    exec('acc1403 = np.append(acc1403, acc1403_' + date + ',axis=1)')
    exec('acc_norm1398 = np.append(acc_norm1398, acc_norm1398_' + date + ',axis=1)')
    exec('acc_norm1403 = np.append(acc_norm1403, acc_norm1403_' + date + ',axis=1)')
    exec('auc1398 = np.append(auc1398, auc1398_' + date + ',axis=1)')
    exec('auc1403 = np.append(auc1403, auc1403_' + date + ',axis=1)')
    exec('epm1398 = np.append(epm1398, epm1398_' + date + ',axis=1)')
    exec('epm1403 = np.append(epm1403, epm1403_' + date + ',axis=1)')
    exec('pkh1398 = np.append(pkh1398, pkh1398_' + date + ',axis=1)')
    exec('pkh1403 = np.append(pkh1403, pkh1403_' + date + ',axis=1)')
    
    exec('data1398_days.append(data1398_days[-1] + acc1398_' + date +'[0].shape[0])')
    exec('data1403_days.append(data1403_days[-1] + acc1403_' + date +'[0].shape[0])') 

#acc_norm1398 = np.append([[np.nan, np.nan, np.nan],[np.nan, np.nan, np.nan]],acc_norm1398, axis=1)    

data1398_PT['Events per minute'] = epm1398[0,:]
data1398_PT['Peak height'] = pkh1398[0,:]
data1398_PT['Area under curve'] = auc1398[0,:]
data1398_PT['Asym.correlation'] = acc1398[0,:]
#data1398_PT['Asym.correlation_norm'] = acc_norm1398[0,:]
data1398_PT['Animal'] = 'Control(ncb1398)'


data1398_IT['Events per minute'] = epm1398[1,:]
data1398_IT['Peak height'] = pkh1398[1,:]
data1398_IT['Area under curve'] = auc1398[1,:]
data1398_IT['Asym.correlation'] = acc1398[1,:]
#data1398_IT['Asym.correlation_norm'] = acc_norm1398[1,:]
data1398_IT['Animal'] = 'Control(ncb1398)'


data1403_PT['Events per minute'] = epm1403[0,:]
data1403_PT['Peak height'] = pkh1403[0,:]
data1403_PT['Area under curve'] = auc1403[0,:]
data1403_PT['Asym.correlation'] = acc1403[0,:]
#data1403_PT['Asym.correlation_norm'] = acc_norm1403[0,:]
data1403_PT['Animal'] = 'mitoPark(ncb1403)'


data1403_IT['Events per minute'] = epm1403[1,:]
data1403_IT['Peak height'] = pkh1403[1,:]
data1403_IT['Area under curve'] = auc1403[1,:]
data1403_IT['Asym.correlation'] = acc1403[1,:]
#data1403_IT['Asym.correlation_norm'] = acc_norm1403[1,:]
data1403_IT['Animal'] = 'mitoPark(ncb1403)'

for k, pair in enumerate(zip(data1398_days[:-1], data1398_days[1:])):
    data1398_PT.loc[pair[0]+1:pair[1], 'Day'] = dates[k]
    data1398_IT.loc[pair[0]+1:pair[1], 'Day'] = dates[k]
    
for k, pair in enumerate(zip(data1403_days[:-1], data1403_days[1:])):
    data1403_PT.loc[pair[0]+1:pair[1], 'Day'] = dates[k]
    data1403_IT.loc[pair[0]+1:pair[1], 'Day'] = dates[k]
    
data_PT = pd.concat([data1398_PT, data1403_PT])
data_IT = pd.concat([data1398_IT, data1403_IT])

#%%
import pandas as pd
import numpy as np

animals = ['1279', '1280', '1398', '1403', '1790', '1796']

data18w_PT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])
data18w_IT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])

acc18w = np.array([[],[]])        #asymmetric correlation

acc_norm18w = np.array([[],[]])   #asymmetric correlation normalized

auc18w = np.array([[],[]])        #area under the curve

epm18w = np.array([[],[]])        #events per minute

pkh18w = np.array([[],[]])        #peak heights

animal_name = np.array([])

for animal in animals:
    exec('acc18w = np.append(acc18w, acc' + animal + '_18 ,axis=1)')
    exec('acc_norm18w = np.append(acc_norm18w, acc_norm' + animal + '_18 ,axis=1)')
    exec('auc18w = np.append(auc18w, auc' + animal + '_18 ,axis=1)')
    exec('epm18w = np.append(epm18w, epm' + animal + '_18 ,axis=1)')
    exec('pkh18w = np.append(pkh18w, pkh' + animal + '_18 ,axis=1)')
    exec('length = len(acc' + animal + '_18[0])')
    
    for k in range(length):
        animal_name = np.append(animal_name, animal)
        

data18w_PT['Events per minute'] = epm18w[0,:]
data18w_PT['Peak height'] = pkh18w[0,:]
data18w_PT['Area under curve'] = auc18w[0,:]
data18w_PT['Asym.correlation'] = acc18w[0,:]
data18w_PT['Asym.correlation_norm'] = acc_norm18w[0,:]
data18w_PT['Animal'] = animal_name


data18w_IT['Events per minute'] = epm18w[1,:]
data18w_IT['Peak height'] = pkh18w[1,:]
data18w_IT['Area under curve'] = auc18w[1,:]
data18w_IT['Asym.correlation'] = acc18w[1,:]
data18w_IT['Asym.correlation_norm'] = acc_norm18w[1,:]
data18w_IT['Animal'] = animal_name

#%%

mito = ['1789', '1855',  '1895', '1899', '1280', '1398', '1403']

data18w_PT['Day'] = '18w'
data18w_PT['Group'] = 'Control'
data18w_PT.loc[data18w_PT['Animal'].isin(mito),'Group'] = 'mitoPark'

data18w_IT['Day'] = '18w'
data18w_IT['Group'] = 'Control'
data18w_IT.loc[data18w_IT['Animal'].isin(mito),'Group'] = 'mitoPark'

#%%
import pandas as pd

data12_16_PT = pd.concat([data12w_PT, data14w_PT, data16w_PT], ignore_index=True)
data12_16_IT = pd.concat([data12w_IT, data14w_IT, data16w_IT], ignore_index=True)

#%%

animals = ['1279', '1280', '1398', '1403']

data17w_PT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])
data17w_IT = pd.DataFrame(columns = ['Animal','Day','Events per minute','Peak height','Area under curve','Asym.correlation',
                                   'Asym.correlation_norm'])

acc17w = np.array([[],[]])        #asymmetric correlation

acc_norm17w = np.array([[],[]])   #asymmetric correlation normalized

auc17w = np.array([[],[]])        #area under the curve

epm17w = np.array([[],[]])        #events per minute

pkh17w = np.array([[],[]])        #peak heights

animal_name = np.array([])

for animal in animals:
    exec('acc17w = np.append(acc17w, acc' + animal + '_17 ,axis=1)')
    exec('acc_norm17w = np.append(acc_norm17w, acc_norm' + animal + '_17 ,axis=1)')
    exec('auc17w = np.append(auc17w, auc' + animal + '_17 ,axis=1)')
    exec('epm17w = np.append(epm17w, epm' + animal + '_17 ,axis=1)')
    exec('pkh17w = np.append(pkh17w, pkh' + animal + '_17 ,axis=1)')
    exec('length = len(acc' + animal + '_17[0])')
    
    for k in range(length):
        animal_name = np.append(animal_name, animal)
        

data17w_PT['Events per minute'] = epm17w[0,:]
data17w_PT['Peak height'] = pkh17w[0,:]
data17w_PT['Area under curve'] = auc17w[0,:]
data17w_PT['Asym.correlation'] = acc17w[0,:]
data17w_PT['Asym.correlation_norm'] = acc_norm17w[0,:]
data17w_PT['Animal'] = animal_name


data17w_IT['Events per minute'] = epm17w[1,:]
data17w_IT['Peak height'] = pkh17w[1,:]
data17w_IT['Area under curve'] = auc17w[1,:]
data17w_IT['Asym.correlation'] = acc17w[1,:]
data17w_IT['Asym.correlation_norm'] = acc_norm17w[1,:]
data17w_IT['Animal'] = animal_name


#%% Calculate the simple statistic on the data

import scipy
import seaborn as sns
import matplotlib.pyplot as plt

pvalues = []
pvalue_asterisks = []
df = per_animal_IT.reset_index()
x_values = df["Day"].unique()
hue1=df["Group"].unique()[0]
hue2=df["Group"].unique()[1]
parameter = 'Events per minute'
label = 'Events per minute'


box = sns.boxplot(x="Day", y=parameter, hue="Group", data=df,
            whis=1.5, width=.6, palette="vlag")

sns.stripplot(x="Day", y=parameter, hue="Group", data=df, dodge=True, legend=None)

box.set_ylabel(parameter, fontdict= { 'fontsize': 14, 'fontweight':'bold'})
box.set_xlabel('')
#box.set_ylim([0,1000])
box.set_ylabel(label)
box.get_legend().remove()




def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"



for x in x_values:
    stat, pvalue = scipy.stats.ttest_ind(
        df.loc[(df["Group"] == hue1) & (df["Day"] == x) , parameter],
        df.loc[(df["Group"] == hue2) & (df["Day"] == x) , parameter]
    )
    
    pvalues.append(pvalue)
    pvalue_asterisks.append(convert_pvalue_to_asterisks(pvalue))
    

y_position = df[parameter].max() * 0.93
for idx, pval in enumerate(pvalue_asterisks):
     plt.text(x=idx, y=y_position, s=pval, fontdict= { 'fontsize': 14, 'fontweight':'bold'} )
    
    
stars=dict(zip(x_values,pvalues))
print(stars)

#%%

data12w_PT.to_csv('PT_Spontaneous_Ca_data_12wo.csv')
data12w_IT.to_csv('IT_Spontaneous_Ca_data_12wo.csv')