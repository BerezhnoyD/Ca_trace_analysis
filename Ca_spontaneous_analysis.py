# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 15:29:33 2023


Pipeline to run the scripts for spontaneous activity analysis

@author: Daniil.Berezhnoi
"""
#%%

directory = 'C:/Users/daniil.berezhnoi.VAI/Desktop/Current Experiments/Ca Imaging/1279-5/Apr14-ncb1279/'
file = '14Apr_ncb1279.csv'
accel_file = '14Apr_ncb1279-accel.csv'
gpio_file = '14Apr_ncb1279-gpio.csv'



#%% Open the traces

Ca_df, Ca_cells, Ca, image = open_traces(directory, file)

Flir_Frames, Insc_Frames, Trigger_Frames, TriggerOff_Frames = open_gpio(directory, gpio_file)
accel,ori= open_accel(directory, accel_file)

#%% align with accelerometer data

Insc_Frames.tail(10)

#%%

accel=stretch1d(accel,26406)
ori=stretch1d(ori,26406)


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
plt.savefig(directory+'Ca_accel_PT')
raw_plot_accel(Ca_IT_max, accel, ori)
plt.savefig(directory+'Ca_accel_IT')
raw_plot_accel(Ca_PTall_max, accel, ori)
plt.savefig(directory+'Ca_accel_PTall')

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

epm1280_0317=plot_scatterbox(PT_events[0],IT_events[0], yaxis='Events per minute')
plt.savefig(directory+'Events_per_minute')
plt.close()
auc1280_0317=plot_scatterbox(PT_events[9],IT_events[9], plot='mean',  yaxis='Area under the curve')
plt.ylim([0,600])
plt.savefig(directory+'Area_under_curve')

#%% calculate the correlations

PT_correl=compute_asymmetric(PT_events[10])
IT_correl=compute_asymmetric(IT_events[10])
#PTall_correl=compute_asymmetric(PTall_events[10])
PT_correl_norm=compute_asymmetric(PT_events[10],'normalized')
IT_correl_norm=compute_asymmetric(IT_events[10],'normalized')

#%% plot the jaccard correlations

plot_correlation(PT_correl, diagonal=False)
plt.savefig(directory+'PT_correl')
plot_correlation(IT_correl, diagonal=False)
plt.savefig(directory+'IT_correl')


#%% compare two groups of neurons

acc1280_0317=plot_scatterbox(PT_correl,IT_correl, plot='mean',  yaxis='Asymmetric correlation')
plt.savefig(directory+'Cell_correlations')


#%%
acc_norm1279_0317=plot_scatterbox(PT_correl_norm,IT_correl_norm, plot='mean',  yaxis='Norm.Asymmetric correlation')
plt.savefig(directory+'Cell_correlations_normalized')

#%%

pkh1280_0317=plot_scatterbox(PT_events[1],IT_events[1], plot='mean', yaxis='Peak height')
plt.savefig(directory+'Peak heights')

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
acc=plot_scatterbox(acc1279_0331[0],acc1280_0331[0], plot='mean',  yaxis='Asymmetric correlation')
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