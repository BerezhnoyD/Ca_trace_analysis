# Ca_trace_analysis
Ca_traces - analysis/visualization pipeline

First version of the repository with all the functions for Ca imaging traces
analysis and visualization. Includes the script files with following functions:

- open and preprocess the data 
- perform event detection and analyze event parameters
- extract the population activity metrics
- plot the data
- perform the reach related activity analysis (needs additional input from behavioral data)
- perform the spontaneous activity analysis based on events detected

Here is the visualisation of the whole pipeline stored in multiple toolboxes:
- Ca_open_data.py 
- Ca_plot_data.py 
- Ca_event_detection.py 
- Ca_population_metrics.py 

## Pipeline overview
Here is the description of methods used in the simplified Ca_Analysis-Spontaneous_Activity Notebook:

![plot](https://github.com/BerezhnoyD/Ca_trace_analysis/blob/main/pictures/Ca_signal_spontaneous.png)

All of the functions for this notebook are contained in the toolboxes above. The notebook should lead you through
opening the data (by now the function supports only the data from Inscopix data processing software as csv table or
the data processed by CaImAn package), plotting the traces for quality control, detecting Ca events and plotting
the basic parameters for the peaks. Ca_population_metrics is not used at this point and contains more advances analysis
functions.

Here are the few results from this pipeline allowing to 
### Visualize the raw signal (Plot the whole recording)

![alt text](https://github.com/BerezhnoyD/Ca_trace_analysis/blob/main/pictures/Raw_signal.png?raw=true)

### Perform the event detection (Detect and analyze the peaks)

![alt text](https://github.com/BerezhnoyD/Ca_trace_analysis/blob/main/pictures/Event_detection.png?raw=true)


### Visualize the basic parameters for the detected Ca events (Plot the peak parameters)

![alt text](https://github.com/BerezhnoyD/Ca_trace_analysis/blob/main/pictures/Basic_parameters.png?raw=true)


### And correlaton matrix between cells (Calculate and plot asymmetric correlation)

![alt text](https://github.com/BerezhnoyD/Ca_trace_analysis/blob/main/pictures/Correlation_visualization.png?raw=true)


## Secondary analyses overview
The whole toolbox allows to perform many more secondary analyses and the documentation on them is yet to come. Most
of the functions for this part are containde in Ca_population_metrics toolbox.

![alt text](https://github.com/BerezhnoyD/Ca_trace_analysis/blob/main/pictures/Ca_signal_overview.png?raw=true)

By now we're providing the extended info on all the functions in the **'Functions summary.txt'** file:
