Files Required: pun.c, init.c, dataplot.c, rand.c, hist.c, fplot.c, allplot.c, sajob.c, xsum.c, sa.config

Compiling and running in pc:
  - set parameter CLUSTER to 0 in pun.c
  - make
  - ./pun sa.config
  sa.config is configuration file where you can change parameters values

Compiling and running in pc:
  - set parameter CLUSTER to 1 in pun.c
  - make
  - ./sajob
  In sajob.c file, you can change parameter values.

Parameters defined pun.c
  Accumulatedpayoff: To use accumulated payoff in simulation, set it to 1 else 0
  PUNISH           : Set it to 1 to use punishment model, else use 0 (no punishment)
  AGGR             : Use aggressiveness trait in model if set to 1
  FORESIGHT        : If set to 1, simulation uses foresight model to punish
  CLUSTER          : To run in cluster set it to 1; if running in pc, set it to 0
  EGALITARIAN      : Set to 1 to use egalitarian model, i.e., all individuals get equal share of benefit and have equal strengths
  DISP_MATRIX      : Set to 1 to display strenth matrix and delta matrix
  ALLDATAFILE      : Set to 1 if you need all data files. If it is set, summary data files are also generated else only whole run files are generated
  GRAPHS           : Set to 1 if you want to save graphs as files.
  SKIP             : Interval in which statistics are calculated
  EVENTS           : Number of events used in simulations
  TRAITS           : Number of traits used in strategy
  NGUT             : Number of groups taking part in 'us vs them' games
  
Output
  Data files
    *e(i).dat      : (i) is index certain run. Data file for efforts evolution in ith run.
    *f(i).dat      : (i) is index certain run. Data file for payoffs evolution in ith run.
    *pi_g(i).dat   : (i) is index certain run. Data file for group payoffs evolution in ith run.
    *dxi(i).dat    : (i) is index certain run. Data file for threshold evolution in ith run.
    *dsi(i).dat    : (i) is index certain run. Data file for aggressiveness evolution in ith run.
    *pun_i(k).dat  : (k) is index certain run. Data file for cost of punishing kth run.
    *pun_j(k).dat  : (k) is index certain run. Data file for punishment kth run.
    *e.dat         : Data file for efforts evolution with time averaged over number of runs.
    *f.dat         : Data file for payoffs evolution with time averaged over number of runs.
    *pi_g.dat      : Data file for group payoffs evolution with time averaged over number of runs.
    *dxi.dat       : Data file for threshold evolution with time averaged over number of runs.
    *dsi.dat       : Data file for aggressiveness evolution with time averaged over number of runs.
    *pun_i.dat     : Data file for cost of punishing others with time averaged over number of runs.
    *pun_j.dat     : Data file for punishment from others with time averaged over number of runs.
    *esum.dat      : Data file for summary of efforts averaged over number of runs.
    *fsum.dat      : Data file for summary of payoffs averaged over number of runs.
    *pi_gsum.dat   : Data file for summary of group payoffs averaged over number of runs.
    *dxisum.dat    : Data file for sumamry of threshold averaged over number of runs.
    *dsisum.dat    : Data file for summary of aggressiveness averaged over number of runs.
    *punisum.dat   : Data file for summary of cost of punishing others averaged over number of runs.
    *punjsum.dat   : Data file for summary of punishment from others averaged over number of runs.    
  
  Graphs
    ef_*.png       : This graph contains plots of evolution of efforts, payoffs, threshold, aggressiveness, cost of punishing and punishment for specific combination of parameters.
    sum_*.png      : This graph contains plots of summary of efforts, payoffs, threshold, aggressiveness, cost of punishing and punishment for specific combination of parameters.
    
  Generating graphs
    - for individual runs and averaged graphs:
      - set parameters in allplot.c file
      - compile changes by 'make' command in terminal
      - run './allplot' command to generate graphs
    
    - for summary graphs
      - set parameters in hist.c file
      - compile changes by 'make' command in terminal
      - run './hist' command to generate graphs
      
  Generic graph
    - it plots summary of runs but more in customized way for either one of parameters (efforts, payoffs, threshold, aggressiveness, cost of punishing, punishment)
    - one can change parameters to display it in clustered form or row stacked form along with options of choosing variation in parameter in plots along x-axis and y-axis
    - one limitation it is not working well for row stacked graphs, if a cluster has only one vertical bar, i.e. VBAR parameter has only one value
    Parameters to change layout of plots
      STACKED     : 1 or 0; if 1, row stacked graphs will be plotted; if 0, clustered normal bar graphs will be plotted
      xAXIS       : parameter for clusters in a single plot; it must be name of array of parameter values
      yAXIS       : trait or value of interest graphs are plotted for
      XAXIS       : parameter for different plots along a row; it must be name of array of parameter values
      YAXIS       : parameter for different plots along a column; it must be name of array of parameter values
      VBAR        : parameter for STACKED option which denote which parameter value should differ along vertical stacked bars in a cluster; if only one option in parameter value, then plot is not working properly at present; it must be name of array of parameter values
      
  Generic graph script xsum.c
    - parameters to be included in graph plots are under the section "Parameters for plotting"
    - parameters that define layout of graph plots are defined under "Parameters to specify type of plots"
    - Pn stores value of parameter currently under multiple loops of parameters for each indexed parameter
    - P array points to reference of parameter used in plotting
    - Pa array stores alphbet strings array used to denote specific parameter in creating data file 
    - Ps array stores no. of values under one parameter indexed by its index.
    - number of parameters are 11 for now, so array size is also fixed at 11. 
    
    Methods in xsum.c
      set_param():
	It sets arrays P, Pa and Ps according to layout for easy manipulation of data files. 
	P is array of pointer to parameter array. 
	Pa is array of alphabet strings representing parameters name.
	Ps is array of number of values in parameter array.

      prep_infile(char *fname, char *appnd):
	It prepares name of file from which data to be read. The method creates string which combination of parameter alphabets (Pa) and current parameter value (Pn) in the loop according to format of data files generated from simulation.
	'appnd' is string to be attached at the end of the string created by the method at the end to specify data file with property of interest like effort or payoff or threshold or punishment.
	The created filename representing data file is assigned to 'fname'.
	It needs to be modified if there is modification in parameters like addition, removal and change in parameter.
	
      prep_outfile(char *file, char *appnd):
	It prepares name of file which stores data ready for plotting graphs. The method creates string which combination of parameter alphabets (Pa) and current parameter value (Pn) in the loop.
	'appnd' is string to be attached at the end of the string created by the method at the end to specify data file with property of interest like effort or payoff or threshold or punishment.
	The created filename representing data file is assigned to 'file'.
	It needs to be modified if there is modification in parameters like addition, removal and change in parameter.
	
      prepare_data_files():
	It prepares all data file needed to plot graphs from data files generated by simulations. It prepares data file for only one of property specified in 'yAXIS' parameter.
	It needs to be modified if there is modification in parameters like addition, removal and change in parameter.
	
      plot_graphs():
        It plots graphs for all parameter combinations.
        It needs to be modified if there is modification in parameters like addition, removal and change in parameter.
        
      write_data(FILE *fp, double *d, char *ba, double bi, int n, int pad):
	It writes data stored in 'd' array to file referenced by file pointer 'fp'. 
	'ba' is alphabet representation of parameter for which data is to be written in a row.
	'bi' is value of parameter for which data is to be written in a row.
	'n' is no. of columns of data to be written or stored in 'd' array.
	'pad' can be 1 or 0. Use 1 if no. of columms of data for each row is different. It pads columns with zero to max possible columns; for example in case, data is collected for different individuals in a group, group size may vary and number of columns of data in different row is different.
	
      calc_data(FILE *fp, double *d, int n, bool toNormal):
	It reads data from file pointed by file pointer 'fp' and stores it in 'd' array.
	'n' is columns to be read from file.
	'toNormal' is option to normalize data; 1 or 0.
	
      hist_sum(FILE *gp, char *edata, int n, char *title, char *outputfile, int wid):
	It plots histogram graphs rowstacked or cluster depending upon parameter value of 'STACKED'.
	'gp' is file pointer to gnuplot. 
	'edata' is initial string of file name of data file to be plotted.
	'n' is number of individuals in a group or group size and also accounts for number of columns of data in a row.
	'title' is title of a graph.
	'outputfile' is file name to which graph is saved as a file.
	'wid' is index used when graphs are generated as window pop ups instead of storing them as file.
	
      
    
