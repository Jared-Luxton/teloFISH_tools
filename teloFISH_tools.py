# enables access to directories/files
import os

# for handling data
import numpy as np
from numpy import array
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile

# graphing
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import seaborn as sns

# for loading telo data column containing individual
# telomere length values
from ast import literal_eval

# statistics, numbers
from scipy import stats
import math



####################################################################################################

# FUNCTIONS FOR EXTRACTING INDIVIDUAL TELOMERE LENGTH MEASUREMENTS TO DATAFRAME

# Functions for extracting individual telomere length measurements from excel files,
# where the individual telomere length measurements are contained in a column with one
# foci measurement per cell, into a pandas dataframe. 
# 
# These functions are intended for longitudinal data.

####################################################################################################



def generate_dictionary_for_telomere_length_data(patharg):
  
    """
    USAGE:
    telomere_data_dict = generate_dictionary_for_telomere_length_data(path/to/telomere_data/directory)

    The directory should contain Excel files containing individual telomere length data in
    a predefined format. This function is written specifically for the Excel file templates
    that our lab uses (provided in this repository) but could be altered for any format.
   
    The function loops through the data files and extracts the column containing individual telomere length
    measurements, then removes missing values & DAPI-intensity values, and outliers (3 std devs from mean of column).
    The sample's filename (which should contain sample ID & timepoint information) is stored in a dictionary as a KEY,
    with it's corresponding VALUE being the telomere data.

    Args:
        patharg (PATH): 
        Path to directory containing telomere length data

    Returns:
        telomere_dict: 
        dictionary containing filenames as KEYs and corresponding telomere data as VALUES

        dict: 
        	  {sample1: SAMPLE1_telomere_data,
          	   sample2: SAMPLE2_telomere_data,
          	   etc.}


    """
    
    # initialize dictionary to hold our data
    telomere_dict = {}
    status = 0

    # loop through directory to grab files
    for file in os.scandir(patharg):
        if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
            print(f'{file.name} telomere data acquisition in progress..') 
            try:
                df = pd.read_excel(file)
                status = 1
            except:
                print(f'{file.name} File issue loading..')
                return -1

            # grabbing individual telomere length data from the file 
            individual_telos_lengths = df.iloc[:, 3]
            
            # these numbers correspond to rows containing information about the DAPI counterstain, NOT telomeres, so we drop
            DAPI_values_to_drop=[5, 192, 379, 566, 753, 940, 1127, 1314, 1501, 1688, 1875, 2062,
                                 2249, 2436, 2623, 2810, 2997, 3184, 3371, 3558, 3745, 3932, 4119, 
                                 4306, 4493, 4680, 4867, 5054, 5241, 5428]
            
            # drop DAPI info, grab only telo measurements in the col
            # enforce the telomere measurements as a numeric data type
            individual_telos_lengths = individual_telos_lengths.drop(labels=DAPI_values_to_drop)
            individual_telos_lengths = individual_telos_lengths.iloc[7:5611]
            telos_str_toNaN = pd.to_numeric(individual_telos_lengths, errors='coerce')
            
            # drop any missing values, make data into a dataframe
            individual_telos_cleaned = telos_str_toNaN.dropna(axis=0, how='any')
            telos_df = individual_telos_cleaned.to_frame(name=None)
            
            # remove telomere measurements/visual artifacts that lie beyond 3 standard deviations of the mean
            # the data is relatively normal in shape, & this process removes about ~10-30 telos from ~5520
            telos_individ_df = telos_df[(np.abs(stats.zscore(telos_df)) < 3).all(axis=1)]
            file_name_trimmed = file.name.replace('.xlsx', '')
            telomere_dict[file_name_trimmed] = telos_individ_df
            
    if status == 0:
        print('Issue finding excel files.. please check your directory')
    elif status == 1:
        print('Done collecting all telomere length excel files')
    
    return telomere_dict



def make_dataframe_from_telomere_data_dict(telo_dict, ordered_timepoint_list, 
										   id_start=3, id_stop=7):
    """
    USAGE:
    telomere_dataframe = make_dataframe_from_telomere_data_dict(telomere_dict, timepoint_list)

	This function makes a dataframe from individual telomere length data and sample info from telomere_dict
	and timepoint_list. 

    Args:
        telomere_dict (dict): 
        The telomere_dict is returned from the generate_dictionary_for_telomere_length_data() 
        function, and is a dictionary which has KEYs as sample filenames, which contain sample ID & timepoint info, 
		and VALUEs which contain individual telomere lenght data. 

		ordered_timepoint_list (list):
		An ordered (first to last) list of all sample timepoints associated with the telomere data in telomere_dict.

		id_start (int):
		Location of ID start position in sample file name KEYs in telomere_dict.

		id_stop (int):
		Location of ID stop position in sample file name KEYs in telomere_dict.

    Returns:
        df (pandas dataframe):
        Dataframe contains:
        columns: 'sample id', 'timepoint', 'telo data', 'telo means', where the sample ID is the unique identifier
        or individual, timepoint is the time sample was person, telo data is a list containing all individual telomere
        length measurements, and telo means is the MEAN of all individual telomere measurements for that sample.
    """
    data = []
    
    for name_key, telo_value in telo_dict.items():
        id_num = name_key[id_start:id_stop]
        time_point = get_timepoint(name_key, ordered_timepoint_list)
        telo_value = standardize_individual_telomere_counts(30, 184, pd.Series(telo_value.values.reshape(-1,)))
        data.append([id_num, time_point, telo_value, np.mean(telo_value.values)])

    df = pd.DataFrame(data, columns = ['sample id', 'timepoint', 'telo data', 'telo means'])

    df['Q1'] = 'telos quartile 1 <0.25'
    df['Q2-3'] = 'telos quartile 2-3 >0.25 & <0.75'
    df['Q4'] = 'telos quartile 4 >0.75'
    df['timepoint'] = df['timepoint'].astype('category')
    df = df.sort_values(['sample id', 'timepoint']).reset_index(drop=True)
    
    return df



def get_timepoint(name_key, ordered_timepoint_list):
    
    for timepoint in ordered_timepoint_list:
        if timepoint in name_key:
            cleaned_timepoint = name_key[-len(timepoint):]
            return cleaned_timepoint.strip()
        
        
        
def standardize_individual_telomere_counts(n_cells, telos_per_cell, telo_data):
	"""
    USAGE:
    df = standardize_individual_telomere_counts(30, 184, telo_data)

    Standardizes the counts of individual telomeres in a sample to that sample's theoretical # of 
    individual telomeres. For samples with fewer telomeres than the theoretical #, the missing values
    are imputed by randomly sampling telomeres from the observed distribution, up to the theoretical #
    of individual telomeres. For samples with more telomeres than the theoretical #, telomeres from 
    the observed distribution are randomly removed.


    Args:
    	n_cells (int):
    	Number of cells counted per sample.

    	telos_per_cell (int):
    	Number of telomeres per cell (species/cell-type specific).

    	telo_data (np.array):
    	Array of individual telomeres.

    Returns:
    	sample_down (pandas dataframe) or 
    	imputed_sample (pandas dataframe),

    	Which are dataframes containing the theoretical # of telomeres for a given number of cells
    	and given number of telomeres per cell, specified by the user, where the theoretical # is achieved
    	by sampling down (sample_down) or sampling up (imputed_sample) using the observed distribution of telomeres
    	for that sample.
    """

    if telo_data.size > 5520:
        sample_down = telo_data.sample(5520)
        return sample_down

    # enables multiple sampling of a given telomere measurement
    if telo_data.size > 25 and telo_data.size <= 2760:
    	# calculating number of missing telomeres & randomly sampling that number
        missing_data_difference = abs((n_cells * telos_per-cell) - telo_data.size)
        imputed_telos = telo_data.sample(missing_data_difference, replace=True, random_state=28)

        # concatenating telomeres imputed by random sampling to original data
        imputed_sample = pd.concat([imputed_telos, telo_data], sort=False)
        np.random.shuffle(imputed_sample.to_numpy())
        imputed_sample.reset_index(drop=True, inplace=True)
        return imputed_sample

    # does not allow multiple sampling of a given measurement
    if telo_data.size > 25 and telo_data.size < 5520:
    	# calculating number of missing telomeres & randomly sampling that number
        missing_data_difference = abs((n_cells * telos_per_cell) - telo_data.size)
        imputed_telos = telo_data.sample(missing_data_difference, random_state=28)

         # concatenating telomeres imputed by random sampling to original data
        imputed_sample = pd.concat([imputed_telos, telo_data], sort=False)
        np.random.shuffle(imputed_sample.to_numpy())
        imputed_sample.reset_index(drop=True, inplace=True)
        return imputed_sample

    else:
        print(f'unable to standardize individual telomere counts.. please check your dataframe for errors')
        return df
    
    
    
####################################################################################################

# FUNCTIONS FOR CONVERTING MANIPULATING INDIVIDUAL TELOMERE LENGTH DATA IN
# A PANDAS DATAFRAME

# Functions for exploding individual telomere length measurements into a rows within a dataframe
# and counting short/medium/long telomeres

####################################################################################################



def feature_engineer_short_long_telomeres(df, ordered_timepoint_list):
	"""
	USAGE:
    df = feature_engineer_short_long_telomeres(df, ordered_timepoint_list)


    Args:


    Returns:
	"""
    
    q1_row, q2_3_row, q4_row = 'Q1', 'Q2-3', 'Q4'

    df['timepoint'] = df['timepoint'].astype('category')
    df['timepoint'].cat.set_categories(ordered_timepoint_list, inplace=True)
    df = df.sort_values(['sample id', 'timepoint']).reset_index(drop=True)
    
    for i, row in df.iterrows():

        if ordered_timepoint_list[0] in row['timepoint']:
            init_timepoint = row['telo data']
            df.at[i, q1_row], df.at[i, q2_3_row], df.at[i, q4_row] = (quartile_cts_rel_to_df1(init_timepoint, init_timepoint))

        elif ordered_timepoint_list[0] not in row['timepoint']:
            df.at[i, q1_row], df.at[i, q2_3_row], df.at[i, q4_row] = (quartile_cts_rel_to_df1(init_timepoint, row['telo data']))

        else:
            print('unknown label in row["timepoint"] of dataframe.. please check timepoint names')
            
    return df



def order_timepoint_col(df, ordered_timepoint_list):
    df['timepoint'] = df['timepoint'].astype('category')
    df['timepoint'].cat.set_categories(ordered_timepoint_list, inplace=True)
    df = df.sort_values(['sample id', 'timepoint']).reset_index(drop=True)
    
    return df



def quartile_cts_rel_to_df1(d_1, d_2):
    d_1 = pd.DataFrame(d_1)
    d_2 = pd.DataFrame(d_2)
    
    quartile_1 = d_2[d_2 <= d_1.quantile(0.25)].count()
    quartile_2_3 = d_2[(d_2 > d_1.quantile(0.25)) & (d_2 < d_1.quantile(0.75))].count()
    quartile_4 = d_2[d_2 >= d_1.quantile(0.75)].count()
    
    return quartile_1.values, quartile_2_3.values, quartile_4.values



def explode_individual_telos(df):
    exploded_telos = (df['telo data'].apply(pd.Series)
        .merge(df, right_index = True, left_index = True)                  
        .drop('telo data', axis=1)
        .melt(id_vars = [col for col in df.columns if col != 'telo data'], value_name = "individual telos") 
        .drop("variable", axis = 1)
        .dropna())
    
    exploded_telos['individual telos'] = exploded_telos['individual telos'].astype('int64')
    
    return exploded_telos



####################################################################################################

# FUNCTIONS FOR GRAPHING INDIVIDUAL TELOMERE LENGTH DATA

####################################################################################################


# patient_ids = list(all_patients_df['patient id'].unique())
# telo_mrp.histogram_plot_groups(x='telo data exploded', data=exploded_telos_all_patients_df, groupby='patient id', iterable=patient_ids)


def histogram_plot_groups(x=None, df=None, sample_id_col=None, groupby=None, 
                          num_samps_per_group=None, ordered_timepoint_list=None, n_bins=40):
    
    all_samples = list(df[sample_id_col].unique())
    group_df = df.groupby(groupby)
    
    for sample in all_samples:
        timepoint_telo_values_dict = {}
        
        plot_df = group_df.get_group(sample).copy()
        plot_df = order_timepoint_col(plot_df, ordered_timepoint_list)

        sample_unique_timepoints = list(plot_df['timepoint'].unique())
        
        for timepoint in sample_unique_timepoints:
            
            sample_ID_complete = str(list(plot_df[sample_id_col].unique())[0]) + ' ' + str(timepoint)
            timepoint_telo_values_dict[sample_ID_complete] = plot_df[plot_df['timepoint'] == timepoint][x]
            
        if len(timepoint_telo_values_dict.keys()) > 1:
            
            n_rows = len(timepoint_telo_values_dict.keys()) / 2
            n_rows = int(math.ceil(n_rows))
  
            n_bins = n_bins
            fig, axes = plt.subplots(n_rows, 2, sharey=True, sharex=True, constrained_layout=True, figsize = (12, 8))
            sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})

            for ax, item in zip(axes.flatten(), timepoint_telo_values_dict.items()):

                name, data = item[0], item[1]
                initial_timepoint = timepoint_telo_values_dict[list(timepoint_telo_values_dict.keys())[0]]
                histogram_stylizer_divyBins_byQuartile(fig, ax, n_bins, data, initial_timepoint, f'{name}')
                
                

def histogram_stylizer_divyBins_byQuartile(fig, ax, n_bins, data, initial_timepoint, name):

    data = data.to_numpy()
    initial_timepoint = initial_timepoint.to_numpy()

    N, bins, patches = ax.hist(data, bins=n_bins, edgecolor='black')

    for a in range(len(patches)):
        if bins[a] <= np.quantile(initial_timepoint, 0.25):
            patches[a].set_facecolor('#fdff38')

        elif np.quantile(initial_timepoint, 0.25) < bins[a] and bins[a] <= np.quantile(initial_timepoint, 0.50):
            patches[a].set_facecolor('#d0fefe')

        elif np.quantile(initial_timepoint, 0.50) < bins[a] and bins[a] <= np.quantile(initial_timepoint, 0.75):
            patches[a].set_facecolor('#d0fefe')

        elif bins[a] > np.quantile(initial_timepoint, 0.75): 
            patches[a].set_facecolor('#ffbacd')
            
    ax.set_title(f"{name}", fontsize=18,)
    ax.tick_params(labelsize=12)
#     ax.xaxis.set_major_locator(plt.MaxNLocator(12))



def order_timepoint_col(df, ordered_timepoint_list):
    df['timepoint'] = df['timepoint'].astype('category')
    df['timepoint'].cat.set_categories(ordered_timepoint_list, inplace=True)
    df = df.sort_values(['sample id', 'timepoint']).reset_index(drop=True).copy()
    
    return df