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
    telomere_data_dict = generate_dictionary_for_telomere_length_data(directory)
    Where the directory contains images of files containing telomere length data in
    a predefined format. This function is written specifically for the Excel file templates
    that I use, and will provide in this repository, but could be altered for any format.
   
    The individual telomere lengths column is extracted, cleansed of missing values & DAPI-intensity 
    values; outliers (3 std devs from mean of column) are removed; and the telomere length values are 
    standardized to each other by use of fluorescent beads which calibrate according to inherent 
    differences between microscope imaging sessions. The individual's ID & timepoint (from filename) (KEY) 
    is associated with its respective individual telomere length data (VALUE) as a KEY:VALUE pair 
    in the dictionary. The dictionary can then be looped over to initialize all timepoint data
    for that individual for analysis, i.e visualizations, statistics, etc.
    
    Example function with types documented in the docstring.

    `PEP 484`_ type annotations are supported. If attribute, parameter, and
    return types are annotated according to `PEP 484`_, they do not need to be
    included in the docstring:

    Args:
        param1 (int): The first parameter.
        param2 (str): The second parameter.

    Returns:
        bool: The return value. True for success, False otherwise.

    .. _PEP 484:
        https://www.python.org/dev/peps/pep-0484/
    """
    
    # initialize dictionary to hold our data
    dict_dfs = {}
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
            DAPI_values_to_drop=[5, 192, 379, 566, 753, 940, 1127, 1314, 1501, 1688,] 
                                 
                                 
#                                  1875, 2062,
#                                  2249, 2436, 2623, 2810, 2997, 3184, 3371, 3558, 3745, 3932, 4119, 
#                                  4306, 4493, 4680, 4867, 5054, 5241, 5428]
            
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
            dict_dfs[file_name_trimmed] = telos_individ_df
            
    if status == 0:
        print('Issue finding excel files.. please check your directory')
    
    elif status == 1:
        print('Done collecting all telomere length excel files')
    
    return dict_dfs



def make_dataframe_from_telomere_data_dict(telo_dict, timepoint_list):
    """
    fxn

    Args:
        big_table: An open Bigtable Table instance.
        keys: A sequence of strings representing the key of each table row
        to fetch .. other_silly_variable: Another optional variable, that has a much
        longer name than the other args, and which does nothing.

    Returns:
        A dict mapping keys to the corresponding data..

                {'Serak': ('Rigel VII', 'Preparer'),
                 'Zim': ('Irk', 'Invader'),
                 'Lrrr': ('Omicron Persei 8', 'Emperor')}

        If a key from the keys argument is missing from the dictionary,
        then that row was not found in the table.

    Raises:
        IOError: An error occurred accessing the bigtable.Table object.
    """
    data = []
    
    for name_key, telo_value in telo_dict.items():
        id_num = name_key[3:7]
        time_point = get_timepoint(name_key, timepoint_list)
        telo_value = gen_missing_values_and_impute_or_randomsampledown(30, 184, pd.Series(telo_value.values.reshape(-1,)))
        data.append([id_num, time_point, telo_value, np.mean(telo_value.values)])

    df = pd.DataFrame(data, columns = ['sample id', 'timepoint', 'telo data', 'telo means'])

    df['timepoint'] = df['timepoint'].astype('category')


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
        
        
        
def gen_missing_values_and_impute_or_randomsampledown(n_cells, telosPercell, df):
    #if wanted to do for max. possible telomeres, just replace the subtraction with max telos
    # print('substracts second astro from first.. equalizing second to first')

    if df.size > 5520:
        df_sampled = df.sample(5520)
        return df_sampled

    if df.size > 25 and df.size <= 2760:
        missing_data_difference = abs( (n_cells * telosPercell) - df.size )
        sampled = df.sample(missing_data_difference, replace=True, random_state=28)
        concat_ed = pd.concat([sampled, df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        concat_ed.reset_index(drop=True, inplace=True)
        return concat_ed

    if df.size > 25 and df.size < 5520:
        missing_data_difference = abs( (n_cells * telosPercell) - df.size )
        sampled = df.sample(missing_data_difference, random_state=28)
        concat_ed = pd.concat([sampled, df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        concat_ed.reset_index(drop=True, inplace=True)
        return concat_ed

    else:
        print(f'unable to standardize individual telomere counts.. please check your dataframe for errors')
        return df
    
    
    
####################################################################################################

# FUNCTIONS FOR CONVERTING MANIPULATING INDIVIDUAL TELOMERE LENGTH DATA IN
# A PANDAS DATAFRAME

# Functions for exploding individual telomere length measurements into a rows within a dataframe
# and counting short/medium/long telomeres

####################################################################################################



def calculate_apply_teloQuartiles_dataframe(df, ordered_timepoint_list):
    
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