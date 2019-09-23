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

# statistics
from scipy import stats

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
            
            # these numbers correspond to rows containing information about the DAPI counterstain, NOT telomeres, so we drop
            DAPI_values_to_drop=[5, 192, 379, 566, 753, 940, 1127, 1314, 1501, 1688, 1875, 2062,
                    2249, 2436, 2623, 2810, 2997, 3184, 3371, 3558, 3745, 3932, 4119, 4306, 4493, 
                    4680, 4867, 5054, 5241, 5428]

            # grabbing individual telomere length data from the file & dropping DAPI info
            individual_telos_lengths = df.iloc[:, 3]
            individual_telos_lengths = individual_telos_lengths.drop(labels=DAPI_values_to_drop)
            
            # first pass at generating synthetic data for github exposition; to initialize actual
            # data, comment out the line below, and uncomment the .iloc[] line
#             individual_telos_lengths = individual_telos_lengths.sample(2500, random_state=1)
            individual_telos_lengths = individual_telos_lengths.iloc[7:5611]

            # ensure the telomere measurements are a numeric data type, drop any missing values, 
            # make data into a dataframe
            telos_str_toNaN = pd.to_numeric(individual_telos_lengths, errors='coerce')
            individual_telos_cleaned = telos_str_toNaN.dropna(axis=0, how='any')
            telos_df = individual_telos_cleaned.to_frame(name=None)
            
            # remove any telomere measurements that lie beyond 3 standard deviations of the mean
            # the data is relatively normal in shape, & this process removes about ~10-20 telos from ~5520
            # modest loss, acceptable to help standardize
            telos_individ_df = telos_df[(np.abs(stats.zscore(telos_df)) < 3).all(axis=1)]
            file_name_trimmed = file.name.replace('.xlsx', '')
            dict_dfs[file_name_trimmed] = telos_individ_df
            
    if status == 0:
        print('Issue finding excel files.. please check your directory')
    
    elif status == 1:
        print('Done collecting all telomere length excel files')
    
    return dict_dfs



def make_dataframe_from_telomere_data_dict(telo_dict):
    data = []
    
    for name_key, telo_value in telo_dict.items():
        id_num = name_key[3:7]
        time_point = get_timepoint(name_key)
        telo_value = gen_missing_values_and_impute_or_randomsampledown(30, 184, pd.Series(telo_value.values.reshape(-1,)), 'rsamp')

        data.append([id_num, time_point, telo_value, np.mean(telo_value.values)])

    df = pd.DataFrame(data, columns = ['sample id', 'timepoint', 'telo data', 'telo means'])

    df['timepoint'] = df['timepoint'].astype('category')


    df['Q1'] = 'telos quartile 1 <0.25'
    df['Q2-3'] = 'telos quartile 2-3 >0.25 & <0.75'
    df['Q4'] = 'telos quartile 4 >0.75'

    df['timepoint'] = df['timepoint'].astype('category')
    df = df.sort_values(['sample id', 'timepoint']).reset_index(drop=True)
    
    return df



def get_timepoint(name_key):
    timepoint_5_char = ['L-270', 'L-180', 'FD140', 'FD260', 'R+105', 'R+180', 'R+270']
    timepoint_4_char = ['L-60', 'FD45', 'FD90', 'R+60']
    timepoint_3_char = ['R+5', 'R+7'] 
    
    for timepoint in timepoint_5_char:
        if timepoint in name_key:
            timepoint = name_key[-5:]
            return timepoint.strip()
    
    for timepoint in timepoint_4_char:
        if timepoint in name_key:
            timepoint = name_key[-4:]
            return timepoint.strip()
            
    for timepoint in timepoint_3_char:
        if timepoint in name_key:
            timepoint = name_key[-3:]
            return timepoint.strip()

        
        
def gen_missing_values_and_impute_or_randomsampledown(n_cells, telosPercell, astro_df, option=None):
    #if wanted to do for max. possible telomeres, just replace the subtraction with max telos
    # print('substracts second astro from first.. equalizing second to first')

    if astro_df.size > 4600:
        astro_dfsampled = astro_df.sample(4600)
        return astro_dfsampled

    if astro_df.size > 25 and astro_df.size <= 2300:
        missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )
        rsampled = astro_df.sample(missing_data_difference, replace=True, random_state=28)
        concat_ed = pd.concat([rsampled, astro_df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        concat_ed.reset_index(drop=True, inplace=True)
        return concat_ed

    if astro_df.size > 25 and astro_df.size < 4600:
        missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )
        if option == 'rsamp':
            rsampled = astro_df.sample(missing_data_difference, random_state=28)
            concat_ed = pd.concat([rsampled, astro_df], sort=False)
            np.random.shuffle(concat_ed.to_numpy())
            concat_ed.reset_index(drop=True, inplace=True)
            return concat_ed
        else:
            return astro_df
    else:
        return astro_df
    
    
    
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