# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 08:17:22 2023

@author: mudoe
"""

#%% imports
import pandas as pd
import numpy as np
import os
import copy
os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")
from metatlas.io import feature_tools as ft
import math
from scipy.signal import argrelextrema

#%%data
os.chdir("C:/Users/mudoe/Desktop/InverSIL/high res/")
df=pd.read_csv("56861213.csv") #Define mzMine2 alignment table
df = df.dropna(axis = 1)
df.columns = pd.Series(['RT','MZ_12C','13'])

#%%nan
df = df.replace(0, np.nan)

#%%drop equal 12 & 13
align_12C_13C = df[df['MZ_12C'].notnull() & df['13'].notnull()].index 
df.drop(align_12C_13C, inplace = True) 
#df.shape
    
#%%Make individual feature lists for each condition
tw = df[df["MZ_12C"].notnull()].loc[:,["RT","MZ_12C"]]
th = df[df["13"].notnull()].loc[:,["RT","13"]]
#%%mz offset params
mz_offset = 4 #Define mass offset [m/z]
mz_tolerance = 1    # Define m/z tolerance [m/z]
rt_tolerance = 0.1  # Define RT tolerance [min]

#%%combo stuff

tw["_temp"] = 0
th["_temp"] = 0
combined = th.merge(tw, 'outer', '_temp', suffixes=['_13C', '_12C'])

#featurelist_12C['_temp'] = 0
#featurelist_13C['_temp'] = 0
#combined = featurelist_13C.merge(featurelist_12C, 'outer', '_temp', suffixes=['_13C', '_12C'])

#%%pare 12C and 13C
combined = combined[
    ((combined['13']  - combined['MZ_12C']) > mz_tolerance) &
    ((combined['13']  - combined['MZ_12C']) <= 51) &
    ((combined['RT_13C'] - combined['RT_12C']).abs() < rt_tolerance)]
combined = combined.drop(['_temp'], axis=1)
combined.sort_values(by=['13'], inplace=True)
combstore = copy.deepcopy(combined)

#%%crunch the pares
trt = combined["RT_12C"].unique()
for i in range(len(combined)):
    for j in range(i+1,len(combined)):
        if combined.iat[i,0] != 0:
            if combined.iat[i,3] == combined.iat[j,3]:
                if abs(combined.iat[i,0] - combined.iat[j,0]) < .5:
                    if combined.iat[j,1] - combined.iat[i,1] < 2:
                        combined.iat[i,0] = 0
combined = combined.replace(0, np.nan)
combined = combined.dropna(axis = 0)

#%%generic
data = {}
data["atlas"] = 0
data["lcmsrun"] = "5686P1.h5" #placeholder for the file location. unsure if I'm gonna make this dynamic
data["file_index"] = 1
data["polarity"] = "positive"
ppm = 3000 #maybe 5k just to be safe?

#%%neutron steps and getting data
import time
start_time = time.time()
os.chdir("C:/Users/mudoe/Desktop/InverSIL/high res/5686")
incorp = []
for i in range(len(combined)):
    steps = math.ceil(combined.iat[i,1] - combined.iat[i,3])
    subst = []
    RT = (combined.iat[i,0] + combined.iat[i,2])/2
    for j in range(-1,steps+1):
        subst.append([combined.iat[i,3] + (j*1.003),RT,RT+.1,RT-.1,ppm,j+2,0,j+3,j+1])
    subst = pd.DataFrame(subst)
    subst = subst.rename(columns={0:"mz",1:"rt_peak",2:"rt_max",3:"rt_min",4:"ppm_tolerance",
                                  5:"group_index",6:"extra_time",7:"label",8:"index"})
    data.update({"atlas": subst})
    d = ft.get_data(data,return_data=True,save_file=False)
    intensity = np.array(d["ms1_summary"]["peak_height"])
    imax = argrelextrema(intensity,np.greater)[0]
    d["ms1_summary"] = d["ms1_summary"].loc[imax,:]
    d["ms1_summary"].sort_values(by=['peak_height'], inplace=True, ascending = False)
    incorp.append([combined.iat[i,3],combined.iat[i,2],combined.iat[i,1],combined.iat[i,0]])
    incorp.append(d["ms1_summary"].iloc[:,4])
    incorp.append(d["ms1_summary"].iloc[:,3])
    
print("--- %s seconds ---" % (time.time() - start_time))
    
#%%set up atlases
#df_standards.drop(columns=['Unnamed: 0'],inplace=True)
df_standards['rt_min'] = df_standards['rt_peak'] - 0.1
df_standards['rt_max'] = df_standards['rt_peak'] + 0.1
df_standards['ppm_tolerance'] = ppm_tolerance
df_standards['extra_time'] = extra_time

df_standards['group_index'] = ft.group_consecutive(df_standards['mz'].values[:],
                                         stepsize=ppm_tolerance,
                                         do_ppm=True)

data_list = []
for i,row in df_files.iterrows():
    data_setup = {}
    data_setup['lcmsrun'] = row['filename']
    # data_setup['ppm_tolerance'] = ppm_tolerance
    # data_setup['extra_time'] = 0
    data_setup['atlas'] = df_standards
    data_setup['file_index'] = int(os.path.basename(row['filename']).split('_')[-1].replace('.h5','').replace('Run',''))
    data_setup['polarity'] = "positive"
    data_list.append(data_setup)
    
test = np.array([1,2,3,4,2,1,2,1])
testtwo = argrelextrema(test,np.greater)[0]

test = [[1,2,8,9],[3,4],[5,6,7]]
test = pd.DataFrame(test)
