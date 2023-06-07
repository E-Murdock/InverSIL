# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 08:28:31 2023

@author: mudoe
"""

#%% imports
import time
start_time = time.time()
import pandas as pd
import numpy as np
import os
import copy
os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")
from metatlas.io import feature_tools as ft
import math
from scipy.signal import argrelextrema
import glob

#%%files list
os.chdir("c:/users/mudoe/desktop/InverSIL/high res/co")
samples = glob.glob(os.path.join('*.csv'))
samples = pd.DataFrame(samples)
samples[0] = samples[0].apply(lambda x: os.path.basename(x).split('_')[0])

#%%combos
combos = {}
for n in range(len(samples)):
    df = pd.read_csv(samples.iat[n,0] + "_1213.csv")
    df = df.dropna(axis = 1)
    df.columns = pd.Series(['RT','MZ_12C','13'])
    df = df.replace(0, np.nan)
    align_12C_13C = df[df['MZ_12C'].notnull() & df['13'].notnull()].index 
    df.drop(align_12C_13C, inplace = True) 
    tw = df[df["MZ_12C"].notnull()].loc[:,["RT","MZ_12C"]]
    th = df[df["13"].notnull()].loc[:,["RT","13"]]
    mz_tolerance = 1    # Define m/z tolerance [m/z]
    rt_tolerance = 0.1  # Define RT tolerance [min]
    tw["_temp"] = 0
    th["_temp"] = 0
    combined = th.merge(tw, 'outer', '_temp', suffixes=['_13C', '_12C'])
    combined = combined[
        ((combined['13']  - combined['MZ_12C']) > mz_tolerance) &
        ((combined['13']  - combined['MZ_12C']) <= 51) &
        ((combined['RT_13C'] - combined['RT_12C']).abs() < rt_tolerance)]
    combined = combined.drop(['_temp'], axis=1)
    combined.sort_values(by=['13'], inplace=True)
    combstore = copy.deepcopy(combined)
    trt = combined["RT_12C"].unique()
    for i in range(len(combined)):
        for j in range(i+1,len(combined)):
            if combined.iat[i,0] != 0:
                if combined.iat[i,3]-(combined.iat[i,3]/200000) < combined.iat[j,3] < combined.iat[i,3]+(combined.iat[i,3]/200000):
                    if abs(combined.iat[i,2] - combined.iat[j,2]) < .05:
                        #if combined.iat[j,1] - combined.iat[i,1] < 2:
                            combined.iat[i,0] = 0
    combined = combined.replace(0, np.nan)
    combined = combined.dropna(axis = 0)
    combined.sort_values(by=['MZ_12C'], inplace=True,ascending=False)
    for i in range(len(combined)):
        for j in range(i+1,len(combined)):
            if combined.iat[i,0] != 0:
                if combined.iat[i,1]-(combined.iat[i,1]/200000) < combined.iat[j,1] < combined.iat[i,1]+(combined.iat[i,1]/200000):
                    if abs(combined.iat[i,0] - combined.iat[j,0]) < .05:
                        #if combined.iat[j,1] - combined.iat[i,1] < 2:
                            combined.iat[i,0] = 0
    combined = combined.replace(0, np.nan)
    combined = combined.dropna(axis = 0)
    combined.sort_values(by=['13'], inplace=True)
    combos[samples.iat[n,0]] = combined

#%%prep
data = {}
data["atlas"] = 0
data["lcmsrun"] = 0 #placeholder for the file location
data["file_index"] = 1
data["polarity"] = "positive"
ppm = 15 #maybe 5k just to be safe?

data12 = {}
data12["atlas"] = 0
data12["lcmsrun"] = 0 #placeholder for the file location
data12["file_index"] = 2
data12["polarity"] = "positive"
ppm2 = 15

data13 = {}
data13["atlas"] = 0
data13["lcmsrun"] = 0
data13["file_index"] = 3
data13["polarity"] = "positive"
ppm2 = 15

#%%loop to check pairs
collection = {}
for n in range(len(samples)):
    key = list(combos.keys())[n]
    combined = copy.deepcopy(combos.get(key))
    data12['lcmsrun'] = key + '_12.h5'
    data13['lcmsrun'] = key + '_13.h5'
    twelvers = []
    for i in range(len(combined)):
        twelvers.append([combined.iat[i,3],combined.iat[i,2],combined.iat[i,2]+.1,combined.iat[i,2]-.1,
                    ppm2,i+1,0,i+2,i+3]) #load in 12
        twelvers.append([combined.iat[i,1],combined.iat[i,0],combined.iat[i,0]+.1,combined.iat[i,0]-.1,
                    ppm2,len(combined)+i+2,0,len(combined)+i+3,len(combined)+i+4]) #load in 13
    twelvers = pd.DataFrame(twelvers)
    twelvers = twelvers.rename(columns={0:"mz",1:"rt_peak",2:"rt_max",3:"rt_min",4:"ppm_tolerance",
                                      5:"group_index",6:"extra_time",7:"label",8:"index"})
    data12.update({"atlas":twelvers})
    data13.update({"atlas":twelvers})
    d12 = ft.get_data(data12,return_data=True,save_file=False)
    d13 = ft.get_data(data13,return_data=True,save_file=False)
    twelvers.sort_values(by=['label'], inplace=True)
    for i in range(int(len(twelvers)/2)):  
        int13 = d13["ms1_summary"].loc[d13["ms1_summary"]["label"] == (twelvers.iat[i,7] + 1 + len(combined)),:] #13 in 13
        int12 = d12["ms1_summary"].loc[d12["ms1_summary"]["label"] == twelvers.iat[i,7],:] #12 in 12
        twinth = d13["ms1_summary"].loc[d13["ms1_summary"]["label"] == twelvers.iat[i,7],:] #12 in 13
        thintw = d12["ms1_summary"].loc[d12["ms1_summary"]["label"] == (twelvers.iat[i,7] + 1 + len(combined)),:] #13 in 12
        x = 0
        y = 0
        z = 0
        if len(twinth) == 1:
            if (twinth.iat[0,3]/int13.iat[0,3]) > .5:
                x = 1
        if len(thintw) == 1:
            if (thintw.iat[0,3]/int12.iat[0,3]) > .5:
                y = 1
        if len(int13) == 1 and len(int12) == 1:
            bigger = max(int13.iat[0,3],int12.iat[0,3])
            smaller = min(int13.iat[0,3],int12.iat[0,3])
            if bigger/smaller > 5:
                z = 1
        else:
            z = 1
        if (x+y+z) > 0:
            combined.iat[i,0] = np.nan
    combined = combined.dropna(axis = 0)
    combos[key] = combined
print("--- %s seconds ---" % (time.time() - start_time))

#%%neutron steps and the posts?
import time
start_time = time.time()
#os.chdir("C:/Users/mudoe/Desktop/InverSIL/high res/OBBP")
specs = {}
choice = {}
for n in range(len(samples)):
    key = list(combos.keys())[n]
    combined = copy.deepcopy(combos.get(key))
    specs[key] = {}
    choice[key] = {}
    data12['lcmsrun'] = key + '_12.h5'
    data13['lcmsrun'] = key + '_13.h5'
    for o in [1,2,3]:
        data['lcmsrun'] = key + '_P' + str(o) + ".h5"
        specs[key]['P' + str(o)] = {}
        choice[key]['P' + str(o)] = {}
        for i in range(len(combined)):
            steps = math.ceil(combined.iat[i,1] - combined.iat[i,3])
            subst = []
            RT = (combined.iat[i,0] + combined.iat[i,2])/2
            for j in range(-1,steps+1):
                subst.append([combined.iat[i,3] + (j*1.003586),RT,RT+.1,RT-.1,ppm,j+2,0,j+3,j+1])
            subst = pd.DataFrame(subst)
            subst = subst.rename(columns={0:"mz",1:"rt_peak",2:"rt_max",3:"rt_min",4:"ppm_tolerance",
                                          5:"group_index",6:"extra_time",7:"label",8:"index"})
            data12.update({"atlas": subst})
            data13.update({"atlas": subst})
            data.update({"atlas": subst})
            d = ft.get_data(data,return_data=True,save_file=False)
            d12 = ft.get_data(data12,return_data=True,save_file=False)
            d13 = ft.get_data(data13,return_data=True,save_file=False)
            intensity = np.array(d["ms1_summary"]["peak_height"])
            imax =  intensity.argmax()##different than old version, reporting global max not locals
            #d["ms1_summary"] = d["ms1_summary"].loc[imax,:]
            #d["ms1_summary"].sort_values(by=['peak_height'], inplace=True, ascending = False)
            d = d['ms1_summary']
            d12 = d12['ms1_summary']
            d13 = d13['ms1_summary']
            #d['12'] = d12['peak_height']
            #d['13'] = d13['peak_height']
            d['12'] = 0
            d['13'] = 0
            for l in range(len(d)):
                h12 = d12.loc[d12['label'] == d.iat[l,0],:]
                if len(h12) > 0:
                    d.iat[l,6] = h12.iat[0,3]
            for l in range(len(d)):
                h13 = d13.loc[d13['label'] == d.iat[l,0],:]
                if len(h13) > 0:
                    d.iat[l,7] = h13.iat[0,3]
            specs[key]['P' + str(o)][str([d.iat[0,5],d.iat[imax,3]])] = d[['label','mz_centroid','12','peak_height','13','rt_peak']]
            worth = 0
            for k in range(2,len(d)-2):
                if d.iat[k,3] > d.iat[k,6]*2 and d.iat[k,3] > d.iat[k,7]*2:
                    worth = 1
            if worth == 1:
                choice[key]['P' + str(o)][str([d.iat[0,5],d.iat[imax,3]])] = d[['label','mz_centroid','12','peak_height','13','rt_peak']]
                #choice[str(max(d['peak_height']))] = d[['label','mz_centroid','12','peak_height','13','rt_peak']]
print("--- %s seconds ---" % (time.time() - start_time))

import winsound
duration = 500  # milliseconds
freq = 660  # Hz
winsound.Beep(freq, duration)


#print("please figure out how to use pickle to save these dictionaries, or write some crazy loop that'll pull the dataframes out and save to csvs. whatever floats your boat")
#%%save
import os
import pickle
os.chdir('c:/users/mudoe/desktop')
with open('InverSIL.pickle', 'wb') as handle:
    pickle.dump(specs, handle, protocol=pickle.HIGHEST_PROTOCOL)

#%%load
import os
import pickle
os.chdir('c:/users/mudoe/desktop')
with open('InverSIL.pickle', 'rb') as handle:
    unserialized_data = pickle.load(handle)