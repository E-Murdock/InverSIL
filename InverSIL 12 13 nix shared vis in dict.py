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
os.chdir("C:/Users/mudoe/Desktop/InverSIL/high res/co")
df=pd.read_csv("AWP114_1213.csv") #Define mzMine2 alignment table
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

#%%crunch the pairs
trt = combined["RT_12C"].unique()
for i in range(len(combined)):
    for j in range(i+1,len(combined)):
        if combined.iat[i,0] != 0:
            if combined.iat[i,3] == combined.iat[j,3]:
                if abs(combined.iat[i,0] - combined.iat[j,0]) < .05: ##try to remember what this does
                    if combined.iat[j,1] - combined.iat[i,1] < 2:
                        combined.iat[i,0] = 0
combined = combined.replace(0, np.nan)
combined = combined.dropna(axis = 0)

#%%generic
data = {}
data["atlas"] = 0
data["lcmsrun"] = "185P2.h5" #placeholder for the file location. unsure if I'm gonna make this dynamic
data["file_index"] = 1
data["polarity"] = "positive"
ppm = 15 #maybe 5k just to be safe?

data12 = {}
data12["atlas"] = 0
data12["lcmsrun"] = "18512.h5" #placeholder for the file location. unsure if I'm gonna make this dynamic
data12["file_index"] = 2
data12["polarity"] = "positive"
ppm2 = 15

data13 = {}
data13["atlas"] = 0
data13["lcmsrun"] = "18513.h5" #placeholder for the file location. unsure if I'm gonna make this dynamic
data13["file_index"] = 3
data13["polarity"] = "positive"
ppm2 = 15
os.chdir("C:/Users/mudoe/Desktop/InverSIL/high res/AWP185")
incorp = []

#%%neutron steps and the posts?
import time
start_time = time.time()
#os.chdir("C:/Users/mudoe/Desktop/InverSIL/high res/OBBP")
incorp = []
specs = {}
choice = {}
for i in range(len(combined)):
    twelvers = []
    twelvers.append([combined.iat[i,3],combined.iat[i,2],combined.iat[i,2]+.1,combined.iat[i,2]-.1,
                ppm2,1,0,2,3]) #load in 12
    twelvers.append([combined.iat[i,1],combined.iat[i,0],combined.iat[i,0]+.1,combined.iat[i,0]-.1,
                ppm2,2,0,3,4]) #load in 13
    twelvers = pd.DataFrame(twelvers)
    twelvers = twelvers.rename(columns={0:"mz",1:"rt_peak",2:"rt_max",3:"rt_min",4:"ppm_tolerance",
                                  5:"group_index",6:"extra_time",7:"label",8:"index"})
    data12.update({"atlas":twelvers})
    data13.update({"atlas":twelvers})
    d12 = ft.get_data(data12,return_data=True,save_file=False)
    d13 = ft.get_data(data13,return_data=True,save_file=False)
    un = list(set(list(d12["ms1_summary"].iloc[:,0])) & set(list(d13["ms1_summary"].iloc[:,0])))
    int13 = d13["ms1_summary"].loc[d13["ms1_summary"]["label"] == 3,:] #13 in 13
    int12 = d12["ms1_summary"].loc[d12["ms1_summary"]["label"] == 2,:] #12 in 12
    twinth = d13["ms1_summary"].loc[d13["ms1_summary"]["label"] == 2,:] #12 in 13
    thintw = d12["ms1_summary"].loc[d12["ms1_summary"]["label"] == 3,:] #13 in 12
    x = 0
    y = 0
    if len(twinth) == 1:
        if (twinth.iat[0,3]/int13.iat[0,3]) > .5:
            x = 1
    if len(thintw) == 1:
        if (thintw.iat[0,3]/int12.iat[0,3]) > .5:
            y = 1
    if (x + y) == 0:
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
        imax = argrelextrema(intensity,np.greater)[0]
        #d["ms1_summary"] = d["ms1_summary"].loc[imax,:]
        #d["ms1_summary"].sort_values(by=['peak_height'], inplace=True, ascending = False)
        d = d['ms1_summary']
        d12 = d12['ms1_summary']
        d13 = d13['ms1_summary']
        d['12'] = d12['peak_height']
        d['13'] = d13['peak_height']
        
        specs[str(max(d['peak_height']))] = d[['label','mz_centroid','12','peak_height','13','rt_peak']]
        worth = 0
        for i in range(2,len(d)-2):
            if d.iat[i,3] > d.iat[i,6]*3 and d.iat[i,3] > d.iat[i,7]*3:
                worth = 1
        if worth == 1:
            choice[str(max(d['peak_height']))] = d[['label','mz_centroid','12','peak_height','13','rt_peak']]
        #for j in range(len(d["ms1_summary"])):
        # thirts = d["ms1_summary"].loc[abs(d["ms1_summary"]["mz_centroid"] - int13.iat[0,4]) < .5,:]
        # if len(thirts) == 0:
        #     if len(d["ms1_summary"]) != 0:
        #         ratio = d["ms1_summary"].iat[0,3]/int13.iat[0,3] ##reports P1max/13
        #     else:
        #         ratio = 0
        # if len(thirts) == 1:
        #     ind = d["ms1_summary"]["label"] == thirts.iat[0,0]
        #     d["ms1_summary"].iloc[ind,0] = np.nan
        #     zac = d["ms1_summary"].dropna(axis = 0)
        #     if len(zac) == 0:
        #         ratio = 0
        #     else:
        #         ratio = (thirts.iat[0,3] + zac.iat[0,3])/int13.iat[0,3]
        # fill = []
        # fill.append(int13.iat[0,5])
        # fill.append(int13.iat[0,4])
        # fill.append(int13.iat[0,3])
        # fill.append(int12.iat[0,4])
        # fill.append(int12.iat[0,3])
        # fill.append(ratio)
        # for j in range(len(d["ms1_summary"])):
        #     fill.append(d["ms1_summary"].iat[j,4])
        #     fill.append(d["ms1_summary"].iat[j,3])
        # incorp.append(fill)
print("--- %s seconds ---" % (time.time() - start_time))

import winsound
duration = 500  # milliseconds
freq = 660  # Hz
winsound.Beep(freq, duration)

# incorpt = pd.DataFrame(incorp)
# incorpt = incorpt.rename(columns={0:"RT",1:"13mz",2:"13i",3:"12mz",4:"12i",
#                                   5:"ratio",6:"Pmz",7:"Pi",8:"etc"})
# reduced = copy.deepcopy(incorpt)
# reduced["run"] = reduced.groupby(["RT"]).RT.transform('size')
# idx = reduced.groupby(["RT"])['13mz'].transform(max) == reduced['13mz']
# reduced = reduced[idx]
# idx = reduced.groupby(["RT"])["12mz"].transform(min) == reduced['12mz']
# reduced = reduced[idx]
# column_to_move = reduced.pop("run")
# reduced.insert(0, "run", column_to_move)
    
    
    
    
#     steps = math.ceil(combined.iat[i,1] - combined.iat[i,3])
#     subst = []
#     RT = (combined.iat[i,0] + combined.iat[i,2])/2
#     for j in range(-1,steps+1):
#         subst.append([combined.iat[i,3] + (j*1.003),RT,RT+.1,RT-.1,ppm,j+2,0,j+3,j+1])
#     subst = pd.DataFrame(subst)
#     subst = subst.rename(columns={0:"mz",1:"rt_peak",2:"rt_max",3:"rt_min",4:"ppm_tolerance",
#                                   5:"group_index",6:"extra_time",7:"label",8:"index"})
#     data.update({"atlas": subst})
#     d = ft.get_data(data,return_data=True,save_file=False)
#     intensity = np.array(d["ms1_summary"]["peak_height"])
#     imax = argrelextrema(intensity,np.greater)[0]
#     d["ms1_summary"] = d["ms1_summary"].loc[imax,:]
#     d["ms1_summary"].sort_values(by=['peak_height'], inplace=True, ascending = False)
#     #for j in range(len(d["ms1_summary"])):
#     thirts = d["ms1_summary"].loc[abs(d["ms1_summary"]["mz_centroid"] - int13.iat[0,4]) < .5,:]
#     if len(thirts) == 0:
#         if len(d["ms1_summary"]) != 0:
#             ratio = d["ms1_summary"].iat[0,3]/int13.iat[0,3] ##reports P1max/13
#         else:
#             ratio = 0
#     if len(thirts) == 1:
#         ind = d["ms1_summary"]["label"] == thirts.iat[0,0]
#         d["ms1_summary"].iloc[ind,0] = np.nan
#         zac = d["ms1_summary"].dropna(axis = 0)
#         if len(zac) == 0:
#             ratio = 0
#         else:
#             ratio = (thirts.iat[0,3] + zac.iat[0,3])/int13.iat[0,3]
#     fill = []
#     fill.append(int13.iat[0,5])
#     fill.append(int13.iat[0,4])
#     fill.append(int13.iat[0,3])
#     fill.append(int12.iat[0,4])
#     fill.append(int12.iat[0,3])
#     fill.append(ratio)
#     for j in range(len(d["ms1_summary"])):
#         fill.append(d["ms1_summary"].iat[j,4])
#         fill.append(d["ms1_summary"].iat[j,3])
#     incorp.append(fill)
# incorpt = pd.DataFrame(incorp)
# incorpt = incorpt.rename(columns={0:"RT",1:"13mz",2:"13i",3:"12mz",4:"12i",
#                                   5:"ratio",6:"Pmz",7:"Pi",8:"etc"})
    


#%%neutron steps and getting data
# import time
# start_time = time.time()
# os.chdir("C:/Users/mudoe/Desktop/InverSIL/high res/OBBP")
# incorp = []
# for i in range(len(combined)):
#     steps = math.ceil(combined.iat[i,1] - combined.iat[i,3])
#     subst = []
#     RT = (combined.iat[i,0] + combined.iat[i,2])/2
#     for j in range(-1,steps+1):
#         subst.append([combined.iat[i,3] + (j*1.003),RT,RT+.1,RT-.1,ppm,j+2,0,j+3,j+1])
#     subst = pd.DataFrame(subst)
#     subst = subst.rename(columns={0:"mz",1:"rt_peak",2:"rt_max",3:"rt_min",4:"ppm_tolerance",
#                                   5:"group_index",6:"extra_time",7:"label",8:"index"})
#     data.update({"atlas": subst})
#     d = ft.get_data(data,return_data=True,save_file=False)
#     intensity = np.array(d["ms1_summary"]["peak_height"])
#     imax = argrelextrema(intensity,np.greater)[0]
#     d["ms1_summary"] = d["ms1_summary"].loc[imax,:]
#     d["ms1_summary"].sort_values(by=['peak_height'], inplace=True, ascending = False)
#     incorp.append([combined.iat[i,3],combined.iat[i,2],combined.iat[i,1],combined.iat[i,0]])
#     incorp.append(d["ms1_summary"].iloc[:,4])
#     incorp.append(d["ms1_summary"].iloc[:,3])
    
# print("--- %s seconds ---" % (time.time() - start_time))
    
#%%set up atlases
# #df_standards.drop(columns=['Unnamed: 0'],inplace=True)
# df_standards['rt_min'] = df_standards['rt_peak'] - 0.1
# df_standards['rt_max'] = df_standards['rt_peak'] + 0.1
# df_standards['ppm_tolerance'] = ppm_tolerance
# df_standards['extra_time'] = extra_time

# df_standards['group_index'] = ft.group_consecutive(df_standards['mz'].values[:],
#                                          stepsize=ppm_tolerance,
#                                          do_ppm=True)

# data_list = []
# for i,row in df_files.iterrows():
#     data_setup = {}
#     data_setup['lcmsrun'] = row['filename']
#     # data_setup['ppm_tolerance'] = ppm_tolerance
#     # data_setup['extra_time'] = 0
#     data_setup['atlas'] = df_standards
#     data_setup['file_index'] = int(os.path.basename(row['filename']).split('_')[-1].replace('.h5','').replace('Run',''))
#     data_setup['polarity'] = "positive"
#     data_list.append(data_setup)
    
# test = np.array([1,2,3,4,2,1,2,1])
# testtwo = argrelextrema(test,np.greater)[0]

# test = [[1,2,8,9],[3,4],[5,6,7]]
# test = pd.DataFrame(test)
