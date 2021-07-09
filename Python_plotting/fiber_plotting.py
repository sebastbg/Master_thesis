#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 10:23:27 2021

@author: sebastiangrorud
"""
# Slip mode plott for different 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
file = ['/Users/sebastiangrorud/GitHub_projects/ProjectWork_TaylorLin/Python_plotting/Results/fiber_data-kopi.csv']

df = pd.read_csv(file[0])
dfs = pd.melt(df, id_vars='Method')

Y =  ['Kesten', 'fibre intensity']
f1 =[1.5929,1.1315,0.8435,1.2828,1.9836,2.3862,2.4426]
f2 = [1.3430,0.9504,0.6920,1.0005,1.85358,3.13210,3.8463]
x = [1.3430,0.9504,0.6920,1.00050,1.85358,3.13210,3.8461]
f3 = [41.8585,30.5797,21.4187,14.5374,17.7305,21.3512,25.1304]
y_ax = [60, 65, 70,75,80,85,90]
# =============================================================================
# 
# =============================================================================
# dfs = pd.DataFrame(data={'mode': Y[1], 
#                           '20': [39.1749,26.6523,21.4436,15.9664,15.0546,11.6945,12.0106], 
#                           '50': [38.9946,27.0838,23.0095,16.5840,17.8644,16.3857,18.9754], 
#                           '1000':[41.8585,30.5797,21.4187,14.5374,17.7305,21.3512,25.1304]})
# dfs1 = pd.melt(dfs, id_vars = "mode")
# print(dfs1)
# =============================================================================
# =============================================================================

sns.lineplot(y_ax,f1)
sns.lineplot(y_ax,f2)
sns.lineplot(y_ax,f3)

#sns.lineplot(data=dfs1, x="mode", y="value", hue="variable")
#dfs1.value 