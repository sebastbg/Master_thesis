#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 00:08:31 2021

@author: sebastiangrorud
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob


A = [14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14]
len(A)
#A = [14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]


path = '/Users/sebastiangrorud/matlab/Resultater/Zeta_scaling_2/MODE_data/'
#path = '/Users/sebastiangrorud/matlab/Resultater/active_slipmodes/'
file = glob.glob(path + '*.dat')
file = sorted(file)

#files = ['/Users/sebastiangrorud/matlab/Til rapport/CHIN_plotting/n_5.dat']
df = pd.read_fwf(file[0], sep='  ',header=None, widths=A)
df = df.apply(pd.to_numeric, errors='coerce')
df.fillna(0,inplace=True)
df.to_numpy
np.transpose(df)

for num in range (0,len(file)):
    test = np.loadtxt(file[num], dtype=float)
    resultat = np.zeros(len(test))
    for i in range (1,len(resultat)):
        teller = 0;
        temp = sum(test[i,:])
        res = np.zeros(96)
        res = (test[i,:])/temp
        for k in range(0,95):
            if (res[k] > 0.05):
                teller = teller +1
        resultat[i] = teller
    
    plt.bar(resultat)
    
