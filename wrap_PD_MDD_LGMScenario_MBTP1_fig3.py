#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is wrapper code for the production-diffusion code, applied to multiple diffusion domains.
Everything the user defines is input into this wrapper.

This code was modified from the original MATLab format created by Marissa Tremblay.
Author: Taylor Bourikas
Contact: tbourika@purdue.edu
Last modified: 2023.11.06

Copyright 2021, Marissa Tremblay
All rights reserved
Developed in part with funding from the National Science Foundation and the American Association for the Advancement of Science.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
"""
'''clear all
close all'''
import numpy as np
import pandas as pd
import matplotlib as mpl
from FormatNormLnD0aa import *
from MDD_retention_EDTScenario import *

"1) Provide Sample-Specific Information"
#Production parameters

dataSpe=pd.read_excel('./Data/SampleSpe_MBTP1.xlsx') #USER INPUT: (ID, lat, long, elev, Pref, shield, thick, rho, fsp, lsp, edot, ztoday)
dataSpe=np.array(dataSpe)

"2) Provide Observational Data"
#3He retention, 3He age, 10Be age

dataObs=pd.read_excel('./Data/SampleObs_MBTP1.xlsx') #USER INPUT: (RetObs, RetObsUnc, HeObs, HeObsUnc, BeObs, BeObsUnc, HeConc and unc, BeConc and unc)
dataObs=np.array(dataObs)

#Grain size on which measurement conducted
r=0.045 #USER INPUT; cosmogenic grain size in cm (radius); og: dataObs.r=0.045

"3) Load Diffusion Kinetics Data"

#Load raw data diffusion kinetics
dataKD=pd.read_excel('./Data/KD_BestFit_MBTP9.xlsx') #USER INPUT: (radius, domain, Ea, LnDoaa, fraction)
#Normalize LnD0aa and format
dataKD=FormatNormLnD0aa(dataKD, r)

"4) Load EDT Scenarios"
d = np.loadtxt('FixedLGMmin15_MBTP1.txt', delimiter= '\t') #importing the txt file

T0 = 0 #if you want to impose a surface temperature amplitude

c=1
TotHestep = [] #for use in the following for loops

fname = np.zeros((len(d), len(d[0]))) #initializing matrix

for j in range(1,len(d[0]),1):
    fname[:, 0]=d[:, 0] #fname is the matrix created by this line from the txt file
    fname[:, 1]=d[:, j]
    results = MDD_retention_EDTScenario(fname, r)
    TotHestep[:, 0]=fname[:, 0] #time
    TotHestep[:, c]=results.TC #IsoEDT
    TotHestep[:, c+1]=results.actTotHestep #concentration variation
    c=c+2
    
# "5) Plot Results vs. Observed Concentration"
# #1. Plot IsoEDT He Concentration
# for i in range(2,len(TotHestep[0]),2):
#     time=np.flip(TotHestep[:, 1])
    
#     #Figure 1
#     mpl.plot(time, TotHestep[:, i])
#     xlabel('Time (yr)')
#     ylabel('Temperature (Celsius)')
#     mpl.show()
    
#     #Figure 2:
#     mpl.plot(time, TotHestep(i+1))
#     xlabel('Time (yr)')
#     ylabel('3He Concentration')
#     mpl.show()

# #2. Plot Observed He Concentrations
# #Figure 3
# expBe=(time(1)-round(data.BeObs, -3))
# mpl.errorbar(expBe, data.HeConc, data.BeObsUnc) #not sure if this is actually plotting the horizontal error bars
# mpl.errorbar(expBe, data.HeConc, data.HeObsUnc) #not sure if this is actually plotting the vertical error bars
# mpl.plot(expBe, data.HeConc, 'o')
# mpl.show()
    
    