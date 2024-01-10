#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Setup to do production and diffusion with MDD model.
Takes user defined inputs from the wrapper code and packages them to be passed to the production-diffusion code.

Written by Marissa Tremblay.
Contact: tremblam@purdue.edu
Translated to Python and edited by Taylor Bourikas
Contact: tbourika@purdue.edu
Last modified: 2023.11.06

Copyright 2021, Marissa Tremblay
All rights reserved
Developed in part with funding from the National Science Foundation and the American Association for the Advancement of Science.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

"""
import pandas as pd
import numpy as np
import math
import statistics as stat
import wrap_PD_MDD_LGMScenario_MBTP1_fig3
from ProdDiff_EDTScenario import *

def MDD_retention_EDTScenario(fname, r):
    #Read SampleSpe
    dataSpe=pd.read_excel('./Data/SampleSpe_MBTP1.xlsx')
    
    #Assign Variables
    lat=np.array(dataSpe.iloc[0,0]) #Sample latitude in decimal degrees
    long=np.array(dataSpe.iloc[0,1]) #Sample longitude in decimal degrees
    elev=np.array(dataSpe.iloc[0,2]) #Sample elevation in meters
    Pref=np.array(dataSpe.iloc[0,3]) #Sea level high latitude producation rate in atoms/g/yr
    shield=np.array(dataSpe.iloc[0,4]) #Correction factor for shielding
    thick=np.array(dataSpe.iloc[0,5]) #Sample thickness in cm
    rho=np.array(dataSpe.iloc[0,6]) #Mineral denisty in g/cm^3
    fsp=np.array(dataSpe.iloc[0,7]) #Fraction of spallogenic cosmogenic nuclide production
    Lsp=np.array(dataSpe.iloc[0,8]) #Effective attenuation length for spallation in g/cm^2
    ee=np.array(dataSpe.iloc[0,9]*(10**(-6))*100*dataSpe.iloc[0,6]) #Erosion rate in g/cm^2/yr
    ztoday=np.array(dataSpe.iloc[0,10]) #if sample is not at surface today, depth in g/cm^2
    #radius = data.r average radius of fragments analyzed, in cm
    
    #Read SampleObs
    dataObs=pd.read_excel('./Data/SampleObs_MBTP1.xlsx')
    
    #Assign variables
    retObs=np.array(dataObs.iloc[0,0]) #observed noble gas retention
    retObsUnc=np.array(dataObs.iloc[0,1]) #1 standard deviation in observed noble gas average retention
    heObs=np.array(dataObs.iloc[0,2]) #apparent 3He exposure age in years
    heObsUnc=np.array(dataObs.iloc[0,3]) #1 standard deviation in apparent 3He exposure age
    beObs=np.array(dataObs.iloc[0,4]) #observed 10Be exposure age in years
    beObsUnc=np.array(dataObs.iloc[0,5]) #1 standard deviation in 10Be exposure age
    
    
    #Assign Variables- this should eventually be able to run multiple rows
    a=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[0]) #radius of proton-irradiated grain analyzed (cm)
    ndom=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[1]) #number of domains from the best-fit MDD model
    allEa=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[2]) #vector of activation energies for each domain
    allLnD0aa=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[3]) #vector of ln(D0/aa) for each domain. Technically unitless (normalized to 1/s)
    fracs=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[4]) #fraction of gas apportioned to each diffusion domain
    
    n=512 #number of mesh points over which production-diffusion calculation takes place
    
    temphist = fname #load in the user-defined time-temperature histories
    times = []
    temps = []
    for inner in temphist:
        times.append(inner[0]) #the first column in this file should correspond to the time steps in years
        temps.append(inner[1])
    maxt = max(times) #total duration of the time-temperature history, in years
    times = np.atleast_2d(times).T
    Tt = np.transpose(times)
    dt = maxt/(len(Tt[0])-1) #size of the time step, in years
    
    totProduced = [] #needed for for loop later
 
    #Run the production-diffusion function for each temperature history the user wants to test
    for i in range(1, len(temps)):
        #temperatures corresponding to each time step for a particular thermal history in degrees C
        temps = np.atleast_2d(temps).T
        TC = np.transpose(temps)
        TC = np.fliplr(TC)
        TZ = np.fliplr(Tt)*ee*ztoday #calculate sample depth at each time step in g/cm^2
        output = ProdDiff_EDTScenario(r, n, maxt, Tt, dt, TC, TZ)
        
        actRtstep = ProdDiff_EDTScenario.actRtstep 
        romRtstep = ProdDiff_EDTScenario.romRtstep 
        
        #results.actHeatomsgtotal(i) = output.actHeatomsgtotal
        totProduced[i] = ProdDiff_EDTScenario.totProducedConcstep[-1] 
        scale = ProdDiff_EDTScenario.P0/Pref
        totProduced = ProdDiff_EDTScenario.totProducedConcstep
        actTotHestep = ProdDiff_EDTScenario.actTotHestep 
        romTotHestep = ProdDiff_EDTScenario.romTotHestep 
        #Careful each results step corresponds to the last EDT-t history scenario, while final R and He values are saved for all the EDT-t scenarios given
        P0 = ProdDiff_EDTScenario.P0 #3He production rate at site