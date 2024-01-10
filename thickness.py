#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the thickness correction for cosmogenic nuclide production by spallation, in a sample of thickness zmax (cm) and density rho (g/cm3), with effective attenuation length Lambda (g/cm2). 

Syntax: correction = thickness(zmax, Lambda, density);

Vector arguments are OK. All arguments must be the same size.

Written by Greg Balco -- UW Cosmogenic Nuclide Lab
balcs@u.washington.edu
First version, Feb. 2001
checked March, 2006
Part of the CRONUS-Earth online calculators: 
      http://hess.ess.washington.edu/math

Copyright 2001-2007, University of Washington
All rights reserved
Developed in part with funding from the National Science Foundation.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).

Translated from MATLab to Python 3 by Taylor Bourikas.
Contact: tbourika@purdue.edu
Last modified: 2023.11.06
"""
import numpy as np
import pandas as pd
###
#Read SampleSpe
dataSpe=pd.read_excel('./Data/SampleSpe_MBTP1.xlsx')
#Read KD_BestFit
dataKD=pd.read_excel('./Data/KD_BestFit_MBTP9.xlsx')
zmax=np.array(dataSpe['thick']) #Sample thickness in cm
L=np.array(dataSpe['Lsp'])
rho=np.array(dataSpe['rho']) #Mineral denisty in g/cm^3

def thickness(zmax, L, rho):
    out= (L/(rho*zmax))*(1-np.e**(((-1*rho*zmax)/L)))
    #Make vector horizontal
    #TAYLOR'S COMMENT: Data should always be 1xn, so there isn't a need to change the shape of the data.
    # if out.shape[0] > out.shape[1]: #If there are more rows than columns....
    #     out = np.transpose(out)
    # else:
    #     out = out
    return out
   
        