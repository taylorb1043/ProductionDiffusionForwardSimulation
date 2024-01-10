#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function transforms values of ln(D0/a^2) from diffusion experiment grain sizes to the size appropriate for cosmogenic noble gas measurements.

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

import numpy as np
import math

def FormatNormLnD0aa(dataKD, r):
    a=dataKD.iloc[:,0] #Grabs data from first column
    ndom=dataKD.iloc[:, 1] #Grabs data from second column
    Ea=dataKD.iloc[:, 2] #Grabs data from third column
    #Normalization to a
    lnD0aa=[]
    for i in range(0, len(dataKD)):
        #print(math.exp(dataKD.iloc[i,3])*(a[i]**2))
        norm = np.log(math.exp(dataKD.iloc[i, 3])*(a[i]**2)/(r**2))
        #print(norm)
        lnD0aa.append(norm)
    fracs=[]
    for j in range(0, len(dataKD)):
        fracs.append(dataKD.iloc[j, 4])
    return a, ndom, Ea, lnD0aa, fracs