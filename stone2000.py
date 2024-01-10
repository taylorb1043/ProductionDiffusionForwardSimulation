#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the geographic scaling factor for cosmogenic-nuclide production as a function of site latitude and atmospheric pressure, according to:

Stone, J., 2000, Air Pressure and Cosmogenic Isotope Production. JGR 105:B10, p. 23753. 

Syntax: scalingfactor = stone2000(latitude,pressure,fsp)

Units: 
latitude in decimal degrees
pressure in hPa
fsp is the fraction (between 0 and 1) of production at sea level and high latitude due to spallation (as opposed to muons). 
This argument is optional and defaults to 0.978, which is the value used by Stone (2000) for Be-10. The corresponding value for Al-26 is 0.974. Note that using 0.844 for Be-10 and 0.826 for Al-26 will closely reproduce the Lal, 1991 scaling factors as long as the standard atmosphere is used to convert sample elevation to atmospheric pressure. Also note that this function will yield the scaling factor for spallation only when fsp=1, and that for muons only when fsp=0.  

Elevation can be converted to pressure with the functions
stdatm.m (general use) and antatm.m (Antarctica). 
 
Vector argments are OK. All arguments must be the same size. 

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

Translated into Python and edited by: Taylor Bourikas
Contact: tbourika@purdue.edu
Last modified: 2023.11.09
"""
import numpy as np

def stone2000(lat, P, Fsp):    
    #Check for obvious errors
    for i in range(0, len(lat)):
        if abs(lat[i]) > 90:
            print("Latitudes below 90, please")
            return
    
    if len(lat) != len(P):
        print("Vectors the same size, please")
        return
        
    #Default Fsp
    if Fsp is None:
        Fsp = 0.978
        
    #Spallogenic production at index latitudes
    
    #Enter constants from table 1
    a = [31.8518, 34.3699, 40.3153, 42.0983, 56.7733, 69.0720, 71.8733]
    b = [250.3193, 258.4759, 308.9894, 512.6857, 649.1343, 832.4566, 863.1927]
    c = [-0.083393, -0.089807, -0.106248, -0.120551, -0.160859, -0.199252, -0.207069]
    d = [7.4260e-5, 7.9457e-5, 9.4508e-5, 1.1752e-4, 1.5463e-4, 1.9391e-4, 2.0127e-4]
    e = [-2.2397e-8, -2.3697e-8, -2.8234e-8, -3.8809e-8, -5.0330e-8, -6.3653e-8, -6.6043e-8]
    
    ilats = [0, 10, 20, 30, 40, 50, 60] #Need to multiply values so that the size of this matches the size of s and m lists
    
    #Calculate index latitudes at given Ps
    lat0 = a[0] + (b[0] * np.exp(P/(-150))) + (c[0]*P) + (d[0]*(P**2)) + (e[0]*(P**3))
    lat10 = a[1] + (b[1] * np.exp(P/(-150))) + (c[1]*P) + (d[1]*(P**2)) + (e[1]*(P**3))
    lat20 = a[2] + (b[2] * np.exp(P/(-150))) + (c[2]*P) + (d[2]*(P**2)) + (e[2]*(P**3))
    lat30 = a[3] + (b[3] * np.exp(P/(-150))) + (c[3]*P) + (d[3]*(P**2)) + (e[3]*(P**3))
    lat40 = a[4] + (b[4] * np.exp(P/(-150))) + (c[4]*P) + (d[4]*(P**2)) + (e[4]*(P**3))
    lat50 = a[5] + (b[5] * np.exp(P/(-150))) + (c[5]*P) + (d[5]*(P**2)) + (e[5]*(P**3))
    lat60 = a[6] + (b[6] * np.exp(P/(-150))) + (c[6]*P) + (d[6]*(P**2)) + (e[6]*(P**3))
   
    #Initialize output
    correction = np.zeros(np.size(P))
    
    #Northernize southern-hemisphere inputs
    lat = np.array(abs(lat))
    
    #Set high lats to 60
    #lat(find(lat>60)) = (zeros(size(find(lat > 60)))+60) might be able to do this with an if statement?? just need to know if lat is an array or an int    
    #Loop
    b = 0
    s = []
    while b < len(lat):
        samples = []
    #for i in range(0, len(lat)):
        samples.append(lat0[b])
        samples.append(lat10[b])
        samples.append(lat20[b])
        samples.append(lat30[b])
        samples.append(lat40[b])
        samples.append(lat50[b])
        samples.append(lat60[b])
        #interpolate for actual elevation
        x = np.interp(lat[b], ilats, samples)
        s.append(x)
        #continue loop
        b = b+1    
    #Production by muons
    #Constants
    mk = [0.587, 0.600, 0.678, 0.833, 0.933, 1.000, 1.000]
    
    #index latitudes at given Ps
    ml0 = mk[0] * np.exp((1013.25 - P)/242)
    ml10 = mk[1] * np.exp((1013.25 - P)/242)
    ml20 = mk[2] * np.exp((1013.25 - P)/242)
    ml30 = mk[3] * np.exp((1013.25 - P)/242)
    ml40 = mk[4] * np.exp((1013.25 - P)/242)
    ml50 = mk[5] * np.exp((1013.25 - P)/242)
    ml60 = mk[6] * np.exp((1013.25 - P)/242)
    
    #Loop
    b = 0
    m = []
    while b < len(lat):
        samples = []
    #for i in range(0, len(lat)):
        samples.append(ml0[b])
        samples.append(ml10[b])
        samples.append(ml20[b])
        samples.append(ml30[b])
        samples.append(ml40[b])
        samples.append(ml50[b])
        samples.append(ml60[b])
        #interpolate for actual elevation
        y = np.interp(lat[b], ilats, samples)
        #interpolate for actual elevation
        m.append(y)
        #continue loop
        b = b+1
        
    #Combine spallogenic and muogenic production; return
    fm = 1-Fsp
    out = ((s * Fsp)+(m * fm))
    
    #Make vectors horizontal
    #TAYLOR'S COMMENT: Data should always be in a certain structure, so there isn't a need to change the shape of the data.
    # if out.shape[0] > out.shape[1]: #If there are more rows than columns....
    #     out = np.transpose(out)
    # else:
    #     out = out
    return out