#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spherical-geometry production-diffusion model

This has spatially homogeneous production that can vary with time. 
Code was originally written by Greg Balco & David Shuster for modeling apatite 4He/3He datasets. Modified here by Marissa Tremblay to model production and diffusion of cosmogenic noble gases, and once again by Taylor Bourikas.
Contact: tbourika@purdue.edu
Last modified: 2024.02.05

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
import MDD_retention_EDTScenario #imports variables from wrapper code
from FormatNormLnD0aa import *
from NCEPatm_2 import *
from thickness import *
from stone2000 import *

#Read SampleSpe
dataSpe=pd.read_excel('./Data/SampleSpe_MBTP1.xlsx')


def ProdDiff_EDTScenario(r, n, maxt, Tt, dt, TC, TZ):
    "1)Assign Variables"
    #Assign Variables
    lat=np.array(dataSpe['lat']) #Sample latitude in decimal degrees
    long=np.array(dataSpe['long']) #Sample longitude in decimal degrees
    elev=np.array(dataSpe['elv']) #Sample elevation in meters
    Pref=np.array(dataSpe['Pref']) #Sea level high latitude producation rate in atoms/g/yr
    shield=np.array(dataSpe['shield']) #Correction factor for shielding
    thick=np.array(dataSpe['thick']) #Sample thickness in cm
    rho=np.array(dataSpe['rho']) #Mineral denisty in g/cm^3
    fsp=np.array(dataSpe['fsp']) #Fraction of spallogenic cosmogenic nuclide production
    Lsp=np.array(dataSpe['Lsp']) #Effective attenuation length for spallation in g/cm^2
    edot=np.array(dataSpe['edot']) #Erosion rate in g/cm^2/yr
 
    #Assign Variables
    a=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[0]) #radius of proton-irradiated grain analyzed (cm)
    ndom=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[1]) #number of domains from the best-fit MDD model
    allEa=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[2]) #vector of activation energies for each domain
    allLnD0aa=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[3]) #vector of ln(D0/aa) for each domain. Technically unitless (normalized to 1/s)
    fracs=np.array(wrap_PD_MDD_LGMScenario_MBTP1_fig3.dataKD[4]) #fraction of gas apportioned to each diffusion domain
    #print(a, ndom, allEa, allLnD0aa, fracs)
    actHeTot = [] #initialize list for total number of atoms
    romHeTot = [] #initialize list for total number of romberg(?) atoms
    totwtAll = [] #initialize list for total weight
    actHeatomsg = [] #initialize list for total grams
    romHeatomsg = [] #initialize list for total grams (romberg)
    totProduced = [] #initialize list of atoms(?) produced
    checkTotal = [] #initialize list for checking number of atoms(?) produced
    actRsave = []
    romRsave = []
    actTotHesave = []
    romTotHesave = []
    totProducedConc = [] #total produced concentration of atoms
    "2) Set Up Production-Diffusion Calculation"
    for dom in range(0,len(ndom),1):
        thisdom = dom
        ea = allEa[thisdom]
        lnD0aa = allLnD0aa[thisdom]
        #Make the calculation mesh
        dx = r/n #mode dimension, according to Cassata units and values are irrelevant
        x = np.vstack(np.arange((dx/2), ((r-dx/2)+dx), dx)) #radial distance of the center of each node
        x = x.flatten() #Make vector (512,)
        xlo = x-(dx/2) #radial distance of the inside of each node
        xhi = x+(dx/2) #radial distance of the outside of each node
        shellvol = (4/3)*math.pi*((xhi**3)-(xlo**3)) #volume of each node
        #Weight of shells
        shellwt = np.multiply(shellvol, rho[:, None]) #g
        shellwt = np.transpose(shellwt)
        totwt = (4/3)*math.pi*(r**3)*rho
        #Convert to Kelvin
        T = TC + 273
        #Deal with diffusion constants
        R = 8.3144 #j/K/mol
        #Convert Ea to J/mol
        ea = ea*(10**3)
        D0 = (60*60*24*365.25)*((math.exp(lnD0aa))*(r**2)) #frequency factor (cm^2/yr)
        #Deal with production rate constants
        pressure = NCEPatm_2(lat, long, elev) #looks up surface pressure and 1000 mb from NCEP reanalysis. Need this for site-specific production rate calc
        thickcorr = thickness(thick, Lsp, rho) #calculate the corrections needed for sample thickness 
        #calculate geographic scaling factor and multiply by the reference production rate, shielding correction factor, and thickness correction factor to determine the local cosmogenic muclide production rate
        stone = stone2000(lat, pressure, fsp) 
        P0 = []
        actRsave = []
        romRsave = []
        checkRlist = []
        for p in range(0, len(Pref)):     
            P0sub = Pref[p]*shield[p]*thickcorr[p]*stone[p]
            P0.append(P0sub)
        "3) Initialize New Variables"
        oldC = np.zeros(np.size(x)) #initialize concentration vector
        tv = np.arange(0, (maxt+dt), dt) #initialize time vector
        Ps = np.zeros((len(Pref), len(Tt[0]))) #initialize production vector
        totHe = [] #initialize vector for total He produced
        "4) Start the Solver Loop"
        #Production in time step 1
        for a in range(0, len(tv)): #Step 1 is all zeros, start at step 2
            #Changed from a in range(2, length(tv)); need production in the first time step
            #Change of old variables
            oldU = []
            for u in range(0, len(oldC)):
                oldU.append(oldC[u]*x[u])
            oldU = np.array(oldU)
            #Obtain K
            Tt = Tt.flatten() #numpy interp function needs 1D arrays, and Tt is techically 2D (1, 181). Same for T and TZ.
            T = T.flatten()
            thisT = np.interp(tv[a], Tt, T) #temperature vector, interpolated between
            D = D0*math.exp(-ea/R*(1/(thisT))) #diffusivity in cm^2/yr
            K = D
            #Solver setup
            beta = 2*(dx**2)/(K*dt)
            A1 = np.diag(np.diag(np.ones((n, n))))*(-2-beta) + np.diag(np.diag(np.ones((n-1, n-1))), -1) + np.diag(np.diag(np.ones((n-1, n-1))), 1)
            A2 = np.diag(np.diag(np.ones((n, n))))*(2-beta) - np.diag(np.diag(np.ones((n-1, n-1))), -1) - np.diag(np.diag(np.ones((n-1, n-1))), 1)
            #No-flux LH boundary
            A1[0][0] = A1[0][0]-1
            A2[0][0] = A2[0][0]+1
            #He source
            #This is due to cosmic ray production
            #Obtain depth
            TZ = TZ.flatten()
            thisz = np.interp(tv[a], Tt, TZ)
            depth = thisz #depth(a) = thisz
            #thisz = (maxt - dt*(a-1)*ee)
            #Not sure which is the preferred way to calculate thisz
            #Compute production at depth
            Pofx = P0*np.exp(-thisz/Lsp) #Atoms/g/yr
            #store for tot-up calculation later
            for val in range(0, len(Pofx)):
                Ps[val][a] = Pofx[val] 
            #Add source to totHe - remember totHe is total He production, not how much He is in grain now. totHe should be exact always
            newhe = [] #atoms
            for g in range(0, len(Pofx)):
                newhesub = []
                for h in range(0, len(shellwt)):
                    newhesub.append((Pofx[g]*shellwt[h][0]*dt))
                newhesub = np.array(newhesub)
                newhe.append(newhesub)
            newhe = np.transpose(newhe)
            if a==0:
                sumhelist = [] #for summing later
                for v in range(0, len(newhe[0])): 
                    sumhe = 0 #for summing together all helium atoms
                    for u in range(0, len(newhe)):
                        sumhe += newhe[u][v] 
                    sumhelist.append(sumhe)
                totHe.append(sumhelist) #atoms
            else:
                newsumlist = []
                for v in range(0, len(newhe[0])):
                    sumhe = 0
                    for u in range(0, len(newhe)):
                        sumhe += newhe[u][v]
                        newsum = sumhe + totHe[a-1][v]
                    newsumlist.append(newsum)
                totHe.append(newsumlist) #atoms
            
            #Build RHS
            mult = [] #atoms
            for g in range(0, len(Pofx)):
                multsub = []
                for h in range(0, len(x)):
                    multsub.append(Pofx[g]*x[h]*beta*dt) 
                multsub = np.array(multsub)
                mult.append(multsub)
            mult = np.array(mult)
            mult = np.transpose(mult)
            b = [] #results for following loop
            for f in range(0, len(mult[0])):
                if a == 0:
                    sub = (np.matmul(A2, oldU)) - mult[:,f]
                    b.append(sub)
                else:
                    sub = (np.matmul(A2, oldU)[:,f]) - mult[:,f]
                    b.append(sub)
            b = np.array(b) #make list into an array
            b = np.transpose(b) #Each sample has its own column
            #Solve
            if np.isnan(np.reciprocal(np.linalg.cond(A1))):
                np.disp("K = ", str(K))
            newU = []
            for sq in range(0, len(b[0])):
                c = b[:,sq]
                newUsub = np.linalg.solve(A1, c)
                newU.append(newUsub)
            newU = np.transpose(newU)
            #Unchange of variables
            newC = np.empty((len(newU), len(newU[0])))
            for d in range(0, len(newU[:,0])):
                newC[d] = (newU[d]/x[d]) #Divide each row in newU by the value in the corresponding row in x
            #Update
            oldC = newC
            
            #Sum up total He left in domain.
            #Use trapezoidal integration
            actHe = sum(oldC*shellwt) #atoms
            
            #This code taken from rombint.m
            #Figure number of iterations
            romHe = []
            for rhocount in range(0, len(rho)):
                decdigs = 1+round(math.log2(n-1))
                rom = [[0 for col in range(decdigs)] for row in range(2)]
                romall = oldC[:,rhocount]*4*math.pi*(x**2)*rho[rhocount]
                romall = np.append(romall, 0)
                h = r
                rom[0][0] = h*(romall[0]+romall[-1])/2
                for j in range(1, decdigs):
                    st=2**(decdigs-j)
                    romallsum = []
                    romallsmall = romall[int((st/2)):(2**(decdigs-1)):st] #Creating an array of the specific values we want to sum over
                    for num in romallsmall:
                        romallsum.append(num)
                    romallsum = np.array(romallsum) 
                    rom[1][0] = (rom[0][0]+h*sum(romallsum))/2 
                    for k in range(2, j+2):
                        rom[1][k-1] = ((4**(k-1))*rom[1][k-2]-rom[0][k-2])/((4**(k-1)-1)) 
                    rom[0][0:j] = rom[1][0:j]
                    h = h/2
                romHesub = rom[0][decdigs-2] #this should also yield atoms.
                romHe.append(romHesub)
            #Number of atoms
            actHeTot.append(actHe)
            romHeTot = list(romHeTot)
            romHeTot.append(romHe)
            totwtAll.append(totwt)
            
            #Compute total produced
            totProducedsub = []
            romRlist = []
            for samp in range(0, len(totwtAll[thisdom])):
                totProducedsub.append(sum(Ps[samp]*dt*totwtAll[thisdom][samp])) #atoms
            totProduced.append(totProducedsub)
            sumlist = []
            for psNum in range(0, len(Ps)):
                sums = sum(Ps[psNum]*dt)
                sumlist.append(sums)
            totProducedConc.append(sumlist) #atoms
            
            romHeTot = np.array(romHeTot) #converts list to array so that later calculations are possible
        checkTotalsub = [] #creates list to calculate checkTotal for every sample in this domain
        actRlist = []
        romRlist = []
        checkRlist = []
        
        for sam in range(0, len(Ps)):
            checkTotalsub.append(stat.mean(Ps[sam])*maxt)
            actR = actHeTot[a][sam]/totProduced[a][sam] 
            actRlist.append(actR)
            romR = romHeTot[dom][sam]/totProduced[dom][sam]
            romRlist.append(romR)
            checkR = actHeTot[a][sam]/checkTotalsub[sam]/totwtAll[a][sam]
            checkRlist.append(checkR)
        actRsave.append(actRlist)
        romRsave.append(romRlist)
        checkTotal.append(checkTotalsub) #append the checkTotalsub amounts for this domain to the masterlist
        romHeTot = list(romHeTot) #converts array back into list so more values can be appended to it
        
        #a this the time step, so here we have R and He for each time step for the given domain

        # actRsave.append(actRlist) does this get returned?
        #romRsave.append(romRsave) does this get returned?
        for sample in range(0, len(actHeTot[thisdom])):
            actHeatomsg.append(actHeTot[thisdom][sample]/totwtAll[thisdom][sample])
        for romSample in range(0, len(romHeTot[thisdom])):
            romHeatomsg.append(romHeTot[thisdom][romSample]/totwtAll[thisdom][romSample])
        #Add Nat
        actTotHesave.append(actHeatomsg)
        romTotHesave.append(romHeatomsg)
        
        for num in range(0, len(actRsave[:][dom])):
            actRsave[:][dom][num] = actRsave[:][dom][num]*fracs[dom]
            romRsave[:][dom][num] = romRsave[:][dom][num]*fracs[dom]
            #Add Nat
            actTotHesave[:][dom][num] = actTotHesave[:][dom][num]*fracs[dom]
            romTotHesave[:][dom][num] = romTotHesave[:][dom][num]*fracs[dom]
            totProducedConc[:][dom][num] = totProducedConc[:][dom][num]*fracs[dom]
        print(actRsave) #only shows very last row of what's supposed to be actRsave.
    
    "5) Sum Up Retention Over All Domains"
    actRtotal = sum(actR*fracs)
    romRtotal = sum(romR*fracs)
    actHeatomsgtotal = sum(actHeatomsg*fracs)
    
    actRtstep = sum(actRsave, 2)
    romRtstep = sum(romRsave, 2)
    actTotHestep = sum(actTotHesave, 2)
    romTotHestep = sum(romTotHesave, 2)
    totProducedConcstep = sum(totProducedConc, 2)
    
    return P0, actHeTot, romHeTot, totwtAll, actHeatomsg, romHeatomsg, tv, depth
    print(tv, depth)