#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:25:42 2023

@author: taylo
"""
import math
import numpy as np
import scipy as sp
import pandas as pd
import scipy.io as sio

# decdigs = 1+round(math.log2(512-1))
# print(decdigs)

# h = 2
# oldC = 4
# rho = 2300
# x = np.array([2, 4, 6])
# rom = np.zeros((2, decdigs))
# print(rom)
# romall = oldC*4*math.pi*(x**2)*rho
# print(romall)
# for j in range(1, decdigs):
#     st = 2**(decdigs-j+1)
#     ''''rom[1,0] = rom[0,0]+h*sum(romall[np.arange((st/2)+1, 2**(decdigs-1), st)])/2'''
#     testfunc = sp.integrate.romb(romall, st)
#     print(testfunc)

# silly = np.array([[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3], [4, 4, 4, 4], [5, 5, 5, 5], [6, 6, 6, 6], [7, 7, 7, 7]])
# print(silly) #a 2D array
# print(silly[:,0])

# one = np.array([1, 1, 1, 1])
# two = np.array([2, 2, 2, 2])
# three = np.array([3, 3, 3, 3])
# four = np.array([4, 4, 4, 4])
# five = np.array([5, 5, 5, 5])
# six = np.array([6, 6, 6, 6])
# seven = np.array([7, 7, 7, 7])
# goofy = []
# goofy.append(one)
# goofy.append(two)
# goofy.append(three)
# goofy.append(four)
# goofy.append(five)
# goofy.append(six)
# goofy.append(seven)
# print(goofy) #a list of 7 1D arrays (which is what we have rn)

# cheese = []
# cheese.append(goofy[0][0])
# cheese.append(goofy[1][0])
# cheese.append(goofy[2][0])
# cheese.append(goofy[3][0])
# cheese.append(goofy[4][0])
# cheese.append(goofy[5][0])
# cheese.append(goofy[6][0])
# print(cheese) #a 1D array made of the first value from each array in goofy

# NCEP2data = sio.loadmat('NCEP2.mat')
# print(NCEP2data['NCEPlon'])

# ones = np.ones((1,512)) #Makes a 1x512 array of ones
# print(np.shape(ones))
# ones2 = np.ones((512, 512))
# print(np.shape(ones2))
# diag1 = np.diag(np.diag(np.ones((512, 512))), -1)
# print(diag1)

# Pofx = [1, 2, 3, 4]
# shellwt = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# dt = 100

# newhe = []
# for i in range(0, len(Pofx)):
#     newhesub = []
#     for j in range(0, len(shellwt)):
#         newhesub.append(Pofx[i]*shellwt[j]*dt)
#     print(newhesub)
#     newhe.append(newhesub)
# print(np.shape(newhe))
# print(newhe)
# print(type(Pofx))
# print(type(shellwt))


# def testing(site_lon, site_lat, site_elev):
#     NCEP2data = sio.loadmat('NCEP2.mat')
#     iP = sp.interpolate.RectBivariateSpline(NCEP2data['NCEPlat'][::-1], NCEP2data['NCEPlon'], NCEP2data['meanslp'])
#     site_slp = iP(site_lat, site_lon, grid=False)
#     iT = sp.interpolate.RectBivariateSpline(NCEP2data['NCEPlat'][::-1], NCEP2data['NCEPlon'], NCEP2data['meant1000'])
#     site_T = iT(site_lat, site_lon, grid=False)
#     site_T_degK = site_T + 273.15
#     gmr = -0.03417
#     dtdz = 0.0065
#     NCEP2 = site_slp*np.exp((gmr/dtdz)*(np.log(site_T_degK)-np.log(site_T_degK-(site_elev*dtdz))))
#     print(NCEP2)

# n = 512
# oldC = np.ones((512, 4))
# rho = [1, 2, 3, 4]
# r = 0.045
# x = np.ones((512, 4))

# romHe = []
# for rhocount in range(0, len(rho)):
#     decdigs = 1+round(math.log2(n-1))
#     rom = np.zeros((2, decdigs))
#     romall = oldC*4*math.pi*(x**2)*rho[rhocount]
#     romall = np.append(romall, 0)
#     h = r
#     rom[0,0] = h*(romall[0]+romall[-1])/2
#     for j in range(1, decdigs):
#         st=2**(decdigs-j)
#         print(st/2)
#         print((2**(decdigs-1)-1))
#         print(st)
#         romallsum = []
#         for num in range(int(st/2), int(2**(decdigs-1)-1), st):
#             number = romall[num]
#             romallsum.append(number)
#         romallsum = np.array(romallsum)    
#         print(romallsum)
#         rom[1,0] = (rom[0,0]+h*sum(romallsum))/2
#         for k in range(0, j-2):
#             rom[1, k] = ((4**k)*rom[1, k-1]-rom[0, k-1])/((4**k)-1)
#         rom[0, 0:j] = rom[1, 0:j]
#         h = h/2
#         romHesub = rom[1, decdigs] #this should also yield atoms
#         romHe.append(romHesub)
# print(romHe)

A2 = [[0 for _ in range(3)] for _ in range(3)]
A2[0][0], A2[0][1], A2[0][2] = 2.9899, -1, 0
A2[1][0], A2[1][1], A2[1][2] = -1, 1.9899, -1
A2[2][0], A2[2][1], A2[2][2] = 0, -1, 0.9899

oldU = [[0 for _ in range(1)] for _ in range(3)]
oldU[0][0] = 0
oldU[1][0] = 0
oldU[2][0] = 0

mult = [[0 for _ in range(2)] for _ in range(3)]
mult[0][0], mult[0][1] = 0.02136, 0.02228
mult[1][0], mult[1][1] = 0.06409, 0.06685
mult[2][0], mult[2][1] = 0.10683, 0.11141

b = []
for f in range(0, len(oldU[0])):
    sub = (np.matmul(A2, oldU)) - mult[:][f]
    b.append(sub)
b = np.array(b) #make list into an array
b = np.transpose(b) #Each sample has its own column
print(b)