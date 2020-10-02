#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 09:57:53 2020

@author: HayleyAgler
"""

#import necessary modules
import numpy as np
from gaussxw import gaussxw 
import matplotlib.pyplot as plt


#Ta avg hourly temperature in C
Ta=np.arange(-60,60,1)
th= (24, 48, 72) #hours
u10=(6, 8, 10) #m/s

#define function we want to integrate
def P(u, u_mean1,stdev):
    return np.exp(-(u_mean1-u)**2/(2*stdev**2))

stdev=np.zeros(len(Ta))
u_mean1=np.zeros((len(th),len(Ta)))

N=100
a=0
b=u10
for i in range(len(Ta)):
    for j in range(len(th)):
        stdev[i]=4.3+0.145*Ta[i] + 0.00196*Ta[i]**2
        u_mean1[j,i]=11.2+0.365*Ta[i] +0.00706*Ta[i]**2 +0.9*np.log(th[j])
   
#gaussian quad calculation
        
x,w = gaussxw(N)
xp1=[0,0,0]
wp1=[0,0,0]
for j in range(len(th)):
    xp1[j] =0.5*(b[j]-a)*x + 0.5*(b[j]+a)
    wp1[j] = 0.5*(b[j]-a)*w

s_gaus1=np.zeros((len(th),len(u10),len(Ta))) 
# Perform the integration
for i in range(len(Ta)):
    for k in range(N):
        for j in range(len(th)): #row,col
            for l in range(len(th)):
                s_gaus1[j,l,0] = 0.0
                s_gaus1[j,l,i] += wp1[j][k]*P(xp1[j][k], u_mean1[l,i], stdev[i])
        

I_gaus1=np.zeros((len(th),len(u10),len(Ta)))
for i in range(len(Ta)):
    for j in range(len(th)): #row,col
        for l in range(len(th)):
            I_gaus1[j,l,i]=s_gaus1[j,l,i]*1/(np.sqrt(2*np.pi)*stdev[i])





plt.xlabel("Temperature (\u2103)")
plt.ylabel("Probability of Blowing Snow")
plt.title("Probability of Blowing Snow for Different Temperatures and Parameters")
plt.plot(Ta, I_gaus1[0,0,:], label="$u_{10}$=6, $t_{h}$=24", color="blue")
plt.plot(Ta, I_gaus1[0,1,:], label="$u_{10}$=6, $t_{h}$=48", linestyle='dashed', color="blue")
plt.plot(Ta, I_gaus1[0,2,:], label="$u_{10}$=6, $t_{h}$=72",linestyle='dashdot', color="blue")
plt.plot(Ta, I_gaus1[1,0,:], label="$u_{10}$=8, $t_{h}$=24", color="green")
plt.plot(Ta, I_gaus1[1,1,:],linestyle='dashed', label="$u_{10}$=8, $t_{h}$=48", color="green")
plt.plot(Ta, I_gaus1[1,2,:], label="$u_{10}$=8, $t_{h}$=72",linestyle="dashdot", color="green")
plt.plot(Ta, I_gaus1[2,0,:], label="$u_{10}$=10, $t_{h}$=24", color="red")
plt.plot(Ta, I_gaus1[2,1,:],linestyle='dashed', label="$u_{10}$=10, $t_{h}$=48", color="red")
plt.plot(Ta, I_gaus1[2,2,:], label="$u_{10}$=10, $t_{h}$=72",linestyle='dashdot', color="red")
plt.legend()

