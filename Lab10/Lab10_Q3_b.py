#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 12:48:03 2020

@author: HayleyAgler
"""

#Q3b
import numpy as np
import matplotlib.pyplot as plt
from scipy import special

N=10000 #sampling points

#define integrand
def f(x):
    return np.exp(-2*np.abs(x-5))

#bounds of integration
a=0
b=10

repeat=100

#mean val method
f_avg=np.zeros(repeat)
f_stored=np.zeros(N)
f_avg_stored=np.zeros(repeat)
for j in range(repeat):
    for i in range(N):
        x=np.random.uniform(a,b) #random points between a and b
        f_stored[i]= f(x)
        f_avg = (b-a)*np.sum(f_stored)/N 
        f_avg_stored[j]=f_avg
print(f_avg_stored) 


#%%
#importance sampling method
#define weighting function
def w(x):
    return 1/(np.sqrt(2*np.pi))*np.exp((-(x-5)**2)/2)

#define mu and sigma for use in the erf
mu=5
sigma=1

#define empty arrays to store values in
fw_stored=np.zeros(N)
fw_avg_stored=np.zeros(repeat)
for j in range(repeat): #loop over 100 times
    for i in range(N): #loop over 10000 times
        xi=np.random.normal(mu,sigma) #compute nonuniformly random points
        fw_stored[i]=f(xi)/w(xi) 
        fw_avg=special.erf(5/np.sqrt(2))*np.sum(fw_stored)/N #compute I
        fw_avg_stored[j]=fw_avg #store for later
print(fw_avg_stored)    

#%%
#make histograms

plt.figure(1)
plt.title("Integral values using importance sampling")
plt.xlabel("Integral value")
plt.ylabel("Number of values")
plt.hist(fw_avg_stored, 10, range=[0.96,1.06])
plt.show()

plt.figure(2)
plt.title("Integral values using mean value method")
plt.hist(f_avg_stored, 10, range=[0.96,1.06])
plt.xlabel("Integral value")
plt.ylabel("Number of values")
plt.show()