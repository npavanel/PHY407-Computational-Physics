#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 17:03:29 2020

@author: HayleyAgler
"""

#from numpy.random import random
import numpy as np
import math

dim=10 #n
N=1000000

#define the functon that returns 1 if the points are inside the hypersphere, 0 otherwise
def f(r):
    r_squared = 0
    for i in range(len(r)):
        r_squared += r[i]**2
    if r_squared > 1:
        return 0
    else:
        return 1 

#compute the volume using the equation in the handout for comparison
V=(1**(dim)*np.pi**(dim/2))/((dim/2)*math.factorial((dim/2)-1))

#%%
#mean value method 
f_avg=0
f_stored=np.zeros(N)
for i in range(N):
    r=np.random.uniform(-1,1, size=(dim,1)) #choses 10 random values between -1 and 1
    for j in range(dim):
        f_stored[i]= f(r) #stores the value of f(r) for later
        f_avg = 2**(dim)*np.sum(f_stored)/N  #computes I
print(f_avg)

#%%
#compute the error
var=(np.sum(f_stored**2))/N-((np.sum(f_stored))/N)**2
sigma=2**(dim)*np.sqrt(var/N)
