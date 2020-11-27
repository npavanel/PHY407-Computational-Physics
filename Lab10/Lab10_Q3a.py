#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:27:45 2020

@author: HayleyAgler
"""

#Q3
import numpy as np
import matplotlib.pyplot as plt

N=10000 #sampling points

#define integrand
def f(x):
    return x**(-1/2)/(1+np.exp(x)) 

#bounds of integration
a=0
b=1

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
#importance sampling

#define function for transformation method of sampling
def rand(z):
    return z**2

#define weighting function
def w(x):
    return x**(-1/2)

#define function g(x)=f(x)/w(x) for convenience
def g(x): 
    return 1/(1+np.exp(x))

#empty arrays to store values in
g_stored=np.zeros(N)
g_avg_stored=np.zeros(repeat)
for j in range(repeat):
    for i in range(N):
        x=np.random.uniform(a,b) #generate random value
        xi=rand(x) #plug random value into rand to generate nonuniform random value
        g_stored[i]=g(xi) #compute g(x) and store for later
        g_avg=2*np.sum(g_stored)/N #factor of 2 comes from integral of w(x)
        g_avg_stored[j]=g_avg
print(g_avg_stored)    
    
#%%
#make histograms

plt.figure(1)
plt.title("Integral values using importance sampling")
plt.xlabel("Integral value")
plt.ylabel("Number of values")
plt.hist(g_avg_stored, 10, range=[0.8, 0.88])
plt.show()

plt.figure(2)
plt.title("Integral values using mean value method")
plt.hist(f_avg_stored, 10, range=[0.8, 0.88])
plt.xlabel("Integral value")
plt.ylabel("Number of values")
plt.show()






