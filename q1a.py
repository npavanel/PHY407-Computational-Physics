#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 20:09:41 2020

@author: HayleyAgler
"""

#import necessary modules
import numpy as np
from scipy import special
from gaussxw import gaussxw 
import matplotlib.pyplot as plt

#taking from my solution in lab2
#define function we want to integrate
def f(x):
    return np.exp(x**2)

 #number of slices
N=8
a=0
b=4


###error estimation
# =============================================================================
# practical estimation
# =============================================================================
#define new N2 and corresponding h2

N1=np.array([2**3,2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11])
h1=np.zeros(len(N1))
h2=np.zeros(len(N1))
N2=np.array([2**3*2,2**4*2,2**5*2,2**6*2,2**7*2,2**8*2,2**9*2,2**10*2,2**11*2])
s2_trap=np.zeros(len(N1))
s_trap=np.zeros(len(N1))
I2_trap=np.zeros(len(N1))
I1_trap=np.zeros(len(N1))
error_trap=np.zeros(len(N1))
s_simp=np.zeros(len(N1))
I1_simp=np.zeros(len(N1))


for i in range(len(N1)):
    h1[i]=(b-a)/N1[i]
    h2[i]=(b-a)/N2[i]

    s_trap[i] = 0.5*f(a) + 0.5*f(b)

    for k in range(1,np.int(N1[i])):
        s_trap[i] += f(a+k*h1[i])

    I1_trap[i]=s_trap[i]*h1[i]*np.exp(-4**2)

    #redo the integral using trapezoidal method with the new value of N2
    s2_trap[i] = 0.5*f(a) + 0.5*f(b)
    for k in range(np.int(N2[i])):
        s2_trap[i] += f(a+k*h2[i])

    I2_trap[i]=s2_trap[i]*h2[i]*np.exp(-4**2)



#redo the integral using simpsons method with the new value of N2
for i in range(len(N1)):
    
    s_simp[i]=f(a)+f(b)
    #odd terms
    for k in range(1,np.int(N1[i]),2):
        s_simp[i] += 4*f(a+k*h1[i])
        #even terms
    for k in range(2,np.int(N1[i]),2):
        s_simp[i] += 2*f(a+k*h1[i])

    #multiply by other part of function
    I1_simp[i]=s_simp[i]*(1/3)*h1[i]*np.exp(-4**2)
    
    
    s2_simp[i]=f(a)+f(b)
    #odd terms
    for k in range(1,np.int(N2[i]),2):
        s2_simp[i] += 4*f(a+k*h2[i])
        #even terms
    for k in range(2,np.int(N2[i]),2):
        s2_simp[i] += 2*f(a+k*h2[i])
    
    I2_simp[i]=s2_simp[i]*(1/3)*h2[i]*np.exp(-4**2)

    #find the error using textbook eqn
    #do i use this equation or the onefrom this ps
    error_simp[i]=1/15*(I2_simp[i]-I1_simp[i])


#gaussian quad
x2=[0,0,0,0,0,0,0,0,0]
w2=[0,0,0,0,0,0,0,0,0]
xp2=[0,0,0,0,0,0,0,0,0]
wp2=[0,0,0,0,0,0,0,0,0]
s2_gaus=np.zeros(len(N1))
I2_gaus=np.zeros(len(N1))
x=[0,0,0,0,0,0,0,0,0]
w=[0,0,0,0,0,0,0,0,0]
xp=[0,0,0,0,0,0,0,0,0]
wp=[0,0,0,0,0,0,0,0,0]
s_gaus=np.zeros(len(N1))
I1_gaus=np.zeros(len(N1))
error_gaus=np.zeros(len(N1))

for i in range(len(N1)):
    x[i],w[i] = gaussxw(N1[i])
    xp[i] =0.5*(b-a)*x[i] + 0.5*(b+a)
    wp[i] = 0.5*(b-a)*w[i]

    # Perform the integration
    s_gaus[i] = 0.0
    for k in range(N1[i]):
        s_gaus[i] += wp[i][k]*f(xp[i][k])

    I1_gaus[i]=s_gaus[i]*np.exp(-4**2)    
    
    
    
    x2[i],w2[i] = gaussxw(N2[i])
    
    xp2[i] =0.5*(b-a)*x2[i] + 0.5*(b+a)
    wp2[i] = 0.5*(b-a)*w2[i]

    # Perform the integration
    s2_gaus[i] = 0.0
    for k in range(N2[i]):
        s2_gaus[i] += wp2[i][k]*f(xp2[i][k])

    I2_gaus[i]=s2_gaus[i]*np.exp(-4**2)

    error_gaus[i]=I2_gaus[i]-I1_gaus[i]
    
#%% plotting relative error
    
#calculate relative error
#error/real value
rel_gaus=np.zeros(len(N1))
for i in range(len(N1)):
    rel_gaus[i]=(np.abs(I1_gaus[i]-special.dawsn(4)))/special.dawsn(4)
    
rel_simp=np.zeros(len(N1))
for i in range(len(N1)):
    rel_simp[i]=np.abs(I1_simp[i]-special.dawsn(4))/special.dawsn(4)
    
rel_trap=np.zeros(len(N1))
for i in range(len(N1)):
    rel_trap[i]=np.abs(I1_trap[i]-special.dawsn(4))/special.dawsn(4)


plt.figure(1)
plt.title("Error as a Function of N")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("N (number of slices)")
plt.ylabel("Error")
#plt.plot(N1, rel_gaus, label="Relative Gaussian")
#plt.plot(N1, rel_simp, label="Relative Simpson's")
#plt.plot(N1, rel_trap, label="Relative Trapezoidal")
plt.plot(N1, error_gaus, label="$I_{2N}-I_N$")
plt.legend()

