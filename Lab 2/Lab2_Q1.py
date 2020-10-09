# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:34:51 2020

@author: Nicholas Pavanel
"""
#import needed modules
import numpy as np
import matplotlib.pyplot as plt

#define forward differnce method
def fwd_dif(x,hs):
    results=[]
    for i in hs:
        results.append((np.exp(-(x+i)**2) - np.exp(-(x)**2))/((x+i) - x))
    return results

#define h values
hs=[10**(-16),10**(-15),10**(-14),10**(-13),10**(-12),10**(-11),10**(-10),10**(-9),10**(-8),10**(-7),10**(-6),10**(-5),10**(-4),10**(-3),10**(-2),10**(-1),10**(0)]

#run forward differnce method at x=0.5
results=fwd_dif(0.5,hs)

#compute analytical answer for comparison
analytical_answer=-2*(0.5)*np.exp(-(0.5)**2)

#find error
derivative_error=np.abs(np.array(results)-analytical_answer)

#plot fwd dif error 
plt.figure(figsize=[8,6])
plt.plot(hs,derivative_error)
plt.title('Derivative Errors and Step Sizes')
plt.xlabel('Step Size')
plt.ylabel('Derivative Error')
plt.yscale('log')
plt.xscale('log')
plt.savefig('1_c')
plt.show()

#define central difference method
def cntrl_dif(x,hs):
    results=[]
    comp=[]
    for i in hs:
        results.append((np.exp(-(x+i)**2) - np.exp(-(x-i)**2))/((x+i) - (x-i)))
        comp.append((np.exp(-(x+i)**2) - np.exp(-(x-i)**2), (x+i) - (x-i)))
    return results, comp

#run central dif method at x=0.5
results_cntrl,comp=cntrl_dif(0.5,hs)

#find error
derivative_error_cntrl=np.abs(np.array(results_cntrl)-analytical_answer)

#plot both fwd and cntrl difference 
plt.figure(figsize=[8,6])
plt.plot(hs,derivative_error,label='Forward Difference Method',color='black')
plt.plot(hs,derivative_error_cntrl,label='Central Difference Method',color='red')
plt.title('Derivative Errors and Step Sizes')
plt.xlabel('Step Size')
plt.ylabel('Derivative Error')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig('1_d')
plt.show()
