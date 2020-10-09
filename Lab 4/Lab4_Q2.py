# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:31:45 2020

@author: Nicholas Pavanel
"""
import numpy as np
import matplotlib.pyplot as plt
from time import time
import astropy.constants as c
import astropy.units as u
from scipy.linalg import eigh

#6.9 b)

#define constants
hbar = 1.0545718 * 10**(-34) # in J * s
L = 5 * 10**(-10) # in m
M = 9.1094 * 10**(-31) # in kg
a = 1.602176634 * 10**(-18) # in J

#define a function that compute H_{mn} for a given m and n
#below is modelled from page 4 of lab handout
def H_mn(m,n):
    if m != n and m % 2 == 0 and n % 2 == 0: 
        return 0
    elif m != n and m % 2 == 1 and n % 2 == 1:
        return 0
    elif m != n and m % 2 == 0 and n % 2 == 1 or m != n and m % 2 == 1 and n % 2 == 0:
        return - (8 * a * m * n) / (np.pi**2 * (m**2 - n**2)**2)
    elif m == n:
        return 0.5 * a + (np.pi**2 * hbar**2 * m**2) / (2 * M * L**2)
    
#6.9 c)

#modify H_mn to return an array of size nxn, with m,n up to n
def H_matrix(size):
    output = np.zeros((size,size))
    for i in range(1, size + 1):
        for j in range(1, size + 1):
            output[i-1, j-1] = H_mn(i,j) #using lab handout method
    return output

#compute H_mn for the 10x10 case
H_10_10 = H_matrix(10)

#compute eigenvalues and eigenvectors in KMS units? for 10x10 case
evals_H1010, evecs_H1010 = eigh(H_10_10)

#compute H_mn for the 100x100 case
H_100_100 = H_matrix(100)

#compute eigenvalues and eigenvectors in KMS units? for 100x100 case
evals_H100100, evecs_H100100 = eigh(H_100_100)
    
