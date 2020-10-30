#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 13:54:33 2020

@author: HayleyAgler
"""

import numpy as np
import matplotlib.pyplot as plt


# Constants
m = 9.1094e-31     # Mass of electron
hbar = 1.0546e-34  # Planck's constant over 2*pi
e = 1.6022e-19     # Electron charge   
N = 1000
a=5.3e-11 # Bohr radius
h = 0.003*a #r_f-r_0/N
r_0=h #R=0 here
r_f=20*a  #R=0 here
ep_not=8.8e-12


# Potential function
def V(r):
    return -(e**2)/(4*np.pi*ep_not*r)

#separated solutions to the ODE
def f(p,r,E,l):
    R = p[0] #psi
    S = p[1] #phi
    fR = S/(r**2)
    fS = ((2*(r**2)*m)/(hbar**2))*(V(r)-E)*R+l*(l+1)*R
    return np.array([fR,fS])

# Calculate the wavefunction for a particular energy
def solve(E,l):
    R = 0.0
    S = 1.0
    p = np.array([R,S],float)
    wavefn=[]
    for r in np.arange(r_0,r_f, h): 
        wavefn.append(p[0])
        #print(p, r, E)
        k1 = h*f(p,r,E,l)
        #print(k1, p+0.5*k1, r+0.5*h)
        k2 = h*f(p+0.5*k1,r+0.5*h,E,l)
        k3 = h*f(p+0.5*k2,r+0.5*h,E,l)
        k4 = h*f(p+k3,r+h,E,l)
        p += (k1+2*k2+2*k3+k4)/6
        #print(p+0.5*k1,k2,k3,k4)

    return p[0], wavefn #returns wavefunctions R


# Main program to find the energy using the secant method
E11 = -15*e/1**2 #l=0,n=1
E21 = -13*e/1**2
R21 = solve(E11,0)[0]

target = e/10 
while abs(E11-E21)>target:
    R11,R21 = R21, solve(E21,0)[0]
    E11,E21 = E21,E21-R21*(E21-E11)/(R21-R11) 
    wavefn1=solve(E21,0)[1]
print("For l=0,n=1, E =",E21/e,"eV")

#l=0,n=2
E12 = -15*e/2**2 #l=0,n=2
E22 = -13*e/2**2
R22 = solve(E12,0)[0]

target = e/10 
while abs(E12-E22)>target:
    R12,R22 = R22, solve(E22,0)[0]
    E12,E22 = E22,E22-R22*(E22-E12)/(R22-R12) 
    wavefn2=solve(E22,0)[1]
print("For l=0,n=2, E =",E22/e,"eV")


#l=1,n=2
E13 = -15*e/2**2 #l=1,n=2
E23 = -13*e/2**2
R23 = solve(E13,1)[0]

target = e/10 
while abs(E13-E23)>target:
    R13,R23 = R23, solve(E23,1)[0]
    E13,E23 = E23,E23-R23*(E23-E13)/(R23-R13) 
    wavefn3=solve(E23,1)[1]
print("For l=1,n=2, E =",E23/e,"eV")




#normalize R(r) using trapezoidal rule
#modified trapezoidal for an array
def trapz_int(N, wavefn):
    f=np.array(wavefn)**2
    # 2.: the beginning and end
    result = 0.5*(f[0] + f[N-1])

    # 3. Then, loop over interior points
    for k in range(1, N):
        result += f[0+k]

    return h/a*result  

#%%
#define explicit solutions 
def explicit01(r): #l=0,n=1
    return 2*np.exp(-r)

#define the wavefunction squared for the normalization
def explicit01_2(r): #l=0,n=1
    return (2*np.exp(-r))**2

def explicit02(r): #l=0,n=2
    return 1/(2*np.sqrt(2))*(2-r)*np.exp(-r/(2))

def explicit02_2(r): #l=0,n=2
    return (1/(2*np.sqrt(2))*(2-r)*np.exp(-r/(2)))**2

def explicit12(r):
    return 1/(2*np.sqrt(6))*r*np.exp(-r/2)

def explicit12_2(r):
    return (1/(2*np.sqrt(6))*r*np.exp(-r/2))**2



#define different trapezoidal integration for a function not an array
def trapz(z_0,z_f, N, f):

    # 1.: from a, b and N, compute h
    h = 0.003

    # 2.: the beginning and end
    result = 0.5*(f(z_0)+ f(z_f))

    # 3. Then, loop over interior points
    for k in range(1, N):
        result += f(z_0 + k*h)

    return h*result  #

rvals=np.arange(r_0,r_f, h)/a

#normalize the explicit solutions
explicit01_norm=explicit01(rvals)/np.sqrt(trapz(r_0,r_f,1000,explicit01_2))
explicit02_norm=explicit02(rvals)/np.sqrt(trapz(r_0,r_f,1000,explicit02_2))
explicit12_norm=explicit12(rvals)/np.sqrt(trapz(r_0,r_f,1000,explicit12_2))

#normalize the wave functions and convert to units of bohr radii
wavefn_norm1=np.array(wavefn1)*a**1.5/np.sqrt(trapz_int(N, np.array(wavefn1)*a**1.5)) #l=0,n=1
wavefn_norm2=np.array(wavefn2)*a**1.5/ np.sqrt(trapz_int(N, np.array(wavefn2)*a**1.5)) #l=0,n=2
wavefn_norm3=np.array(wavefn3)*a**1.5/ np.sqrt(trapz_int(N, np.array(wavefn3)*a**1.5)) #l=1,n=2
rvals=np.arange(r_0,r_f, h)/a

#plot the wavefunctions
plt.figure(1)
plt.title("l=0,n=1")
plt.plot(rvals[0:2500], wavefn_norm1[0:2500], label="shooting method")
plt.plot(rvals[0:2500], explicit01_norm[0:2500], label="explicit solution", linestyle="dashed")
plt.xlabel("Bohr radii [5.3e-11]")
plt.ylabel("Normalized Wavefunction R(r)")
plt.legend()


plt.figure(2)
plt.title("l=0,n=2")
plt.plot(rvals, wavefn_norm2, label="shooting method")
plt.plot(rvals, explicit02_norm, label="explicit solution", linestyle="dashed")
plt.xlabel("Bohr radii [5.3e-11]")
plt.ylabel("Normalized Wavefunction R(r)")
plt.legend()


plt.figure(3)
plt.title("l=1,n=2")
plt.plot(rvals, wavefn_norm3, label="shooting method")
plt.plot(rvals, explicit12_norm, label="explicit solution", linestyle="dashed")
plt.xlabel("Bohr radii [5.3e-11]")
plt.ylabel("Normalized Wavefunction R(r)")
plt.legend()