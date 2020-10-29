# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 15:11:33 2020

@author: Nicholas Pavanel
"""

# i modify the code from the text in example 8.7

import numpy as np
import matplotlib.pyplot as plt
from math import sin,pi
from numpy import empty,array,arange,sqrt
import matplotlib.pyplot as plt

# define constants in SI units - kg, s, m
G = 6.6738 * 10**(-11)
M = 1.9891 * 10**(30)
Peri = 1.4710 * 10**(11)
Vperi = 3.0287 * 10**(4)
delta = 1000     # Required position accuracy per unit time - 1km in m

#define time domain
H = 604800      # Size of "big steps" - 1 week in s
N = 5           # do 5 revolutions around the sun
yr = 365.25*24*3600.   # seconds in a yr
T = N * yr         # time domain in s
Nsteps = int(T/H)

# define a function to compute equations of motion
def x_eom(x,y):
    return - G * M * x / (sqrt(x**2 + y**2))**3
def y_eom(x,y):
    return - G * M * y / (sqrt(x**2 + y**2))**3

def f(ics):
    x = ics[0]
    vx = ics[1]
    y = ics[2]
    vy = ics[3]
    r = sqrt(x**2 + y**2)
    return array([vx,x_eom(x,y),vy,y_eom(x,y)],float)

tpoints = arange(0,T,H)
xpoints = []
ypoints = []
r = array([Peri,0,0,Vperi],float)  #ics in [x,vx,y,vy]

# Do the "big steps" of size H
for t in tpoints:

    xpoints.append(r[0])
    ypoints.append(r[2])

    # Do one modified midpoint step to get things started
    n = 1
    r1 = r + 0.5*H*f(r)
    r2 = r + H*f(r1)

    # The array R1 stores the first row of the
    # extrapolation table, which contains only the single
    # modified midpoint estimate of the solution at the
    # end of the interval
    R1 = empty([1,4],float)
    R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

    # Now increase n until the required accuracy is reached
    error = 2*H*delta
    while error>H*delta:

        n += 1
        h = H/n

        # Modified midpoint method
        r1 = r + 0.5*h*f(r)
        r2 = r + h*f(r1)
        for i in range(n-1):
            r1 += h*f(r2)
            r2 += h*f(r1)

        # Calculate extrapolation estimates.  Arrays R1 and R2
        # hold the two most recent lines of the table
        R2 = R1
        R1 = empty([n,4],float)
        R1[0] = 0.5*(r1 + r2 + 0.5*h*f(r2))
        for m in range(1,n):
            epsilon = (R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1)
            R1[m] = R1[m-1] + epsilon
        error = abs(epsilon[0])

    # Set r equal to the most accurate estimate we have,
    # before moving on to the next big step
    r = R1[n-1]
    
# change xpoints and ypoints to AU
x_au = array(xpoints[:53])*6.68459e-12
y_au = array(ypoints[:53])*6.68459e-12
# Plot the results
plt.figure(figsize=[9,9])
plt.plot(x_au,y_au,'b.')
plt.xticks(size=16)
plt.yticks(size=16)
plt.xlabel('$x$ (AU)',size=18)
plt.ylabel('$y$ (AU)',size=18)
plt.savefig('Q2a')
plt.show()

#define ics for pluto
Peri = 4.4368 * 10**(12)
Vperi = 6.1218 * 10**(3)
delta = 1000     # Required position accuracy per unit time - 1km in m

#define time domain for pluto
H = 604800 * 248      # Size of "big steps" - 1 week in s - Pluto's year is 248 times larger than Earth's
N = 2           # do 5 revolutions around the sun
yr = 365.25*24*3600. * 248   # seconds in a yr - again with Pluto's year
T = N * yr         # time domain in s
Nsteps = int(T/H)

tpoints = arange(0,T,H)
xpoints = []
ypoints = []
r = array([Peri,0,0,Vperi],float)  #ics in [x,vx,y,vy]

# Do the "big steps" of size H
for t in tpoints:

    xpoints.append(r[0])
    ypoints.append(r[2])

    # Do one modified midpoint step to get things started
    n = 1
    r1 = r + 0.5*H*f(r)
    r2 = r + H*f(r1)

    # The array R1 stores the first row of the
    # extrapolation table, which contains only the single
    # modified midpoint estimate of the solution at the
    # end of the interval
    R1 = empty([1,4],float)
    R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

    # Now increase n until the required accuracy is reached
    error = 2*H*delta
    while error>H*delta:

        n += 1
        h = H/n

        # Modified midpoint method
        r1 = r + 0.5*h*f(r)
        r2 = r + h*f(r1)
        for i in range(n-1):
            r1 += h*f(r2)
            r2 += h*f(r1)

        # Calculate extrapolation estimates.  Arrays R1 and R2
        # hold the two most recent lines of the table
        R2 = R1
        R1 = empty([n,4],float)
        R1[0] = 0.5*(r1 + r2 + 0.5*h*f(r2))
        for m in range(1,n):
            epsilon = (R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1)
            R1[m] = R1[m-1] + epsilon
        error = abs(epsilon[0])

    # Set r equal to the most accurate estimate we have,
    # before moving on to the next big step
    r = R1[n-1]
    
# change xpoints and ypoints to AU
x_au = array(xpoints[:53])*6.68459e-12
y_au = array(ypoints[:53])*6.68459e-12
# Plot the results
plt.figure(figsize=[8,8])
plt.plot(x_au,y_au,'b.')
plt.xticks(size=16)
plt.yticks(size=16)
plt.xlabel('$x$ (AU)',size=18)
plt.ylabel('$y$ (AU)',size=18)
plt.xlim([-55,55])
#plt.ylim([-55,55])
plt.savefig('Q2b')
plt.show()

