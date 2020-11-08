#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 21:36:38 2020

@author: HayleyAgler
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
M = 100         # Grid squares on a side
V = 1.0         # Voltage at the plates
target = 10e-6 # Target accuracy

# Create arrays to hold potential values
phi = np.zeros([M+1,M+1],float) #voltage is 0 on all four walls
phi[20:80,20] = V #want column 20, rows 20 to 80 to have V #rows,cols
phi[20:80,80] = -V #want column 80, rows 20 to 80 to have -V


omega=0.5
count=0
# Main loop
delta = 1.0
while delta>target:
    count +=1
    temp_phi = phi.copy() #copy the previous potential    
    for i in range(M): # Calculate new values of the potential
        for j in range(M):
            if i==0 or i==M or j==0 or j==M: 
                phi[i,j] = 0  #voltage never changes on boundaries
            elif 20<=i<=80 and j==20:
                phi[i,j] = 1
            elif 20<=i<=80 and j==80:
                phi[i,j] = -1
            else:
                phi[i,j] = (1+omega)*(phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4 - omega*phi[i,j] 
                        #uses the four points around the current point
    delta=max(abs(phi - temp_phi)) # Calculate maximum difference from old values
    
    
print(count)



#contour plot
plt.figure(1)
plt.title("Potential of Capacitor")
plt.imshow(phi, vmax=1, vmin=-1)
plt.colorbar(label="V")
plt.xlabel("x [cm]")
plt.ylabel("y [cm]")
plt.savefig("potential", dpi=300)

#make arrays for the streamplot
x=np.arange(0,101,1)
y=np.arange(0,101,1)
X,Y= np.meshgrid(x, y)
Ey, Ex = np.gradient(-phi, y, x)


fig = plt.figure()
strm=plt.streamplot(X,Y, Ex, Ey, color=phi, cmap='autumn')
cbar = fig.colorbar(strm.lines)
cbar.set_label('Potential $V$')
plt.title('Electric field lines')
plt.xlabel('$x$ [cm]')
plt.ylabel('$y$ [cm]')
plt.show()
plt.savefig("efl", dpi=300)


