# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 16:51:39 2020

@author: Nicholas Pavanel
"""
#a)

#Pseudocode

#creating w(x,y)
#open file
#initialize an 2d array of zeros size 1201x1201
#create a double for loop over the length of the tile: 1201
#for each iteration of the loop: 
    #read a line of the file
    #unpacked the line 
    #set the correct location in the zeros array to the unpacked line value
        #this will be arcoss the first row of the array, then the second, and so on
        #the double for loop travels across the top row, then the second, and so on
    
#calculating the gradient of w
#create a double for loop that will run over the rows of w(x,y)
#initialize an list that will hold the derivatives at each point
#for each iteration (point in w(x,y)), append to the initialized list a tuple: (dw/dx,dw/dy)
#at each iteration of the outside loop, a new list is initalized to mimic w(x,y)
#for the the first point, use a forward dif scheme, middle points central, end backwards

#create plots of w and I
#imshow the array for w
#compute I from the eqn at the bottom of page 212 in the text (assume light shining horizontally w unit intensity)
#imshow I

#b)

#import needed modules
import struct
import numpy as np
import matplotlib.pyplot as plt
from time import time

#initialize w to store the elevation data
w = np.zeros((1201,1201))
whis = np.zeros((1201,1201))
haw1 = np.zeros((1201,1201))
haw2 = np.zeros((1201,1201))

#double loop through the import data
f = open('N46E006.hgt','rb')
for i in range(1201): 
    for j in range(1201): 
        buf = f.read(2)
        val = struct.unpack('>h',buf)[0]
        if val > 25000 or val < - 25000:
            w[i][j] = 0
        else:
            w[i][j] = val #set w[0][0] northwesternmost point, w[-1][-1] southeasternmost point
            
#set x labels
x = np.array([0,200,400,600,800,1000,1200])
x = x*420/1000
#plot elevation data w
plt.figure(figsize=[12,10])
plt.xticks([0,200,400,600,800,1000,1200],x,fontsize=14)
plt.yticks([1200,1000,800,600,400,200,0],x,fontsize=14)
plt.ylabel('x (km)',fontsize=18)
plt.xlabel('x (km)',fontsize=18)
plt.imshow(w)
plt.colorbar(label='y (m)')
#plt.savefig('Q3_b_1')

"""#get hawaii1 data
f = open('N21W159.hgt','rb')
for i in range(1201): 
    for j in range(1201): 
        buf = f.read(2)
        val = struct.unpack('>h',buf)[0]
        if val > 25000 or val < - 25000:
            haw1[i][j] = 0
        else:
            haw1[i][j] = val #set w[0][0] northwesternmost point, w[-1][-1] southeasternmost point
            
#get hawaii2 data
f = open('N21W158.hgt','rb')
for i in range(1201): 
    for j in range(1201): 
        buf = f.read(2)
        val = struct.unpack('>h',buf)[0]
        if val > 25000 or val < - 25000:
            haw2[i][j] = 0
        else:
            haw2[i][j] = val #set w[0][0] northwesternmost point, w[-1][-1] southeasternmost point"""
            
"""#subplot hawaii
fig = plt.figure(figsize=[12,10])
ax1 = fig.add_axes([0,0.5,0.5,0.5])
ax2 = fig.add_axes([0.42,0.5,0.5,0.5],yticklabels=[])
ax1.imshow(haw1)
ax2.imshow(haw2)"""

#define derivative methods
def fwd_dif(fx_h,fx,h):
    result=0
    result += (fx_h - fx)/(h)
    return result

def bwd_dif(fx_h,fx,h):
    result=0
    result += (fx - fx_h)/(h)
    return result

#calculate the derivatives of w 
#set conditions
h = 420
phi = (-5 * np.pi) / 6
#initialize the list that will hold the derivatives
dw = []
for i in range(1201):
    hold=[]
    for j in range(1201):
        if i == 0 and j == 0:
            dwdx = fwd_dif(w[i][j+1], w[i][j], h)
            dwdy = fwd_dif(w[i+1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        elif i == 0 and j == 1200:
            dwdx = bwd_dif(w[i][j-1], w[i][j], h)
            dwdy = fwd_dif(w[i+1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        elif i == 1200 and j == 0:
            dwdx = fwd_dif(w[i][j+1], w[i][j], h)
            dwdy = bwd_dif(w[i-1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        elif i == 1200 and j == 1200:
            dwdx = bwd_dif(w[i][j-1], w[i][j], h)
            dwdy = bwd_dif(w[i-1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        elif i == 0:
            dwdx = fwd_dif(w[i][j+1], w[i][j], h)
            dwdy = fwd_dif(w[i+1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        elif j == 0:
            dwdx = fwd_dif(w[i][j+1], w[i][j], h)
            dwdy = fwd_dif(w[i+1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        elif i == 1200:
            dwdx = fwd_dif(w[i][j+1], w[i][j], h)
            dwdy = bwd_dif(w[i-1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        elif j == 1200:
            dwdx = bwd_dif(w[i][j-1], w[i][j], h)
            dwdy = fwd_dif(w[i+1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
        else:
            dwdx = fwd_dif(w[i][j+1], w[i][j], h)
            dwdy = fwd_dif(w[i+1][j], w[i][j], h)
            partials = (dwdx,dwdy)
            hold.append(partials)
    dw.append(hold)
    
#note that dw has dimensions below:
#dw[i][j][k] - i determines the row, j the position in the row (column), k dwdx[0] or dwdy[1]

#calculate I
I = np.zeros((1201,1201))
for i in range(1201):
    for j in range(1201):
        I[i][j] += - (np.cos(phi) * dw[i][j][0] + np.sin(phi) * dw[i][j][1]) / np.sqrt(dw[i][j][0]**2 + dw[i][j][1]**2 + 1)
        
# fix broken values
for i in range(1201):
    for j in range(1201):
        if I[i][j] > 0.025 or I[i][j] < -0.025:
            I[i][j] = 0
        else:
            continue
        
plt.figure(figsize=[12,10])
plt.imshow(I)
plt.xticks([0,200,400,600,800,1000,1200],x,fontsize=14)
plt.yticks([1200,1000,800,600,400,200,0],x,fontsize=14)
plt.ylabel('x (km)',fontsize=18)
plt.xlabel('x (km)',fontsize=18)
plt.colorbar(label='I(x)')
#plt.savefig('Q3_b_2')
plt.show()
plt.close()

"""#set y labels
y = [1200,1000,800,600,400,200,0,0]
y = np.array(y)*420/1000
#set x labels
x1 = np.array([0,0,200,400,600,800,1000])
x1 = x1*420/1000
x2 = np.array([1200,1200,1400,1600,1800,2000,2200,2400])
x2 = x2*420/1000"""

"""#subplot hawaii
fig = plt.figure(figsize=[12,10])
ax1 = fig.add_axes([0,0.5,0.5,0.5])
ax2 = fig.add_axes([0.42,0.5,0.5,0.5],yticklabels=[])

ax1.imshow(I_haw1)
ax2.imshow(I_haw2) 
ax1.set_xticklabels(x1,fontsize=14)
ax1.set_yticklabels(y,fontsize=14)
ax2.set_xticklabels(x2,fontsize=14)
ax1.set(ylabel= 'x (km)')
ax1.set(xlabel="                                                                                                                      x (km)")
plt.savefig('Q3_haw')"""
