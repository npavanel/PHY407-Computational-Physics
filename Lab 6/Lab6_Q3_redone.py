#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 13:03:35 2020

@author: HayleyAgler
"""

import numpy as np
import matplotlib.pyplot as plt


dt = 0.01
T = 1000
t1 = 0
t2 = 0 + dt*T
h = (t2-t1)/T
t = np.arange(t1,t2,h) 

N = 16
Lx = 4.0
Ly = 4.0
dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)
x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)
xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()



#for each paticle i, calculate the distance rij between particle i and every particle j
#calculate the acceleration of particle i due to each particle j then sum up to get total acc for particle i
#then use verlet to advance positions and velocities 

def radius(x_initial, y_initial):
    x_ij=np.zeros([16,16])
    y_ij=np.zeros([16,16])
    for i in range(len(x_initial)):
        for j in range(len(x_initial)):
             x_ij[j][i] =x_initial[j]-x_initial[i]
             y_ij[j][i]=y_initial[j]-y_initial[i]
    return x_ij, y_ij
        

def acceleration(x_ij, y_ij):
    acc_x=np.zeros([16,16])
    acc_y=np.zeros([16,16])
    acc_x_tot=np.zeros(16)
    acc_y_tot=np.zeros(16)
    for i in range(16):
        for j in range(16):
            if i==j:
                acc_x[j][i]=0
                acc_y[j][i]=0
            else:
                acc_x[j][i]=(-12 * (x_ij[j][i]**2) * ((1/(((x_ij[j][i])**2+(y_ij[j][i])**2)**4))-2/(((x_ij[j][i])**2+(y_ij[j][i])**2)**7)))
                acc_y[j][i]=(-12 * (y_ij[j][i]) * ((1/(((x_ij[j][i])**2+(y_ij[j][i])**2)**4))-2/(((x_ij[j][i])**2+(y_ij[j][i])**2)**7)))
                acc_x_tot[j]+=acc_x[j][i]
                acc_y_tot[j]+=acc_y[j][i]
    return acc_x_tot,acc_y_tot
#make this be valid for j
x_ij, y_ij=radius(x_initial,y_initial)

print(acceleration(x_ij,y_ij))


def potential(x_ij,y_ij): #where x_ij, y_ij are distances between two particles like in acceleration
    pot_x=np.zeros([16,16])
    pot_y=np.zeros([16,16])
    pot_x_tot=np.zeros(16)
    pot_y_tot=np.zeros(16)
    for i in range(16):
        for j in range(16):
            if i==j:
                pot_x[j][i]=0
                pot_y[j][i]=0
            else:
                pot_x[j][i]=4*(1/(x_ij[j][i]**12)-1/(x_ij[j][i]**6))
                pot_x[j][i]=4*(1/(y_ij[j][i]**12)-1/(y_ij[j][i]**6))
            pot_x_tot[j]=np.sum(pot_x[j])
            pot_y_tot[j]=np.sum(pot_y[j])
    return pot_x_tot, pot_y_tot
            
                

                
x_pos=np.zeros([16,1000])
y_pos=np.zeros([16,1000])
x_pos[:,0]=x_initial #sets initial x positions as first column
y_pos[:,0]=y_initial
kx=np.zeros([16,1000])
ky=np.zeros([16,1000])
vx=np.zeros([16,1000])
vy=np.zeros([16,1000])
vx_half=np.zeros([16,1000])
vy_half=np.zeros([16,1000])
pot_ex=np.zeros([16,1000])
pot_ey=np.zeros([16,1000])
kin_ex=np.zeros([16,1000])
kin_ey=np.zeros([16,1000])
kin_ex_v=np.zeros([16,1000])
kin_ey_v=np.zeros([16,1000])


#for the first step only:
#set first element of every row of vx_half to 
vx_half[:,0]=0+1/2*h*acceleration(x_ij,y_ij)[0]
vy_half[:,0]=0+1/2*h*acceleration(x_ij,y_ij)[1]

#now do verlet
for i in range(1000-1):
    for j in range(16):
        x_pos[j,i+1]=x_pos[j,i]+h*vx_half[j,i]
        y_pos[j,i+1]=y_pos[j,i]+h*vy_half[j,i]
        
        kx[:,i+1]=h*np.array(acceleration(radius(x_pos[:,i],y_pos[:,i])[0], radius(x_pos[:,i],y_pos[:,i])[1]))[0]
        ky[:,i+1]=h*np.array(acceleration(radius(x_pos[:,i],y_pos[:,i])[0], radius(x_pos[:,i],y_pos[:,i])[1]))[1]

#verify these steps are doing what wwe want
        vx[j,i+1]=vx_half[j,i]+0.5*kx[j,i+1]
        vy[j,i+1]=vy_half[j,i]+0.5*ky[j,i+1]
        
        vx_half[j,i+1]=vx_half[j,i]+kx[j,i+1]
        vy_half[j,i+1]=vy_half[j,i]+ky[j,i+1]
        
        pot_ex[:,i]=potential(radius(x_pos[:,i],y_pos[:,i])[0], radius(x_pos[:,i],y_pos[:,i])[1])[0]
        pot_ey[:,i]=potential(radius(x_pos[:,i],y_pos[:,i])[0], radius(x_pos[:,i],y_pos[:,i])[1])[1]
        
        kin_ex_v[j,i]=i+h*vx[j,i]
        kin_ey_v[j,i]=i+h*vy[j,i]
        kin_ex[j,i]=0.5*1* kin_ex_v[j,i]**2
        kin_ey[j,i]=0.5*1* kin_ey_v[j,i]**2

#%%
total_kin_j=kin_ex+kin_ey
total_pot_j=pot_ex+pot_ey

total_kin=np.sum(total_kin_j,axis=0)
total_pot=np.sum(total_pot_j,axis=0)
#%%
for i in range(16):
    plt.plot(x_pos[i],y_pos[i])
    

#%%
#plt.plot(x_pos[1],y_pos[1])
plt.plot(t,total_kin)
plt.plot(t,total_pot)
plt.plot(t,total_kin+total_pot)


#%%





