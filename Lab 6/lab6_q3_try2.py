#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 17:28:30 2020

@author: HayleyAgler
"""
#write pseudocode


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


vx_grid=np.zeros(4)
vy_grid=np.zeros(4)
vxx_grid, vyy_grid = np.meshgrid(vx_grid, vy_grid)
vx_initial, vy_initial=vxx_grid.flatten(), vyy_grid.flatten()

x_ij=np.zeros([16,16])
y_ij=np.zeros([16,16]) #np.meshgrid(np.zeros(16), np.zeros(16))

accx_i=np.zeros([16,16])
accy_i=np.zeros([16,16])

accx_=np.zeros(16)
accy_=np.zeros(16)

def accx(x,y):
    return (-12 * (x) * ((1/(((x)**2+(y)**2)**4))-2/(((x)**2+(y)**2)**7)))

def accy(x,y):
    return (-12 * (y) * ((1/(((x)**2+(y)**2)**4))-2/(((x)**2+(y)**2)**7)))

#finds the initial accelerations for all particles
#by finding the acceleration on particle i due to all other partiles j and then summing 
#to find total acceleration of particle i
for i in range(16):
    for j in range(16):
        x_ij[j][i] =x_initial[j]-x_initial[i]
        y_ij[j][i]=y_initial[j]-y_initial[i]
#for particle 1, acceleration due to the other particles is:
        if i==j:
            accx_i[i][j]=0
            accy_i[i][j]=0
            
        else:
            accx_i[i][j]=accx(x_ij[j][i],y_ij[j][i])
            accy_i[i][j]=accy(x_ij[j][i],y_ij[j][i])
            
        accx_[i]=np.sum(accx_i[i])
        accy_[i]=np.sum(accy_i[i])




r1=np.zeros([2,N,1000]) #x matrix is 16x1000, y matrix is 16x1000
acc_x_matrix=np.zeros([N,N,1000])
acc_y_matrix=np.zeros([N,N,1000])
for i in range(N):
    r1[0,i,0]=x_initial[i] #sets first element in every row of x matrix to the x0
    r1[1,i,0]=y_initial[i]
    k1=np.zeros([2,N,1000])
    v1=np.zeros([2,N,1000]) #vx,vy
    v1_half=np.zeros([2,N,1000]) 
    for i in range(N):
        v1_half[0,i,0], v1_half[1,i,0]=np.array([0,0])+1/2*h*np.array(accx_[i], accy_[i])#sets inital vx0,vy0
        acc_x_sum=np.zeros([16,1000])
        acc_y_sum=np.zeros([16,1000])
        acc_x_sum[i,0]=accx_[i]
        acc_y_sum[i,0]=accy_[i]
        
    for i in range(len(t)-1):
        for j in range(N):
            for k in range(N):
                r1[0,j,i+1], r1[1,j,i+1] = r1[:,j,i] + h * np.array(v1_half[:,j,i])

                if k==j:
                    acc_x_matrix[k][j][i]=0
                    acc_y_matrix[k][j][i]=0
                else:
                    acc_x_matrix[k][j][i]=accx(r1[0,k,i]-r1[0,j,i],r1[1,k,i]-r1[1,j,i])
                    acc_y_matrix[k][j][i]=accy(r1[0,k,i]-r1[0,j,i],r1[1,k,i]-r1[1,j,i])
                acc_x_sum[k,i+1]=np.sum(acc_x_matrix[k,j,i+1])
                acc_y_sum[k,i+1]=np.sum(acc_y_matrix[k,j,i+1])
                    
                
                k1[:,j,i+1] = h * np.array([acc_x_sum[j][i+1], acc_y_sum[j][i+1]])
                v1[:,j,i+1] = v1_half[:,j,i] + 0.5 * k1[:,j,i+1]
                v1_half[:,j,i+1] = v1_half[:,j,i] + k1[:,j,i+1]
#%%




rx=r1[0]
ry=r1[1]
                
                
for i in range(16):
    plt.scatter(rx[i],ry[i], s=1)