#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:26:08 2020

@author: HayleyAgler
"""
import numpy as np
from numpy import empty,copy
import matplotlib.pyplot as plt

#Lab_4_Q1_c
#define constants
R1=R3=R5=1000 # resistance in ohms
R2=R4=R6=2000 # ohms
C1=1e-6 # capacitance in F
C2=0.5e-6 # F
x_p=3 #volts
omega=1000 # rad/s
L=R6/omega
#define complex arrays
v=np.array([x_p/R1,x_p/R2,x_p/R3], complex)
#[x1,x2,x3]
A=np.array([[1/R1+1/R4+1j*omega*C1, -1j*omega*C1, 0],
         [-1j*omega*C1, 1/R2+1/R5+(1j*omega*C1)+1j*omega*C2, -1j*omega*C2],
            [0, -1j*omega*C2, 1/R3+1/(R6)+1j*omega*C2]], complex)


##from part 1

N=len(v)

def partialpivot(A,v):
    #want it to find the row with the largest 1st/diagonal element and switch that row with that one
    #for first row, first element, compare first element of each row, biggest one switches with first row
    #then the gaussian happens again
    #then second row first element should be 0, so find largest diagonal element in each row and swap with second row
    for i in range(N): #!= means not equal
        max_index=A[i:,i].argmax()+i #finds index of biggest element in ith column
        #if that index isnt the index of the diagonal element in the row i, switch the rows
        if max_index !=i: #if not equal
            #row i, max row = #max row, row i
            A[i, :], A[max_index, :] = copy(A[max_index, :]), copy(A[i, :])
            #A[[i,max_index]] =  A[[max_index, i]]
            #switch rows for the rhs too
            #v[[i,max_index]] =  v[[max_index, i]]
            v[i], v[max_index] = copy(v[max_index]), copy(v[i])



# Gaussian elimination
for m in range(N):
# Divide by the diagonal element
    
    partialpivot(A,v)
    div = A[m,m] #sets div=to the diagonal element of the matrix
    A [m, : ] /= div #divides each element in the mth row by the diagonal element
    v[m] /= div  #divides the mth element of the right hand side by the diagnoal element also
    
    # Now subtract from the lower rows
    for i in range(m+1,N):
        mult = A [i ,m] #mth element of next row down from mth row (row,col)
        #subtract second row by mult*first row
        A[i,:] -= mult*A[m,:] #c-=a -> c=c-a #the next row from m is equal to taht row minus the first element*row above
        v[i] -= mult*v[m] #same as above but for rhs


# Backsubstitution
x = empty(N,complex)
for m in range(N-1,-1,-1):
    x[m] = v[m]
    for i in range(m+1,N):
        x [m] -= A [m, i] *x [i]
print(x) 


#%%

#calculate voltage amplitudes
t=0
V1=np.abs(x[0]*np.exp(1j*omega*t)) #Volts
V2=np.abs(x[1]*np.exp(1j*omega*t))
V3=np.abs(x[2]*np.exp(1j*omega*t))

ph1=np.degrees(np.arctan(x[0].imag/x[0].real))#degrees
ph2=np.degrees(np.arctan(x[1].imag/x[1].real))
ph3=np.degrees(np.arctan(x[2].imag/x[2].real))

print("Amplitudes", V1,V2,V3, "Phases",ph1,ph2,ph3)

#%%
t2=np.linspace(0,0.02,1000)
#[0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.0011,0.0012,0.0013,0.0014,0.0015]

Vr1=[0]*len(t2)
Vr2=[0]*len(t2)
Vr3=[0]*len(t2)

for i in range(len(t2)):
    Vr1[i]=x[0]*np.exp(1j*(omega)*t2[i])
    Vr2[i]=x[1]*np.exp(1j*(omega)*t2[i])
    Vr3[i]=x[2]*np.exp(1j*(omega)*t2[i])

plt.figure(1)
plt.title("Voltage vs. Time with Inductor")
plt.plot(t2, np.real(Vr1), label='V1')
plt.plot(t2, np.real(Vr2), label='V2')
plt.plot(t2, np.real(Vr3), label='V3')
plt.xlabel("Time (s)")
plt.ylabel("Real Voltage")
plt.xlim(0,0.02)
plt.legend()





    

