# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 09:48:43 2020

@author: Nicholas Pavanel
"""

#import needed modules
import numpy as np
import matplotlib.pyplot as plt
from time import time

#choose values of n for the first computation
ns=[2,4,6,8,10,15,20,25,30,45,50,60,70,80,90,100,150,200,250,300,400,500,700,900,1100,1300,1400,1500]

#create lists that will hold time information
times_mm=[]
times_nd=[]
times_mc=[]

#loop through all chosen values of n
for i in ns:
    start1=time() #set a start time for the matrix creation (unneeded)
    A = np.ones([i,i],float)*7
    B = np.ones([i,i],float)*4 #set up matricies
    matrix_creation=time()-start1 #compute unused matrix creation time
    start2=time() #set time that will compute matrix multiplication
    C = A*B #multiply matricies 
    matrix_multiplication=time()-start2 #compute time of matrix multiplication
    times_mc.append(matrix_creation) #append time for value of n
    times_mm.append(matrix_multiplication) #append time for value of n
    start3=time() #set time that will compute numpy dot calculation
    D = np.dot(A,B) #compute numpy dot calculation
    matrix_dot=time()-start3 #compute time of numpy dot calculation
    times_nd.append(matrix_dot) #append time for value of n
    
#the same process but for n^3
#values of n^3 above 30 require too much disk space
#the total computation time should be approximately 5minutes
ns3=[2**3,4**3,6**3,8**3,10**3,15**3,20**3,25**3,30**3]
times_mm3=[]
times_nd3=[]
times_mc3=[]
for i in ns3:
    start4=time()
    A3 = np.ones([i,i],float)*7
    B3 = np.ones([i,i],float)*4
    matrix_creation3=time()-start4
    start5=time()
    C3 = A3*B3
    matrix_multiplication3=time()-start5
    times_mc3.append(matrix_creation3)
    times_mm3.append(matrix_multiplication3)
    start6=time()
    D3 = np.dot(A3,B3)
    matrix_dot3=time()-start6
    times_nd3.append(matrix_dot3)
    
#plot a comparison between numpy matrix multiplication and numpy dot to n=1500
plt.plot(ns,times_mm,color='orange',label='n - Matrix Multiplication')
plt.plot(ns,times_nd,color='blue',label='n - Numpy Dot Product')
plt.xlabel('Dimension of Square Matrix [n]')
plt.ylabel('Time Taken For Computation [s]')
plt.title('Matrix Operation Time Trials')
plt.legend()
plt.show()
plt.close()
#plot a comparison between numpy matrix multiplication and numpy dot to n=27000
plt.plot(ns3,times_mm3,color='red',label='$n^{3}$ - Matrix Multiplication')
plt.plot(ns3,times_nd3,color='black',label='$n^{3}$ - Numpy Dot Product')
plt.xlabel('Dimension of Square Matrix [n]')
plt.ylabel('Time Taken For Computation [s]')
plt.title('Matrix Operation Time Trials')
plt.legend()

"""
Clearly numpy array multiplication is faster than the equivalent matrix 
multiplication using numpy dot. Although the increased speed of computation is 
practically negligible for square matricies with 3000 columns or less 
(approximately half a second), the increased speed of compuatation is great 
for matricies with 3000 columns or more. At the extreme, multiplying two 
square matricies with 27000 columns is 25 times faster using numpy arrays than 
it is with numpy dot, a difference of 279 seconds, or 4 and a half minutes; 
numpy array calculations take approximately 11 seconds, whereas numpy dot 
calculations take almost 5 minutes.
"""
    
