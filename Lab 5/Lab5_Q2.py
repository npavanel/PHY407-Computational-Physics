#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:13:19 2020

@author: HayleyAgler
"""
#import necessart modules
from numpy.fft import rfft2, irfft2
import numpy as np
import matplotlib.pyplot as plt

blur=np.loadtxt("/Users/owner/407/blur.txt")

#define values for point spread function
sigma=25
rows, cols = blur.shape
gauss=np.empty([1024,1024], float)

#routine from the lab to create a 2d array of point spread values
for i in range(rows):
    ip = i
    if ip > rows/2:
        ip -= rows # bottom half of rows moved to negative values
    for j in range(cols):
        jp = j
        if jp > cols/2:
            jp -= cols # right half of columns moved to negative values
            
        gauss[i, j] = np.exp(-(ip**2 + jp**2)/(2.*sigma**2)) # compute gaussian
        
#fourier transform the blurry photo and point spread function     
blur_trans=rfft2(blur)
gauss_trans=rfft2(gauss)

#define empty array for fourier trans of the division to be stored in
img_trans=np.empty(gauss_trans.shape, complex)

#divide one by the other, except when gauss_trans is small
epsilon=1e-3
for i in range(len(gauss_trans)):
    for j in range(len(gauss_trans.T)):
        if np.abs(gauss_trans[i,j])>epsilon:
            img_trans[i,j]=blur_trans[i,j]/gauss_trans[i,j]
        else:
            img_trans[i,j]=blur_trans[i,j]


#inverse transform the transform of pure image
img=irfft2(img_trans)

#%% plotting 
plt.figure(1)
plt.imshow(blur, cmap="Greys_r")
plt.figure(2)
plt.imshow(gauss, cmap="Greys_r", vmax=0.1)
plt.figure(3)
plt.imshow(img, cmap="Greys_r")
