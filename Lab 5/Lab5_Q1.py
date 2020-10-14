# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 12:04:34 2020

@author: Nicholas Pavanel
"""

#needed imports
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# a)

#load in the data
sun_data = np.loadtxt('PHY407_Lab5_sunspots.txt')

#plot the data
plt.figure(figsize=[8,6])
plt.plot(sun_data[:,0],sun_data[:,1])
plt.xlabel('t (months)',size=18)
plt.ylabel('n (# of sunspots)',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
#plt.savefig('1_a')
plt.show()
plt.close()

#import the discrete fourier transform function from page 5 of class notes
def dft(y):
    N = len(y)
    c = np.zeros(N//2+1, complex) # we only need the first 1/2 of the points
    for k in range(N//2+1): # let's do loops, pedagogy > speed today
        for n in range(N): # trapezoidal integration with y0=yend
            c[k] += y[n]*np.exp(-2j*np.pi*k*n/N)
    return c

#fourier transform the sunspot data with dft
fft_sun_data = dft(sun_data[:,1])

#get data ready to plot
sun_fft = np.abs(fft_sun_data[1:])**2

#plot the figure
plt.figure(figsize=[8,6])
plt.plot(sun_fft)
plt.xlabel('K Value',size=18)
plt.ylabel('Amplitude',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
#plt.xlim([-5,100])
plt.xscale('log')
#plt.savefig('1_b1')
#plt.yscale('log')
plt.show()
plt.close()

#to find the approximate value of k, we may argsort
peak_k = np.argsort(sun_fft)[-1]
print('K-value of peak:',peak_k)

#find the period of the sine wave with k=23
sun_sine = sun_fft[23] * np.sin( (2 * np.pi * peak_k * sun_data[:,0]) / (len(sun_data[:,0])) )

#plot sign curve
plt.figure(figsize=[10,6])
plt.plot(sun_sine)
plt.xlim([-2,200])
plt.xticks(size=16)
plt.yticks(size=16)
plt.xlabel('k=23 Sine Wave',size=18)
plt.ylabel('Amplitude',size=18)
plt.savefig('1_b2')
plt.show()
plt.close()

# b)

#load in the data
dow_data = np.loadtxt('PHY407_Lab5_dow.txt')

#define a function that smooths a data set by taking its Fourier transformation,
#choosing to keep a subset of the Fourier coefficients (co_per = what percentage to keep),
#and inverse Fourier transforming back with those coefficients
def smooth(data, co_per):
    fft_data = np.fft.rfft(data) # take a Fourier transformation of the data
    cutoff_index = int(co_per * len(data)) # determine cutoff index from percentage that we want t0 keep
    fft_data_subset = np.array(fft_data)
    fft_data_subset[cutoff_index:] = 0 # set all indicies past cutoff_index equal to zero
    ifft_data = np.fft.irfft(fft_data_subset) # take inverse Fourier transform with a subset of coefficients
    return ifft_data

#smooth dow data with smooth
smooth_dowdata_2per = smooth(dow_data,0.02)
smooth_dowdata_10per = smooth(dow_data,0.1)

#plot dow data w ifft'd data
plt.figure(figsize=[10,6])
plt.plot(dow_data,label='Dow Data',color='orange')
plt.plot(smooth_dowdata_10per,label='10% Non-zero FFT-IFFT Recreation',color='red')
plt.plot(smooth_dowdata_2per,label='2% Non-zero FFT-IFFT Recreation',color='blue')
plt.xlabel('t (days)',size=18)
plt.ylabel('Dow Closing Value',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
plt.legend()
plt.savefig('1_b3')
plt.show()
plt.close()


# c)

#load in the data
dow2_data = np.loadtxt('PHY407_Lab5_dow2.txt')

#smooth dow data 2 with smooth
smooth_dowdata2_2per = smooth(dow2_data,0.02)

#define a function that smooths a data set by taking its Fourier transformation,
#choosing to keep a subset of the Fourier coefficients (co_per = what percentage to keep),
#and inverse Fourier transforming back with those coefficients
def smooth_dct(data, co_per):
    dct_data = sp.fft.dct(data) # take a Fourier transformation of the data
    cutoff_index = int(co_per * len(data)) # determine cutoff index from percentage that we want t0 keep
    dct_data_subset = np.array(dct_data)
    dct_data_subset[cutoff_index:] = 0 # set all indicies past cutoff_index equal to zero
    idct_data = sp.fft.idct(dct_data_subset) # take inverse Fourier transform with a subset of coefficients
    return idct_data

#smooth dow data 2 with smooth_dct
smooth_dowdata2_2per_dct = smooth_dct(dow2_data,0.02)

#plot dow data w ifft'd data
plt.figure(figsize=[10,6])
plt.plot(dow2_data,label='Dow Data',color='orange')
plt.plot(smooth_dowdata2_2per,label='2% Non-zero FFT-IFFT Recreation',color='red')
plt.plot(smooth_dowdata2_2per_dct,label='2% Non-zero DCT-IDCT Recreation',color='blue')
#plt.plot(ifft_dow200_sp,label='20% Non-zero DCT-IDCT Recreation',color='black')
plt.xlabel('t (days)',size=18)
plt.ylabel('Dow2 Closing Value',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
plt.legend()
#plt.savefig('1_c')
plt.show()
plt.close()

# d)

#load data
trumpet_data = np.loadtxt('PHY407_Lab5_trumpet.txt')
piano_data = np.loadtxt('PHY407_Lab5_piano.txt')

#plot raw data
plt.figure(figsize=[11,6])
plt.plot(trumpet_data,label='Trumpet')
plt.plot(piano_data,label='Piano')
plt.legend()
plt.xticks(size=16)
plt.yticks(size=16)
plt.xlabel('Note Number',size=18)
plt.ylabel('Note Intensity',size=18)
#plt.savefig('1_d1')
plt.show()
plt.close()

#calculate the Fourier coefficients with numpy rfft
fft_trumpet = np.fft.rfft(trumpet_data)
fft_piano = np.fft.rfft(piano_data)

#plot fft data
plt.figure(figsize=[8,6])
plt.plot(fft_trumpet,label='Trumpet')
plt.plot(fft_piano,label='Piano')
plt.legend()
plt.xticks(size=16)
plt.yticks(size=16)
plt.xlabel('Coefficent Number',size=18)
plt.ylabel('Coefficent Magnitude',size=18)
#plt.xlim([750,2000])
#plt.savefig('1_d2')
plt.show()
plt.close()

#determine what note is being played
n =  np.argsort(fft_piano)[-1] * (44100) / 100000
print('Note being played is', n, 'Hz')
