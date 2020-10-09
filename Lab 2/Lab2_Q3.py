# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 20:04:20 2020

@author: Nicholas Pavanel
"""
#import needed modules
import numpy as np
import matplotlib.pyplot as plt

#3c

#define alpha
omega=10*20 #total width of the grating in um
alpha=np.pi/(20) #slit separation of 20um 
grating=np.linspace(-omega/2,omega/2,300)
screen=np.linspace(-0.05,0.05,500)
#define the intensity transmission function of the diffraction grating at distance u from the central axis
def q(u):
    return (np.sin(alpha*u))**2
plt.figure(figsize=[8,6])
plt.plot(grating,q(grating),color='black')
plt.title('Transmission Intensity Profile 1',fontsize=18)
plt.ylabel(r'$\rm u(x) $',fontsize=18)
plt.xlabel(r'$\rm x \ (um)$',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('3c_grating')
plt.show()

#define a function that computes (*)
def f(x,u):
    return np.sqrt(q(u))*np.exp((1j*2*np.pi*x*u)/(wl*fl))

#define needed constants - all units in micrometers
wl=0.5 #wavelength: 500nm in um
fl=1 #focal length: 1m in m

def simpson(xs, u, a, b, n):
    intensity=[]
    f_a=0
    f_b=0
    for i in u:
        f_a+=f(a,i)
        f_b+=f(b,i)
    
    for k in xs:
        h=(b-a)/n
        odd_sum=0.0
        even_sum=0.0
        
        step=a + h
        for i in range(1,n,2):
            odd_sum += f(k,step)
            step += 2*h

        step = a + 2*h
        for i in range(2,n,2):
            even_sum += f(k,step)
            step += 2*h
        
        simp_sum = (h/3)*(f_a+f_b+4*odd_sum+2*even_sum)
        intensity.append(np.abs(simp_sum)**2)
            
    return np.array(intensity) * 10**(-12)

intensity=simpson(screen,grating,-omega/2,omega/2,250)

plt.figure(figsize=[8,6])
plt.plot(screen,intensity,color='green')
plt.ylabel(r'$\rm I (x) $',fontsize=18)
plt.xlabel(r'$\rm x \ (m)$',fontsize=18)
plt.title('Diffraction Intensity 1',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('q3_c_1')

#3d

vis=np.tile(intensity,(400,1))

plt.figure(figsize=[8,5])
plt.imshow(vis)
plt.xlabel(r'$\rm x \ (m)$',fontsize=18)
plt.yticks([0],'.')
plt.xticks([83,166,250,333,416],labels=['-0.04','-0.02','0','0.02','0.04'])
plt.title('Diffraction Visualization 1',fontsize=18)
plt.colorbar(label=r'$\rm I(x) $')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('q3_d')

#3ei

#define alpha
omega=10*20 #total width of the grating in m
alpha=np.pi/(20) #slit separation of 20um in m
beta=alpha/2
grating=np.linspace(-omega/2,omega/2,300)
screen=np.linspace(-0.05,0.05,500)
#define the intensity transmission function of the diffraction grating at distance u from the central axis
def q_2(u):
    return np.sin(alpha*u)**2 * np.cos(beta*u)**2
plt.figure(figsize=[8,6])
plt.plot(grating,q_2(grating),color='black')
plt.title('Transmission Intensity Profile 2',fontsize=18)
plt.ylabel(r'$\rm u(x) $',fontsize=18)
plt.xlabel(r'$\rm x \ (um)$',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('3ei_grating')
plt.show()

#define a function that computes (*)
def f_2(x,u):
    return np.sqrt(q_2(u))*np.exp((1j*2*np.pi*x*u)/(wl*fl))

def simpson_2(xs, u, a, b, n):
    intensity=[]
    f_a=0
    f_b=0
    for i in u:
        f_a+=f_2(a,i)
        f_b+=f_2(b,i)
    
    for k in xs:
        h=(b-a)/n
        odd_sum=0.0
        even_sum=0.0
        
        step=a + h
        for i in range(1,n,2):
            odd_sum += f_2(k,step)
            step += 2*h

        step = a + 2*h
        for i in range(2,n,2):
            even_sum += f_2(k,step)
            step += 2*h
        
        simp_sum = (h/3)*(f_a+f_b+4*odd_sum+2*even_sum)
        intensity.append(np.abs(simp_sum)**2)
            
    return np.array(intensity) * 10**(-12)

intensity2=simpson_2(screen,grating,-omega/2,omega/2,250)
vis2=np.tile(intensity2,(400,1))

plt.figure(figsize=[8,6])
plt.plot(screen,intensity2,color='green')
plt.ylabel(r'$\rm I (x) $',fontsize=18)
plt.xlabel(r'$\rm x \ (m)$',fontsize=18)
plt.title('Diffraction Intensity 2',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('q3_ei_1')
plt.show()
plt.close()
plt.figure(figsize=[8,5])
plt.imshow(vis2)
plt.xlabel(r'$\rm x \ (m)$',fontsize=18)
plt.yticks([0],'.')
plt.xticks([83,166,250,333,416],labels=['-0.04','-0.02','0','0.02','0.04'],fontsize=14)
plt.title('Diffraction Visualization 2',fontsize=18)
plt.colorbar(label=r'$\rm I(x) $')
plt.savefig('q3_ei_2')

#3eii

#define alpha
omega=10*20 #total width of the grating in m
grating=np.linspace(-omega/2,omega/2,200)
screen=np.linspace(-0.05,0.05,500)
#to get an idea, i make an example grating
u3=np.zeros(len(grating)) #create an array of zeros the length of the grating
u3[55:65]=1 #55 um of no transmission up to the first slit of 10um with transmission set = 1
u3[65:125]=0 #60um of transmission between slits
u3[125:145]=1 #second slit of width 20um
#define transmission function
def q_3(a):
    if a > grating[55] and a < grating[65]:
        return 1
    elif a > grating[125] and a < grating[145]:
        return 1
    else:
        return 0

#define alpha
omega=10*20 #total width of the grating in m
grating=np.linspace(-omega/2,omega/2,200)
screen=np.linspace(-0.05,0.05,500)
#to get an idea, i make an example grating
u3=np.zeros(len(grating)) #create an array of zeros the length of the grating
u3[55:65]=1 #55 um of no transmission up to the first slit of 10um with transmission set = 1
u3[65:125]=0 #60um of transmission between slits
u3[125:145]=1 #second slit of width 20um
#define transmission function
def q_3(a):
    if a > grating[55] and a < grating[65]:
        return 1
    elif a > grating[125] and a < grating[145]:
        return 1
    else:
        return 0

plt.figure(figsize=[8,6])
plt.plot(grating,u3,color='black')
plt.title('Transmission Intensity Profile 3',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel(r'$\rm u(x) $',fontsize=18)
plt.xlabel(r'$\rm x \ (um)$',fontsize=18)
plt.savefig('3eii_grating')
plt.show()

#define a function that computes (*)
def f_3(x,u):
    return np.sqrt(q_3(u))*np.exp((1j*2*np.pi*x*u)/(wl*fl))

def simpson_3(xs, u, a, b, n):
    intensity=[]
    f_a=0
    f_b=0
    for i in u:
        f_a+=f_3(a,i)
        f_b+=f_3(b,i)
    
    for k in xs:
        h=(b-a)/n
        odd_sum=0.0
        even_sum=0.0
        
        step=a + h
        for i in range(1,n,2):
            odd_sum += f_3(k,step)
            step += 2*h

        step = a + 2*h
        for i in range(2,n,2):
            even_sum += f_3(k,step)
            step += 2*h
        
        simp_sum = (h/3)*(f_a+f_b+4*odd_sum+2*even_sum)
        intensity.append(np.abs(simp_sum)**2)
            
    return np.array(intensity) * 10**(-12)

intensity3=simpson_3(screen,grating,-omega/2,omega/2,250)
vis3=np.tile(intensity3,(400,1))

plt.figure(figsize=[8,6])
plt.plot(screen,intensity3,color='green')
plt.xlabel(r'$\rm x \ (m)$',fontsize=18)
plt.ylabel(r'$\rm I(x) $',fontsize=18)
plt.title('Diffraction Intensity 3',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('q3_eii_1')
plt.show()
plt.close()
plt.figure(figsize=[8,5])
plt.imshow(vis3)
plt.xlabel(r'$\rm x \ (m)$',fontsize=18)
plt.yticks([0],'.')
plt.xticks([83,166,250,333,416],labels=['-0.04','-0.02','0','0.02','0.04'],fontsize=14)
plt.title('Diffraction Visualization 3',fontsize=18)
plt.colorbar(label=r'$\rm I(x) $')
plt.savefig('q3_eii_2')
