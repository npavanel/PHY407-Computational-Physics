
import numpy as np
import matplotlib.pyplot as plt
from time import time
import astropy.constants as c
import astropy.units as u
from scipy.linalg import eigh

#6.9 b)

#define constants
hbar = 1.0545718 * 10**(-34) # in J * s
L = 5 * 10**(-10) # in m
M = 9.1094 * 10**(-31) # in kg
a = 1.602176634 * 10**(-18) # in J

#define a function that compute H_{mn} for a given m and n
#below is modelled from page 4 of lab handout
def H_mn(m,n):
    if m != n and m % 2 == 0 and n % 2 == 0: 
        return 0
    elif m != n and m % 2 == 1 and n % 2 == 1:
        return 0
    elif m != n and m % 2 == 0 and n % 2 == 1 or m != n and m % 2 == 1 and n % 2 == 0:
        return - (8 * a * m * n) / (np.pi**2 * (m**2 - n**2)**2)
    elif m == n:
        return 0.5 * a + (np.pi**2 * hbar**2 * m**2) / (2 * M * L**2)
    
#6.9 c)

#modify H_mn to return an array of size nxn, with m,n up to n
def H_matrix(size):
    output = np.zeros((size,size))
    for i in range(1, size + 1):
        for j in range(1, size + 1):
            output[i-1, j-1] = H_mn(i,j) #using lab handout method
    return output

#compute H_mn for the 10x10 case
H_10_10 = H_matrix(10)

#compute eigenvalues and eigenvectors in KMS units? for 10x10 case
evals_H1010, evecs_H1010 = eigh(H_10_10)
evals_H1010=evals_H1010* 6.242e18

#6.9 d)
#compute H_mn for the 100x100 case
H_100_100 = H_matrix(100)

#compute eigenvalues and eigenvectors in KMS units? for 100x100 case
evals_H100100, evecs_H100100 = eigh(H_100_100)
#evals_H100100=evals_H100100#* 6.242e18
#evecs_H100100=evecs_H100100#* 6.242e18

#6.9 e

x=np.linspace(0, L, 100)
psi_0=[0]*100#np.zeros(100)
psi_1=np.zeros(100)
psi_2=np.zeros(100)
psi_n=np.zeros((100,100))

#%%

def psi_n(n, x):
    psi = 0
    for m in range(100):
        psi +=  evecs_H100100[n][m] * np.sin(np.pi * (m+1) * x / L)
    return psi



def psi0(x):
    return psi_n(0, x)

def psi1(x):
    return psi_n(1, x)

def psi2(x):
    return psi_n(2, x)


#%%
#normalizing 
def psi0_2(x):
    return psi0(x)**2

def psi1_2(x):
    return psi1(x)**2
def psi2_2(x):
    return psi2(x)**2

N1=8
a=0
b=L

h1=(b-a)/N1

s_trap1 = 0.5*psi0_2(a) + 0.5*psi0_2(b)
s_trap2 = 0.5*psi1_2(a) + 0.5*psi1_2(b)
s_trap3 = 0.5*psi2_2(a) + 0.5*psi2_2(b)

for k in range(1,N1):
    s_trap1 += psi0_2(a+k*h1)
    s_trap2 += psi1_2(a+k*h1)
    s_trap3 += psi2_2(a+k*h1)

    I1_trap1=s_trap1*h1
    I1_trap2=s_trap2*h1
    I1_trap3=s_trap3*h1
    
    
#%% Plotting
    
plt.title("Wave Functions")
plt.plot(x, psi0_2(x)/(np.sqrt(I1_trap1)), label='Ground')
plt.plot(x, psi1_2(x)/(np.sqrt(I1_trap1)), label='1st')
plt.plot(x, psi2_2(x)/(np.sqrt(I1_trap1)), label='2nd')
plt.xlabel("x values")
plt.ylabel('Probability')
plt.legend()
