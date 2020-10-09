
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
evals_H100100=evals_H100100* 6.242e18
evecs_H100100=evecs_H100100* 6.242e18

#6.9 e
x=np.linspace(0, 2*L, 100)
psi_0=np.zeros(100)
psi_1=np.zeros(100)
psi_2=np.zeros(100)
psi_n=np.zeros((100,100))

def psi_n(n,x):
    for i in range(len(evecs_H100100)):
        return np.sum(evecs_H100100[n][i]*np.sin((np.pi*n*x)/L))

for i in range(len(x)):
    psi_0[i]=psi_n(1,x[i])
    psi_1[i]=psi_n(2,x[i])
    psi_2[i]=psi_n(3,x[i])


#integral to find norm
def psi_int(n,x):
    return (x/2-(L*np.sin((2*np.pi*n*x)/L))/(4*np.pi*n))
#-L/(np.pi*n)*np.cos((np.pi*x*n)/L)

psi_0_norm=psi_int(1,L)-psi_int(1,0)
psi_0_normed=np.array(psi_0)/psi_0_norm

psi_1_norm=psi_int(2,L)-psi_int(2,0)
psi_1_normed=np.array(psi_1)/psi_1_norm

plt.plot(x, np.abs(psi_0_normed)**2)
plt.plot(x, (psi_1_normed)**2)




'''
def psi_n(n, x):
    psi0 = 0
    for m in range(100):
        psi0 += np.sqrt(2/L)*evecs_H100100[n][m]*np.sin(np.pi*(m+1)*x/L)
    return psi0

psi_0=psi_n(0,x)'''
