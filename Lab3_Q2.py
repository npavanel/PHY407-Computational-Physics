# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:13:14 2020

@author: Nicholas Pavanel
"""
#import needed modules
import numpy as np
import matplotlib.pyplot as plt

#to check function is working properly, define function for first 5 Hermite polynomials
#h_0 = 1
def h_1(x):
    return 2 * x
def h_2(x):
    return 4 * x**2 - 2
def h_3(x):
    return 8 * x**3 - 12 * x
def h_4(x):
    return 16 * x**4 * - 48 * x**2 + 12
def h_5(x):
    return 32 * x**5 - 160 * x**3 + 120 * x

#2a)

#define a recursive function to compute any Hermite polynomial at x
#for understanding: substitute n=n-1 into eqn 6 from lab handout
def h_poly(n,x):
    if n < 0:
        return -1
    elif n == 0:
        return 1
    elif n == 1:
        return 2 * x
    else:
        return 2 * x * h_poly(n-1,x) - 2 * (n-1) * h_poly(n-2,x)
    
#check that h_poly(n,x) works
print('1st Hermite Polynomial @ x=2 -','h_poly(1,2):',h_poly(1,2),'Answer:',h_1(2))
print('2nd Hermite Polynomial @ x=3 -','h_poly(2,3):',h_poly(2,3),'Answer:',h_2(3))
print('3rd Hermite Polynomial @ x=5 -','h_poly(3,5):',h_poly(3,5),'Answer:',h_3(5))
print('5th Hermite Polynomial @ x=15 -','h_poly(5,15):',h_poly(5,15),'Answer:',h_5(15))

#define a function to compute the harmonic oscillator wavefunction psi
def psi(n,x):
    return (1/np.sqrt(2**n * np.math.factorial(n) * np.sqrt(np.pi))) * np.exp((-x**2)/(2)) * h_poly(n,x) 

#plot the harmonic oscillator wavefunction for a range of ns on -4<=x<=4
x_a=np.linspace(-4,4,100)
plt.figure(figsize=[8,6])
plt.plot(x_a,psi(0,x_a),label='n=0')
plt.plot(x_a,psi(1,x_a),label='n=1')
plt.plot(x_a,psi(2,x_a),label='n=2')
plt.plot(x_a,psi(3,x_a),label='n=3')
plt.legend()
plt.xlabel(r'$\rm x $',fontsize=18)
plt.ylabel(r'$\rm \psi_{n}(x) $',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('1_a')
plt.show()

#2b)

#plot the harmonic oscillator wavefunction for a range of ns on -4<=x<=4
x_b=np.linspace(-10,10,500)
plt.figure(figsize=[8,6])
plt.plot(x_b,psi(30,x_b),color='black')
plt.xlabel(r'$\rm x $',fontsize=18)
plt.ylabel(r'$\rm \psi_{30}(x) $',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('1_b')
plt.show()

#2c)

#from lecture notes, define gaussxw (for -1->1)
def gaussxw(N):
    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))
    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))
    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
    return x, w

#define the integrand for 1c_x
def integrand_1c_x(n,x,dx):
    return x**2 * np.abs(psi(n,x))**2 * dx


# for both integrals i use change of variables 5.73 from the text

#compute <x^2>
#set integrand bounds
a=-1
b=1
#calculate the sample points and weights
N=100
x,w = gaussxw(N)
ns=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#compute the integral using change of variables from the text 5.73
sums_1c_x=[]
for n in ns:
    sum_1c_x=0
    for i in range(N):
        sum_1c_x += w[i] * integrand_1c_x(n, x[i]/(1 - x[i]**2), (1 + x[i]**2)/(1 - x[i]**2)**2)
    sums_1c_x.append(sum_1c_x)
#check the integration has been successful by printing n=5
print(' Sqrt of n=5 should equal approx. 2.35: n =',str(ns[5]),'=',str(np.sqrt(sums_1c_x[5])))

#compute <p^2>
#the derivative of the wave function uses a similar but different formula
def h_poly_d(n,x):
    if n < 0:
        return -1
    elif n == 0:
        return 1
    elif n == 1:
        return 2 * x
    else:
        return -x * h_poly(n,x) + 2 * n * h_poly(n-1,x)
    
#define a function to compute the harmonic oscillator wavefunction psi
def d_psi(n,x):
    return (1/np.sqrt(2**n * np.math.factorial(n) * np.sqrt(np.pi))) * np.exp((-x**2)/(2)) * h_poly_d(n,x)

#define the integrand for 1c_p
def integrand_1c_p(n,x,dx):
    return np.abs(d_psi(n,x))**2 * dx

#compute the integral using change of variables from text 5.73
sums_1c_p=[]
for n in ns:
    sum_1c_p=0
    for i in range(N):
        sum_1c_p += w[i] * integrand_1c_p(n, x[i]/(1 - x[i]**2), (1 + x[i]**2)/(1 - x[i]**2)**2)
    sums_1c_p.append(sum_1c_p)
    
#compute E = 1/2 (<x^2> + <p^2>)
E = 0.5 * (np.array(sums_1c_x) + np.array(sums_1c_p))

