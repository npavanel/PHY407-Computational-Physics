# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 09:18:21 2020

@author: Nicholas Pavanel
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#pre-a) test an example from the book for understanding

#from page 251 of the text 
def test(x):
    return 2 - np.exp(-x)

#define relaxation as per pg 251 of the text
def relaxation_test(x0,steps):
    for k in range(steps):
        x0 = test(x0)
    print(x0)
    
#test with in text example
relaxation_test(1,50)

#a) - 6.10 a)

#define function for 6.10
def f(c,x):
    return 1 - np.exp(-c * x)

#modify relaxation_test to relax with f(c,x)
def relaxation_1(c,x0,steps):
    for k in range(steps):
        x0 = f(c,x0)
        print(np.round(x0,6))
        
#run relaxation with printed values
relaxation_1(2,2,25)

#a) - 6.10 b)

#modify relaxation_1 to relax over a list of c values
def relaxation_2(c,x0,steps):
    solns = []
    for i in c:
        x = x0
        for k in range(steps):
            x = f(i,x)
        solns.append(x)
    return solns

#define the domain of c values
cs = np.arange(0,3.01,0.01)
#check solutions
solns = relaxation_2(cs,1,10)

#plot the relaxed solutions vs the cs
plt.figure(figsize=[8,6])
plt.plot(cs,solns)
plt.xticks(size=16)
plt.yticks(size=16)
plt.xlabel('C Value',size=18)
plt.ylabel('Relaxation Solution',size=18)
#plt.savefig('3')

#b) - 6.11 b)

#modify relaxation_1 so that it checks how many steps to get to the soln
def relaxation_3(c,x0,steps):
    num_steps = 0
    for k in range(steps):
        if num_steps != 0 and np.round(x0,6) == np.round(f(c,x0),6):
            print('It took',num_steps,'steps to get a soln accurate to 10^(-6)')
            break
        else:
            num_steps += 1
            x0 = f(c,x0)
            
#run relaxation 3
relaxation_3(2,2,25)

#b) - 6.11 c)

#modify relaxation_3 for to employ overrelaxation
def over_relaxation(c, x0, w, steps):
    num_steps = 0
    for k in range(steps):
        if num_steps != 0 and np.round(x0,6) == np.round(f(c,x0),6):
            print('It took',num_steps,'steps to get a soln accurate to 10^(-6)')
            break
        else:
            num_steps += 1
            x0 = x0 + (1 + w) * (f(c,x0) - x0)
            print(np.round(x0,6))
            
#run over_relaxation
over_relaxation(2,2,0.5,25)

#b) - 6.11 d)

#find answer in lab hand in sheet

#c) - 6.13 b)

#Relaxation

#define a function to be used for relaxation
def g(x):
    return 5 - 5 * np.exp(-x)

#modify relaxation_3 so that it uses g(x)
def relaxation_4(x0,steps):
    num_steps = 0
    for k in range(steps):
        if num_steps != 0 and np.round(x0,6) == np.round(g(x0),6):
            print('It took relaxation',num_steps,'steps to get a soln accurate to 10^(-6):',np.round(x0,6))
            break
        else:
            num_steps += 1
            x0 = g(x0)
            
#modify over_relaxation so that it uses g(x)
def over_relaxation_2(x0, w, steps):
    num_steps = 0
    for k in range(steps):
        if num_steps != 0 and np.round(x0,6) == np.round(g(x0),6):
            print('It took over-relaxation',num_steps,'steps to get a soln accurate to 10^(-6):',np.round(x0,6))
            break
        else:
            num_steps += 1
            x0 = x0 + (1 + w) * (g(x0) - x0)
            #print(np.round(x0,6))
            
#compute relaxation and over-relaxation methods 
relaxation_4(2,50)
over_relaxation_2(2,0.035,50)

#Binary

#define a function to be used for binary
def h(x):
    return 5 * np.exp(-x) + x - 5

#define binary search method
def binary(x1,x2,accuracy):
    num_steps = 0 #define varible to count 
    if h(x1) < 0 and h(x2) > 0:
        while abs(x1 - x2) > accuracy: 
            x_mid = 0.5 * (x1 + x2)
            f_at_mid = h(x_mid)
            if f_at_mid < 0:
                x1 = x_mid
                num_steps += 1
            else:
                x2 = x_mid
                num_steps += 1
        print('It took binary',num_steps,'steps to get a soln accurate to 10^(-6):',np.round(0.5 * (x1 + x2),6)) 
    elif h(x1) > 0 and h(x2) < 0:
        while abs(x1 - x2) > accuracy: 
            x_mid = 0.5 * (x1 + x2)
            f_at_mid = h(x_mid)
            if f_at_mid > 0:
                x1 = x_mid
                num_steps += 1
            else:
                x2 = x_mid
                num_steps += 1
        print('It took binary',num_steps,'steps to get a soln accurate to 10^(-6):',np.round(0.5 * (x1 + x2),6)) 
    else:
        print('Initial points do not enclose a zero.')
        
#run binary method
binary(2, 6, 10**(-6))

#Newton's Method

#note that newton uses the same root method as binary, so h(x) will work here. 
#define h prime
def h_prime(x):
    return - 5 * np.exp(-x) + 1

#define a function with both h_prime and h
def newton_h(x):
    return x - h(x) / h_prime(x)

#define a function for the newton method
def newton(x0,accuracy):
    num_steps = 0
    while abs(x0 - newton_h(x0)) > accuracy:
        x0 = newton_h(x0)
        num_steps += 1
    print('It took Newton',num_steps,'steps to get a soln accurate to 10^(-6):',np.round(newton_h(x0),6))
    
#run newton's method
newton(2,10**(-6))

#c) - 6.13 c)

#set initial conditions
x = 4.965114
lam = 5.02 * 10**(-7)
kb = 1.380649 * 10**(-23)
c = 2.9979246 * 10**8

#calculate temperature
T = (h * c) / (kb * x * lam)

#print temp
T

        
