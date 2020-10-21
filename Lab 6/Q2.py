# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 12:45:54 2020

@author: Nicholas Pavanel
"""
#import sympy to differentiate for me
import sympy as sp
import numpy as np

#define variables
x = sp.Symbol('x')
y = sp.Symbol('y')

#diff. wrt x
sp.diff(2 * ( (1)/(sp.sqrt(x**2 + y**2))**12 - (1)/(sp.sqrt(x**2 + y**2))**6 ),x)

#diff. wrt y
sp.diff(2 * ( (1)/(sp.sqrt(x**2 + y**2))**12 - (1)/(sp.sqrt(x**2 + y**2))**6 ),y)

#pseudocode
#define time step of dt=0.01
#define time domain of 100steps
#define initial velocitites to zero
#initialize first velocities, eqn 7 from lab handout
#define RHS equation for similar use as Q1
#define function for use in Verlet algorithm f(r,v0,t), eqn 6 from lab handout
#loop over eqns 8, 9, 10, 11 from lab handout
    #in the above loop compute the force due to the 1 other particle in domain
    #AND the force due to the 8N other fictious particles because of periodic boundary conditions
#at the end of above loop, ensure boundary conditions with np.mod(x,lx) from handout

#define timestep, time domain
dt = 0.01
N = 100
t1 = 0
t2 = 0 + 0.01*100
h = (t2-t1)/N

#define sets of initial conditions, [x,vx,y,vy]
ics_p1_1, ics_p2_1 = [4,0,4,0], [5.2,0,4,0]
ics_p1_2, ics_p2_2 = [4.5,0,4,0], [5.2,0,4,0]
ics_p1_3, ics_p2_3 = [2,0,3,0], [3.5,0,4.4,0]

def RHS_q2(x_y,r):
    return 12 * x_y * ( (1)/(r**2)**4 - (2)/(r**2)**7 )

#we basically have to do Q1 (before the RK4 loop) twice: need the positions of TWO particles

#define function to use in the Verlet algorithm
def f(r_p1, r_p2, t):
    x_p1, vx_p1, y_p1, vy_p1 = r_p1[0], r_p1[1], r_p1[2], r_p1[3] #set initial conditions
    x_p2, vx_p2, y_p2, vy_p2 = r_p2[0], r_p2[1], r_p2[2], r_p2[3]
    r = np.sqrt(x**2 + y**2)
    RHSx_p1, RHSy_p1 = RHS_q2()
    RHSx_p2, RHSy_p2
    
    