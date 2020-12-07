# imports 
from math import sqrt,exp
from numpy import empty
from random import random,randrange
from random import seed
import matplotlib.pyplot as plt
import numpy as np
from time import time

# initial conditions
N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e4

# Function to calculate the magnitude of a vector
def mag(x):
    return sqrt(x[0]**2+x[1]**2)

# Function to calculate the total length of the tour
def distance():
    s = 0.0
    for i in range(N):
        s += mag(r[i+1]-r[i])
    return s

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    seed(i)
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

# Main loop
t = 0
T = Tmax
#seed(3)
while T>Tmin:

    # Cooling
    t += 1
    T = Tmax*exp(-t/tau)

    # Choose two cities to swap and make sure they are distinct
    i,j = randrange(1,N),randrange(1,N)
    while i==j:
        i,j = randrange(1,N),randrange(1,N)

    # Swap them and calculate the change in distance
    oldD = D
    r[i,0],r[j,0] = r[j,0],r[i,0]
    r[i,1],r[j,1] = r[j,1],r[i,1]
    D = distance()
    deltaD = D - oldD

    # If the move is rejected, swap them back again
    if random()>exp(-deltaD/T):
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = oldD
        
# parse through the route info to create a list of x and y coordinates
xs = []
ys = []
for i in range(len(r)):
    xs.append(r[i][0])
    ys.append(r[i][1])
    
# plot the route over the cities
plt.figure(figsize=[8,6])
plt.plot(xs,ys,color='black',label='Route')
plt.scatter(xs,ys,color='blue',label='Cities')
plt.legend()
plt.show()

# varying the route seed

t0 = time()
# initalize lists to hold all trajectoies and distances
Xs_1 = []
Ys_1 = []
Ds_1 = []

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    seed(i+100)
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

for k in range(4):
    # Main loop
    t = 0
    T = Tmax
    seed(k)
    while T>Tmin:

        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Choose two cities to swap and make sure they are distinct
        i,j = randrange(1,N),randrange(1,N)
        while i==j:
            i,j = randrange(1,N),randrange(1,N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = distance()
        deltaD = D - oldD

        # If the move is rejected, swap them back again
        if random()>exp(-deltaD/T):
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,1],r[j,1] = r[j,1],r[i,1]
            D = oldD

    #parse through the route info to create a list of x and y coordinates
    xs = []
    ys = []
    for l in range(len(r)):
        xs.append(r[l][0])
        ys.append(r[l][1])

    Xs_1.append(xs)
    Ys_1.append(ys) 
    Ds_1.append(D)
    
print('Time taken: ', time() - t0)

fig, axs = plt.subplots(2,2,figsize=(10,10),sharex=True,sharey=True)
axs[0,0].scatter(Xs_1[0],Ys_1[0],label='Cities',color='Black')
axs[0,0].plot(Xs_1[0],Ys_1[0],label='Route 1',color='Blue')
axs[0,1].scatter(Xs_1[0],Ys_1[0],label='Cities',color='Black')
axs[0,1].plot(Xs_1[1],Ys_1[1],label='Route 2',color='Green')
axs[1,0].scatter(Xs_1[0],Ys_1[0],label='Cities',color='Black')
axs[1,0].plot(Xs_1[2],Ys_1[2],label='Route 3',color='Orange')
axs[1,1].scatter(Xs_1[0],Ys_1[0],label='Cities',color='Black')
axs[1,1].plot(Xs_1[3],Ys_1[3],label='Route 4',color='Red')

axs[0,0].legend(loc=(0.4,0.3),fontsize=14)
axs[0,1].legend(loc=(0.4,0.3),fontsize=14)
axs[1,0].legend(loc=(0.4,0.3),fontsize=14)
axs[1,1].legend(loc=(0.4,0.3),fontsize=14)

axs[0,0].set_yticklabels([0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],fontsize=16)
axs[1,0].set_yticklabels([0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],fontsize=16)

axs[1,1].set_xticklabels([0,0,0.2,0.4,0.6,0.8,1.0],fontsize=16)
axs[1,0].set_xticklabels([0,0,0.2,0.4,0.6,0.8,1.0],fontsize=16)

plt.tight_layout()

#plt.savefig('Lab11_Q1a')

# varying the time constant

tau = 1e3

t0 = time()

# initalize lists to hold all trajectoies and distances
Xs_2 = []
Ys_2 = []
Ds_2 = []

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    seed(i+100)
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

for k in range(4):
    # Main loop
    t = 0
    T = Tmax
    seed(k)
    while T>Tmin:

        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Choose two cities to swap and make sure they are distinct
        i,j = randrange(1,N),randrange(1,N)
        while i==j:
            i,j = randrange(1,N),randrange(1,N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = distance()
        deltaD = D - oldD

        # If the move is rejected, swap them back again
        if random()>exp(-deltaD/T):
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,1],r[j,1] = r[j,1],r[i,1]
            D = oldD

    #parse through the route info to create a list of x and y coordinates
    xs = []
    ys = []
    for l in range(len(r)):
        xs.append(r[l][0])
        ys.append(r[l][1])

    Xs_2.append(xs)
    Ys_2.append(ys) 
    Ds_2.append(D)
    
print('Time taken: ', time() - t0)

tau = 5e3

t0 = time()

# initalize lists to hold all trajectoies and distances
Xs_4 = []
Ys_4 = []
Ds_4 = []

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    seed(i+100)
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

for k in range(4):
    # Main loop
    t = 0
    T = Tmax
    seed(k)
    while T>Tmin:

        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Choose two cities to swap and make sure they are distinct
        i,j = randrange(1,N),randrange(1,N)
        while i==j:
            i,j = randrange(1,N),randrange(1,N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = distance()
        deltaD = D - oldD

        # If the move is rejected, swap them back again
        if random()>exp(-deltaD/T):
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,1],r[j,1] = r[j,1],r[i,1]
            D = oldD

    #parse through the route info to create a list of x and y coordinates
    xs = []
    ys = []
    for l in range(len(r)):
        xs.append(r[l][0])
        ys.append(r[l][1])

    Xs_4.append(xs)
    Ys_4.append(ys) 
    Ds_4.append(D)
    
print('Time taken: ', time() - t0)

tau = 1e5

t0 = time()

# initalize lists to hold all trajectoies and distances
Xs_3 = []
Ys_3 = []
Ds_3 = []

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    seed(i+100)
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

for k in range(4):
    # Main loop
    t = 0
    T = Tmax
    seed(k)
    while T>Tmin:

        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Choose two cities to swap and make sure they are distinct
        i,j = randrange(1,N),randrange(1,N)
        while i==j:
            i,j = randrange(1,N),randrange(1,N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = distance()
        deltaD = D - oldD

        # If the move is rejected, swap them back again
        if random()>exp(-deltaD/T):
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,1],r[j,1] = r[j,1],r[i,1]
            D = oldD

    #parse through the route info to create a list of x and y coordinates
    xs = []
    ys = []
    for l in range(len(r)):
        xs.append(r[l][0])
        ys.append(r[l][1])

    Xs_3.append(xs)
    Ys_3.append(ys) 
    Ds_3.append(D)
    
print('Time taken: ', time() - t0)

tau = 5e4

t0 = time()

# initalize lists to hold all trajectoies and distances
Xs_5 = []
Ys_5 = []
Ds_5 = []

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    seed(i+100)
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

for k in range(4):
    # Main loop
    t = 0
    T = Tmax
    seed(k)
    while T>Tmin:

        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Choose two cities to swap and make sure they are distinct
        i,j = randrange(1,N),randrange(1,N)
        while i==j:
            i,j = randrange(1,N),randrange(1,N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = distance()
        deltaD = D - oldD

        # If the move is rejected, swap them back again
        if random()>exp(-deltaD/T):
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,1],r[j,1] = r[j,1],r[i,1]
            D = oldD

    #parse through the route info to create a list of x and y coordinates
    xs = []
    ys = []
    for l in range(len(r)):
        xs.append(r[l][0])
        ys.append(r[l][1])

    Xs_5.append(xs)
    Ys_5.append(ys) 
    Ds_5.append(D)
    
print('Time taken: ', time() - t0)

plt.figure(figsize=[8,6])
plt.plot(Ds_2,label='$\\tau = 10^{3}$')
plt.plot(Ds_4,label='$\\tau = 5 \\times 10^{3}$')
plt.plot(Ds_1,label='$\\tau = 10^{4}$')
plt.plot(Ds_5,label='$\\tau = 5 \\times 10^{4}$')
plt.plot(Ds_3,label='$\\tau = 10^{5}$')
plt.yticks(size='16')
plt.xticks([0,1,2,3],[1,2,3,4],size='16')
plt.xlabel('Route #',size='18')
plt.ylabel('Route Distance',size='18')
plt.legend(fontsize=14)
#plt.savefig('Lab11_Q1a_2')
plt.show()

# b)

# define a function to compute eqn (9) from handout, function to be analyzed
def eqn_9(x,y):
    return x**2 - np.cos(4 * np.pi * x) + (y - 1)**2

# inital conditions
Tmax = 1
Tmin = 1e-3
tau = 1e4
x0 = 2
y0 = 2

# define lists to hold time evolution of choices num time steps
Xs = []
Ys = []
times = []

# set initial conditions 
T = Tmax
t = 0
x = x0
y = y0

# main loop
while T > Tmin:
    
    # cooling 
    t += 1
    T = Tmax*exp(-t / tau)
    
    # store old points and function value
    oldx = x
    oldy = y
    oldfxy = eqn_9(x,y)
    
    # generate new steps from a random dist
    dx = np.random.standard_normal()
    dy = np.random.standard_normal()
    
    # determine new points and func value
    x += dx
    y += dy
    fxy = eqn_9(x,y)
    
    # compute difference between old and new func values
    delta_fxy = fxy - oldfxy
    
    # reject values under normal conditions or if they are out of domain
    if random() > exp(-delta_fxy / T):
        x = oldx
        y = oldy
        fxy = oldfxy
        
    Xs.append(x)
    Ys.append(y)
    times.append(t)
        
print('(x,y): ', (x,y))
print('f(x,y): ', eqn_9(x,y))

plt.figure(figsize=[8,6])
plt.scatter(times,Xs,color='black')
plt.xlabel('Timestep (arb)',size=16)
plt.ylabel('x (arb)',size=16)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.savefig('Lab11_Q1b_i_1')
plt.show()
plt.close()

plt.figure(figsize=[8,6])
plt.scatter(times,Ys,color='black')
plt.xlabel('Timestep (arb)',size=16)
plt.ylabel('y (arb)',size=16)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.savefig('Lab11_Q1b_i_2')
plt.show()

# now define a funciton to compute eqn 11
def eqn_11(x,y):
    return np.cos(x) + np.cos(np.sqrt(2) * x) + np.cos(np.sqrt(3) * x) + (y-1)**2

# inital conditions
Tmax = 1
Tmin = 1e-3
tau = 1e6
x0 = 2
y0 = 2

# define lists to hold time evolution of choices num time steps
Xs = []
Ys = []
times = []

# set initial conditions 
T = Tmax
t = 0
x = x0
y = y0

# main loop
while T > Tmin:
    
    # cooling 
    t += 1
    T = Tmax * exp(-t / tau)
    
    # store old points and function value
    oldx = x
    oldy = y
    oldfxy = eqn_11(x,y)
    
    # generate new steps from a random dist
    dx = np.random.standard_normal()
    dy = np.random.standard_normal()
    
    # determine new points and func value
    x += dx
    y += dy
    fxy = eqn_11(x,y)
    
    # compute difference between old and new func values
    delta_fxy = fxy - oldfxy
    
    # reject values under normal conditions or if they are out of domain 
    if random() > exp(-delta_fxy / T) or x > 50 or x < 0 or y > 20 or y < -20:
        x = oldx
        y = oldy
        fxy = oldfxy
        
    Xs.append(x)
    Ys.append(y)
    times.append(t)
        
print('(x,y): ', (x,y))
print('f(x,y): ', eqn_11(x,y))

plt.figure(figsize=[8,6])
plt.scatter(times,Xs,color='black')
plt.xlabel('Timestep (arb)',size=16)
plt.ylabel('x (arb)',size=16)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.savefig('Lab11_Q1b_ii_1')
plt.show()
plt.close()

plt.figure(figsize=[8,6])
plt.scatter(times,Ys,color='black')
plt.xlabel('Timestep (arb)',size=16)
plt.ylabel('y (arb)',size=16)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.savefig('Lab11_Q1b_ii_2')
plt.show()

