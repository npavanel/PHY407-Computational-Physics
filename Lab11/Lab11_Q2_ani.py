# This program calculates the total energy and magnetization
# for a 1D Ising model with N dipoles
# Author: Nico Grisouard, University of Toronto
# Date: 24 November 2020

# import modules
import numpy as np
from random import random, randrange
import matplotlib.pyplot as plt
import matplotlib.animation as ani

def energyfunction(J_, dipoles):
    """ function to calculate energy """
    sum1=np.sum(dipoles[:-1,:]*dipoles[1:,:])
    sum2=np.sum(dipoles[:,:-1]*dipoles[:,1:])
    energy = -J_*(sum1+sum2)
    return energy


def acceptance(Enew, Eold):
    """ Function for acceptance probability; to be completed """

    if dE>0:
        if random()<np.exp(-beta*dE):
            result = True #flip accepted
        else:
            result = False #flip denied
    else:
        result =  True #flip always accepted if new energy is less than old one
    return result  # result is True of False

# define constants
kB = 1.0
T = 3.0
J = 1.0
num_dipoles = 100 
N = 100 #length of chain
beta=1/(kB*T)
steps=100000#1000000

# generate array of dipoles and initialize diagnostic quantities
dipoles = (2*np.random.randint(0, 2, size=400))-1 #initialize array of dipoles
dipoles=dipoles.reshape(20,20)
energy = []  # empty list; to add to it, use energy.append(value)
magnet = []  # empty list; to add to it, use magnet.append(value)

E = energyfunction(J, dipoles)
energy.append(E)
magnet.append(np.sum(dipoles))
print(dipoles)
dipoles_array=np.zeros((steps,20,20))
dipoles_array[0]=dipoles

for i in range(steps):
    Eold=energyfunction(J, dipoles)
    picked1 = randrange(20)  # choose a victim
    picked2 = randrange(20) 
    dipoles[picked1, picked2] *= -1  # propose to flip the victim
    Enew = energyfunction(J, dipoles)  # compute Energy of proposed new state
    
    dE=Enew-Eold
    
    if dE>0: 
         if random()<np.exp(-beta*dE): #flip accepted
            energy.append(Enew)
            magnet.append(np.sum(dipoles))
            dipoles_array[i]=dipoles
         else:#flip denied
            energy.append(Eold)
            dipoles[picked1, picked2] *= -1  #flip it back
            magnet.append(np.sum(dipoles))
            dipoles_array[i]=dipoles
    else:#flip accepted
        energy.append(Enew)
        magnet.append(np.sum(dipoles))
        dipoles_array[i]=dipoles
        
#%%

dip_ani=np.zeros((int(steps/500), 20, 20)) #so that we plot every 5th 
for i in range(int(steps/500)):
    dip_ani[i]=dipoles_array[i*500]

fig = plt.figure( figsize=(8,8) )
im = plt.imshow(dipoles_array[0], interpolation='none', aspect='auto', vmin=-1, vmax=1) #how to label?
plt.colorbar()
def animate_func(i):
    if i % steps == 0:
        print( '.', end ='' )

    im.set_array(dip_ani[i])
    return [im]

anim = ani.FuncAnimation(
                               fig, 
                               animate_func, 
                               frames = int(steps/500),
                               interval = 1, 
                               )
anim.save('/Users/owner/407/lab11/animation.mp4', dpi=100, writer="ffmpeg")

