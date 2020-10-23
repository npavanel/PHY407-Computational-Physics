
import numpy as np
import matplotlib.pyplot as plt


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
t2 = 0 + dt*N

h = (t2-t1)/N
t = np.arange(t1,t2,h) 


#gives acceleration
def acc(x1,x2,y1,y2):
    return (-12 * (x1-x2) * ((1/(((x1-x2)**2+(y1-y2)**2)**4))-2/(((x1-x2)**2+(y1-y2)**2)**7))), (-12 * (y1-y2) * ((1/(((x1-x2)**2+(y1-y2)**2)**4))-2/(((x1-x2)**2+(y1-y2)**2)**7)))


#verlet

def verly(x10,y10, x20, y20):
    r1=np.zeros([2,100]) 
    r1[0,0], r1[1,0]=(x10,y10) #sets initial x1 and y1
    k1=np.zeros([2,100])
    v1=np.zeros([2,100]) #x1,y1
    v1_half=np.zeros([2,100])
    r2=np.zeros([2,100]) 
    r2[0,0], r2[1,0]=(x20,y20) #sets initial x2 and y2
    k2=np.zeros([2,100])
    v2=np.zeros([2,100])
    v2_half=np.zeros([2,100])
    v1_half[0,0], v1_half[1,0]=np.array([0,0])+1/2*h*np.array(acc(x10, x20, y10, y20)) #sets inital vx0,vy0
    v2_half[0,0], v2_half[1,0]=np.array([0,0])+1/2*h*np.array((acc(x20, x10, y20, y10))) #sets inital vx0,vy0
    for i in range(len(t)-1):
        r1[0,i+1], r1[1,i+1]=r1[:,i]+h*np.array(v1_half[:,i]) #x[i+1],y[i+1] #cols, rows , so : means both columns, ith row
        r2[0,i+1], r2[1,i+1]=r2[:,i]+h*np.array(v2_half[:,i]) #x[i+1],y[i+1] #cols, rows , so : means both columns, ith row
        k1[:,i+1]=h*np.array(acc(r1[0,i+1], r2[0,i+1], r1[1,i+1], r2[1,i+1]))
        k2[:,i+1]=h*np.array(acc(r2[0,i+1], r1[0,i+1], r2[1,i+1], r1[1,i+1]))
        v1[:,i+1]=v1_half[:,i]+0.5*k1[:,i+1]
        v1_half[:,i+1]=v1_half[:,i]+k1[:,i+1]
        v2[:,i+1]=v2_half[:,i]+0.5*k2[:,i+1]
        v2_half[:,i+1]=v2_half[:,i]+k2[:,i+1]
    return (r1,r2)





plt.figure(1)
plt.scatter(verly(4,4,5.2,4)[0][0], verly(4,4,5.2,4)[0][1], label="p1", marker=".", s=2)
plt.scatter(verly(4,4,5.2,4)[1][0], verly(4,4,5.2,4)[1][1] , label="p2", marker=".", s=2)
plt.ylim(3.5,4.5)
plt.xlabel("X position")
plt.ylabel("Y position")

plt.figure(2)
plt.scatter(t,verly(4,4,5.2,4)[0][0], label="p1", marker=".", s=2)
plt.scatter(t,verly(4,4,5.2,4)[1][0], label="p2", marker=".", s=2)
plt.ylabel("X position")
plt.xlabel("time (s)")

plt.figure(3)
plt.scatter(verly(4.5,4,5.2,4)[0][0], verly(4.5,4,5.2,4)[0][1], label="p1", marker=".", s=2)
plt.scatter(verly(4.5,4,5.2,4)[1][0], verly(4.5,4,5.2,4)[1][1] , label="p2", marker=".", s=2)
plt.xlabel("X position")
plt.ylim(3.5,4.5)
plt.ylabel("Y position")

plt.figure(4)
plt.scatter(t,verly(4.5,4,5.2,4)[0][0], label="p1", marker=".", s=2)
plt.scatter(t,verly(4.5,4,5.2,4)[1][0], label="p2", marker=".", s=2)
plt.ylabel("X position")
plt.xlabel("time (s)")

plt.figure(5)
plt.scatter(verly(2,3,3.5,4.4)[0][0], verly(2,3,3.5,4.4)[0][1], label="p1", marker=".", s=2)
plt.scatter(verly(2,3,3.5,4.4)[1][0], verly(2,3,3.5,4.4)[1][1] , label="p2", marker=".", s=2)
plt.xlabel("X position")
plt.ylabel("Y position")

plt.figure(6)
plt.scatter(t,verly(2,3,3.5,4.4)[0][0], label="p1", marker=".", s=2)
plt.scatter(t,verly(2,3,3.5,4.4)[1][0], label="p2", marker=".", s=2)
plt.ylabel("X position")
plt.xlabel("time (s)")




