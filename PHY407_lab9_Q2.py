import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#define 1D sine transformation
def dst(y):    # Type-I discrete sine transform (DST) of real data y    
    N = len(y)    
    y2 = np.zeros(2*N)    
    y2[0] = y2[N] = 0.0    
    y2[1:N] = y[1:]    
    y2[:N:-1] = -y[1:]    
    a = -np.imag(np.fft.rfft(y2))[:N]    
    a[0] = 0.0    
    return a

#define 1D inverse-sine transformation
def idst(a):    # Type-I inverse DST of a    
    N = len(a)    
    c = np.zeros(N+1,dtype=complex)    
    c[0] = c[N] = 0.0    
    c[1:N] = -1j*a[1:]    
    y = np.fft.irfft(c)[:N]    
    y[0] = 0.0    
    return y

def dct(y):    
    N = len(y)    
    y2 = np.zeros(2*N,float)    
    y2[:N] = y[:]    
    y2[N:] = y[::-1]    
    c = np.fft.rfft(y2)    
    phi = np.exp(-1j*np.pi*np.arange(N)/(2*N))    
    return np.real(phi*c[:N])

#define 1D inverse-cosine transformation
def idct(a):    
    N = len(a)    
    c = np.zeros(N+1,complex)    
    phi = np.exp(1j*np.pi*np.arange(N)/(2*N))    
    c[:N] = phi*a    
    c[N] = 0.0    
    return np.fft.irfft(c)[:N]

# define 2d sine transform
def dst2(y):    
    M = y.shape[0]    
    N = y.shape[1]    
    a = np.zeros([M,N],float)    
    b = np.zeros([M,N],float)    
    for i in range(N):        
        a[:,i] = dst(y[:,i])    
    for j in range(M):        
        b[j,:] = dst(a[j,:])    
    return b

# define 2d sine transform
def idst2(b):    
    M = b.shape[0]    
    N = b.shape[1]    
    a = np.zeros([M,N],float)    
    y = np.zeros([M,N],float)    
    for i in range(M):        
        a[i,:] = idst(b[i,:])    
    for j in range(N):        
        y[:,j] = idst(a[:,j])    
    return y

#define a function that does hx
def hx_2d(y):
    M = y.shape[0]    
    N = y.shape[1]    
    a = np.zeros([M,N],float)    
    b = np.zeros([M,N],float)    
    for i in range(N):        
        a[:,i] = dct(y[:,i])    
    for j in range(M):        
        b[j,:] = dst(a[j,:])    
    return b

def ihx_2d(b):    
    M = b.shape[0]    
    N = b.shape[1]    
    a = np.zeros([M,N],float)    
    y = np.zeros([M,N],float)    
    for i in range(M):        
        a[i,:] = idst(b[i,:])    
    for j in range(N):        
        y[:,j] = idct(a[:,j])    
    return y

#define a function that does hy
def hy_2d(y):
    M = y.shape[0]    
    N = y.shape[1]    
    a = np.zeros([M,N],float)    
    b = np.zeros([M,N],float)    
    for i in range(N):        
        a[:,i] = dst(y[:,i])    
    for j in range(M):        
        b[j,:] = dct(a[j,:])    
    return b

def ihy_2d(b):    
    M = b.shape[0]    
    N = b.shape[1]    
    a = np.zeros([M,N],float)    
    y = np.zeros([M,N],float)    
    for i in range(M):        
        a[i,:] = idct(b[i,:])    
    for j in range(N):        
        y[:,j] = idst(a[:,j])    
    return y

# discretizing time
tau = 0.01 # time step
N = 2000 # total number of steps
T = N * tau # end time
n = np.arange(1,N,1) # individual steps
t = n * tau # time domain
# discretizing space
P = 32
ax = 1/32
ay = 1/32
x = np.arange(1,P) * ax
y = np.arange(1,P) * ay
X, Y = np.meshgrid(x,y)
p, q = X / ax, Y / ay 

# setting constants
Lx = Ly = j0 = m = n = c = 1
w = 3.75

Dx = np.pi * c * tau / (2 * Lx)
Dy = np.pi * c * tau / (2 * Ly)

#set output lists
Ez = []
Hx = []
Hy = []
Jz = []

#set initial conditions
ez = np.zeros((31,31))
hx = np.zeros((31,31))
hy = np.zeros((31,31))
jz = np.zeros((31,31))

for tn in t:
    jz = j0 * np.sin((m * np.pi * X)/Lx) * np.sin((n * np.pi * Y)/Ly) * np.sin(w * tn) # calc current pattern - eqn 8 

    # get the fourier coefficients - eqn 10
    Cjz = dst2(jz)
    Cez = dst2(ez)
    Chx = hx_2d(hx)
    Chy = hy_2d(hy)

    # evolve fourier coefficients - eqn 11
    Cez_next = ((1 - p**2 * Dx**2 - Q**2 * Dy**2) * Cez + (2 * Q * Dy * Chx) + (2 * p * Dx * Chy) + (tau * Cjz)) / (1 + p**2 * Dx**2 + Q**2 + Dy**2)
    Chx_next = Chx - Q * Dy * (Cez_next + Cez)
    Chy_next = Chy - p * Dx * (Cez_next + Cez)
    
    # inverse fourier back to regular
    ez = idst2(Cez_next)
    hx = ihx_2d(Chx_next)
    hy = ihy_2d(Chy_next) 
            
    Ez.append(ez[14][14])
    Hx.append(hx[30][14])
    Hy.append(hy[14][30])
    
plt.figure(figsize=[10,8])
plt.plot(t,Ez,label='Ez',color='blue')
plt.plot(t,Hx,label='Hx')
plt.plot(t,Hy,label='Hy')
plt.xlabel('t (arb)',size=18)
plt.ylabel('amp (arb)',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
plt.legend(fontsize=14)
plt.savefig('Q2_c')

# define a function to integrate the cavity for any value of w
def res_cavity2d_partd(ws):
    maxAs = []
    for w in ws:
        #set output lists
        Ez = []
        Hx = []
        Hy = []
        Jz = []

        #set initial conditions
        ez = np.zeros((31,31))
        hx = np.zeros((31,31))
        hy = np.zeros((31,31))
        jz = np.zeros((31,31))

        for tn in t:
            jz = j0 * np.sin((m * np.pi * X)/Lx) * np.sin((n * np.pi * Y)/Ly) * np.sin(w * tn) # calc current pattern - eqn 8 

            # get the fourier coefficients - eqn 10
            Cjz = dst2(jz)
            Cez = dst2(ez)
            Chx = hx_2d(hx)
            Chy = hy_2d(hy)

            # evolve fourier coefficients - eqn 11
            Cez_next = ((1 - p**2 * Dx**2 - q**2 * Dy**2) * Cez + (2 * q * Dy * Chx) + (2 * p * Dx * Chy) + (tau * Cjz)) / (1 + p**2 * Dx**2 + q**2 + Dy**2)
            Chx_next = Chx - q * Dy * (Cez_next + Cez)
            Chy_next = Chy - p * Dx * (Cez_next + Cez)

            # inverse fourier back to regular
            ez = idst2(Cez_next)
            hx = ihx_2d(Chx_next)
            hy = ihy_2d(Chy_next) 

            Ez.append(ez[15][15])
            
        maxAs.append(max(Ez))
        
    return maxAs

# define a uniformly separated list of ws
ws = np.linspace(0,9,21)

maxAs = res_cavity2d_partd(ws)

plt.figure(figsize=[10,8])
plt.plot(ws, maxAs)
plt.xlabel('$\omega$ (arb)',size=18)
plt.ylabel('max amp (arb)',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
plt.savefig('Q2_d')

norm_freq = np.pi * c * np.sqrt((n * Lx)**(-2) + (m * Ly)**(-2))

norm_freq

# define a function to integrate the cavity for any value of w
def res_cavity2d(w):
    #set output lists
    Ez = []
    Hx = []
    Hy = []
    Jz = []

    #set initial conditions
    ez = np.zeros((31,31))
    hx = np.zeros((31,31))
    hy = np.zeros((31,31))
    jz = np.zeros((31,31))
    
    for tn in t:
        jz = j0 * np.sin((m * np.pi * X)/Lx) * np.sin((n * np.pi * Y)/Ly) * np.sin(w * tn) # calc current pattern - eqn 8 

        # get the fourier coefficients - eqn 10
        Cjz = dst2(jz)
        Cez = dst2(ez)
        Chx = hx_2d(hx)
        Chy = hy_2d(hy)

        # evolve fourier coefficients - eqn 11
        Cez_next = ((1 - p**2 * Dx**2 - q**2 * Dy**2) * Cez + (2 * q * Dy * Chx) + (2 * p * Dx * Chy) + (tau * Cjz)) / (1 + p**2 * Dx**2 + q**2 + Dy**2)
        Chx_next = Chx - q * Dy * (Cez_next + Cez)
        Chy_next = Chy - p * Dx * (Cez_next + Cez)

        # inverse fourier back to regular
        ez = idst2(Cez_next)
        hx = ihx_2d(Chx_next)
        hy = ihy_2d(Chy_next) 

        Ez.append(ez[14][14])
        Hx.append(hx[30][14])
        Hy.append(hy[14][30])
        
    return Ez, Hx, Hy

Ez_freq, Hx_freq, Hy_freq = res_cavity2d(norm_freq)

plt.figure(figsize=[10,8])
plt.plot(t,Ez,label='Ez - w = 3.75')
plt.plot(t,Ez_freq,label='Ez - w = 4.44')
plt.xlabel('t (arb)',size=18)
plt.ylabel('amp (arb)',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
plt.legend()
plt.savefig('Q2_e1')
plt.show()
plt.close()

plt.figure(figsize=[8,6])
plt.plot(t,Hx,label='Hx - w = 3.75')
plt.plot(t,Hx_freq,label='Hx - w = 4.44')
plt.xlabel('t (arb)',size=18)
plt.ylabel('amp (arb)',size=18)
plt.legend()
plt.savefig('Q2_e2')
plt.show()
plt.close()

plt.figure(figsize=[8,6])
plt.plot(t,Hy,label='Hy - w = 3.75')
plt.plot(t,Hy_freq,label='Hy - w = 4.44')
plt.xlabel('t (arb)',size=18)
plt.ylabel('amp (arb)',size=18)
plt.legend()
plt.savefig('Q2_e3')
plt.show()
plt.close()