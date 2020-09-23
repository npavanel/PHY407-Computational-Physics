# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 20:04:20 2020

@author: Nicholas Pavanel
"""
#import needed modules
import numpy as np

#b)

#define alpha
alpha=20
#define the intensity transmission function of the diffraction grating at distance u from the central axis
def q(u):
    return (np.sin(alpha*u))**2

#c)

"""
Rephrasing the question: I am supposed to numerically integrate the integral $I(x) = \lvert \int^{\frac{\omega}{2}}_{-\frac{\omega}{2}} \sqrt{q(u)} e^{\frac{i 2 \pi x u}{\lambda f}} du \rvert^{2}$ representing the intensity of a diffraction pattern produced by a grating of 10 slits, from light with wavelength $\lambda = 500nm$, on a screen that has a total width $\omega = 10cm$, that is focused with a lens of focal length $f = 1m$, and has an intensity transmission function given by q(u) that was defined in part b).
​
When considering which integration method to use, I mainly thought about what a diffraction pattern looks like, and what our freedom of choice for number of steps meant. I think that I will use Simpson's integration as it is likely to give a higher degreen of accuracy than the trapezoid integration. Although the text states that trouble can arise when using Simpson's integration if the integrand is noisy or not smooth, I believe that diffraction patterns are relatively smooth functinos. Diffraction patterns smoothly increase and decrease in intensity across a given area, and thus Simpson's integration should work beautifully. 
​
I chose not to use Romberg integration as the text states that it is best applied to smooth functions whose form can be determined accurately from only a small number of sample points. I believe that a good amount of sample points are needed to accurately display a diffraction pattern, as there needs to be enough points not only to identify where each peak and trough is, but also the rapidity at which peaks transition into troughs. 
​
I chose not to use Gaussian quadrature because the text states that the integration points need to be unequally spaced. I believe that equally spaced integration points will be able to determine a diffraction pattern better than unequally spaced integration points because diffraction patterns demonstrate a good amount of symmetry.
​
As for the number of steps that will be used, I shall repeat the process of section 5.3 in the text until I get an accuracy of 3 decimal places.
​
From the text, equation 5.9, Simpson's Rule reads: $I(a,b) = \frac{1}{3} h  [f(a) + f(b) + 4 \Sigma_{k odd} f(a + kh) + 2 \Sigma_{k even} f(a + kh)]$ . I recognize that f in Simpson's rule will be $\sqrt{q(u)} e^{\frac{i 2 \pi x u}{\lambda f}}$ as written in the first paragraph of the text. a and b in Simpson's Rule will be the bounds of the integral, $\frac{\omega}{2}$ and $-\frac{\omega}{2}$. Finally, to get the correct answer, the absolute square of the answer of Simpson's rule will be taken as with the equation in paragraph 1.
"""

#define needed constants - all units in centimeters
wl=5*np.exp(-5) #wavelength
wi=10 #width
fl=100 #focal length
ns=10 #number of slits

#pseudocode
#define a function that computes (*) from above
#the function will take u and x as arguments, and return the value of (*)
#loop through for values of u across the width of the screen?
#pseudocode
#define a function that takes a position x, and integration bounds a and b as arguments to perform simposon's integration
#use the loops from the text to set up the sums: odd-for k in range(1,N,2) and even-for k in range(2,N,2)
#within the loops, reference (*) with (a+kh) as their x values
#multiple
