import numpy as np
import matplotlib as mpl
import math #can not remeber if this is how we do this ...
#import scipy.fftpack ... this might be a little awkward

#this is for the solver part ...

#in the function header we need to have the timestep and duration of simulation specified

def cont_solver(dt,DT):
    #add doc string here

    #need to add some assert statements ...

    Ns = DT/dt #need to then round this off to the nearest integer ...

    vF = 1#then need to specify the value of the fermi velocity in terms of the tight binding model

    #k_x and k_y need to be extracted from the gaussian ... maybe have all of this information in a file to start with!

    XIT = dt*vF*(k_x**2 + k_y**2)**0.5

    CXIT = math.cos(XIT)
    SXIT = math.sin(XIT)/XIT
    XIP = k_x + k_y# do these need to be the modulus?  also need the imaginary part
    XIN = k_x - k_y# do these need to be the modulus?  also need the imaginary part

    #might want to define some arrays for other outputs ...

    #need to determine if there is a potential or not ... have some kind of if statment that checks if all of the elements
    #of the array from the potential function is zero or not

    #then need to have a for loop ...

    #for when there is no potential
    #we only need to perfrom one matrix multiplication
    #this needs to be done in reciprocal space though ...
    #since we do not need to flick between the two, we can do the fourier transform before and after the loop ... ?

    #column vector form of the spinor ...
    phi_A = c1*wavefunction
    phi_B = c2*wavefunction

    #fourier transfrom the wavefunction ...
    FT_A = np.fft.fft2(phi_A)
    FT_B = np.fft.fft2(phi_B)#might be able to get rid of this one ... it is only important for magnetic fields ... ?

    for i in range(Ns):

        #perform matrix multiplication in reciprocal space ... use some simplified way of splitting it up

        #need to get i ...
        FT_A = CXIT*FT_A - SXIT*XIN*FT_B
        FT_B = CXIT*FT_B - SXIT*XIP*FT_A

        #when we are plotting we will need to think of how to get specific frames out

    #when the loop has finished we need to take the inverse transform

    PHI_A = np.fft.ifft2(FT_A)
    PHI_B = np.fft.ifft2(FT_B)

    #then we need to determine the probability density
    #this is a simple matrix multiplication ...



    #for when there is a potential

    #there are different parts of this in reciprocal space and others in real monospace
    #need to perform the fourier transform before and after multiplication ...

    for i in range(Ns):

        #add simple matrix multiplication for the potential ...

        #fourier transform the function to reciprocal space
        #perform matrix multiplication in reciprocal space
        #inverse fourier transform to real space again

        #simple matric multiplication in real monospace


    #then need to determine the probability density ...
