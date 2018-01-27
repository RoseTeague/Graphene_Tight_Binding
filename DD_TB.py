import math #can not remeber if this is how we do this ...
import numpy as np
import matplotlib as mpl
#import scipy.fftpack ... this might be a little awkward

#this is for the solver part ...

#in the function header we need to have the timestep and duration of simulation specified

def TB_solver(dt,DT):
    #add doc string here

    #need to add some assert statements ...

    #probably want to have a file with all of the details ...
    Ns = DT/dt #need to then round this off to the nearest integer ...

    HI = -2.7#Hopping integral in terms of eV

    #can also work out the matrix element that we will multiply by to solve the
    #propogation of the wave packet ... this is math.i*math.e/math.hbar ...

    #should call the wavefunction ...
    #need to specify how it is constructed so we can make the tridagonly matrices

    #construct the first tri diagonal matrix ...
    #two matrices are required here ...

    #construct the second tri diagonal matrix ...
    #two matrices are required here ...

    for i in range(Ns):

        #matrix multiply the wavefunctino by the + version of the first tri diagonal

        #solve the linear equation with the - version of the first tri diagonal and
        #the vector just calcualted

        #Now we need to shift the result around such that another tri diagonal
        #matrix equation can be solved

        #After we have reordered the result, we multiply with the + version of the
        #second tri diagonal

        #then solbe the linear equation with the - version of he second tri diagonal
        #and the vector just calculated

        #The result needs to be shifted around again so we can solve with the
        #first tri diagonal matrix again

        #repeat the first few steps ...

        #now we have our wavefunction at dt
