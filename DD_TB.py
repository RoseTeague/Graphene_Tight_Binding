"""Module for solver of with TB split operator
"""
import math #can not remeber if this is how we do this ...
import numpy as np
import matplotlib as mpl
from scipy import constants
import scipy.linalg
#import scipy.fftpack ... this might be a little awkward

#this is for the solver part ...

#in the function header we need to have the timestep and duration of simulation specified

def TB_solver(dt,DT):
    #add doc string here

    #need to add some assert statements ...

    #probably want to have a file with all of the details ...
    Ns = DT/dt #need to then round this off to the nearest integer ...
    N = 100
    n = 10

    #should call the wavefunction ...
    #need to specify how it is constructed so we can make the tridagonly matrices
    wvf = np.ones((N,1))#wavefunction ... just put some matrix there for now ...

    #then need to call the hamiltonian ...

    #for i in range(Ns):

        #matrix multiply the wavefunctino by the + version of the first tri diagonal

        #https://stackoverflow.com/questions/44388358/python-numpy-matrix-multiplication-with-one-diagonal-matrix
        #c = (a*b.T).T ... fastest way and we have it in this form ...

        #solve the linear equation with the - version of the first tri diagonal and
        #the vector just calcualted

        # scipy.linalg.solve_banded((1,1), m_packed, b) ... not sure what the first tupel is

        #Now we need to shift the result around such that another tri diagonal
        #matrix equation can be solved
        #might be able to do this with the range thing from above ...
        #should be able to do this with something similar to DD_TBH

        #After we have reordered the result, we multiply with the + version of the
        #second tri diagonal

        #then solbe the linear equation with the - version of he second tri diagonal
        #and the vector just calculated

        #The result needs to be shifted around again so we can solve with the
        #first tri diagonal matrix again

        #repeat the first few steps ...

        #now we have our wavefunction at dt
