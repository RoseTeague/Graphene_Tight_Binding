#!/usr/bin/env python3
"""
Module for 1D correlated Gaussian disorder potential on graphene
"""

import numpy as np

def oneDdisorderpotential(m,n,lc,pos):
    """
    ============================================================================
    Create one-dimensional disorder potential on graphene lattice
    ============================================================================
    
    This function creates an one-dimensional spatially-correlated Gaussian 
    disorder potential for graphene. The method used is the same as that in the 
    paper: Choi, SangKook, Cheol-Hwan Park, and Steven G. Louie. 
    "Electron supercollimation in graphene and Dirac Fermion materials 
    using one-dimensional disorder potentials." 
    Physical review letters 113.2 (2014): 026802. 
    
    
    To be clearer, the one-dimensional spatially-correlated Gaussian disorder 
    potential can be in the form of a random vector having the two-point 
    spatial correlation property. Hence, firstly, a random vector consisting of 
    spatially-uncorrelated Gaussian-random variables is composed. Next, using the 
    positions of atoms taken as input parameters, the two-point spatial correlation
    matrix is created and Cholesky decomposition method is used to obtain the 
    matrix with desired spatial correlation property. Finally, the final vector 
    is the dot product of the random vector and matrix with the required spatial 
    correlation property.
    
    
    Inputs
    ----------
    m : integer
        Number of atoms along the x-direction

    n : integer
        Number of atoms along the y-direction

    lc : float
        correlation length
    
    pos: float, tuple
        A tuple containing position information of atoms

    Returns
    -------
    Wfinal: float, array
         The final potential at each x position
    """
    
    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert type(lc) is float or int, "The correlation length must be numeric"
    assert type(pos) is tuple, "The pos must be a tuple"
    
    #Exctracting each unique x position.
    X = pos[0]
    x_1 = X[0:n*m:n,0]
    x_2 = X[1:n*m:n,0]

    #Parameters for the strength of the disorder potential in units of eV.
    Delta = 0.3

    #Generating sample of random numbers of Guassian distribution.
    V = np.random.normal(0,1,m)

    #Generate the two-point matrix for two rows along x axis
    X1 = np.tile(x_1,(m,1))
    X2 = X1.T
    X3 = np.tile(x_2,(m,1))
    X4 = X3.T

    #Generate the two-point spatial correlation matrix for two rows
    C1 = Delta**2*np.exp(-abs(X1-X2)/lc)
    C2 = Delta**2*np.exp(-abs(X3-X4)/lc)

    #Cholesky decomposition of the two-point correlation matrix and generate the final random vector
    L1 = np.linalg.cholesky(C1)
    W1 = np.dot(L1,V)
    L2 = np.linalg.cholesky(C2)
    W2 = np.dot(L2,V)

    #Reshaping for further calculations.
    Wf=np.zeros((n,m))
    Wf[0:n:2,0:m]=W1
    Wf[1:n:2,0:m]=W2
    Wfinal=Wf.T.reshape((n*m,1))
    
    return Wfinal

if __name__ == "__main__":
    from DD_WP_G import Crystal

    n = 10
    m = 10

    pos = Crystal(m,n)

    oneDdisorderpotential(m,n,10,pos)
