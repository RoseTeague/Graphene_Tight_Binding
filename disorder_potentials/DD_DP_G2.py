#!/usr/bin/env python3
"""
Module for 2D correlated Gaussian disorder potential on graphene
"""

import numpy as np

def twoDdisorderpotential(m,n,lc,pos):
    """
    ============================================================================
    Create two-dimensional disorder potential on graphene lattice
    ============================================================================
    
    This function creates an two-dimensional spatially-correlated Gaussian 
    disorder potential for graphene. The method used is the same as that in the 
    paper: Choi, SangKook, Cheol-Hwan Park, and Steven G. Louie. 
    "Electron supercollimation in graphene and Dirac Fermion materials 
    using one-dimensional disorder potentials." 
    Physical review letters 113.2 (2014): 026802. 
    
    
    To be clearer, the two-dimensional spatially-correlated Gaussian disorder 
    potential can be the sum of two one-dimensional potentials respectively along x
    and y axis. The one-dimensional potential can be in the form of a random vector 
    having the two-point spatial correlation property. Hence, firstly, a random vector 
    consisting of spatially-uncorrelated Gaussian-random variables is composed. 
    Next, using the positions of atoms taken as input parameters, the two-point
    spatial correlation matrix is created and Cholesky decomposition method is 
    used to obtain the matrix with desired spatial correlation property. Finally, 
    the final vector is the dot product of the random vector and matrix with the 
    required spatial correlation property.
    
    
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
         The final potential at each atom
    """
    
    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert type(lc) is float or int, "The correlation length must be numeric"
    assert type(pos) is tuple, "The pos must be a tuple"


    #Exctracting each unique x position.
    X = pos[0]
    x_1 = X[0:n*m:n,0]
    x_2 = X[1:n*m:n,0]

    #Exctracting each unique y position.
    Y = pos[1]
    y = Y[0:n,0]

    #Parameters for the strength of the disorder potential in units of eV.
    Delta = 0.3

    #Generating sample of random numbers of Guassian distribution along x axis
    Vx = np.random.normal(0,1,m)

    #Generating sample of random numbers of Guassian distribution along y axis
    Vy = np.random.normal(0,1,n)

    #Generate the two-point matrix for two rows along x axis
    X1 = np.tile(x_1,(m,1))
    X2 = X1.T
    X3 = np.tile(x_2,(m,1))
    X4 = X3.T

    #Generate the two-point matrix for y axis
    X5 = np.tile(y,(n,1))
    X6 = X5.T

    #Generate the two-point spatial correlation matrix for x (two rows) and y axis
    Cx1 = Delta**2*np.exp(-abs(X1-X2)/lc)
    Cx2 = Delta**2*np.exp(-abs(X3-X4)/lc)
    Cy = Delta**2*np.exp(-abs(X5-X6)/lc)

    #Cholesky decomposition of the two-point correlation matrix and generate the final random vector
    Lx1 = np.linalg.cholesky(Cx1)
    Wx1 = np.dot(Lx1,Vx)
    Lx2 = np.linalg.cholesky(Cx2)
    Wx2 = np.dot(Lx2,Vx)
    Ly = np.linalg.cholesky(Cy)
    Wy = np.dot(Ly,Vy)

    #Put the potentials on all atoms and reshape them to the requried output format
    Wxf=np.zeros((n,m))
    Wyf=np.zeros((n,m))

    Wxf[0:n:2,0:m]= Wx1
    Wxf[1:n:2,0:m]= Wx2
    Wxfinal=Wxf.T.reshape((n*m,1))

    Wyf[:,0:m]=Wy.reshape((n,1))
    Wyfinal=Wyf.T.reshape((n*m,1))

    #The final potential for each atom is the sum of the x- and y- potential component
    Wfinal=Wxfinal+Wyfinal

    return Wfinal

