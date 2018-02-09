#!/usr/bin/env python3
"""
Module for correlated Gaussian disorder potential on square lattice
"""

import numpy as np

def oneDdisorderpotential(m,n,lc,pos):
    """
    ============================================================================
    Create one-dimensional disorder potentials on square/rectangular lattice
    ============================================================================

    Inputs
    ----------
    m : integer
        Number of atoms along the x-direction

    n : integer
        Number of atoms along the y-direction

    lc : float
        correlation length
        
    pos: float, list
        A list containing position information of atoms

    Returns
    -------
    Wfinal: float, array
         The final potential at each x position

    """
    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert type(lc) is float or int, "The correlation length must be numeric"
    assert type(pos) is list, "The pos must be a list"
    
    #Exctracting each unique x position.
    X = pos[0]
    x = X[0:n*m:n,0]

    #Create One-dimensional disorder potential

    #disorder magnitude in the unit of eV
    Delta = 0.3

    #Generate a vector containing spatially-uncorrelated Gaussian-random variables
    V=np.random.normal(0,1,m)

    #Generate a two-point position matrix for further calculations
    X1=np.tile(x,(m,1))
    X2=X1.T

    #Generate a two-point spatial correlation matrix
    C=Delta**2*np.exp(-abs(X1-X2)/lc)

    #Cholesky decomposition of the correlation matrix
    L=np.linalg.cholesky(C)

    #Get the final random vector having desired two-point correlation matrix.
    W=np.dot(L,V)

    #Reshape the vector into a column for further calculations.
    Wf=np.zeros((n,m))
    Wf[:,0:m] = W
    Wfinal = Wf.T.reshape((n*m,1))

    return Wfinal

if __name__ == "__main__":

    from DD_WP_S import Crystal

    n = 10
    m = 10

    pos = Crystal(m,n)

    oneDdisorderpotential(m,n,10,pos)
