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

    This function creates an one-dimensional spatially-correlated Gaussian
    disorder potential for a square/rectangular lattice. The method used is
    the same as that in the paper: Choi, SangKook, Cheol-Hwan Park, and
    Steven G. Louie. "Electron supercollimation in graphene and Dirac Fermion
    materials using one-dimensional disorder potentials." Physical review
    letters 113.2 (2014): 026802.

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
    m : integer,
        Number of atoms along the x-direction

    n : integer,
        Number of atoms along the y-direction

    lc : float,
        correlation length

    pos: float, tuple,
        A tuple containing position information of atoms

    Returns
    -------
    Wfinal: float, array,
         The final potential at each x position

    """
    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert type(lc) is float or int, "The correlation length must be numeric"
    assert type(pos) is tuple, "The pos must be a tuple"

    #Exctracting each unique x position.
    X = pos[0]
    x = X[0:n*m:n,0]

    #disorder magnitude in the unit of eV
    Delta = 0.3

    #Generate a vector containing spatially-uncorrelated Gaussian-random variables
    V = np.random.normal(0,1,m)

    #Generate a two-point position matrix for further calculations
    X1 = np.tile(x,(m,1))
    X2 = X1.T

    #Generate a two-point spatial correlation matrix
    C = Delta**2*np.exp(-abs(X1-X2)/lc)

    #Cholesky decomposition of the correlation matrix
    L = np.linalg.cholesky(C)

    #Get the final random vector having desired two-point correlation matrix.
    W = np.dot(L,V)

    #Reshape the vector into a column for further calculations.
    Wf = np.zeros((n,m))
    Wf[:,0:m] = W
    Wfinal = Wf.T.reshape((n*m,1))

    return Wfinal
