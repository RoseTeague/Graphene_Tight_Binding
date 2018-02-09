#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for 1D correlated Gaussian disorder potential on graphene
"""

import numpy as np

def oneDdisorderpotential(m,n,lc,pos):
    """
    ============================================================================
    Create one-dimensional disorder potential on graphene lattice
    ============================================================================

    Inputs
    ----------
    m : integer
        Number of atoms along the x-direction

    n : integer
        Number of atoms along the y-direction

    lc : float
        correlation length

    Returns
    -------
    Wfinal: float, array
         The final potential at each x position
    """

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
