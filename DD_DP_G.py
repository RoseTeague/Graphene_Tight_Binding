#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for 1D correlated Gaussian disorder potential on graphene
"""

import numpy as np
import DD_WP_G
from DD_WP_G import *

def oneDdisorderpotential(m,n,lc):
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

    #Calling Crystal to obtain the positions of all carbon atoms.
    points = Crystal(m, n)

    #Exctracting each unique x position.
    X = points[0]
    x_1 = X[0:n*m:n,0]
    x_2 = X[1:n*m:n,0]


    #Parameters for the strength of the disorder potential in units of eV.
    Delta = 0.3

    #Generating sample of random numbers of Guassian distribution. 
    V = np.random.normal(0,1,m)

    #Generate the two-point matrix for the first row of atoms
    X1 = np.tile(x_1,(m,1))
    X2 = X1.T
    #Generate the two-point matrix for the second row of atoms
    X3 = np.tile(x_2,(m,1))
    X4 = X3.T

    #Generate the two-point spatial correlation matrix for two rows
    C1 = Delta**2*np.exp(-abs(X1-X2)/lc)
    C2 = Delta**2*np.exp(-abs(X3-X4)/lc)

    #Cholesky decomposition of the two-point correlation matrix and generate the final random vector for each row
    L1 = np.linalg.cholesky(C1)
    W1 = np.dot(L1,V)
    L2 = np.linalg.cholesky(C2)
    W2 = np.dot(L2,V)

    #Reshaping for further calculation.
    Wfinal = np.sort(W1.tolist()+W2.tolist())

    return Wfinal

if __name__ == "__main__":
    oneDdisorderpotential(10,10)
