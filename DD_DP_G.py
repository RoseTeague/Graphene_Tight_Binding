#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for 1D correlated Gaussian disorder potential on graphene
"""

import numpy as np
import DD_WP_G
from DD_WP_G import *

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

    #Generate all x positions in an ascending order
    x = np.sort(x_1.tolist()+x_2.tolist())

    #Parameters for the strength of the disorder potential in units of eV.
    Delta = 0.3

    #Generating sample of random numbers of Guassian distribution.
    V = np.random.normal(0,1,2*m)

    #Generate the two-point matrix
    X1 = np.tile(x,(2*m,1))
    X2 = X1.T

    #Generate the two-point spatial correlation matrix
    C = Delta**2*np.exp(-abs(X1-X2)/lc)


    #Cholesky decomposition of the two-point correlation matrix and generate the final random vector
    L = np.linalg.cholesky(C)
    W = np.dot(L,V)

    #Reshaping for further calculation.
    Wfinal = W.reshape((2*m,1))

    return Wfinal

if __name__ == "__main__":
    DP = oneDdisorderpotential(10,10,10,pos)
    plt.plot(DP)
    plt.show()
