#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for correlated Gaussian disorder potential on graphene
"""

import numpy as np
import DD_WP_G
from DD_WP_G import *

def oneDdisorderpotential(m,n):
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

    lc : real number
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

    #Sorting positions into ascending order.
    x = np.sort(x_1.tolist()+x_2.tolist())

    #Parameters for the strength of the disorder potential in units of eV.
    Delta = 0.3
    lc = 200

    #Generating sample of random numbers.
    V = np.random.normal(0,1,2*m)

    #
    X1 = np.tile(x,(2*m,1))
    X2 = X1.T

    #
    C = Delta**2*np.exp(-abs(X1-X2)/lc)

    #
    L = np.linalg.cholesky(C)
    W = np.dot(L,V)

    #Reshaping for further calculation.
    Wfinal = W.reshape((2*m,1))

    return Wfinal

if __name__ == "__main__":
    oneDdisorderpotential(10,10)
