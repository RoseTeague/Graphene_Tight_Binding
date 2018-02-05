#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for 2D correlated Gaussian disorder potential on graphene
"""

import numpy as np
import DD_WP_G
from DD_WP_G import *

def twoDdisorderpotential(m,n,lc):
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
    Wfinalx: float, array
         The final potential at each x position
    Wfinaly: float, array
         The final potential at each y position
    """

    #Calling Crystal to obtain the positions of all carbon atoms.
    points = Crystal(m, n)

    #Exctracting each unique x position.
    X = points[0]
    x_1 = X[0:n*m:n,0]
    x_2 = X[1:n*m:n,0]

    #Exctracting each unique y position.
    Y = points[1]
    y = Y[0:n,0]

    #Generate all x positions in an ascending order
    x=np.sort(x_1.tolist()+x_2.tolist())

    #Parameters for the strength of the disorder potential in units of eV.
    Delta = 0.3

    #Generating sample of random numbers of Guassian distribution along x axis
    Vx = np.random.normal(0,1,2*m)

    #Generating sample of random numbers of Guassian distribution along y axis
    Vy = np.random.normal(0,1,n)

    #Generate the two-point matrix for x axis
    X1 = np.tile(x,(2*m,1))
    X2 = X1.T

    #Generate the two-point matrix for y axis
    X3 = np.tile(y,(n,1))
    X4 = X3.T

    #Generate the two-point spatial correlation matrix for x and y axis
    C1 = Delta**2*np.exp(-abs(X1-X2)/lc)
    C2 = Delta**2*np.exp(-abs(X3-X4)/lc)

    #Cholesky decomposition of the two-point correlation matrix and generate the final random vector
    L1 = np.linalg.cholesky(C1)
    W1 = np.dot(L1,Vx)
    L2 = np.linalg.cholesky(C2)
    W2 = np.dot(L2,Vy)

    #Reshaping for further calculation.
    Wfinalx = W1.reshape((2*m,1))
    Wfinaly = W2.reshape((n,1))

    return Wfinalx, Wfinaly

if __name__ == "__main__":
    DP2 = twoDdisorderpotential(10,10,10)
