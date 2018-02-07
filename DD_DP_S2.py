#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for 2D correlated Gaussian disorder potential on square lattice
"""

import numpy as np

def twoDdisorderpotential(m,n,lc,pos):
    """
    ============================================================================
    Create two-dimensional disorder potential on square lattice
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
        positions of all atoms

    Returns
    -------
    Wfinal: float, array
         The final potential at each atom
    """

    #Exctracting each unique x position.
    X = pos[0]
    x_1 = X[0:n*m:n,0]

    #Exctracting each unique y position.
    Y = pos[1]
    y = Y[0:n,0]

    #Parameters for the strength of the disorder potential in units of eV.
    Delta = 0.3

    #Generating sample of random numbers of Guassian distribution along x axis
    Vx = np.random.normal(0,1,m)

    #Generating sample of random numbers of Guassian distribution along y axis
    Vy = np.random.normal(0,1,n)

    #Generate the two-point matrix for x axis
    X1 = np.tile(x_1,(m,1))
    X2 = X1.T

    #Generate the two-point matrix for y axis
    X3 = np.tile(y,(n,1))
    X4 = X3.T

    #Generate the two-point spatial correlation matrix for x and y axis
    Cx = Delta**2*np.exp(-abs(X1-X2)/lc)
    Cy = Delta**2*np.exp(-abs(X3-X4)/lc)
  

    #Cholesky decomposition of the two-point correlation matrix and generate the final random vector
    Lx = np.linalg.cholesky(Cx)
    Wx = np.dot(Lx,Vx)
    Ly = np.linalg.cholesky(Cy)
    Wy = np.dot(Ly,Vy)

    #Put the potentials on all atoms and reshape them to the requried output format
    Wxf=np.zeros((n,m))
    Wyf=np.zeros((n,m))

    Wxf[:,0:m]= Wx1
    Wxfinal=Wxf.T.reshape((n*m,1))

    Wyf[:,0:m]=Wy.reshape((n,1))
    Wyfinal=Wyf.T.reshape((n*m,1))

    #The final potential for each atom is the sum of the x- and y- potential component
    Wfinal=Wxfinal+Wyfinal

    return Wfinal

#if __name__ == "__main__":
    #Need to import and call crystal here ... 
    #twoDdisorderpotential(5,5,10,pos)
