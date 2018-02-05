#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:03:11 2018

@author: amywang
"""

import numpy as np
import matplotlib.pyplot as plt
import DD_WP_S#this does not work with the graphene ... ?!
from DD_WP_S import *#this does not work with the graphene ... ?!

def oneDdisorderpotential(m,n):
    """
    Create one-dimensional disorder potentials

    Parameters
    ----------
    x : float, array
        The x-position of carbon atoms
    y : float, array
        The y-position of carbon atoms
    lc: float
        The disorder correlation length
    v0: float
        The band velocity of the electrons
    a: float
        Lattice size parameter

    Returns
    -------
    Wfinal: float, array
         The final potential at each atom
    """
    #this could just be replaced by n and m ...

    points = Crystal(m, n)

    X = points[0]
    x = X[0:n*m:n,0]
    #print(np.shape(x))
    lc = 20
    Delta = 0.3
    V=np.random.normal(0,1,m)

    X1=np.tile(x,(m,1))
    X2=X1.T
    C=Delta**2*np.exp(-abs(X1-X2)/lc)
    L=np.linalg.cholesky(C)  
    W=np.dot(L,V)
    Wfinal=W.reshape((m,1))

    #print(W)
    return Wfinal

if __name__ == "__main__":
    oneDdisorderpotential(100,100)