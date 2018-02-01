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

    C=np.empty((m,m))
    for j in range(m):
        C[:,j]=Delta**2*np.exp(-abs(x-x[j])/lc)

    Wfinal=np.empty((n,m))
    for z in range(n):

        L=np.linalg.cholesky(C)
        W=np.dot(L,V)
        Wfinal[z]=W

    #Wfinal_mn = Wfinal.reshape((n,m))
    Wfinal_T = Wfinal.T
    Wfinal = Wfinal_T.reshape((n*m,1)).T
    #print(Wfinal)
    W = Wfinal[0,0:n*m:n]
    W = W.reshape((m,1))

    #print(W)
    return W

if __name__ == "__main__":
    oneDdisorderpotential(100,100)
