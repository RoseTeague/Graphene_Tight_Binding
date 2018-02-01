#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:03:11 2018

@author: amywang
"""

import numpy as np
import matplotlib.pyplot as plt

def oneDdisorderpotential(x,y,lc,v0,a):
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
    Nx=x.size() #get the number of atoms along x axis
    Ny=y.size() #get the number of atoms along y axis
    hbar=1.054*1E-34
    delta=4*np.pi*hbar*v0/lc
    V=np.random.normal(0,1,Nx) #Create a Gaussian-random distribution array
    Wfinal=np.empty((Ny,Nx)) #Define the final potential matrix
    sign=1
    for z in range(Ny):
        C=np.empty((Nx,Nx))
        for i in range(Nx):
            for j in range(Nx):
                C[i,j]=delta**2*np.exp(-abs(x[i]-x[j])/lc)
        L=np.linalg.cholesky(C)   
        W=np.dot(L,V)
        Wfinal[z]=W
        sign=sign*(-1)
        x = x+sign*3*a/2
    return Wfinal


