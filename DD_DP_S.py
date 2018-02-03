#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for correlated Gaussian disorder potential on square lattice
"""

import numpy as np
import DD_WP_S
from DD_WP_S import *

def oneDdisorderpotential(m,n):
    """
    ============================================================================
    Create one-dimensional disorder potentials on square/rectangular lattice
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

    points = Crystal(m, n)

    X = points[0]
    x = X[0:n*m:n,0]

    Delta = 0.3
    V=np.random.normal(0,1,m)
    X1=np.tile(x,(m,1))
    X2=X1.T
    C=Delta**2*np.exp(-abs(X1-X2)/lc)
    L=np.linalg.cholesky(C)
    W=np.dot(L,V)
    Wfinal=W.reshape((m,1))

    return Wfinal

if __name__ == "__main__":
    oneDdisorderpotential(10,10)
