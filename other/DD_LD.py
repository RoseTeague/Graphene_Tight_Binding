#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for Tight-Binding Hamiltonian for the graphene with line defects
"""

def linedefect(p,q,m,n,psi,H):
    """
    ============================================================================
                Function for creating a line defect in the system
    ============================================================================

    Inputs
    ------
    p : integer,
        index which controls which column of atoms has  defect

    q : integer,
        index which controls the starting position of the line defect in a column

    m : integer,
        Number of atoms along the x-direction

    n : integer,
        Number of atoms along the y-direction

    psi : arrary,
        imported wave packet from module

    H : tuple of arrays,
        imported matrices for split operator technique

    Returns
    -------
    psi - array,
        wave packet with line defect

    TH1P - array (3,n*m),
        tridiagonal matrix hamiltonian with line defect

    TH1N - array (3,n*m),
        tridiagonal matrix hamiltonian with line defect

    TH2P - array (3,n*m),
        tridiagonal matrix hamiltonian with line defect

    TH2N - array (3,n*m),
        tridiagonal matrix hamiltonian with line defect

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert type(p) is int, "Index for line defect must be an integer"
    assert type(q) is int, "Index for line defect must be an integer"

    #Generate a new wave-packet with the values of the pth column to be zero.
    psi[p*n + q:(p + 1)*n] = 0

    #Splitting the Hamiltonian output into corresponding matrices
    TH1P = H[0]
    TH1N = H[1]
    TH2P = H[2]
    TH2N = H[3]

    #Generate a new set of Hamiltonian with the values of the pth column to be zero.
    TH1P[0:3,p*n + q:(p + 1)*n] = 0
    TH1N[0:3,p*n + q:(p + 1)*n] = 0

    TH2P[0:3,p + q*m:p + m*(n - 1):m] = 0
    TH2N[0:3,p + q*m:p + m*(n - 1):m] = 0

    return psi, TH1P, TH1N, TH2P, TH2N
