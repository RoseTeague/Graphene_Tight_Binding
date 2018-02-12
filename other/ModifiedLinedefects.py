#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module for Tight-Binding Hamiltonian for the graphene with line defects
"""


def linedefect(p,m,n,psi,TH1P,TH1N,TH2P,TH2N):
    
    """
    Function for creating a line defect in the system
    
    Inputs
    -----------
    p : integer
        The index of the column of the line defect

    m : integer
        Number of atoms along the x-direction
    
    n : integer
        Number of atoms along the y-direction

    Returns
    ----------
    LISTH : list
        A list containing the four Hamitonians (with a line defect) needed in further calculations.

    """
    
    
    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert n % 2 == 0, "The Hamiltonian can only be constructed for an even number of rows"
    assert m % 2 == 0, "The Hamiltonian can only be constructed for an even number of columns"
    
    
    
    #Generate a new wave-packet with the values of the pth column to be zero.  
    psinew=psi
    psinew[(p-1)*n:p*n]=0
    
    #Generate a new set of Hamiltonian with the values of the pth column to be zero.
    ListH=[TH1P,TH1N,TH2P,TH2N]
    for i in ListH:
        i[0:3,(p-1)*n:p*n]=0
        
    return psinew, ListH