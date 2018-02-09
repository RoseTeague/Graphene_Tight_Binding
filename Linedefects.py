#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module for Tight-Binding Hamiltonian for the graphene with line defects
"""


def linedefect(p,m,n,psi,LISTH):
    
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
        
    LISTH: list
        A list containing four original Hamiltonians

    Returns
    ----------
    LISTHnew : list
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
    ListHnew=LISTH
    for i in ListHnew:
        i[0:3,(p-1)*n:p*n]=0
        
    return psinew, ListHnew
    
if __name__ == "__main__":
    from  DD_GH import TBH
    from  DD_WP_G import Crystal, Psi
    pos = Crystal(5,5)
    psi = Psi(14, 0.1, 0.1, 5, 5, pos)
    LISTH=TBH(5,5,dt=0.1e-15,V='one dimensional')
    psinew, ListHnew=linedefect(3,5,5,psi,LISTH)
    