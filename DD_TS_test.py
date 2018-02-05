#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:05:02 2018

@author: georgios
"""
import numpy as np

import scipy.linalg

from DD_2D import TB_solver_2D
from DD_TB_S import *
from DD_WP_S import *



def Comparison(n, m, dt, DT):
        
    points = Crystal(m, n)
    pos = Crystal(100,100)

    wvf2D = TB_solver_2D(n, m, dt, DT)
    wvf_S = TB_solver_S(n, m, dt, DT)
    differnce = wvf2D - wvf_S

    
    Plotting
    plt.contourf(points[0].reshape((n,m)),points[1].reshape((n,m)),difference)
    plt.show()
    return wvf2D, wvf_S, differnce


if __name__ == "__main__":
    #Crystal(5,5)
    #TBH(5,5,dt=0.1e-15)
    n = 100
    m = 100
    dt=0.1e-15
    DT = 1e-15
    Comparison(n, m, dt, DT)