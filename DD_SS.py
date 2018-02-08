"""Module for 2 D problem
"""

import numpy as np
from scipy.sparse.linalg import expm_multiply
import matplotlib.pyplot as plt

def TB_ss(n, m, pos, wfc, H, DT, dt,video=False):
    """Split operator technique for propogation of wave packets with square
        lattice tight-binding hamiltonian.

    """

    Ns = round(DT/dt) + 1

    for i in range(Ns):

        wfc = expm_multiply(H, wfc)

        if video:

            wfc_c = np.conjugate(wfc)

            pd = np.multiply(wfc_c, wfc)
            pd = np.reshape(pd,(n,m))
            plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)),pd, 100, cmap = 'gnuplot')
            plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
            plt.savefig('Images/'+str(i))

    wfc_c = np.conjugate(wfc)

    pd = np.multiply(wfc_c, wfc)
    pd = np.reshape(pd,(n,m))

    return pd
