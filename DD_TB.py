"""Module for solver of with TB split operator
"""
import math #can not remeber if this is how we do this ...
import numpy as np
import matplotlib as mpl
from scipy import constants
from scipy import sparse
import scipy.linalg
from DD_TBH import TBH

#import scipy.fftpack ... this might be a little awkward

#this is for the solver part ...

#in the function header we need to have the timestep and duration of simulation specified

def TB_solver(dt,DT):
    #add doc string here

    #need to add some assert statements ...

    #probably want to have a file with all of the details ...
    Ns = DT/dt #need to then round this off to the nearest integer ...
    N = 100
    n = 10
    m=10
    H = TBH(n=10,m=10,dt=0.1e-15,V=False)
    TH1P = H[0]
    TH1N = H[1]
    TH2P = H[2]
    TH2N = H[3]    
    #should call the wavefunction ...
    #need to specify how it is constructed so we can make the tridagonly matrices
    wvf = np.ones((N,1))#wavefunction ... just put some matrix there for now ...

    #import from TBH

    #then need to call the hamiltonian ...

    for i in range(Ns):

        #matrix multiply the wavefunctino by the + version of the first tri diagonal
        

        #https://stackoverflow.com/questions/44388358/python-numpy-matrix-multiplication-with-one-diagonal-matrix
        #c = (a*b.T).T ... fastest way and we have it in this form ...

        #solve the linear equation with the - version of the first tri diagonal and
        #the vector just calcualted
#       NX1 matrix
        
        
        
        hnm = sparse.csr_matrix.dot(TH1N,wvf)
        A1 = scipy.linalg.solve_banded((1,1), TH1P, hnm) 
        Anm = A1.reshape((n,m))
        A2 = Anm.reshape((N,1))
        
        f=sparse.csr_matrix.dot(TH2N,A2)
        B = scipy.linalg.solve_banded((1,1), TH2P, f) 
        Bnm = B.reshape((m,n))
        B2 = Bnm.reshape((N,1))
        
        
        g = sparse.csr_matrix.dot(TH1N,B2)
        wvf = scipy.linalg.solve_banded((1,1), TH1P, g) 
#    Calculation of Conjugate wave function
    wvf_conj=np.conjugate(wvf)

#   Probability density function by element wise multiplication
    pd = np.multiply(wvf_conj,wvf)
    
    
    
    
    
#contour plot
def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(1),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)




def visualize(Nx=50,Ny=50,yl=-2,xl=-2,yu=2,xu=2,noise=False,na = 0,display=True):
    """Display cost function with and without noise on an Ny x Nx grid, which is
    determined by the input parameters. Calls the vcost2d subrouine in cost.f90

    Input
    ----------------------------------------------------------------------------
    Nx - number of points on the x1 axis; integer
    Ny - number of points on the x2 axis; integer
    xl - lower value of x1; real number; smaller than xu
    xu - upper value of x1; real number; larger than xl
    yl - lower value of x2; real number; smaller than yu
    yu - upper value of x2; real number; larger than yl
    noise - switches noise on and off; boolean
    na - noise amplitude; real number
    display - determines if a figure will be displayed; boolean

    Output
    ----------------------------------------------------------------------------
    Displays figure if display == True
    x1 - meshgrid of x1
    x2 - meshgrid of x2
    z - cost of j in meshgrid
    """

    assert type(Nx) is int, "Number of points in x1 must be an integer"
    assert type(Ny) is int, "Number of points in x2 must be an integer"
    assert xu > xl, "xu must be larger than xl"
    assert yu > yl, "yu must be larger than yl"
    assert type(noise) is bool, "noise must be a boolean"
    assert type(display) is bool, "display must be a boolean"

    #setting the noise values for global parameters in cost module
    cost.c_noise = noise
    cost.c_noise_amp = na

    #sets up mesh grid for plotting and arrays for vcost2d subrouine
    x1_g = np.arange(xl, xu, ((xu - xl)/Nx))
    x2_g = np.arange(yl, yu, ((yu - yl)/Ny))
    x1, x2 = np.meshgrid(x1_g, x2_g)

    #Calling the subroutine
    z = cost.vcost2d(x1_g, x2_g)

    #plotting
    if display:

        plt.figure()
        plt.contourf(x1, x2, z, locator=ticker.LogLocator(),cmap='RdGy')
        cb = plt.colorbar()
        cb.set_label('Z')
        plt.scatter(1,1,c='black',marker='o')
        plt.xlabel(r'$x_{1}$', fontsize = 15)
        plt.ylabel(r'$x_{2}$', fontsize = 15)
        plt.locator_params('x',tight=True, nbins =8)
        plt.locator_params('y',tight=True, nbins =8)
        plt.title('Zachary Goodwin - visualize')
        plt.tick_params(axis='both',which='major',labelsize =15,direction='in',right=True,top=True)
        plt.tight_layout()
        #plt.savefig('hw212.pdf')
        #plt.savefig('hw212.png')
        plt.show()

    return x1,x2,z

#End of plotting stuff


 
        #Now we need to shift the result around such that another tri diagonal
        #matrix equation can be solved
        #might be able to do this with the range thing from above ...
        #should be able to do this with something similar to DD_TBH

        #After we have reordered the result, we multiply with the + version of the
        #second tri diagonal

        #then solbe the linear equation with the - version of he second tri diagonal
        #and the vector just calculated

        #The result needs to be shifted around again so we can solve with the
        #first tri diagonal matrix again

        #repeat the first few steps ...

        #now we have our wavefunction at dt
