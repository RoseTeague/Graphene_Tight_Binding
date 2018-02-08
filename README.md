# Tight Binding Solver
#### Version 1.0  -   13/02/2018

A Tight binding solver to study the evolution of Gaussian wavepackets in 2D lattices.
Uses a split Hamiltonian approach to produce fast results and is capable of creating
.mp4 videos of the propagation.

## Setup

This script is written in python 3 and uses only pre-installed modules
**(numpy, matplotlib, scipy)**

All modules for the **crystal structure**, **initial wavepacket**,
**Tight Binding Hamiltonian**, and **solvers** have been included for
**square (S)** and **hexagonal (G)** lattices. Additionally, an **animate.py**
file in included to produce .mp4 files and, most importantly, the **DD_Main.py**
file, which should be the only file you need to modify for simulate different
systems.

Once you have downlaoded these, install cv2 using the command below which will be
necessary for producing the video files.

```ruby
pip install opencv-python
```

**Create an empty folder called Images into the work directory**


## Modules
If you wish to use your own input modules for different structures, they must have
the following form:

### Crystal
Describes the structure of your 2D lattice. This module defines the set of
coordinates in real space which are occupied by atoms. The coordinates should move
down each column, then to the top of the next, until all lattice points have been
described.

*Inputs*
m - number of atoms along the x direction
n - number of atoms along the y direction

*Returns*
X   - complete list of atomic x-positions
Y   - complete list of atomic y-positions
X_0 - initial x-position of the centre of the wavepacket
Y_0 - initial y-position of the centre of the wavepacket
    *N.B - X_0 and Y_0 __must__ be atomic positions*

#### Example: Square Lattice; m=4, n=4
.  .  .  .  
.  .  .  .
.  .  .  .  
**.**  .  .  .  

X  = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]
Y  = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
X0 = 0
Y0 = 3

### Psi
Describes the initial wavepacket to be propagated through the system. This
module determines the value of a gaussian centred at (X_0, Y_0) at each atomic
position on the lattice, again moving down each column from left to right.

*Inputs*
s     - width of gaussian Wavepacket
kx,ky - wavenumbers along x and y directions
m     - number of atoms along the x direction
n     - number of atoms along the y direction
pos   - output from the Crystal function

*Returns*
Psi   - an (n*m)x1 matrix with the value of the wavefunction defined at each atom
          elements moving from top left to bottom right

#### Example : Square Lattice; m=4, n=4
.  .  .  .  
.  .  .  .
.  .  .  .  
**.**  .  .  .  

Psi = [[0,0,1,5,  0,0,0,1,  0,0,0,0,  0,0,0,0]]

### Potential
Desribes the external potential that the wavepacket is within. This module can
be defined as a **1D** or **2D** potential, where the value of V is determined
at each x-position **or** at each atomic location respectively.

*Inputs*
m      - number of atoms along the x direction
n      - number of atoms along the y direction
pos    - output from the Crystal function
------
lc     - (optional), the disorder length

*Returns*
Wfinal - value of V at each: x-position
                  **OR**   : atomic location

#### Example : Square Lattice; m=4, n=4
