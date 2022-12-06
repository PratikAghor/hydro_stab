import numpy as np  # Import NumPy
from numpy import pi  # Import pi from numpy
import math
# Import plotting functions:
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

"""
Code to solve Rayleigh's equation using shooting method

Related to Rayleigh's inflextion pt theorem:
phi'' - alpha^2 phi - U''/(U-c) phi = 0

Decompose into two first order odes
dphi0dt = phi1
dphi1dt = alpha^2 phi0 + U''/(U-c) phi0
ref:
    1. `Hydrodynamic stability' by Drazin and Reid.
    2. bvp.py written by Jonathan Senning, Gordon College

Author: Pratik Aghor
"""

##########################################
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}
#matplotlib.rc('font', **font)
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Lucida Grande']})
# for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})
##########################################

def getf_1z(f, z):
    """
    first derivative of f
    1. assume periodicity in z i.e., ex. assume z = [0, 2 pi), excluded last point
    2. assume equi-spaced points, i.e., dz = const 
    """
    n = len(z)
    f_1z = np.zeros(n)
    dz = z[1] - z[0];
    D = np.zeros((n, n)) # differentiation matrix

    # periodic BC
    D[0, 1] = 1; D[0, n-1] = -1
    for j in range(1, n-1):
        D[j, j+1] = 1
        D[j, j-1] = -1
    D[n-1, 0] = 1; D[n-1, n-2] = -1

    D = (1/dz)*D

    f_1z = np.matmul(D, f)

    return f_1z

#-----------------------------------------------------------------------------
def Velocity(phi, z, alpha, c, U, UdoublePrime):
    """
    Velocity function for the M equation

    Inputs:
    ssp: State space vector. dx1 NumPy array: y=[y0, y1]
    t: Time. Has no effect on the function, we have it as an input so that our
       ODE would be compatible for use with generic integrators from
       scipy.integrate

    Outputs:
    vel: velocity at phi. dx1 NumPy array: vel = [dphi0/dt, dphi1/dt]
    here U and UdoublePrime are values of U, U'' at a particular z 

    for RK4, this would vary with the time variable z
    """
    phi0, phi1 = phi  # Read state space points
    

    # Rayleigh equation:
    dphi0dt = phi1
    dphi1dt = alpha**2*phi0 + (UdoublePrime/(U-c))*phi0
    # Collect equations in a single NumPy array:
    vel = np.array([dphi0dt, dphi1dt], float)  # Velocity vector
    return vel
#-----------------------------------------------------------------------------
def RK4(velocityFunction, initialCondition, timeArray, alpha, c, U, UdoublePrime):
    """
    Runge-Kutta 4 Integrator.
    Inputs:
    VelocityFunction: Function name to integrate
                      this function must have two inputs namely state space
                      vector and time. For example: velocity(ssp, t)
    InitialCondition: Initial condition, 1xd NumPy array, where d is the
                      dimension of the state space
    TimeArray: 1 x Nt NumPy array which contains instances for the solution
               to be returned.
    Outputs:
    SolutionArray: d x Nt NumPy array which contains numerical solution of the
                   ODE.
    """
    # Generate the solution array to fill in:
    SolutionArray = np.zeros((np.size(timeArray, 0),
                              np.size(initialCondition, 0)))
    #Assign the initial condition to the first element:
    SolutionArray[0, :] = initialCondition

    for i in range(0, np.size(timeArray) - 1):
        # Read time element:
        dt = timeArray[i + 1] - timeArray[i]
        # Runge Kutta k's:
        k1 = dt * velocityFunction(SolutionArray[i, :], timeArray[i], alpha, c, U[i], UdoublePrime[i])
        k2 = dt * velocityFunction(SolutionArray[i, :]+ 0.5 * k1, timeArray[i], alpha, c, U[i], UdoublePrime[i] )
        k3 = dt * velocityFunction(SolutionArray[i, :]+ 0.5 * k2, timeArray[i], alpha, c, U[i], UdoublePrime[i] )
        k4 = dt * velocityFunction(SolutionArray[i, :]+ k3, timeArray[i], alpha, c, U[i], UdoublePrime[i] )
        # Next integration step:
        SolutionArray[i + 1] = SolutionArray[i] + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0
    return SolutionArray

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
"""
parameters for shooting:
a = value of phi at z = 0 (left BC)
b = value of phi at z = Lz = 2 pi (here) (right BC), to be satisfied within tolerance 'tol'

U = base state profile 
z = axis, treated as time
phiPrimeGuess1, phiPrimeGuess2 are two initial guesses for phi1 to solve 
Rayleigh equation

ssp = phi = [phi0, phi1], where phi1 = phi0'
ssp stands for state space point phi
"""
def shooting(a, b, z, alpha, c, U, UdoublePrime, phiPrimeGuess1, phiPrimeGuess2, tol ):

    max_iter = 25   # Maximum number of shooting iterations
    nz = len ( z ) # length of the z-array
    ssp1 = np.array([a, phiPrimeGuess1])

    # Solve the first initial value problem (IVP) with phi'(z[0]) = phiPrimeGuess1.
    phi = RK4( Velocity, ssp1, z, alpha, c, U, UdoublePrime )

    # Store only the end value as phiRightBoundary1 for comparison with the right BC at z[nz-1]
    phiRightBoundary1 = phi[nz-1, 0]
    print(0, "phiPrimeGuess1 = ", phiPrimeGuess1, ", error = ",  b - phiRightBoundary1)
    # Now, we start the actual shooting.
    for i in range(0, max_iter ):
        # Solve one more IVP with phi'(r[0]) = phiPrimeGuess2
        # Now, we keep the full solution, because, if the end value is good enough,
        # we can return the current solution as is.
        ssp2 = np.array([a, phiPrimeGuess2])
        
        phi = RK4( Velocity, ssp2, z, alpha, c, U, UdoublePrime )
        # Store only the end value as phiRightBoundary2 for comparison with the right BC at z[nz-1]
        phiRightBoundary2 = phi[nz-1, 0]
        print(0, "phiPrimeGuess2 = ", phiPrimeGuess2, ", error = ",  b - phiRightBoundary2)

        # Check if the solution obtained is good enough by checking the error
        if abs( b - phiRightBoundary2 ) < tol:
            break

        # Update phiPrimeGuesses.

        phiPrimeGuess1, phiPrimeGuess2 = \
        ( phiPrimeGuess2, phiPrimeGuess2 + ( phiPrimeGuess2 - phiPrimeGuess1 ) / ( phiRightBoundary2 - phiRightBoundary1 ) * ( b - phiRightBoundary2 ) )
        phiRightBoundary1 = phiRightBoundary2

    print("done! Returning solution...")

    return phi[:, 0]

#-----------------------------------------------------------------------------

if __name__ == "__main__":
    # This block will be evaluated if this script is called as the main routine
    # and will be ignored if this file is imported from another script.

    # parameters
    c = 0.1 # wave-speed
    alpha = 0.1 # wavenumber

    nz = 32
    z = np.zeros(nz); 
    z_min = 0; z_max = 2*np.pi
    dz = (z_max - z_min)/(nz+1)
    z =  np.arange(z_min, z_max, dz)
    print("z = ", z)
    # define the base state:
    U = np.sin(z) # Tollmein's counterexample
    UPrime = getf_1z(U, z)
    UdoublePrime = getf_1z(UPrime, z)

    phiPrimeGuess1 = 1
    phiPrimeGuess2 = -1
    tol = 1e-10
    
    phi0 = shooting(1, 1, z, alpha, c, U, UdoublePrime, phiPrimeGuess1, phiPrimeGuess2, tol)
    
    np.savetxt("phi0.txt", phi0)
    print("phi0 = \n", phi0)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    fig = plt.figure()  # Create a figure instance
    ax = fig.gca()  # Get current axes
    ax.plot(z, phi0)  # Plot the solution
    ax.grid();
    ax.set_xlabel('z')  # Set x label
    ax.set_ylabel(r'$\phi_0$')  # Set y label
    plt.tight_layout()
    ax.grid(linestyle='--')
    fig.savefig('phi0.png')

    #-----------------------------------------------------------------------------

