"""
Usage:
  rb_stability.py [--Pr=<prandtl> --RaT=<rayleighT> --RaT_range=<RaT_range> --kx_range=<kx_range> \
  --nRaT=<nRaT> --nkx=<nkx> --nz=<nz>]

Options:
  --Pr=<prandtl>            Prandtl number [default: 1]
  --RaT=<rayleighT>         Rayleigh number for thermal convection [default: 1708]
  --RaT_range=<RaT_range>   range of RaT [default: None]
  --kx_range=<kx_range>     range of kx [default: None]
  --nRaT=<nRaT>             search-grid size of RaT-axis for crit finder [default: 20]
  --nkx=<nkx>               search-grid size of kx-axis for crit finder [default: 20]
  --nz=<nz>                 Chebyshev discretization in z [default: 32]
"""

"""
Neutral stability curves for 2d Rayleigh-Benard convection 
Ref: 
1. Hydrodynamic and hydromagnetic stability, Chandrasekhar
2. Hydrodynamic stability, Drazin and Ried

Author: Pratik Aghor
"""
import numpy as np
from numpy import sqrt
import h5py

from dedalus import public as de
from eigentools import Eigenproblem, CriticalFinder
import matplotlib.pyplot as plt
from docopt import docopt
import os
from mpi4py import MPI

import logging
logger = logging.getLogger(__name__.split('.')[-1])


###########################
args=docopt(__doc__)
###########################
Pr = float(args['--Pr'])
RaT = float(args['--RaT'])
nRaT = int(args['--nRaT'])
nkx = int(args['--nkx'])
nz = int(args['--nz'])
kx=3.117 # initialize the parameter kx

if args['--RaT_range'] == 'None':
    RaT_range = (2400, 7200)
else:
    RaT_range = [float(i) for i in args['--RaT_range'].split(',')]
if args['--kx_range'] == 'None':
    kx_range = (0.1, 8)
else:
    kx_range = [float(i) for i in args['--kx_range'].split(',')]

if MPI.COMM_WORLD.rank == 0:
    logger.info("Computing Stability for Pr={:.3e}, nz = {:d}".format(Pr, nz))

variables = ['psi', 'psi_z', 'omega', 'omega_z', 'T', 'T_z']

# domain
r_basis = de.Chebyshev('z', nz, interval=[0, 1])

bases = [r_basis]
domain = de.Domain(bases, comm=MPI.COMM_SELF)

root_name = "rb_neutral_Pr_{:.3e}_".format(Pr)

# problem
problem = de.EVP(domain, eigenvalue='sigma', variables=variables)

# params into equations
problem.parameters['RaT']=RaT
problem.parameters['Pr']=Pr
problem.parameters['pi']=np.pi        
problem.parameters['RaT'] = RaT
problem.parameters['kx'] = kx

# substitutions
problem.substitutions['dx(f)'] = '1j*kx*f'
problem.substitutions['dt(f)'] = 'sigma*f'

"""
lap --> scalar laplacian
"""
problem.substitutions['lap(f, f_z)'] = "(dx(dx(f)) + dz(f_z))"
problem.substitutions['poission_bracket(f, f_z, g, g_z)'] = "dx(f)*g_z - f_z*dx(g)" # -1*convection in terms of f=psi

# base state if z \in [-1, 1]
# problem.substitutions['T0'] = "0.5*(1-z)"
# problem.substitutions['T0_z'] = "-0.5"

# base state if z \in [0, 1]
problem.substitutions['T0'] = "(1-z)"
problem.substitutions['T0_z'] = "-1"

# vorticity equation
problem.add_equation("omega + lap(psi, psi_z) = 0")

problem.add_equation("dt(omega) - (RaT/Pr)*dx(T) - lap(omega, omega_z) = 0")

# energy equation
problem.add_equation("Pr*dt(T) - Pr*T0_z*dx(psi) - lap(T, T_z) = 0")

# auxillilary equations
problem.add_equation("psi_z - dz(psi) = 0")
problem.add_equation("omega_z - dz(omega) = 0")
problem.add_equation("T_z - dz(T) = 0")

# boundary Conditions
problem.add_bc("left(psi) = 0")
problem.add_bc("right(psi) = 0")
problem.add_bc("left(psi_z) = 0")
problem.add_bc("right(psi_z) = 0")
problem.add_bc("left(T) = 0")
problem.add_bc("right(T) = 0")

ep = Eigenproblem(problem)

cf = CriticalFinder(ep, ("kx","RaT"), find_freq=True)
nx = nkx
ny = nRaT
xpoints = np.linspace(kx_range[0], kx_range[1], nx)
ypoints = np.linspace(RaT_range[0], RaT_range[1], ny)

grid_file = 'results/'+root_name+'grid'
if os.path.exists(grid_file+'.h5'):
    if MPI.COMM_WORLD.rank == 0:
        cf.load_grid(grid_file+'.h5')
else:
    cf.grid_generator((xpoints,ypoints))
    cf.save_grid(grid_file)
if MPI.COMM_WORLD.rank == 0:
    try:
        crits = cf.crit_finder(maxiter=300)
        print("crits = ", crits)
        logger.info("Pr = {:3e}: Critical RaT = {:5.2f}, kx = {:7.5f}".format(Pr, crits[1],crits[0]))
    except ValueError:

        crits = None
        print("Critical finder failed.")
    pax,cax = cf.plot_crit()

    """
    copied from eigentools CriticalFinder.plot_crit()
    find neutral stab points and save them for plotting later
    """
    x = cf.parameter_grids[0][0,:]
    y = cf.roots[:]
    y, x = y[np.isfinite(y)], x[np.isfinite(y)]
    # print("Re_n = ", y)
    # print("kz_n = ", x)
    hf = h5py.File('results/'+root_name+'neutral.h5', 'w')
    hf.create_dataset('kx_n', data=x)
    hf.create_dataset('RaT_n', data=y)

    # pax.collections[0].set_clim(-0.03,-0.08)

    cax.xaxis.set_ticks_position('top')
    cax.xaxis.set_label_position('top')
    pax.contour(cf.parameter_grids[0], cf.parameter_grids[1],cf.evalue_grid.real, levels=(0,), colors='white')
    pax.figure.savefig('figs/'+root_name+'growth_rates.png',dpi=300)
