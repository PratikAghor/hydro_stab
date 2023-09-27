"""
Usage:
  rb_sigma_vs_k.py [--Pr=<prandtl> --RaT_range=<RaT_range> --kx_range=<kx_range> --nRaT=<nRaT> --nkx=<nkx> --nz=<nz>]

Options:
  --Pr=<prandtl>            Prandtl number [default: 1]
  --RaT_range=<RaT_range>   range of RaT [default: None]
  --kx_range=<kx_range>     range of kx [default: None]
  --nRaT=<nRaT>             search-grid size of RaT-axis for crit finder [default: 20]
  --nkx=<nkx>               search-grid size of kx-axis for crit finder [default: 20]
  --nz=<nz>                 Chebyshev discretization in z [default: 32]
"""

"""
Dispersion curves for 2d Rayleigh-Benard convection 
Ref: 
1. Hydrodynamic and hydromagnetic stability, Chandrasekhar
2. Hydrodynamic stability, Drazin and Ried

Author: Pratik Aghor
"""

import numpy as np
from numpy import sqrt
from dedalus import public as de
import logging
import h5py

from docopt import docopt
from eigentools import Eigenproblem
logger = logging.getLogger(__name__)

###########################
args=docopt(__doc__)
###########################
Pr = float(args['--Pr'])
nRaT = int(args['--nRaT'])
nkx = int(args['--nkx'])
nz = int(args['--nz'])

if args['--RaT_range'] == 'None':
    RaT_range = (2400, 7200)
else:
    RaT_range = [float(i) for i in args['--RaT_range'].split(',')]
if args['--kx_range'] == 'None':
    kx_range = (0.1, 8)
else:
    kx_range = [float(i) for i in args['--kx_range'].split(',')]

logger.info("Computing Stability for Pr={:.3e}, nz = {:d}".format(Pr, nz))

variables = ['psi', 'psi_z', 'omega', 'omega_z', 'T', 'T_z']

# domain
r_basis = de.Chebyshev('z', nz, interval=[0, 1])

bases = [r_basis]
domain = de.Domain(bases)


# find growth rate, frequency for each (kx, Ta) and store it

kpoints = np.linspace(kx_range[0], kx_range[1], nkx)
RaTpoints = np.linspace(RaT_range[0], RaT_range[1], nRaT)
sigma_vs_k = np.zeros((nRaT, nkx, 3)) # array to store growth rate and freq for each kx


root_name = "rb_sigma_vs_kx_Pr_{:.3e}".format(Pr)
# save (kz, sigma_r, sigma_i) for each Ta


for i in range(0, nRaT):
    # define a new problem
    for j in range(0, nkx):
        RaT = RaTpoints[i]
        kx = kpoints[j]

        # problem
        problem = de.EVP(domain, eigenvalue='sigma', variables=variables)

        # params into equations
        problem.parameters['Pr']=Pr
        problem.parameters['pi']=np.pi        
        problem.parameters['RaT'] = RaT
        problem.parameters['kx'] = kx
        sigma_vs_k[i, :, 0] = kpoints

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
        problem.add_equation("Pr*dt(T) - Pr*dx(psi)*T0_z - lap(T, T_z) = 0")

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

        growth, index, freq = ep.growth_rate({})
        sigma_vs_k[i, j, 1] = growth
        sigma_vs_k[i, j, 2] = freq

        logger.info("RaT = {:3e}; kx = {:3e}; Growth rate = {:16.15e}; frequency = {:16.15e}".format(RaT, kx, growth, freq))

# write sigma_vs_k

hf = h5py.File('./results/'+root_name+'.h5', 'w')
hf.create_dataset('RaTpoints', data=RaTpoints)
hf.create_dataset('sigma_vs_k', data=sigma_vs_k)


