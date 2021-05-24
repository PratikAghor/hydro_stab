import numpy as np
from scipy.linalg import eig
from numpy.linalg import cholesky, inv, pinv, cond

import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
from numpy import pi, cos, sin, eye, sqrt, arange, real

from cheb import *
from rbc_2d_eqns import *
##############################
"Define the inputs"
suffix = "_rbc_2d"
# Ra = 1709.
Pr = 1
# N = 10
Nz = 32

kxmin = 1
kxmax = 3.3
dkx = 0.1

kxvec = arange(kxmin, kxmax + dkx, dkx)
nz = Nz+1

er = 1
zc, D = cheb(Nz)
"""
zc = az + b
zc in [-1, 1]
z in [0, 1]
a = 2, b = -1
d/dz = a*(d/dzc) ...chain rule
d/dz = a*D
"""
a = 2; b = -1
z = (zc - b)/a
D = a*D
D2 = np.matmul(D, D)
Ra_vs_kx = zeros((len(kxvec), 2))

rmin = 1
#############################
for i in range(0, len(kxvec)):
    kx = kxvec[i]
    L, M = rbc_2d_marginal_eqns(Pr, Nz, z, D, D2, kx, er)
    #############################
    eigvals, Uos = eig(L, M, check_finite=True)
    idx = np.argsort(abs(eigvals))
    eigvals = eigvals[idx]

    # Uos = Uos[:, idx]
    Ra_vs_kx[i, 0] = kx
    Ra_vs_kx[i, 1] = real(eigvals[0])

    print("kx =", kx, " Ra_cr =", Ra_vs_kx[i, 1])

    print("eigvals = \n", eigvals)

    #############################

fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes

# Plot  G(t) vs t
ax.plot(Ra_vs_kx[:, 0], Ra_vs_kx[:, 1], linewidth=2, color='k')

ax.set_xlabel(r'$k_{x}$', fontsize=20)  # Set x label
ax.set_ylabel(r'$Ra$', fontsize=20)  # Set y label
# ax.legend(loc=1)
# ax.grid()
fig.savefig('marginal_stab' + suffix + '.png')

#############################
#############################
