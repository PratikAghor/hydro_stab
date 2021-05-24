import numpy as np
from scipy.linalg import eig

import matplotlib.pyplot as plt
from numpy import pi, cos, sin, eye, zeros

from cheb import *
##############################
"""
OS eqn: Lq = i \omega M q
get the OS operator for plane Poisuille flow (ppf) and plane Couette flow (pcf)
Author: Pratik Aghor
"""
##############################
def rbc_2d_eqns(Ra, Pr, Nz, z, D, D2, kx, er):
    """
    define inputs
    """
    n = Nz+1

    ksq = kx*kx
    Idn = eye(n)
    """
    unknowns are [\hat{u}, \hat{w}, \hat{p}, \hat{\theta}]^{T}
    equations order
    0: x-mom
    1: z-mom
    2: continuity
    3: energy
    """
    # x-mom
    L00 = Pr*(D2 - ksq*Idn)
    L01 = zeros((n, n))
    L02 = -1j*kx*Idn
    L03 = zeros((n, n))

    # z-mom
    L10 = zeros((n, n))
    L11 = Pr*(D2 - ksq*Idn)
    L12 = -D
    L13 = Ra*Pr*Idn

    # continuity
    L20 = 1j*kx*Idn
    L21 = D
    L22 = zeros((n, n))
    L23 = zeros((n, n))

    # energy
    L30 = zeros((n, n))
    L31 = Idn
    L32 = zeros((n, n))
    L33 = (D2 - ksq*Idn)

    L = np.block([[L00, L01, L02, L03], \
    [L10, L11, L12, L13], \
    [L20, L21, L22, L23], \
    [L30, L31, L32, L33]
    ])

    """
    RHS: \sigma M [[\hat{u}, \hat{w}, \hat{p}, \hat{\theta}]^{T}]
    """
    M00 = Idn; M01 = zeros((n, n)); M02 = zeros((n, n)); M03 = zeros((n, n));
    M10 = zeros((n, n)); M11 = Idn; M12 = zeros((n, n)); M13 = zeros((n, n));
    M20 = zeros((n, n)); M21 = zeros((n, n)); M22 = zeros((n, n)); M23 = zeros((n, n));
    M30 = zeros((n, n)); M31 = zeros((n, n)); M32 = zeros((n, n)); M33 = Idn;


    M = np.block([[M00, M01, M02, M03], \
    [M10, M11, M12, M13], \
    [M20, M21, M22, M23], \
    [M30, M31, M32, M33]
    ])

    """
    BC for 2d RBC:
    u = v = w = theta = 0 at walls

    following Schmid and Henningson
    implenment:
    200i v = v and 200i Dv = Dv

    This helps identify "0" eigenvalues
    obtained from the BCs (while sorting,
    we'll have omega = 200i now)
    and find the true 0-eigvals corresponding to transitions
    """
    # er=-200*1j;
    """
    u = 0
    """
    L00[0, :] = 0.; L01[0, :]  = 0.; L02[0, :]  = 0.; L03[0, :] = 0.;
    L00[0, 0] = 1.*er;

    L00[n-1, :] = 0.; L01[n-1, :]  = 0.; L02[n-1, :]  = 0.; L03[n-1, :] = 0.;
    L00[n-1, n-1] = 1.*er;

    """
    w = 0
    """
    L10[0, :] = 0.; L11[0, :]  = 0.; L12[0, :]  = 0.; L13[0, :] = 0.;
    L11[0, 0] = 1.*er;

    L10[n-1, :] = 0.; L11[n-1, :]  = 0.; L12[n-1, :]  = 0.; L13[n-1, :] = 0.;
    L11[n-1, n-1] = 1.*er;

    """
    continuity
    """
    # L20[0, :] = 0.; L21[0, :]  = 0.; L22[0, :]  = 0.; L23[0, :] = 0.;
    #
    # L20[n-1, :] = 0.; L21[n-1, :]  = 0.; L22[n-1, :]  = 0.; L23[n-1, :] = 0.;

    """
    theta = 0
    """
    L30[0, :] = 0.; L31[0, :]  = 0.; L32[0, :]  = 0.; L33[0, :] = 0.;
    L33[0, 0] = 1.*er;

    L30[n-1, :] = 0.; L31[n-1, :]  = 0.; L32[n-1, :]  = 0.; L33[n-1, :] = 0.;
    L33[n-1, n-1] = 1.*er;

    """
    u = 0
    """
    M00[0, :] = 0.; M01[0, :]  = 0.; M02[0, :]  = 0.; M03[0, :] = 0.;
    M00[0, 0] = 1.;

    M00[n-1, :] = 0.; M01[n-1, :]  = 0.; M02[n-1, :]  = 0.; M03[n-1, :] = 0.;
    M00[n-1, n-1] = 1.;

    """
    w = 0
    """
    M10[0, :] = 0.; M11[0, :]  = 0.; M12[0, :]  = 0.; M13[0, :] = 0.;
    M11[0, 0] = 1.;

    M10[n-1, :] = 0.; M11[n-1, :]  = 0.; M12[n-1, :]  = 0.; M13[n-1, :] = 0.;
    M11[n-1, n-1] = 1.;

    """
    continuity
    """
    # M20[0, :] = 0.; M21[0, :]  = 0.; M22[0, :]  = 0.; M23[0, :] = 0.;
    #
    # M20[n-1, :] = 0.; M21[n-1, :]  = 0.; M22[n-1, :]  = 0.; M23[n-1, :] = 0.;

    """
    theta = 0
    """
    M30[0, :] = 0.; M31[0, :]  = 0.; M32[0, :]  = 0.; M33[0, :] = 0.;
    M33[0, 0] = 1.;

    M30[n-1, :] = 0.; M31[n-1, :]  = 0.; M32[n-1, :]  = 0.; M33[n-1, :] = 0.;
    M33[n-1, n-1] = 1.;


    return L,M
#############################
#############################
def rbc_2d_marginal_eqns(Pr, Nz, z, D, D2, kx, er):
    """
    formulate the marginal eigval problem with Ra
    as the eigenvalue by putting \sigma = 0
    change L10, L11, L12, L13 from above,
    take Ra \hat{\theta} to the RHS
    """
    n = Nz+1

    ksq = kx*kx
    Idn = eye(n)
    """
    unknowns are [\hat{u}, \hat{w}, \hat{p}, \hat{\theta}]^{T}
    equations order
    0: x-mom
    1: z-mom
    2: continuity
    3: energy
    """
    # x-mom
    L00 = Pr*(D2 - ksq*Idn)
    L01 = zeros((n, n))
    L02 = -1j*kx*Idn
    L03 = zeros((n, n))

    # z-mom
    L10 = zeros((n, n))
    L11 =-(D2 - ksq*Idn)
    L12 = (1.0/Pr)*D
    L13 = zeros((n, n))

    # continuity
    L20 = 1j*kx*Idn
    L21 = D
    L22 = zeros((n, n))
    L23 = zeros((n, n))

    # energy
    L30 = zeros((n, n))
    L31 = Idn
    L32 = zeros((n, n))
    L33 = (D2 - ksq*Idn)

    L = np.block([[L00, L01, L02, L03], \
    [L10, L11, L12, L13], \
    [L20, L21, L22, L23], \
    [L30, L31, L32, L33]
    ])

    """
    RHS: \sigma M [\hat{u}, \hat{w}, \hat{p}, \hat{\theta}]^{T}
    eigvals, Uos = eig(L, M, check_finite=True)
    idx = np.argsort(abs(eigvals))
    eigvals = eigvals[idx]
    """
    M00 = zeros((n, n)); M01 = zeros((n, n)); M02 = zeros((n, n)); M03 = zeros((n, n));
    M10 = zeros((n, n)); M11 = zeros((n, n)); M12 = zeros((n, n)); M13 = Idn;
    M20 = zeros((n, n)); M21 = zeros((n, n)); M22 = zeros((n, n)); M23 = zeros((n, n));
    M30 = zeros((n, n)); M31 = zeros((n, n)); M32 = zeros((n, n)); M33 = zeros((n, n));


    M = np.block([[M00, M01, M02, M03], \
    [M10, M11, M12, M13], \
    [M20, M21, M22, M23], \
    [M30, M31, M32, M33]
    ])

    """
    BC for 2d RBC:
    u = v = w = theta = 0 at walls

    following Schmid and Henningson
    implenment:
    200 v = v and 200 Dv = Dv

    This helps identify "0" eigenvalues
    obtained from the BCs (while sorting,
    we'll have omega = 200i now)
    and find the true 0-eigvals corresponding to transitions
    """
    # er=-200*1j;
    """
    u = 0
    """
    L00[0, :] = 0.; L01[0, :]  = 0.; L02[0, :]  = 0.; L03[0, :] = 0.;
    L00[0, 0] = 1.*er;

    L00[n-1, :] = 0.; L01[n-1, :]  = 0.; L02[n-1, :]  = 0.; L03[n-1, :] = 0.;
    L00[n-1, n-1] = 1.*er;

    """
    w = 0
    """
    L10[0, :] = 0.; L11[0, :]  = 0.; L12[0, :]  = 0.; L13[0, :] = 0.;
    L11[0, 0] = 1.*er;

    L10[n-1, :] = 0.; L11[n-1, :]  = 0.; L12[n-1, :]  = 0.; L13[n-1, :] = 0.;
    L11[n-1, n-1] = 1.*er;

    """
    continuity
    """
    L20[0, :] = 0.; L21[0, :]  = 0.; L22[0, :]  = 0.; L23[0, :] = 0.;
    L22[0, :]  = er*D[0, :];

    L20[n-1, :] = 0.; L21[n-1, :]  = 0.; L22[n-1, :]  = 0.; L23[n-1, :] = 0.;
    L22[n-1, :]  = er*D[n-1, :];


    """
    theta = 0
    """
    L30[0, :] = 0.; L31[0, :]  = 0.; L32[0, :]  = 0.; L33[0, :] = 0.;
    L33[0, 0] = 1.*er;

    L30[n-1, :] = 0.; L31[n-1, :]  = 0.; L32[n-1, :]  = 0.; L33[n-1, :] = 0.;
    L33[n-1, n-1] = 1.*er;

    """
    u = 0
    """
    M00[0, :] = 0.; M01[0, :]  = 0.; M02[0, :]  = 0.; M03[0, :] = 0.;
    M00[0, 0] = 0.; # 1.;

    M00[n-1, :] = 0.; M01[n-1, :]  = 0.; M02[n-1, :]  = 0.; M03[n-1, :] = 0.;
    M00[n-1, n-1] = 0.; # 1.;

    """
    w = 0
    """
    M10[0, :] = 0.; M11[0, :]  = 0.; M12[0, :]  = 0.; M13[0, :] = 0.;
    M11[0, 0] = 0.; # 1.;

    M10[n-1, :] = 0.; M11[n-1, :]  = 0.; M12[n-1, :]  = 0.; M13[n-1, :] = 0.;
    M11[n-1, n-1] = 0.; # 1.;

    """
    continuity
    """
    M20[0, :] = 0.; M21[0, :]  = 0.; M22[0, :]  = 0.; M23[0, :] = 0.;
    M22[0, :] = 0.; # 1.;

    M20[n-1, :] = 0.; M21[n-1, :]  = 0.; M22[n-1, :]  = 0.; M23[n-1, :] = 0.;
    M22[n-1, :] = 0.; # 1.;


    """
    theta = 0
    """
    M30[0, :] = 0.; M31[0, :]  = 0.; M32[0, :]  = 0.; M33[0, :] = 0.;
    M33[0, 0] = 0.; # 1.;

    M30[n-1, :] = 0.; M31[n-1, :]  = 0.; M32[n-1, :]  = 0.; M33[n-1, :] = 0.;
    M33[n-1, n-1] = 0.; # 1.;

    return L,M
#############################
