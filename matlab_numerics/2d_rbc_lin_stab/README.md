# Linear Stability of 2d Rayleigh-Benard: 
## Author: Pratik Aghor
### Governing equations:
* Can be found in Drazin and Ried (Eqns. 8.6, 8.7, 8.8) with no-slip BCs
* Want to reproduce the marginal stability curve and other linear stability results using `primitive variables'.
* Derived artificial BCs for pressure by taking the ```z``` component of the momentum equations and evaluating it at the boundaries.
* ```dp/dz = 0``` at both rigid boundaries. 
### Files and usage: 

* ```cheb.py``` gives the Chebyshev-Gauss-Lobatto grid and the Cheb differentiation matrix
* ```rbc_2d_eqns.py``` contains the dimensionless linear equations for the growth-rate problem (sigma != 0) and the marginal stability problem (sigma = 0, Ra as the eigenvalue) 
* ```test_rbc_2d_marginal_stab.py``` contains the test for the marginal stability curve. 

### May 24, 2021:
* The codes don't work yet. 




