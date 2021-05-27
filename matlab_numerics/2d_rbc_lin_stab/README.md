# Linear Stability of 2d Rayleigh-Benard Convection: 
## Author: Adhithya Sivakumar
### Governing equations:
* Can be found in Drazin and Ried (Eqns. 8.6, 8.7, 8.8) with no-slip BCs
* Want to reproduce the marginal stability curve and other linear stability results using `primitive variables'.
* Derived artificial BCs for pressure by taking the ```z``` component of the momentum equations and evaluating it at the boundaries.
* ```dp/dz = 0``` at both rigid boundaries. 
* Replaced incompressibility condition with Poisson equation for pressure. 

### Files and usage: 

* ```cheb.m``` produces the Chebyshev-Gauss-Lobatto grid and the Cheb differentiation matrix
* ```get_linop.m``` produces the matrices A and B for the growth rate problem, for a given N, k, R, Pr. For the same parameters, ```get_maxeig.m``` computes the maximum growth rate. For a range of k (with the other parameters fixed), ```grate_study.m``` is a script to plot the growth rate curve and the dominant eigenmodes. 
* ```get_mstab_linop.m``` produces the matrices A and B for the marginal stability problem, for a given N, k, R. For the same parameters, ```get_mstab_mineig.m``` computes the marginal Rayleigh number. For a range of k (with the other parameters fixed), ```mstab_study.m``` is a script to plot the marginal stability curve. 
* ```help_plotV``` plots the eigenmodes corresponding to a particular wavenumber in the selected range.
* ```testcode``` and ```testcode3``` are tests, as the name suggests. 

### May 27, 2021:
* The curves look alright in form, but numerical agreement hasn't been attained yet. 




