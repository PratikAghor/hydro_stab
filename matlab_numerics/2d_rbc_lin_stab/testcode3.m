clear
close all
clc

%% PRELIMINARY PARAMETERS

N = 64; 
k  = 2*pi;
R = 3500;
Pr = 1;

%% GENERATE CGL GRID

[~, z] = cheb(N-1);

%% CONSTRUCT LINEAR OPERATOR

[A, B, z] = get_linop(N, k, R, Pr);

%% SOLVE GENERALIZED EIGENVALUE PROBLEM

[V, D] = eig(A, B);
s = diag(D);

%% SORT EIGENVALUES

[ssr, si] = sort((real(s)),'descend');
ssi = imag(s(si));

%% MAXIMUM FINITE EIGENVALUE AND ASSOCIATED EIGENVECTOR

mind = find(isfinite(ssr), 1);
while abs(ssi(mind)) == Inf
    mind = mind+1;
end
smax = s(si(mind));
Vmax = V(:, si(mind));

