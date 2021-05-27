clear
close all
clc

%% PRELIMINARY PARAMETERS

N = 64;
k = 1; 
Pr = 1;

%% GENERATE CGL GRID

[D1, z] = cheb(N-1);

%% CONSTRUCT LINEAR OPERATOR

[A, B, z] = get_mstab_linop(N, k, Pr);

%% SOLVE GENERALIZED EIGENVALUE PROBLEM

[V, D] = eig(A, B);

%% THROW AWAY NEGATIVE REAL PART ILAs

s = diag(D);

nind = find(real(s)<=0);
s(nind) = [];
V(:, nind) = [];

%% SORT EIGENVALUES (IN THE ASCENDING ORDER)

[ssr, si] = sort(real(s),'ascend');
ssi = imag(s(si));

%% MINIMUM NONZERO EIGENVALUE AND ASSOCIATED EIGENVECTOR

mind = 1;
% while abs(ssi(mind-1)) == Inf
%     mind = mind+1;
% end
smax = s(si(mind));
Vmax = V(:, si(mind));