function [smin, Vmin, z] = get_mstab_mineig(N, k, Pr)

% Finds the minimal 1/La for some specified y-wavenumber and its
% associated eigenvector. 

%% CONSTRUCT LINEAR OPERATOR

[A, B, z] = get_mstab_linop(N, k, Pr);

%% SOLVE GENERALIZED EIGENVALUE PROBLEM

[V, D] = eig(A, B);

%% THROW AWAY NEGATIVE REAL PART ILAs

s = diag(D);

nind = find(real(s)<=0);
s(nind) = [];
V(:, nind) = [];

%% SORT EIGENVALUES IN ASCENDING ORDER

[ssr, si] = sort(real(s),'ascend');

%% MINIMUM NONZERO EIGENVALUE AND ASSOCIATED EIGENVECTOR

mind = 1;
smin = s(si(mind));
Vmin = V(:, si(mind));

end