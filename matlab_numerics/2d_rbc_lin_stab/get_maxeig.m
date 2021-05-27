function [smax, Vmax, z] = get_maxeig(N, k, R, Pr)

% Finds the maximal eigenvalue for some specified y-wavenumber, La, and its
% associated eigenvector. 

%% CONSTRUCT LINEAR OPERATOR

[A, B, z] = get_linop(N, k, R, Pr);

%% SOLVE GENERALIZED EIGENVALUE PROBLEM

[V, D] = eig(A, B);

%% SORT EIGENVALUES

s = diag(D);
[ssr, si] = sort((real(s)),'descend');
ssi = imag(s(si));

%% MAXIMUM FINITE EIGENVALUE AND ASSOCIATED EIGENVECTOR

mind = find(isfinite(ssr), 1);    % Find maximum finite real part.

while abs(ssi(mind)) == Inf       % Make sure the imaginary part of this 
    mind = mind+1;                % eigenvalue isn't infinite.
end

smax = s(si(mind));
Vmax = V(:, si(mind));

end