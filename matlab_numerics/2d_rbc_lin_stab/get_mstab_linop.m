function [A, B, z] = get_mstab_linop(N, k, Pr)

%   Given an x-wavenumber k = 2*pi*K, where K is an integer, the
%   number of CGL gridpoints N, and an appropriate Rayleigh Number R and a
%   Prandtl Number Pr, this code produces the linear operators A and B to
%   solve the marginal stability problem.

%% GET CGL GRID-POINTS AND CHEBYSHEV DIFF'N MATRICES

[D1, z] = cheb(N-1);  D2 = D1^2;  %   Note, z is from -1 to 1.

%% DEFINE SOLUTION ARRAYS

A = zeros(4*N, 4*N);    B = zeros(4*N, 4*N);

%% DEFINE CONVENIENCE VARIABLES

Zr = zeros(1, N);

I  = eye(N);
Lp = 4*D2 - I*(k^2);
Zm = zeros(N, N);

%% DEFINE A and B, PRELIMINARILY

A = [Pr*Lp, Zm, -1i*k*I, Zm;...
     Zm, Lp, -(2/Pr)*D1, Zm;...
     Zm, Zm, (1/Pr)*Lp, Zm;...
     Zm, I, Zm, Lp];
 
B = [Zm, Zm, Zm, Zm;...
     Zm, Zm, Zm, -I;...
     Zm, Zm, Zm, 2*D1;...
     Zm, Zm, Zm, Zm];

%% APPLY BCS

% Homogeneous Neumann BCs on p

A(2*N+1, 1:N) = Zr;                 A(3*N, 1:N) = Zr;
A(2*N+1, N+1:2*N) = Zr;             A(3*N, N+1:2*N) = Zr;
A(2*N+1, 2*N+1:3*N) = D1(1, :);     A(3*N, 2*N+1:3*N) = D1(N, :);
A(2*N+1, 3*N+1:4*N) = Zr;           A(3*N, 3*N+1:4*N) = Zr;
B(2*N+1, 3*N+1:4*N) = Zr;           B(3*N, 3*N+1:4*N) = Zr;

% Homogeneous Dirichlet BCs on u

A(1, :) = [];         A(:, 1) = [];
A(N-1, :) = [];       A(:, N-1) = [];
B(1, :) = [];         B(:, 1) = [];
B(N-1, :) = [];       B(:, N-1) = [];

% Homogeneous Dirichlet BCs on w

A(N-1, :) = [];         A(:, N-1) = [];
A((2*N)-3, :) = [];     A(:, (2*N)-3) = [];
B(N-1, :) = [];         B(:, N-1) = [];
B((2*N)-3, :) = [];     B(:, (2*N)-3) = [];

% Homogeneous Dirichlet BCs on theta

A((3*N)-3, :) = [];     A(:, (3*N)-3) = [];
A((4*N)-5, :) = [];     A(:, (4*N)-5) = [];
B((3*N)-3, :) = [];     B(:, (3*N)-3) = [];
B((4*N)-5, :) = [];     B(:, (4*N)-5) = [];

end