%=======================================================================%
%   Study marginal stability 1/La for various cross-wind modes.         %
%=======================================================================%
clear
close all
clc

%% PARAMETER SETUP

N  = 64; 
k  = linspace(1, 4, 100); 
Pr = 1;

%% GENERATE CGL GRID

[~, z] = cheb(N-1);

%% CREATE SOLUTION ARRAYS

Rmin = zeros(size(k));
Vmin = zeros(4*N-6, length(k));

%% SOLVE EIGENVALUE PROBLEMS

for i = 1:length(k)
   [Rmin(i), Vmin(:, i), ~] = get_mstab_mineig(N, k(i), Pr);
end

%% DECOMPOSE EIGENVALUES

Rr = real(Rmin);
Ri = imag(Rmin);

%% PLOT MARGINAL STABILITY CURVE

figure(1)
hold on

plot(k, Rr, 'r-', 'linewidth', 3, 'markersize', 3)
xlabel('$k$', 'interpreter', 'latex')
ylabel('$R$', 'interpreter', 'latex')
set(gca, 'fontsize', 20)
axis tight
axis square
grid on
box on

%% PLOT SPECIFIC MARGINALLY STABLE MODES

% za = 0.5*(z+1);
% f = help_plotV(Vmin, k, 4*pi, N, za);