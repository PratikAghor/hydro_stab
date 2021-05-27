%=======================================================================%
%   Study growth rates of various streamwise modes at various Rayleigh  %
%   numbers.                                                            %
%=======================================================================%
clear
close all
clc

%% PARAMETER SETUP

N  = 64; 
k  = 2*pi*(0.1:0.1:3);
R = 3000;
Pr = 1;

%% GENERATE CGL GRID

[~, z] = cheb(N-1);

%% CREATE SOLUTION ARRAYS

Smax = zeros(size(k));
Vmax = zeros(4*N-6, length(k));

%% SOLVE EIGENVALUE PROBLEMS

for i = 1:length(k)
   [Smax(i), Vmax(:, i), ~] = get_maxeig(N, k(i), R, Pr); 
end

%% DECOMPOSE EIGENVALUES

Sr = real(Smax);
Si = imag(Smax);

%% PLOT GROWTH-RATE CURVE

figure(1)
hold on

plot(k, zeros(size(k)), 'g--', 'linewidth', 3)
plot(k, Sr, 'k-o', 'linewidth', 3, 'markersize', 4)
xlabel('$k$', 'interpreter', 'latex')
ylabel('$Re(\sigma)$', 'interpreter', 'latex')
set(gca, 'fontsize', 20)
axis tight
axis square
%xlim([0, 80])
%ylim([-7, 2])
grid on
box on

%% PLOT DOMINANT MODES

za = 0.5*(z+1);
f1 = help_plotV(Vmax, k, 2*pi, N, za);
%f2 = help_plotV(Vmax, l, 6*pi, N, z);
