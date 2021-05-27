function f = help_plotV(V, l, lselect, N, z)

%   Plot the eigenmodes u, psi, xi, whether marginal or not, for a
%   particular l, given by lselect.

%% FIND l INDEX

lind = find(l==lselect);

%% EXTRACT V

vind = V(:, lind);

%% EXTRACT u, w, p, theta

uind = vind(1:N-2);                   uind = [0; uind; 0];    
wind = vind(N-1:(2*N)-4);             wind = [0; wind; 0];  
pind = vind((2*N)-3:(3*N)-4); 
tind = vind((3*N)-3:(4*N)-6);         tind = [0; tind; 0];

%% PLOT MODES

f = figure;

subplot(1,4,1)
hold on
plot(real(uind), z, 'b-o', 'linewidth', 2, 'markersize', 3)
plot(imag(uind), z, 'c-o', 'linewidth', 2, 'markersize', 3)
xlabel('$\hat{u}(z;k)$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
set(gca, 'fontsize', 20)
legend('Real', 'Imag', 'Location', 'Best')
axis tight
%xlim([-1, 1])
grid on
box on

subplot(1,4,2)
hold on
plot(real(wind), z, 'b-o', 'linewidth', 2, 'markersize', 3)
plot(imag(wind), z, 'c-o', 'linewidth', 2, 'markersize', 3)
xlabel('$\hat{w}(z;k)$', 'interpreter', 'latex')
set(gca, 'fontsize', 20)
axis tight
%xlim([-0.022, 0.015])
grid on
box on

subplot(1,4,3)
hold on
plot(real(pind), z, 'b-o', 'linewidth', 2, 'markersize', 3)
plot(imag(pind), z, 'c-o', 'linewidth', 2, 'markersize', 3)
xlabel('$\hat{p}(z;k)$', 'interpreter', 'latex')
set(gca, 'fontsize', 20)
axis tight
%xlim([-0.1, 1.1])
grid on
box on

subplot(1,4,4)
hold on
plot(real(tind), z, 'b-o', 'linewidth', 2, 'markersize', 3)
plot(imag(tind), z, 'c-o', 'linewidth', 2, 'markersize', 3)
xlabel('$\hat{\theta}(z;k)$', 'interpreter', 'latex')
set(gca, 'fontsize', 20)
axis tight
%xlim([-0.1, 1.1])
grid on
box on

end