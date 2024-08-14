% EXAMPLE 7 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

%% Solve

K = @(s,t) exp((t-s));
n = 14; dom = [0 1]; a = 1/2;
t = chebpts(n, dom); I = eye(n); z = zeros(1,n-1);
D = diffmat(n, dom);            % Differentiation matrix (www.chebfun.org)
V = K(t,t').*cumsummat(n, dom); % Volterra operator (www.chebfun.org)
A = D + a*I - V ; rhs = 0*t;    % Discrete operator and RHS
A(1,:) = [1, z]; rhs(1) = 1;    % Boundary conditions
y = A\rhs;
b = sqrt(a^2-2*a+5)/2;
sol = @(t) exp(-(a+1)/2*t).*(cosh(b*t) + .5*(1-a)/b*sinh(b*t));

nexttile
plot(t, sol(t), '-k', t, y, '.', 'markersize', 25), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

err = norm(y - sol(t), inf)

%% Convergence

nn = 1:20; err = []; CN = [];
for n = nn
    t = chebpts(n, dom); I = eye(n); z = zeros(1,n-1);
    D = diffmat(n, dom);            % Differentiation matrix (www.chebfun.org)
    V = K(t,t').*cumsummat(n, dom); % Volterra operator (www.chebfun.org)
    A = D + a*I - V ; rhs = 0*t;    % Discrete operator and RHS
    A(1,:) = [1, z]; rhs(1) = 1;    % Boundary conditions
    y = A\rhs;
    
    err(n) = norm(y - sol(t), inf);
    CN(n) = cond(A);
end

nexttile
semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example7

%% Chebfun solution:
dom = [0 1];
N = chebop(@(t,u) diff(u) + a*u - volt(K, u), dom);
N.lbc = 1;
y = N\0

figure(2)
plot(t, sol(t), '-', y, ':'), shg
err = norm(sol(t) - y(t), inf)