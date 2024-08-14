% EXAMPLE 1 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

%% Solution 

set(0, 'defaultlinelinewidth', 3)
set(0, 'DefaultAxesFontSize', 16);

dom = [0 1];                 % Solution interval
n = 12;                      % Discretisation size
[t,~,w] = chebpts(n, [0 1]); % Chebyshev grid
D = Diffmat(t, w);           % Differentiation matrix
I = eye(n);                  % Identity matrix
z = zeros(1, n-1);           % Zero vector
A = D + I; rhs = 0*t;        % Discretised operator and RHS
% Enforce boundary condition via boundary bordering:
A(1,:) = [1, z]; rhs(1) = 1;
y = A\rhs;                   % Solve
sol = exp(-t);               % Exact solution
err = norm(y - sol)          % Error

nexttile
plot(t, sol, '-k', t, y, '.', 'markersize', 30)
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

%% Convergence

nn = 1:20; err = []; CN = [];
for n = nn
    [t,~,w] = chebpts(n, [0 1]); % Chebyshev grid
    D = Diffmat(t, w);           % Differentiation matrix
    I = eye(n);                  % Identity matrix
    z = zeros(1, n-1);           % Zero vector
    A = D + I; rhs = 0*t;        % Discretised operator and RHS
    % Enforce boundary condition via boundary bordering:
    A(1,:) = [1, z]; rhs(1) = 1;
    y = A\rhs;                   % Solve
    sol = exp(-t);               % Exact solution
    err(n) = norm(y - sol);      % Error
    CN(n) = cond(A);              % Condition number
end

nexttile
semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
% hold on, plot(nn, eps*nn.^2.5./(nn(end)^2.5/CN(nn(end))), ':', 'color', .6*[1 0 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example1

%% Chebfun implementation

N = chebop(@(t,y) diff(y) + y(t), [0 1]);
N.lbc = 1;
y = N\0;
err = norm(y - chebfun(@(t) exp(-t), [0 1]), inf)

figure(2)
plot(y), xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)