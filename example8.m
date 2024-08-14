% EXAMPLE 8 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

%% Solution

dom = [0 1];                 % Solution interval
n = 16;                      % Discretisation size
[t,~,w] = chebpts(n, dom);   % Chebyshev grid
D = Diffmat(t,w);            % Differentiation matrix
tau = 1-t.^2;                % Discrete delay
P = Barymat(tau, t, w);      % Barycentric interpolation matrix
I = eye(n);                  % Identity matrix
z = zeros(1,n-1);            % Zero vector
% Discretised operator and RHS:
A = D + I + P; rhs = exp(t.^2-1); 
% Enforce boundary condition via boundary bordering:
A(1,:) = [1, z]; rhs(1) = 1;
y = A\rhs;                   % Solve
sol = exp(-t);                % Exact solution
err = norm(y - sol)          % Error

nexttile
plot(t, sol, '-k', t, y, '.', 'markersize', 30), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

sol = chebfun(y, dom);

%% Convergence

nn = 1:20; err = []; CN = [];
for n = nn
    [t,~,w] = chebpts(n, [0 1]); % Chebyshev grid
    D = Diffmat(t, w);           % Differentiation matrix
    I = eye(n);                  % Identity matrix
    tau = 1-t;                   % Discrete delay
    P = Barymat(tau, t, w);      % Barycentric interpolation matrix    
    z = zeros(1,n-1);            % Zero vector
    % Discretised operator and RHS:
    A = D + I + P; rhs = exp(t-1);
    % Enforce boundary condition via boundary bordering:
    A(1,:) = [1, z]; rhs(1) = 1;
    y = A\rhs;                   % Solve
    err(n) = norm(y - sol(t));   % Error
    CN(n) = cond(A);
end

nexttile
semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example8

%% Chebfun implementation

t = chebfun(t, [0 1]);
N = chebop(@(t,y) diff(y) + y(t) + y(1-t) - exp(t-1), [0 1]);
N.lbc = 1;
y = N\0;
err = norm(y - chebfun(@(t) exp(-t), [0 1]), inf)
figure(2)
plot(y), xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)