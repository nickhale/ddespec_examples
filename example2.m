% EXAMPLE 2 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

%% Solution

dom = [0 1];                 % Solution interval
n = 12;                      % Discretisation size
[t,~,w] = chebpts(n, [0 1]); % Chebyshev grid
D = Diffmat(t, w);           % Differentiation matrix
tau = t/2;                   % Pantograph-type delay
P = Barymat(tau, t, w);      % Barycentric interpolation matrix
I = eye(n);                  % Identity matrix
z = zeros(1,n-1);            % Zero vector
% Discretised operator and RHS:
A = D + I + P; rhs = exp(-t/2); 
% Enforce boundary condition via boundary bordering:
A(1,:) = [1, z]; rhs(1) = 1;
y = A\rhs;                   % Solve
sol = exp(-t);                % Exact solution
err = norm(y - sol)          % Error

nexttile
plot(t, sol, '-k', t, y, '.', 'markersize', 30), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

%% Convergence

nn = 1:20; err = []; CN = [];
for n = nn
    [t,~,w] = chebpts(n, [0 1]); % Chebyshev grid
    D = Diffmat(t, w);           % Differentiation matrix
    I = eye(n);                  % Identity matrix
    tau = t/2;                   % Pantograph-type delay
    P = barymat(tau, t, w);       % Barycentric interpolation matrix    
    z = zeros(1,n-1);            % Zero vector
    % Discretised operator and RHS:
    A = D + I + P; rhs = exp(-t/2);
    % Enforce boundary condition via boundary bordering:
    A(1,:) = [1, z]; rhs(1) = 1;
    y = A\rhs;                   % Solve
    sol = exp(-t);               % Exact solution
    err(n) = norm(y - sol);      % Error
    CN(n) = cond(A);
end

nexttile
semilogy(nn, err(nn))
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
% hold on, plot(nn, eps*nn.^2.5./(nn(end)^2.5/CN(nn(end))), ':', 'color', .6*[1 0 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
set(gca, 'ytick', 10.^(-15:5:0))
set(gca, 'fontsize', 16)

% print -depsc2 ../paper/figures/example2

%% Chebfun implementation

N = chebop(@(t,y) diff(y) + y(t) + y(t/2) - exp(-t/2), [0 1]);
N.lbc = 1;
y = N\0;
err = norm(y - chebfun(@(t) exp(-t), [0 1]), inf)
figure(2)
plot(y), xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)