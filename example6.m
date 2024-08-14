% EXAMPLE 6 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

%% Solve

% Problem parameters:
n = 12; dom = [0 1]; y0 = 0;
% Constant operators:
[t,~,w] = chebpts(n, dom);
D = Diffmat(t, w);
% Initial guess:
y = t;
f = cos(t)+sin(sin(t));
% Newton iteration:
disp('    iter      norm(F)   norm(dy)')
for j = 1:10
    % Linear operators
    P = Barymat(y, t, w);
    F = D*y + P*y - f;
    J = D + P + diag(P*D*y);
    % Boundary conditionsed
    J(1,:) = [1, zeros(1,n-1)]; F(1) = y(1) - y0;  
    % Newton update:
    dy = J\F;
    y = y - dy;
    ndy = norm(dy);
    disp([j-1 norm(F, inf), ndy])
    if ( ndy < 1e-10), break, end
end
% Plot the solution
sol = sin(t);

nexttile
plot(t, sol, '-k', t, y, '.', 'markersize', 30), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

err = norm(y - sol, inf)

%% Convergence

nn = 1:20; err = []; CN = [];
for n = nn 

    % Problem parameters:
    dom = [0 1]; y0 = 0;
    % Constant operators:
    [t,~,w] = chebpts(n, dom);
    D = Diffmat(t, w);
    % Initial guess:
    y = t;
    f = cos(t)+sin(sin(t));
    % Newton iteration:
    for j = 1:10
        % Linear operators
        P = Barymat(y, t, w);
        F = D*y + P*y - f;
        J = D + P + diag(D*P*y);
        % Boundary conditions
        J(1,:) = [1, zeros(1,n-1)]; F(1) = y(1) - y0;  
        % Newton update:
        dy = J\F;
        y = y - dy;
        ndy = norm(dy);
        if ( ndy < 1e-13), break, end
    end

    sol = sin(t);
    err(n) = norm(y - sol, inf)./norm(sol, inf);
    CN(n) = cond(J);
end

nexttile
semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('relative error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example6

%% Chebfun solution:

dom = [0 1];
t = chebfun('t', dom);
f = cos(t)+sin(sin(t));
N = chebop(@(t,y) diff(y) + y(y) - f, dom);
N.lbc = 0;
N.init = t;
y = N\0;
sol = sin(t);

figure(2)
plot(y, '-', sol, ':'), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)
err = norm(y - sol, inf)
