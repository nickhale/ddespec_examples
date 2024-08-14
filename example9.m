% EXAMPLE 9 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

%% Solve

% Problem parameters:
n = 12; dom = [0 1]; y0 = 1;
sol = sin(t);
% Constant operators:
[t, ~, w] = chebpts(n, dom);
D = Diffmat(t, w);
% Initial guess:
y = 1+0*t;

% Newton iteration:
disp('    iter      norm(F)   norm(dy)')
for j = 1:10
    % Linear operators
    E = Barymat(y, t, w);
    F = D*y + E*y;
    J = D + E + diag(E*D*y);
    % Boundary conditions
    J(1,:) = [1, zeros(1,n-1)]; F(1) = y(1) - y0;  
    % Newton update:
    dy = J\F;
    y = y - dy;
    ndy = norm(dy);
    disp([j-1 norm(F, inf), ndy])
    if ( ndy < 1e-10), break, end
end

% Plot the solution
nexttile
plot(t, y, '-k', t, y, '.', 'markersize', 25), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)


%% Convergence

err = []; sol = []; CN = [];
nn = [30 1:20];
for n = nn 
    % Problem parameters:
    dom = [0 1]; y0 = 1;
    % Constant operators:
    [t,~,w] = chebpts(n, dom);
    D = Diffmat(t,w);
    % Initial guess:
    y = 1+0*t;
    % Newton iteration:
    for j = 1:10
        % Linear operators
        E = Barymat(y, t, w);
        F = D*y + E*y;
        J = D + E + diag(E*D*y);
        % Boundary conditions
        J(1,:) = [1, zeros(1,n-1)]; F(1) = y(1) - y0;  
        % Newton update:
        if ( n == 2 ), J = J + 1e-12*randn(2); end % Fix for n = 2.
        dy = J\F;
        y = y - dy;
        ndy = norm(dy);
        if ( ndy < 1e-10), break, end
    end
    if ( isempty(sol) )
        sol = chebfun(y, dom);
    end
    err(n) = norm(y - sol(t), inf);
    CN(n) = cond(J);
end

nexttile
nn(1) = []; semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example9

%% Chebfun solution:

N = chebop(@(t,y) diff(y) + y(y), [0 1]);
N.lbc = 1;
y = N\0;
figure(3)
plot(y, '-k'); hold on, plot(y, '.', 'markersize', 25), hold off, shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)
