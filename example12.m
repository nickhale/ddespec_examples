% EXAMPLE 12 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024


close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

lam = .5; f = @(t) lam*sin(t);          % Parameters
dom = [0 pi];
n = 50; [t, ~, w] = chebpts(n, dom); % Discretisation
D = Diffmat(t,w); P = Barymat(f(t), t, w); I = eye(n); 
A = P - lam*I;                          % Linear operator

A = [D(1,:) ; A(2:end,:)];          % Add Neumann condition
rhs = [1 ; zeros(n-1,1)];
u = A\rhs;                          % Solve

w = t; 
dw = inf;  k = 0;                   % Initialise
while ( dw > eps)                   % iterate
    k = k+1;
    fw = f(w);
    dw = norm(lam*w-fw, inf)/lam^k;
    w = fw;
end
w = w./lam^k;                       % Normalise

nexttile
plot(t, u, t, w, ':')               % Plot
xlim([0, pi])
xlabel('$t$', 'interp', 'latex'); ylabel('$u$', 'interp', 'latex'), grid on
set(gca, 'fontsize', 16)
legend('interpolation', 'iteration')
% title('lambda = 2')
set(gca, 'fontsize', 16)

%%

nn = 1:1:50;
err = zeros(max(nn),1); CN = [];

for n = nn
    lam = .5; f = @(t) lam*sin(t); fp = @(t) lam*cos(t);
    [t, ~, w] = chebpts(n, dom);
    D = Diffmat(t, w); P = Barymat(f(t), t, w); I = eye(n); 
    A = P - diag(fp(0))*I;              % Linear operator

    A = [D(1,:) ; A(2:end,:)];          % Add Neumann condition
    rhs = [1 ; zeros(n-1,1)];
    u = A\rhs;                          % Solve

    w = t; 
    dw = inf;  k = 0;                   % Initialise
    while ( dw > eps)                   % iterate
        k = k+1;
        fw = f(w);
        dw = norm(lam*w-fw, inf)/lam^k;
        w = fw;
    end
    w = w./lam^k;                       % Normalise

    err(n) = norm(u - w, inf)./norm(u, inf);
    CN(n) = cond(A);
    xlim([0, pi])
end

nexttile
semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error', 'ver', 'mid'), grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example12


