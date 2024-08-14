% EXAMPLE 5 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024


%% Determine the exact solution:

% Explicit expression for 0 < t < sqrt(3)/2
dom = [0 .5 sqrt(3)/2 1];
sol = chebfun({@(t) exp(-t), @(t) 1/2*exp(-t)*(exp(1/2)*sqrt(pi)*erf(1/2 - t) + 2)}, dom(1:end-1));
t = chebfun('x', [sqrt(3)/2 1]);
% Integrate to find solution on sqrt(3)/2 < t < 1
f = sum(exp(-t*(t - 1))*(exp(1/4) - (pi^(1/2)*exp(3/4)*erf(t^2 - 3/4))/2), 3^(1/2)/2, t);
f = - exp(-t)*((pi^(1/2)*exp(1/2)*erf(3^(1/2)/2 - 1/2))/2 - 1) - exp(-t)*f;
sol = join(sol, f);
figure(1)
plot(sol, '-'), hold on, shg
set(gca, 'fontsize', 16)
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');

%% Chebfun solution:

dom = [0 .5 sqrt(3)/2 1];
N = chebop(@(t,u) diff(u) + u + u(max(t^2-1/4, 0)), dom);
N.lbc = 1;
t = chebfun('t', dom);
f = t < 1/2; f(0.5) = 0;
y = N\f, plot(y, ':'), shg
hold off
err = norm(y - sol, inf)

%% Solve

% Problem parameters
x0 = 1; n = 12;

% Left discretisation        Middle disc                  Right discretisation
domL = [0, 1/2];             domM = [1/2, sqrt(3)/2];     domR = [sqrt(3)/2 1];
[tL,~,wL] = chebpts(n,domL); [tM,~,wM] = chebpts(n,domM); [tR,~,wR] = chebpts(n,domR);
t = [tL ; tM ; tR];
DL = Diffmat(tL, wL);        DM = Diffmat(tM,wM);         DR = Diffmat(tR,wR);
I  = eye(n);                 Z  = zeros(n,n);

% Boundary and continuity conditions:
z = zeros(1,n-1);
B = [1, z,  z, 0, z, 0];
C = [z, 1, -1, z, z, 0 ;
     z, 0, z, 1, -1, z ];

tauL = tL.^2 - 1/4;          tauM = tM.^2 - 1/4;          tauR = tR.^2 - 1/4;  
PL = Barymat(tauL,tL,wL);    PM = Barymat(tauM,tL,wL);    PR = Barymat(tauR,tM,wM);
PL(tauL<=0,:) = 0;

% Assemble:
AL  = [DL + I, Z, Z];        AM = [PM, DM + I, Z];        AR = [Z, PR, DR + I];
rhs = 0*t; rhs(1) = 1;
A   = [B  ; AL(2:n,:) ; C(1,:) ; AM(2:n,:) ; C(2,:) ; AR(2:n,:)]; 

% Solve and plot:
y = A\rhs;

% Plotting
close all, figure(1)
set(gcf, 'position', [2158 385 1680 420])
tiles = tiledlayout(1,3,'TileSpacing','loose');

nexttile(2)
plot(sol, '-k');
xlabel('$t$', 'interp', 'latex'); 
ylabel('$y$', 'interp', 'latex');

U = reshape(y, n ,3);
T = reshape(t, n, 3);
hold on
plot(T, U, '.', 'markersize', 25), shg
hold off
set(gca, 'fontsize', 16)

nexttile(1)
spy(A)
set(gca, 'fontsize', 16)

err = norm(y - sol(t))

%%

nn = 1:20; err = []; CN = [];
for n = nn
    % Left discretisation         Middle discretisation         Right discretisation
    domL = [0, 1/2];              domM = [1/2, sqrt(3)/2];      domR = [sqrt(3)/2 1];
    [tL,~,wL] = chebpts(n, domL); [tM,~,wM] = chebpts(n, domM); [tR,~,wR] = chebpts(n, domR);
    t = [tL ; tM ; tR];
    DL = Diffmat(tL, wL);         DM = Diffmat(tM, wM);         DR = Diffmat(tR, wR);
    I  = eye(n);                  Z  = zeros(n,n);
    % Boundary and continuity conditions:
    z = zeros(1,n-1);
    B = [1, z,  z, 0, z, 0];
    C = [z, 1, -1, z, z, 0 ;
         z, 0, z, 1, -1, z ];
    tauL = tL.^2 - 1/4;           tauM = tM.^2 - 1/4;           tauR = tR.^2 - 1/4;  
    PL = Barymat(tauL,tL,wL);     PM = Barymat(tauM,tL,wL);     PR = Barymat(tauR,tM,wM);
    PL(tauL<=0,:) = 0;
    % Assemble:
    AL = [DL + I, Z, Z];          AM = [PM, DM + I, Z];         AR = [Z, PR, DR + I];
    rhs = 0*t; rhs(1) = 1;
    A   = [B  ; AL(2:n,:) ; C(1,:) ; AM(2:n,:) ; C(2,:) ; AR(2:n,:)]; 
    % Solve and plot:
    y = A\rhs;
    err(n) = norm(y - sol(t));
    CN(n) = cond(A);
end

n3 = nexttile(3);
semilogy(nn, err(nn));
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
ylim([1e-16 1e0])
set(gca, 'ytick', [1e-15 1e-10 1e-5 1e0])
set(gca, 'fontsize', 16)

% print -depsc2 ../paper/figures/example5



