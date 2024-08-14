% EXAMPLE 3 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

%% Solve

p = 0.5;

% Left discretisation           Right discretisation
domL = [0, p];                  domR = [p, 2*p];
n = 10; nL = n;                 nR = n+1;
[tL,~,wL] = chebpts(nL, domL);  [tR,~,wR] = chebpts(nR, domR); 
t = [tL ; tR];
DL = Diffmat(tL, wL);           DR = Diffmat(tR,wR);
IL  = eye(nL);                  IR  = eye(nR);             
ZL  = zeros(nL,nR);             ZR  = zeros(nR,nL); 
rhsL = 0*tL;                    rhsR = 0*tR;
tauR = tR - 1/2;
PR = Barymat(tauR, tL, wL);
AL = [DL + IL, ZL];             AR = [PR, DR + IR];

% Boundary and continuity conditions:
zL = zeros(1,nL-1);             zR = zeros(1,nR-1);
B = [1, zL,  zR, 0];
C = [zL, 1, -1, zR];
tauR = tR - 1/2;

% Assemble:
A   = [B  ; AL(2:nL,:) ; C ; AR(2:nR,:)]; 
rhs = [1  ; rhsL(2:nL) ; 0 ; rhsR(2:nR)];

% Solve and plot:
y = A\rhs;
sol = [exp(-tL) ; exp(-tR+1/2).*(1/2 - tR + exp(-1/2))];
err = norm(y - sol, inf)

nexttile
plot(t, sol, '-k', t, y, '.', 'markersize', 30), shg
ylim([0 1]), xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

%% Convergence

nn = 1:20; err = []; CN = [];
for n = nn
    % Left discretisation           Right discretisation
    domL = [0, p];                  domR = [p, 2*p];
    nL = n;                         nR = n+1;
    [tL,~,wL] = chebpts(nL, domL);  [tR,~,wR] = chebpts(nR, domR); 
    t = [tL ; tR];
    DL = Diffmat(tL, wL);           DR = Diffmat(tR, wR);
    IL  = eye(nL);                  IR  = eye(nR);             
    ZL  = zeros(nL,nR);             ZR  = zeros(nR,nL); 
    rhsL = 0*tL;                    rhsR = 0*tR;
    tauR = tR - 1/2;
    PR = Barymat(tauR, tL, wL);
    AL = [DL + IL, ZL];             AR = [PR, DR + IR];

    % Boundary and continuity conditions:
    zL = zeros(1,nL-1);             zR = zeros(1,nR-1);
    B = [1, zL,  zR, 0];
    C = [zL, 1, -1, zR];
    tauR = tR - 1/2;

    % Assemble:
    A   = [B  ; AL(2:nL,:) ; C ; AR(2:nR,:)]; 
    rhs = [1  ; rhsL(2:nL) ; 0 ; rhsR(2:nR)];

    % Solve:
    y = A\rhs;
    sol = [exp(-tL) ; exp(1/2-tR).*(1/2 - tR + exp(-1/2))];
    err(n) = norm(y - sol);
    CN(n) = cond(A); 
end

nexttile
semilogy(nn, err(nn))
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

%% Single discretisation for comparison

for n = nn
    
    dom = [0, 1];
    [t,~,w] = chebpts(n, dom);
    D = Diffmat(t, w);
    tau = t - p;
    P = Barymat(tau, t, w);
    P(tau < 0,:) = 0;
    I  = eye(n);
    A = D + I + P;
    rhs = 0*t;
    % Boundary and continuity conditions:
    z = zeros(1,n-1);
    B = [1, z];
    % Assemble:
    A   = [B  ; A(2:n,:)]; 
    rhs = [1  ; rhs(2:n)];
    % Solve:
    y = A\rhs;
    solL = exp(-t);
    solR = exp(p-t).*(p - t + exp(-p));
    sol = 0*y; idx = t < p;
    sol(idx) = solL(idx); sol(~idx) = solR(~idx);
    err(n) = norm(y - sol);

end

col = colororder;
hold on
loglog(nn, err(nn), '-.', 'color', col(2,:)), shg
hold off
set(gca, 'fontsize', 16)

% print -depsc2 ../paper/figures/example3

%% Chebfun implementation

N = chebop(@(t,y) diff(y) + y(t) + y(max(t-.5,0)), [0 1]);
N.lbc = 1;
rhs = chebfun({1,0}, [0 .5 1]); % Chebfun assumes y(t<0) = y(0).
y = N\rhs;
tL = chebfun('t', domL); tR = chebfun('t', domR);
sol = join(exp(-tL), exp(-tR+1/2).*(1/2 - tR + exp(-1/2)));
err = norm(y - sol, inf)

figure(2)
plot(y, '-', sol, ':'), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

