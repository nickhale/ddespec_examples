% EXAMPLE 9 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1680 420])
tiles = tiledlayout(1,3,'TileSpacing','loose');

%% Solution

T = 2*pi;
n = 32;
t = linspace(0, T, n+1)'; t(end) = [];
D1 = diffmat(n, 1, 'periodic', [0 T]);
D2 = diffmat(n, 2, 'periodic', [0 T]);

P1 = barymattrig(t-pi/sqrt(2),t);
P2 = barymattrig(t-pi/2,t);
A = D2+diag(sin(t))*P1*D1...
      +diag(cos(t))*P2;
u = A\ones(n,1);

cols = colororder;
nexttile
h = plot(t, u);
xlim([0 T])
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

uu = @(z) barymattrig(z,t)*u;
c = fftshift(fft(u))/(2*n);
nexttile(3)
semilogy((-n/2:n/2-1)+mod(n,2)/2, abs(c), '.','markersize', 25)
xlabel('$k$', 'interp', 'latex'); ylabel('$|c_k|$', 'interp', 'latex');
ylim([1e-16 1e1]), shg, grid on
xlim([-16 16])
set(gca, 'ytick', [1e-15 1e-10 1e-5 1e0])
set(gca, 'fontsize', 16)


%% Convergence

nn = 3:25;
for n = nn 
    t = linspace(0, T, n+1)'; t(end) = [];
    D1 = diffmat(n, 1, 'periodic', [0 T]);
    D2 = diffmat(n, 2, 'periodic', [0 T]);
    
    P1 = barymattrig(t-pi/sqrt(2),t);
    P2 = barymattrig(t-pi/2,t);
    A = D2+diag(sin(2*pi*t/T))*P1*D1...
          +diag(cos(2*pi*t/T))*P2;
    u = A\ones(n,1);

    err(n) = norm(u - uu(t), inf)./norm(u, inf);
    CN(n) = cond(A);

end

nexttile(2)
semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
ylim([1e-16 1e1]), shg, grid on
set(gca, 'ytick', [1e-15 1e-10 1e-5 1e0])
xlim([0 25])
xlabel('$n$', 'interp', 'latex'); ylabel('relative error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example10

return

%% Chebyshev version

n = 42;
[t,~,w] = chebpts(n, [0, T]);
I = eye(n);
D1 = Diffmat(t, w);
D2 = D1^2;

P1 = barymat(mod(t-pi/sqrt(2), T),t);
P2 = barymat(mod(t-pi/2, T),t);
A = D2+diag(sin(t))*P1*D1...
      +diag(cos(t))*P2;
rhs = ones(n,1);
A(1,:) = I(1,:) - I(end,:); rhs(1) = 0;
A(end,:) = D1(1,:) - D1(end,:); rhs(end) = 0;
u = A\rhs;
uu = chebfun(u, [0 T])
figure(2), plot(t, u, 'k', t, u, '.r', 'markersize', 25), shg
hold on, plot(t, uu(t),':'); hold off
set(gca, 'fontsize', 16)

%% Convergence (Chebyshev version)

nn = 5:50;
for n = nn 
    
    [t,~,w] = chebpts(n, [0, T]);
    I = eye(n);
    D1 = Diffmat(t, w);
    D2 = D1^2;
    
    P1 = barymat(mod(t-pi/sqrt(2), T),t);
    P2 = barymat(mod(t-pi/2, T),t);
    A = D2+diag(sin(t*2*pi/T))*P1*D1...
          +diag(cos(t*2*pi/T))*P2;
    rhs = ones(n,1);
    A(1,:) = I(1,:) - I(end,:); rhs(1) = 0;
    A(end,:) = D1(1,:) - D1(end,:); rhs(end) = 0;
    u = A\rhs;

    err(n) = norm(u - uu(t), inf)./norm(u, inf);
    CN(n) = cond(A);

end

figure(3)
semilogy(nn, err(nn)), shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
ylim([1e-16 1e1]), shg, grid on
set(gca, 'ytick', [1e-15 1e-10 1e-5 1e0])
% xlim([0 25])
xlabel('$n$', 'interp', 'latex'); ylabel('relative error'); grid on
set(gca, 'fontsize', 16)
title('Chebyshev version')

%% Chebfun implementation I (Chebyshev discretiation)

N = chebop(@(t,u) diff(u,2) + sin(2*pi*t/T)*feval(diff(u), mod(t-pi/sqrt(2),T)) + cos(2*pi*t/T)*u(mod(t-pi/2,T)), [0 T]);
N.bc = @(t,u) [u(0) - u(T) ; deriv(u,0) - deriv(u,T)];
figure(2)
hold on, plot(N\1, ':g'), hold off

return

%% Chebfun implementation II (Fourier discretisation)
% Note: requires feature-periodic-DDE branch

N = chebop(@(t,u) diff(u,2) + sin(2*pi*t/T)*feval(diff(u), t-pi/sqrt(2)) + cos(2*pi*t/T)*u(t-pi/2), [0 T]);
N.bc = 'periodic';
figure(2)
hold on, plot(N\1, ':g'), hold off

%%

function P = barymattrig(z, t)
    T = t(end)-2*t(1)+t(2);
    l = length(t);
    if ( mod(l,2) )
        P = csc((z-t.')*pi/T);
    else
        P = cot((z-t.')*pi/T);
    end
    P(:,2:2:end) = -P(:,2:2:end);
    P = P./sum(P,2);
    P(isnan(P)) = 1;
end