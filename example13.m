% EXAMPLE 13 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

n = 19;
tol = 1e-12;
lam = 1.7;

T0 = 4; % Estimate of the period
u0 = 0.5; % Estimate for a point on the limit cycle

clc, close all, figure(1)
set(gcf, 'position', [2158 385 1680 420])
tiles = tiledlayout(1,3,'TileSpacing','loose');

for disc = {'cheb', 'trig'}
    disc = disc{:};

    switch disc
        case 'trig'
            nproxy = 40;
            nn = 5:33;
        case 'cheb'
            nproxy = 70;
            nn = 10:1:54;
    end
    
    uproxy = [];
    err = []; resnorm = [];
    
    for n = [nproxy nn]
    
        T = T0;
        dom = [0 3*T];
        dom2 = [1*T 3*T];
        
        % Solve using dde23 to get an approximate solution
        f = @(t,u,z) (lam - z).*u;
        sol = dde23(f, 1, u0, dom);
        
        % Approximate the solution over a single period
        uu = chebfun(@(x) deval(sol, x), dom2);
        [~, idx0] = min(uu, 'local');
        T = idx0(end)-idx0(end-1);
        idx = idx0(end-1);
        
        uu = restrict(uu, [idx, idx+T]); % Single period solution
        uu = newDomain(uu, [0 1]);       % Change of variables to [0,1]
    
        switch disc
            case 'trig'
                [s, w] = trigpts(n, [0 1]);
                D = diffmat(n, 1, 'periodic', [0 1]);
                P = Barymattrig(s-1/T,s);
            case 'cheb'
                [s, w] = chebpts(n, [0 1]);
                D = diffmat(n, 1, [0 1]);
                P = barymat(s-1/T,s);
        end
        I = eye(size(D));
        u = uu(s);
        
        % Newton iteration
        for ell = 1:10
            uT = [u ; T];

            switch disc
                case 'trig'
                    P = Barymattrig(s-1/T,s);
                case 'cheb'
                    P = barymat(mod(s-1/T, 1),s);
            end
            
            F = D*u - T*(lam-P*u).*u;                   
            Ju = D - T*(lam*I - diag(P*u) - diag(u)*P);
            JT = - (lam-P*u).*u + (1/T*P*D*u).*u;
            J = [Ju JT];
        
            % Boundary conditions
            switch disc
                case 'trig'
                    % Enforce u'(0) to remove translation invariance.
                    F(2*n+1) = D(1,:)*u;
                    J(2*n+1,:) = [D(1,:), 0];
                case 'cheb'
                    F(1) = (I(1,:)-I(end,:))*u;
                    J(1,:) = [I(1,:)-I(end,:), 0];
                    F(2*n+1) = D(1,:)*u;
                    J(2*n+1,:) = [D(1,:), 0];
            end
            
            duT = -(J\F);
            % disp([norm(duT, inf), norm(F, inf), T])
            u = u + duT(1:end-1);
            T = T + duT(end);
        
            if ( norm(duT, inf) < tol ), break, end
        end
        
        % Compute errors
        switch disc
            case 'trig'
                uu = chebfun(u, [0 T], 'periodic');
            case 'cheb'
                uu = chebfun(u, [0 abs(T)]);
        end
        tt = linspace(0, T, 10001).';
        if ( n == nproxy )
            uproxy = uu;
            Tproxy = T;
        else
            err(n) = norm(uproxy(tt)-uu(tt), 2);
            res = deriv(uu, tt) - (lam-uu(mod(tt-1,T))).*uu(tt);
        end
    
    end

    nexttile(3)
    semilogy(nn, err(nn));
    ylim([1e-12, 1e2])
    grid on, hold on


end

%%

nexttile(3)
h1 = legend('Cheb', 'Trig', 'fontsize', 16, 'interpreter', 'latex');
xlabel('$n$','interp', 'latex', 'fontsize', 16);
ylabel('error','interp', 'latex', 'fontsize', 16, 'ver', 'mid');
hold on, plot([0 60], 5e-9*[1 1], ':k', 'linewidth', 2); hold off
grid on
set(h1, 'String', {'Cheb', 'Trig', '$5e-9$'})
ylim([1e-11, 1e1]), grid on, grid minor
set(gca, 'xtick', 0:10:60)
set(gca, 'ytick', 10.^(-15:5:0))
drawnow, shg
set(gca, 'fontsize', 16)


nexttile(1)
plot(uu)
title(['period: T = ' num2str(T,10)])
xlabel('$t$','interp', 'latex', 'fontsize', 16); 
ylabel('$u(t)$','interp', 'latex', 'fontsize', 16)
drawnow, shg
set(gca, 'fontsize', 16)
ylim([0 3.5])


nexttile(2)
y = chebfun(@(t) deval(sol, t), [sol.x(1), sol.x(end)]);
cols = colororder;
h2 = plot(y, diff(y), ':k', uu, diff(uu), '-');
set(h2(2), 'color', cols(1,:))
axis([0.4 3.6 -3.5 2.5])
xlabel('$u(t)$','interp', 'latex', 'fontsize', 16); 
ylabel('$u''(t)$','interp', 'latex', 'fontsize', 16)
grid on, drawnow, shg
set(gca, 'fontsize', 16)

% print -depsc ../paper/figures/example13

%%

function P = Barymattrig(z, t)
    T = t(end)-2*t(1)+t(2);
    n = length(t);
    if ( mod(n, 2) )
        P = csc((z-t.')*pi/T);
    else
        P = cot((z-t.')*pi/T);
    end
    P(:,2:2:end) = -P(:,2:2:end);
    P = P./sum(P,2);
    P(isnan(P)) = 1;
end
