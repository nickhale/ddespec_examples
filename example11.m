% EXAMPLE 11 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024
%
% Accurately compute the limit cycle of a non-dimensionalised 
% Lotka-Volterra Predator-Prey model with a delayed Holling type-II 
% response using a Fourier spectral method.
%
% Summary:
% * Obtain an initial guess for the limit cycle and period using dde23.
% * Solve for x(t) and y(t) using a Fourier spectral method. Since the period 
% T is unknown, we rescale time to [0 1] and solve for T as part of the problem. 
% To square the system we impose that y'(0) = 0. 

clf, close all


%% 
% *Model parameters* (Can be changed)

gam = 2/15;     % Predator death rate
del = 1;        % Predation rate
K = 7/5;        % Prey carrying capacity
s = 1;          % Delay

% Known values/points:
alp = gam/del;
xe = alp/(1-alp); ye = (1-xe/K)*(1+xe); % (critical point)

% Discretisation parameters:
N = 129;        % Size of Fourier discretisation.

%% 
% *DDE23 solve*
% 
% Solver settings:

opts_ode = odeset('abstol', 1e-6, 'reltol', 1e-6);
x0 = xe; 
y0 = 1.5;   % Estimate for limit cycle height above xs. (xs, yt) is used 
            % as the initial condition. Getting this accurately is not very
            % important, it just saves time in integrating too long to get
            % to the limit cycle.
Tmax = 100; % How long are we willing to integrate in order to get to  
            % the limit cycle.

%% 
% Define the DDE and solve:

if ( s > 0 )
    f = @(t, x, z) [x(1) - x(1)^2/K - x(1)*z(2)/(1+x(1)) ;
                    -gam*x(2) +  del*x(1)*z(2)/(1+x(1))] ;
    sol = dde23(f, s, [x0, y0], [0 Tmax], opts_ode);
else
    f = @(t, x) [x(1) - x(1)^2/K - x(1)*x(2)/(1+x(1)) ;
                 -gam*x(2) +  del*x(1)*x(2)/(1+x(1))] ;
    sol = ode15s(f, [0 Tmax], [x0, y0], opts_ode);
end


%% 
% Plot the DDE solution:

figure(1)
plot(sol.x, sol.y)
figure(2)
plot(sol.y(1,:), sol.y(2,:))
hold on, plot(sol.y(1,1), sol.y(2,1), 'og', 'markerfacecolor', 'g')
plot(sol.y(1,end), sol.y(2,end), 'or', 'markerfacecolor', 'r'), hold off

%% 
% *ESTIMATE A SINGLE PERIOD:*

lm = find(islocalmin(-sol.y(2,:)));
dom = [sol.x(lm(end-2)), sol.x(lm(end-1))];
T = diff(dom);
f = 1/T;

%% 
% *FOURIER DISCRETISATION:*

t = linspace(dom(1), dom(2), N+1)'; t(end) = [];
xy = deval(sol, t).';

%% 
% Shift time to [0, 1] for convenience. We must solve for the period, T.

t = linspace(0, 1, N+1)'; t(end) = [];
D = diffmattrig(N, [0 1]);
I = eye(N);

%% 
% *NEWTON ITERATION:*

xy = [xy(:) ; f];           % Initial guess
% disp('Newton iteration:')   % Let's roll
for k = 1:10

    % Extract solution components:
    x = xy(1:N);
    y = xy((N+1):2*N);
    f = xy(end);                 % Time scaling. a = 1/T.

    P = barymattrig(t-f*s, t);   % Update the interpolation operator
    % (Since we scale time, we must also scale s.)

    z = P*y; Dx = D*x; Dy = D*y; Dz = D*z; % Define for convenience/speed

    % RHS:
    F1 = f*Dx - x + 1/K*x.^2 + x.*z./(1+x);
    F2 = f*Dy + gam*y - del*x.*z./(1+x);
    F = [F1 ; F2];

    % Jacobian
    J1x = f*D - I + diag(2*x/K + z./(1+x) - x.*z./(1+x).^2);
    J1y = diag(x./(1+x))*P;
    J2x = -del*diag(z./(1+x) - x.*z./(1+x).^2);
    J2y = f*D + gam*I - del*diag(x./(1+x))*P;
    J1a = Dx - s*x.*Dz./(1+x); 
    J2a = Dy + del*s*x.*Dz./(1+x);
    J = [J1x J1y J1a ; J2x J2y J2a];

    % % Additional boundary constraint to square the system:
    F(2*N+1) = D(1,:)*y;
    J(2*N+1,:) = [zeros(1,N), D(1,:), 0];

    % Newton solve:
    dxy = -J\F;
    xy = xy + dxy;

    fprintf('Iteration %d errors = %8.8g, %8.8g.\n', k, norm(F, inf), norm(dxy, inf))
    if ( norm(dxy, inf) < 1e-12 ), break, end

end

%% 
% Extract solutions:

x = xy(1:N);
y = xy((N+1):2*N);
f = xy(end);
T = 1/f
t = t/f;

%% 
% *PLOTTING:*
close all, figure(1)
set(gcf, 'position', [2158 385 1680 420])
tiles = tiledlayout(1,3,'TileSpacing','loose');

nexttile(3)
FourierCoeffs = conj(fftshift(ifft([x y]), 1));
N = length(FourierCoeffs);
% semilogy(-N/2+1/2:N/2, abs(FourierCoeffs), '.', 'markersize', 16)
semilogy(-N/2+1/2:N/2, abs(FourierCoeffs(:,1)), 'o', 'markersize', 7, 'linewidth', 2); hold on
semilogy(-N/2+1/2:N/2, abs(FourierCoeffs(:,2)), 'x', 'markersize', 7, 'linewidth', 2); hold off
% title('Magnitude of Fourier modes', 'fontsize', 14)
xlim([-65 65]),  ylim([1e-16, 1e0])
xlabel('$k$','interp', 'latex', 'fontsize', 16);
ylabel('$|c|_k$','interp', 'latex', 'fontsize', 16);
legend('$x$','$y$','interpreter', 'latex', 'fontsize', 16);
grid on
drawnow, shg
set(gca, 'fontsize', 16)
pause(eps)

nexttile
plot(t, x, t, y)
% title('Limit cycle solution', 'fontsize', 14)
legend('$x$','$y$','interpreter', 'latex', 'fontsize', 16);
xlabel('$t$','interp', 'latex', 'fontsize', 16);
% grid on
set(gca, 'fontsize', 16)


nexttile
plot([x;x(1)], [y;y(1)]), 
% title('Phase plane' , 'fontsize', 14)
xlabel('$x$','interp', 'latex', 'fontsize', 16);
ylabel('$y$','interp', 'latex', 'fontsize', 16);
xlim([-.1 1.1]), ylim([0 2.4]), grid on

% Add the DDE23 phase plane to the one we computed for visual verification:
nexttile(2)
hold on
plot(sol.y(1,:), sol.y(2,:), ':k'); 
% plot(sol.y(1,1), sol.y(2,1), 'og'); 
% plot(xe, ye,'.b'), 
cols = colororder; blue = cols(1,:); 
plot(xe, ye, 'color', blue);
set(gca, 'fontsize', 16)
hold off

% print -depsc ../paper/figures/example11

%% 
% *SANITY CHECK*

% Check residual using Chebfun
xc = chebfun(x, [0 T], 'trig'); 
yc = chebfun(y, [0 T], 'trig');
tc = chebfun('t', xc.domain);
f1 = diff(xc) - xc + xc.^2/K + xc.*yc(mod(tc-s, T))./(1+xc);
f2 = diff(yc) + gam*yc - del*xc.*yc(mod(tc-s, T))./(1+xc);
res = norm([f1 ; f2], inf)

% 
%% 
% *BARYMAT AND DIFFMAT FUNCTIONS*

function P = barymattrig(z, t)
    T = t(end)-2*t(1)+t(2);
    N = length(t);
    if ( mod(N, 2) )
        P = csc((z-t.')*pi/T);
    else
        P = cot((z-t.')*pi/T);
    end
    P(:,2:2:end) = -P(:,2:2:end);
    P = P./sum(P,2);
    P(isnan(P)) = 1;
end

function D = diffmattrig(N, dom) % From Chebfun.
    h = 2*pi/N;
    if ( mod(N, 2) ) % N is odd.
        column = [0, .5*csc((1:N-1)*h/2)]';
    else % N is even.
        column = [0, .5*cot((1:N-1)*h/2)]';
    end
    column(2:2:end) = -column(2:2:end);
    row = column([1 N:-1:2]);
    D = toeplitz(column, row);
    D = 2*pi*D/diff(dom);
end

%% 
% 

% CHEB  compute D = differentiation matrix, x = Chebyshev grid

function [D,x] = cheb(N)
    if N==0, D=0; x=1; return, end
    x = cos(pi*(0:N)/N)'; 
    c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X = repmat(x,1,N+1);
    dX = X-X';                  
    D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
    D  = D - diag(sum(D'));                 % diagonal entries
end