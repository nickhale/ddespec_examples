% EXAMPLE 4 from "A spectral collocation method for functional and delay
% differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

% Solution
p = 0.5;

% Discretisation sizes:
n = 10:13;
T = 0:.5:2;

% Initialise the block operators:
AA = {};
for j = 1:4
    for k = 1:4
        AA{j,k} = zeros(n(j),n(k));
    end
end

% Form the operator on each subdomain:
t = {}; w = {};
for k = 1:4
    domk = T(k:k+1);
    nk = n(k);
    [t{k}, ~, w{k}] = chebpts(nk, domk);
    D = Diffmat(t{k}, w{k});  
    I = eye(nk);
    A = D + I;
    A(1,:) = I(1,:); % Boundary condition
    AA{k,k} = A;
    if ( k > 1 )     % Interpolation condition
        P = Barymat(t{k}-1/2, t{k-1}, w{k-1});
        P(1,:) = -fliplr(eye(1,n(k-1)));
        AA{k,k-1} = P;
    end
end

% Collate:
AA = cell2mat(AA);
b = zeros(size(AA,2),1);
b(1) = 1;
% Solve:
y = AA\b;

% The global solution:
sol = chebfun({@(t) exp(-t),  
    @(t) (exp(-t)*(exp(1/2) - 2*t*exp(1/2) + 2))/2, 
    @(t) exp(-t)*(exp(1)/2 + exp(1/2)/2 + 1) + exp(-t)*((exp(1)*t^2)/2 + (- exp(1) - exp(1/2))*t),
    @(t) exp(-t)*(exp(1)/2 + exp(1/2)/2 + (9*exp(3/2))/16 + 1) - exp(-t)*((exp(3/2)*t^3)/6 + (- exp(1)/2 - (3*exp(3/2))/4)*t^2 + (exp(1) + exp(1/2) + (9*exp(3/2))/8)*t)}, 0:.5:2);

% Plot:
nexttile
plot(sol, 'k'); hold on
N = [0 cumsum(n)];
Y = {};
for k = 1:4
    Y{k} = y(N(k)+1:N(k+1));
    plot(t{k}, Y{k}, '.', 'markersize', 25); hold on
end
hold off
xlim([0 2]), xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

err = norm(y - sol(vertcat(t{:})), inf)

%% Convergence

nn = 1:20;
err = []; CN = [];
T = 0:.5:2;
for n1 = nn 
    % Discretisation sizes:
    n = n1:n1+3;
    % Initialise the block operators:
    AA = {};
    for j = 1:4
        for k = 1:4
            AA{j,k} = zeros(n(j),n(k));
        end
    end
    % Form the operator on each subdomain:
    t = {};
    for k = 1:4
        domk = T(k:k+1);
        nk = n(k);
        [t{k}, ~, w{k}] = chebpts(nk, domk);
        D = Diffmat(t{k}, w{k});  
        I = eye(nk);
        A = D + I;
        A(1,:) = I(1,:);
        AA{k,k} = A;
        if ( k > 1 )
            P = Barymat(t{k}-1/2, t{k-1}, w{k-1});
            P(1,:) = -fliplr(eye(1,n(k-1)));
            AA{k,k-1} = P;
        end
    end
    % Collate:
    AA = cell2mat(AA);
    b = zeros(size(AA,2),1);
    b(1) = 1;
    % Solve:
    y = AA\b;

    err(n1) = norm(y - sol(vertcat(t{:})));
    CN(n1) = cond(AA);

end

nexttile
semilogy(nn, err(nn)); shg
hold on, plot(nn, eps*CN(nn), ':', 'color', .6*[1 1 1]); hold off
xlabel('$n$', 'interp', 'latex'); ylabel('error'); grid on
set(gca, 'fontsize', 16)
set(gca, 'ytick', 10.^(-15:5:0))

% print -depsc2 ../paper/figures/example4

%% Chebfun solution:
 
N = chebop(@(t, y) diff(y) + y + y(max(t-1/2,0)), 0:.5:2);
N.lbc = 1;
t = chebfun('t', 0:.5:2);
f = (t < .5); % Chebfun assumes y(t<0) = y(0). Adjust RHS accordingly.
y = N\f;

figure(2)
plot(y, 'b', sol, ':r'), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)
err = norm(y - sol, inf)

