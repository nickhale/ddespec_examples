% CHEBFUN EXAMPLE 2 from "A spectral collocation method for functional and
% delay differential equations". 
% Nick Hale - August 2024

close all, figure(1)
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

for s = [2, -lambertw(-2/exp(2))]
    dom = [0, .1];
    N = chebop(@(t,y) diff(y) - (2*cos(2*t)*y(t/2)^(2*cos(t)) + ...
                       log(deriv(y, t/2)) - log(2*cos(t)) - sin(t)), dom);
    N.lbc = 1;
    N.init = chebfun(@(t) 1+s*t, dom);
    y = N\0;
    nexttile, plot(y), title(['s = ' num2str(s, 14)])
end

%%

figure(1), clf
set(gcf, 'position', [2158 385 1168 420])
tiles = tiledlayout(1,2,'TileSpacing','loose');

for s = [2, -lambertw(-2/exp(2))]  % Two consistent initial slopes
    dom = [0, .1];
    N = chebop(@(t,y) diff(y) - (2*cos(2*t)*y(t/2)^(2*cos(t)) + log(deriv(y, t/2)) - log(2*cos(t)) - sin(t)), dom);
    N.lbc = 1;
    N.init = chebfun(@(t) 1+s*t, dom);
    y = N\0;

    nexttile
    plot(y), hold on
    sol = myddex5(s); plot(sol.x, sol.y, ':'); hold off
    title(['s = ' num2str(s, 14)]), shg
    legend('Chebfun', 'ddensd', 'location', 'nw')
    title(['$s$ = ' num2str(s, 14)], 'interp', 'latex'), shg
    xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
    set(gca, 'fontsize', 16)

    if ( s == 2 )
        err1 = norm(y.values-exp(sin(2*y.points)), inf)./norm(exp(sin(2*y.points)), inf)
        err2 = norm(sol.y-exp(sin(2*sol.x)), inf)./norm(exp(sin(2*sol.x)), inf)
    else
        err1 = norm(N(y), inf)
        sol = chebfun(@(t) deval(sol, t), dom);
        err2 = norm(N(sol), inf)
    end 
end

% print -depsc2 ../paper/figures/chebxample2

function sol = myddex5(s)
    % Anonymous function calculating the delay
    delay = @(t,y) t/2;
    y0 = 1;
    t0 = 0;
    tf = 0.1;
    tspan = [t0, tf];
    % Solve as initial-value neutral DDE.
    opts = ddeset('abstol', 1e-5, 'reltol', 1e-5);
    sol = ddensd(@ddex5de,delay,delay,{y0,s},tspan, opts);

end % ddex5
function dydt = ddex5de(t,y,ydel,ypdel)
% Differential equation function for DDEX5.
dydt = 2*cos(2*t)*ydel^(2*cos(t)) + log(ypdel) - log(2*cos(t)) - sin(t);
end % ddex5de