% CHEBFUN EXAMPLE 4 from "A spectral collocation method for functional and
% delay differential equations". 
% Nick Hale - August 2024

y0 = 1;
y1 = .25;

N = chebop(@(t,y,p) diff(y) + y + y(p*t) - exp(-t/2), [0 1]);
N.lbc = @(y,p) y - y0;
N.rbc = @(y,p) y - y1;
yp = N\0;

close all, plot(yp{1}), title(['p = ' num2str(yp{2}, 10)]), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

% print -depsc2 ../paper/figures/chebxample4
