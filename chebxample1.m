% CHEBFUN EXAMPLE 1 from "A spectral collocation method for functional and
% delay differential equations". 
% Nick Hale - August 2024

dom = [0 20];
q = 1/2;
N = chebop(@(x,y) diff(y) - 1/100*(q*x-x-10)*y(q*x) - 1/100*(x+20)*exp(-1) - ...
   1/100*cumsum(y) - 1/1000*feval(volt(@(t,s)(t/q-s), y), q*x), dom);
N.lbc = exp(-1);
y = N\0;
sol = chebfun(@(x) exp(x/10-1), dom);
err = norm(sol-y, inf)./norm(sol, inf)

close all
plot(sol, '-k'); hold on
plot(y, '.', 'markersize', 30); hold off
title(['relative error = ' num2str(err)])
% legend('Exact solution', 'Chebfun solution', 'location', 'nw')
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16);

% print -depsc2 ../paper/figures/chebxample1