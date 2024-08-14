% CHEBFUN EXAMPLE 1 from "A spectral collocation method for functional and
% delay differential equations". 
% Nick Hale - August 2024

N = chebop(@(t,y) diff(y) + y(y), [0 1]);
N.lbc = 1;
u = N\0;

close all, plot(u), shg
title(['$s$ = ' num2str(s, 14)], 'interp', 'latex'), shg
xlabel('$t$', 'interp', 'latex'); ylabel('$u$', 'interp', 'latex');
set(gca, 'fontsize', 16)

% print -depsc2 ../paper/figures/chebxample3 

return

%% Variant 1

N = chebop(@(t,y) diff(y) + y(y^y), [0 1]);
N.lbc = 1;
N.init = chebfun('1-t/2', [0 1]);
u = N\0;
plot(u), shg

%% Variant 2

N = chebop(@(t,y) diff(y, 2) + y(y), [0 1]);
N.lbc = @(y) y - 1;
N.rbc = @(y) y - .1;
t = chebfun('t', [0 1]);
N.init = chebfun(@(t) 1-t/2, [0 1])
u = N\0;
plot(u), shg
