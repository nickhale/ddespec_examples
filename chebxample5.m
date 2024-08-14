% CHEBFUN EXAMPLE 5 from "A spectral collocation method for functional and
% delay differential equations". 
% Nick Hale - August 2024

A = chebop(@(y) diff(y,2), [0 1]);
A.bc = 'dirichlet';
B = chebop(@(t,y) -y(t/2), [0 1]);
[V, D] = eigs(A, B);
fprintf('%10.10g\n', diag(D))

close all, plot(V), 
xlabel('$t$', 'interp', 'latex'); ylabel('$y$', 'interp', 'latex');
set(gca, 'fontsize', 16)

% print -depsc2 ../paper/figures/chebxample5