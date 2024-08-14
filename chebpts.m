function [t,w,v] = chebpts(n, dom)

if ( n == 1 ) % Special case (single point)
    t = 0; w = 2; v = 1;
else          % General case
    
    % Chebyshev points:
    m = n - 1;
    t = sin(pi*(-m:2:m)/(2*m)).';  % (Use of sine enforces symmetry.)
    % Quadrature weights
    c = 2./[1, 1-(2:2:(n-1)).^2];  % Exact integrals of T_k (even)
    c = [c, c(floor(n/2):-1:2)];   % Mirror for DCT via FFT
    w = ifft(c);                   % Interior weights
    w([1,n]) = w(1)/2;             % Boundary weights
    % Barycentric weights
    v = [ones(n-1,1) ; .5];        % Note v(end) is positive.
    v(end-1:-2:1) = -1; 
    v(1) = .5*v(1);
    
end

if ( nargin == 2 )
    t = dom(2)*(t + 1)/2 + dom(1)*(1 - t)/2;
    w = (diff(dom)/2)*w;
end

end
