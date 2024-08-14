function D = diffmat(t, w)  
    D = w.'./(w.*(t-t.'));
    n = length(t); ii = 1:n+1:n^2;       
    D(ii) = 0;  D(ii) = -sum(D,2); 
end