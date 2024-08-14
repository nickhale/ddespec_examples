function P = barymat(tau, t, w)
    P = w.'./(tau-t.');
    P = P./sum(P, 2);
    P(isnan(P)) =  1;
end