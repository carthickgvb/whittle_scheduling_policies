function out = mu(k,l)
    %out = 1/(M+1);
    %e = 2.71828;
    %lam = floor(M/2);
    %out = (e^-lam)*(lam^k)/factorial(k);
    out = poisspdf(k,l);
end