function [A] = condlplq(x,alpha,beta,eta,p,q,ro)
%%% This function computes the metric matrix for the variable metric
%%% Forward-Backward algorithm
    lp = Lpsmooth(x,alpha,p);
    Xpq = (q-1)/((eta^q + ro^q)^(2/q));
    A = Xpq + (1/(lp^p+beta^p))*(x.^2+alpha^2).^(p/2-1);
end