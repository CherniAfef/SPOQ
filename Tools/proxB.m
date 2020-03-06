function [p] = proxB(B,x,xhat,teta)
%%% This function computes the proximity operators of 
%%% f(x) = (teta/2) * ||y-x||_B^2
    p = (x+teta*(B.*xhat))./(1+teta*B);
    p(p<0) = 0;
end