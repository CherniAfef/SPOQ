function [lpx] = Lpsmooth(x,alpha,p)
    % This function computes the smooth Lp norm of the vector x
    lpx = (sum((x.^2 + alpha.^2).^(p/2) - alpha.^p))^(1/p);

end