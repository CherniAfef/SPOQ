function [lqx] = Lqsmooth(x,mu,q)
    % This function computes the smooth Lq norm of the vector x
    lqx = (mu.^q + sum(abs(x).^q))^(1/q);
end