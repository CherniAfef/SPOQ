function [fcost]=Fcost(x,alpha,beta,mu,p,q)
    lp = (sum((x.^2 + alpha.^2).^(p/2)) - alpha.^p).^(1/q);
    lq = (mu.^q + sum(abs(x).^q)).^(1/q);
    fcost = log(((lp.^p + beta.^p).^(1/q)) / lq);
end
