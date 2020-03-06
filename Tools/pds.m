function [xk, refspec]=pds(K,y,eta,nbiter)
    [M,N] = size(K);
    normK = norm2(K, N);
    tau = 1/normK;
    sigma = 0.9/(tau*normK^2);
    ro = 1.;
    refspec = zeros(nbiter,1);
    xk_old = ones(N,1);
    uk_old = K*xk_old;
    prec = 1e-6;
    for i=1:nbiter
        xxk = proxl1(xk_old - tau*K'*uk_old, tau);
        zk = uk_old + sigma*K*(2*xxk-xk_old);
        uuk = zk - sigma*proxl2(zk/sigma, y, eta);
        xk = xk_old + ro*(xxk - xk_old);
        uk = uk_old + ro*(uuk - uk_old);
        ex = norm(xk-xk_old,2)^2 / norm(xk,2)^2;
        eu = norm(uk-uk_old,2)^2 / norm(uk,2)^2;
        if ex < prec && eu < prec
            break;
        end
        refspec(i) = ex;
        xk_old = xk;
        uk_old = uk;
    end
end