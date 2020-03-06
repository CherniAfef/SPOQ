function [zk,j] = proxPPXAplus(D,B,x,y,eta,J,prec)
%%% This function computes the proximity operator 
%%% using the PPXA+ algorithm
    [~,N] = size(D);
    x1k_old = x;
    x2k_old = D*x1k_old;
    A = inv(eye(N) + D'*D);
    zk_old = A*(x1k_old + D'*x2k_old);
    teta = 1.9;
    for j=1:J
        y1k_old = proxB(B,x1k_old,x,teta);        
        y2k_old = proxl2(x2k_old,y,eta);
        vk_old = A*(y1k_old + D'*y2k_old);
        x1k = x1k_old + 2*vk_old - zk_old - y1k_old;
        x2k = x2k_old + D*(2*vk_old - zk_old) - y2k_old;
        zk = vk_old;
        error = norm(zk-zk_old,2)^2;
        if error < prec
            disp(['PPXA stops at j = ',num2str(j)])
            break;
        end
        x1k_old = x1k;
        x2k_old = x2k;
        zk_old = zk;      
    end
end