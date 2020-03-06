function [xk,fcost,Bwhile,Time,mysnr]=FB_PPXALpLq(K,y,p,q,metric,alpha,beta,eta,xi,nbiter,xtrue)
%%% This function defines the Trust region algorihtm based on
%%% Forward-Backward algorithm
    
    %Initialization
    [~,N] = size(K);
    [xk_old, ~] = pds(K,y,xi,10);
    mysnr(1) = -10*log10(sum((xk_old-xtrue).^2) / sum(xtrue.^2));
    fcost(1) = Fcost(xk_old,alpha,beta,eta,p,q);
    gamma = 1;
    prec = 1e-12;
    Bwhile = zeros(nbiter,1);
    fcost = zeros(nbiter,1);
    J = 5000; %ppxa max iterations
    %bmax = 10;% maximum TR iterates
    %metric 0: Lip constant, 1: FBVM without TR, 2: FBVM-TR
    L = ComputeLipschitz(alpha,beta,eta,p,q,N);
    %Algorithm   
    for k=1:nbiter
        if (rem(k,100)==0)
            disp(['it=',num2str(k),': fcost = ',num2str(fcost(k-1))])
        end
        tic;
        switch metric
            case 0         
                A = L*ones(N,1);
                B = A./gamma;
                xxk = xk_old - (1./B).*gradlplq(xk_old,alpha,beta,eta,p,q);
                xk = proxPPXAplus(K,B,xxk,y,xi,J,prec);
                
            case 1  
                A = condlplq(xk_old,alpha,beta,eta,p,q,0);               
                B = A./gamma;
                xxk = xk_old - (1./B).*gradlplq(xk_old,alpha,beta,eta,p,q);
                xk = proxPPXAplus(K,B,xxk,y,xi,J,prec);
                
            case 2
                ro = sum(abs(xk_old.^q))^(1/q); 
                bwhile = 0;
                while 1 
                    A = condlplq(xk_old,alpha,beta,eta,p,q,ro);
                    B = A/gamma;
                    xxk = xk_old - (1./B).*gradlplq(xk_old,alpha,beta,eta,p,q);                                       
                    xk = proxPPXAplus(K,B,xxk,y,xi,J,prec);

                    if sum(abs(xk).^q)^(1/q) < ro 
                        ro = ro/2;
                        bwhile = bwhile + 1;
                    else
                        break;
                    end
                end
            Bwhile(k)= bwhile;
        end
        
        Time(k) = toc;
        error = norm(xk-xk_old)^2 / norm(xk_old)^2;
        mysnr(k+1) = -10*log10(sum((xk-xtrue).^2)/sum(xtrue.^2));
        fcost(k+1) = Fcost(xk,alpha,beta,eta,p,q);

        if error < prec
            break;
        end
        xk_old = xk;
     end
end