function [norm2K] = norm2(K,N)
    nbiter = 50;
    b = rand(N,1);
    for i=1:nbiter
        tmp = K'*(K*b);
        %tmpnorm = sqrt(dot(tmp,tmp));
        tmpnorm = norm(tmp,2);
        b = tmp./tmpnorm;
    end
    norm2K = norm(K*b,2);
end

