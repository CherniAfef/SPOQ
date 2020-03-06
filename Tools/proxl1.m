function [p] = proxl1(x,w)
    % proximity operator of l1 norm: Thresholding
    %y = max(abs(x)-w,0).*sign(x);
    p = zeros(size(x));
    pos = find(x>w);
    p(pos) = x(pos) - w;
    neg = find(x<-w);
    p(neg) = x(neg) +w;
    p(p<0)=0;
end